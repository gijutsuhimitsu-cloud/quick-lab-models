# SPDX-License-Identifier: MIT
"""基礎データ（base_data/）を、配布元から自分で取得した元のファイルから作り直す（任意。約 10 分。SMRT の 8 万分子の特徴量の計算が大半）。
  python train/scripts/build_from_raw.py --raw <元のファイルのフォルダ> [--out train/work/base_from_raw]
  そのあと: python train/scripts/compare_base_data.py train/work/base_from_raw   （同梱の base_data/ と比べる）
元のファイルの置き方（取得の手順は README.md の「元のデータから作り直す」）:
  <raw>/smrt_figshare/SMRT_dataset.csv・SMRT_dataset.sdf
  <raw>/hsm2_mendeley/HSM2_Database.xlsx
  <raw>/weber_isocratic_mendeley/RetentionData
  <raw>/mcart_msd_zenodo/enamine_480.csv
  <raw>/wiczling_osf/<zqtj7・qbv7j・2mcnw・y6s8p・evuj9>/（OSF の各記録の zip を展開した物。1-*.csv・4-compounds-names.csv・5-smiles.csv）
  <raw>/zhang2024_scidb/Dataset NN.xlsx（学習に使うのは下の ZHANG_COLUMN にある番号。Dataset 10 は評価専用なので読まない）
カラムは、製品名ではなく weights/columns.json の id（HSM2 データベースの番号）で対応づけてある（下の表）。
同梱の表（こちらで作った物）はそのまま使う: names_manual.tsv・names_weber.tsv（名前 → 構造）、names_exclude.tsv、retention_form.json。
"""
import os, sys, re, json, glob, argparse, shutil
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from lib import common as C
from lib.table import feature_names

# ---- カラムの対応（weights/columns.json の id。match: exact＝同じ製品、approx＝同じ固定相の粒径違いなどの値を使った）
SMRT_COLUMN = ("hsm2_5", "exact")
WEBER_COLUMN = ("hsm2_404", "approx")
WICZLING_COLUMN = {"zqtj7": ("hsm2_406", "exact"), "qbv7j": ("hsm2_303", "exact"), "2mcnw": ("hsm2_304", "exact"),
                   "y6s8p": ("hsm2_484", "exact"), "evuj9": ("hsm2_405", "exact")}           # この順に読む（学習の並び）
_Z = {"hsm2_521": "exact", "hsm2_267": "exact", "hsm2_93": "approx", "hsm2_404": "approx"}
ZHANG_COLUMN = {n: (c, _Z[c]) for n, c in {
    "01": "hsm2_521", "02": "hsm2_93", "04": "hsm2_267", "05": "hsm2_404", "06": "hsm2_521", "07": "hsm2_93", "08": "hsm2_93", "09": "hsm2_267",
    "11": "hsm2_404", "12": "hsm2_93", "14": "hsm2_521", "15": "hsm2_93", "16": "hsm2_521", "17": "hsm2_267", "18": "hsm2_267", "19": "hsm2_93",
    "20": "hsm2_521", "21": "hsm2_267", "23": "hsm2_521", "24": "hsm2_93", "25": "hsm2_521", "26": "hsm2_521", "27": "hsm2_521", "28": "hsm2_521",
    "29": "hsm2_93", "30": "hsm2_93"}.items()}
# MCaRT の 4 カラムは銘柄の記載が無いので、HSM2 の同じ固定相の中央値（列の名前の先頭 → 固定相、死時間）
MCART_COLUMN = {"C18": ("C18", 0.63, "C18（HSM2 の C18 の中央値。銘柄は配布元に記載なし）"),
                "CN": ("CN", 0.73, "シアノ（HSM2 のシアノの中央値。銘柄は配布元に記載なし）"),
                "Phen": ("Phenyl", 0.83, "フェニル（HSM2 のフェニルの中央値。銘柄は配布元に記載なし）"),
                "AQ": ("C18", 0.83, "C18 の水系対応型（HSM2 の C18 の中央値。銘柄は配布元に記載なし）")}
HSM_COLS = ["H", "A", "B", "C", "V", "D"]
NOTE = {
    "smrt_figshare": "METLIN SMRT。カラムの種類（C18、50 × 2.1 mm、1.8 µm）・勾配・流量・0.1% ギ酸は、SMRT の論文（Domingo-Almenara ほか, Nat. Commun. 10, 5811 (2019), doi:10.1038/s41467-019-13680-7）に書かれた測定条件（事実の記載。値は、その条件を記録した公開のメタデータから写した。論文の本文との 1 項目ずつの照合はしていない）。"
                     "t0 = 1.35 分は保持されない化合物の山（1.3〜1.7 分）から、装置の遅れ 3.8 分は溶出の始まり（約 8.2 分）からのこちらの見積もり。pH 3.0 は添加剤からの見積もり",
    "hsm2_mendeley": "HSM2 データベースのカラムのパラメータと試験溶質のパラメータから再構成した log k（独立の実測ではない）。条件は HSM の標準条件（50% アセトニトリル／リン酸緩衝液 pH 2.8、35 °C）とみなした",
    "weber_isocratic_mendeley": "アセトニトリル／10 mM リン酸（pH 2.1 は見積もり）。割合は体積分率の列を使った。同じ条件の繰り返しは log k の平均",
    "mcart_msd_zenodo": "5 分の 5〜85% アセトニトリルの勾配（配布元の README）。t0 は「調整前 − 調整後」の差。pH は緩衝液からの見積もり。印（*）の付いた値は外した",
    "wiczling_osf": "勾配の遅れ td・死時間 to は配布元の記載の値（遅れは勾配の表に入れてある）。保持時間は RT − te（カラム外の時間）。pH は公称値",
}
# SMRT の測定条件（論文 doi:10.1038/s41467-019-13680-7 に書かれた事実。値はその条件を記録した公開のメタデータから写し、論文の本文との照合はしていない。カラムは C18 50 × 2.1 mm・1.8 µm、0.1% ギ酸の水／アセトニトリル、
# 勾配 0〜3 分 5%、5 分 50%、20 分 85%、23 分 85%。死時間 1.35 分と装置の遅れ 3.8 分は、こちらがデータの分布から見積もった値）
SMRT_METHOD = dict(length_mm=50, id_mm=2.1, particle_um=1.8, temp_C=np.nan, flow_mL_min=0.1, t0_min=1.35, pH_weak=3.0, frac_acn=1.0,
                   add_formic=1, add_acetic=0, add_tfa=0, add_ammonium=0, add_basic=0, add_phosphate=0,
                   gradient=json.dumps([[0, 0.05], [3, 0.05], [5, 0.5], [20, 0.85], [23, 0.85]]), isocratic_phi=np.nan)


def phase_norm(p):
    p = str(p).strip()
    if p in ("C18", "C30", "C28"):
        return "C18"
    if p in ("C8", "CN", "EP", "F"):
        return p
    if p.lower() == "phenyl":
        return "Phenyl"
    if p in ("C4", "C5", "C3", "C1"):
        return "short"
    return "other"


def featurize(name, rows, work, sdf=None):
    """rows: [(id, smiles)] → {id: 特徴量}"""
    out = os.path.join(work, name + ".jsonl")
    if sdf:
        C.node("featurize.cjs", "sdf", sdf, out, "PUBCHEM_COMPOUND_CID")
    else:
        p = os.path.join(work, name + ".tsv")
        with open(p, "w", encoding="utf-8", newline="\n") as f:
            for i, s in rows:
                f.write("%s\t%s\n" % (i, s))
        C.node("featurize.cjs", "tsv", p, out)
    res = {}
    for line in open(out, encoding="utf-8"):
        r = json.loads(line)
        if r.get("ok") and r.get("inchikey"):
            res[r["id"]] = r
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", default=os.path.join(C.WORK, "base_from_raw"))
    A = ap.parse_args()
    RAW, OUT = os.path.abspath(A.raw), os.path.abspath(A.out)
    work = os.path.join(OUT, "_features")
    os.makedirs(work, exist_ok=True)
    rp = lambda *a: os.path.join(RAW, *a)
    desc, grp, ion = feature_names()

    # ---- 構造の特徴量（Node + RDKit.js）
    man = pd.read_csv(os.path.join(C.BASE, "names_manual.tsv"), sep="\t", dtype=str, keep_default_na=False)
    web = pd.read_csv(os.path.join(C.BASE, "names_weber.tsv"), sep="\t", dtype=str, keep_default_na=False)
    mc = pd.read_csv(rp("mcart_msd_zenodo", "enamine_480.csv"))
    wz_dirs = sorted(d for d in glob.glob(rp("wiczling_osf", "*")) if os.path.exists(os.path.join(d, "5-smiles.csv")))
    wsm = pd.read_csv(os.path.join(wz_dirs[0], "5-smiles.csv"))
    zfiles = sorted(p for p in glob.glob(rp("zhang2024_scidb", "Dataset *.xlsx")) if os.path.basename(p)[8:10] in ZHANG_COLUMN)
    import openpyxl
    zsmi = []
    for p in zfiles:
        R = list(openpyxl.load_workbook(p, read_only=True, data_only=True)["RT"].iter_rows(values_only=True))
        hd = [("" if c is None else str(c).replace("\n", " ")) for c in R[0]]
        i_smi = next(i for i, h in enumerate(hd) if "SMILES" in h)
        zsmi += [str(r[i_smi]).replace("\xa0", "").strip() for r in R[1:] if r[i_smi]]
    zsmi = list(dict.fromkeys(zsmi))
    C.log("構造の特徴量を計算する（SMRT の SD ファイルが大半）")
    f_smrt = featurize("smrt", None, work, sdf=rp("smrt_figshare", "SMRT_dataset.sdf"))
    f_man = featurize("manual", list(zip(man["name"], man["smiles"])), work)
    f_web = featurize("weber", list(zip(web["id"], web["smiles"])), work)
    f_mc = featurize("mcart", [("mcart_%d" % i, s) for i, s in enumerate(mc["Smiles"])], work)
    f_wz = featurize("wiczling", [(str(int(a)), str(b).strip()) for a, b in zip(wsm["METID"], wsm["Smiles"])], work)
    f_zh = featurize("zhang", [(s, s) for s in zsmi], work)

    comp = {}

    def reg(r):
        ck = r["inchikey"][:14]
        if ck not in comp:
            comp[ck] = r
        return ck

    methods, rows = [], []
    hsm = pd.read_excel(rp("hsm2_mendeley", "HSM2_Database.xlsx"), sheet_name="Column Parameters").iloc[:, :15]
    hsm["Name"] = hsm["Name"].astype(str).str.strip()
    hsm["phase"] = hsm["Phase"].map(phase_norm)
    for c in HSM_COLS + ["kEB"]:
        hsm[c] = pd.to_numeric(hsm[c], errors="coerce")
    first = hsm.drop_duplicates("Name").set_index("Name")              # 同じ名前の行が 2 つある時は、先の行の値（学習の時と同じ）
    byid = {"hsm2_%s" % i: n for i, n in zip(hsm["ID"], hsm["Name"])}
    means = hsm.groupby("phase")[HSM_COLS + ["kEB"]].median()
    cols = {c["id"]: c["name"] for c in json.load(open(os.path.join(C.WEIGHTS, "columns.json"), encoding="utf-8"))["columns"]}

    def colrec(cid, match):
        row = first.loc[byid[cid]]
        return dict(column_id=cid, column_name=cols[cid], column_match=match, phase=row["phase"], **{c: float(row[c]) for c in HSM_COLS}, logkEB=float(np.log10(row["kEB"])))

    def add_method(src, mid, col, note, **kw):
        base = dict(method_id=mid, source=src, **col)
        base.update(kw)
        base.setdefault("fixed_cal", 0)
        base["exclude"] = 0
        base["note"] = note
        methods.append(base)

    # 1. SMRT（Apache 2.0）
    smrt = pd.read_csv(rp("smrt_figshare", "SMRT_dataset.csv"), sep=";")
    add_method("smrt_figshare", "smrt", colrec(*SMRT_COLUMN), NOTE["smrt_figshare"], **SMRT_METHOD)
    for cid, rt in zip(smrt["pubchem"].astype(str), smrt["rt"]):
        r = f_smrt.get(cid)
        if r is not None:
            rows.append(("smrt", reg(r), rt / 60.0, np.nan))
    # 2. HSM2（CC BY 4.0）: 551 カラム × 試験溶質 16、等組成 50% アセトニトリル（パラメータからの再構成）
    name2ck = {nm.lower(): reg(r) for nm, r in f_man.items()}
    sol = pd.read_excel(rp("hsm2_mendeley", "HSM2_Database.xlsx"), sheet_name="Solute Parameters")
    sol = sol.rename(columns={sol.columns[0]: "name"})
    sol_rows = [(name2ck[str(s["name"]).lower()], [float(s[c]) for c in ["h", "b", "a", "k", "v", "d"]]) for _, s in sol.iterrows() if str(s["name"]).lower() in name2ck]
    eb = name2ck["ethylbenzene"]
    for _, c in hsm.iterrows():
        if not np.isfinite(c["kEB"]) or c[HSM_COLS].isna().any():
            continue
        mid = "hsm2_%s" % c["ID"]
        col = colrec("hsm2_%s" % first.loc[c["Name"], "ID"], "exact")
        note = NOTE["hsm2_mendeley"]
        if col["column_id"] != mid:
            note += "。注意：このカラムは HSM2 の表で同じ名前の行が 2 つあり、学習ではカラムの性質の値（H〜D・logkEB）に、同じ名前の先の行（%s）の値が入っている（再現のためそのまま。保持の値はこのカラム自身の値から再構成）" % col["column_id"]
        col.update(column_id=mid, column_name=cols[mid], column_match="exact" if col["column_id"] == mid else "exact_same_name_first")
        add_method("hsm2_mendeley", mid, col, note, length_mm=np.nan, id_mm=np.nan, particle_um=np.nan, temp_C=35.0, flow_mL_min=np.nan, t0_min=np.nan,
                   pH_weak=2.8, frac_acn=1.0, add_formic=0, add_acetic=0, add_tfa=0, add_ammonium=0, add_basic=0, add_phosphate=1, gradient="", isocratic_phi=0.5)
        lk = float(np.log10(c["kEB"]))
        rows.append((mid, eb, np.nan, lk))
        colp = [float(c[x]) for x in HSM_COLS]
        for ck, sp in sol_rows:
            rows.append((mid, ck, np.nan, lk + float(np.dot(sp, colp))))
    # 2b. Weber（CC BY 4.0）: 105 化合物 × 割合 × 温度
    wd = pd.read_csv(rp("weber_isocratic_mendeley", "RetentionData"), header=None, names=["cid", "w", "phi", "temp", "k"])
    wd = wd[wd["k"] > 0]
    wd["logk"] = np.log10(wd["k"])
    wg = wd.groupby(["cid", "phi", "temp"], as_index=False)["logk"].mean()
    wcol = colrec(*WEBER_COLUMN)
    for (phi, temp), s in wg.groupby(["phi", "temp"]):
        mid = "weber_p%05.3f_T%d" % (phi, temp)
        add_method("weber_isocratic_mendeley", mid, dict(wcol), NOTE["weber_isocratic_mendeley"] + "。カラムの性質は、HSM2 の表で同じ固定相の粒径違いの製品の値（%s）を使った" % wcol["column_name"],
                   length_mm=50, id_mm=2.1, particle_um=1.7, temp_C=float(temp), flow_mL_min=0.35, t0_min=np.nan, pH_weak=2.1, frac_acn=1.0,
                   add_formic=0, add_acetic=0, add_tfa=0, add_ammonium=0, add_basic=0, add_phosphate=1, gradient="", isocratic_phi=float(phi))
        for cid, lk in zip(s["cid"], s["logk"]):
            r = f_web.get(str(int(cid)))
            if r is not None:
                rows.append((mid, reg(r), np.nan, float(lk)))
    # 2c. MCaRT（CC BY 4.0）: 480 化合物 × 4 カラム × 2 緩衝液
    for colname in mc.columns[1:]:
        c, buf = colname.split("-")
        ph, t0, cname = MCART_COLUMN[c]
        mid = "mcart_%s_%s" % (c, buf)
        col = dict(column_id="", column_name=cname, column_match="none", phase=ph, **{k: float(means.loc[ph, k]) for k in HSM_COLS}, logkEB=float(np.log10(means.loc[ph, "kEB"])))
        add_method("mcart_msd_zenodo", mid, col, NOTE["mcart_msd_zenodo"], length_mm=np.nan, id_mm=np.nan, particle_um=np.nan, temp_C=40.0, flow_mL_min=0.4, t0_min=t0,
                   pH_weak=2.0 if buf == "TFA" else 6.8, frac_acn=1.0, add_formic=0, add_acetic=0, add_tfa=int(buf == "TFA"), add_ammonium=int(buf == "NH4"),
                   add_basic=0, add_phosphate=0, gradient=json.dumps([[0, 0.05], [5, 0.85]]), isocratic_phi=np.nan)
        for i, v in enumerate(mc[colname].astype(str)):
            if "*" in v:
                continue
            try:
                rt = float(v)
            except ValueError:
                continue
            r = f_mc.get("mcart_%d" % i)
            if r is not None and np.isfinite(rt):
                rows.append((mid, reg(r), rt, np.nan))
    # 2d. Wiczling グループ（OSF、CC BY 4.0）: 300 化合物 × 5 カラム × 84 実験
    for oid, (cid_, match) in WICZLING_COLUMN.items():
        fs = sorted(glob.glob(rp("wiczling_osf", oid, "1-*.csv")))
        if not fs:
            C.log("Wiczling %s が無い（飛ばす）" % oid); continue
        d = pd.read_csv(fs[0])
        d = d[np.isfinite(d["RT"]) & np.isfinite(d["Area"])]
        d = d.sort_values(["METID", "EXPID", "Area"], ascending=[True, True, False])
        top = d.groupby(["METID", "EXPID"]).head(2).copy()
        top["rank"] = top.groupby(["METID", "EXPID"]).cumcount()
        a1 = top[top["rank"] == 0].set_index(["METID", "EXPID"])
        a2 = top[top["rank"] == 1].set_index(["METID", "EXPID"])[["RT", "Area"]].rename(columns={"RT": "RT2", "Area": "Area2"})
        a1 = a1.join(a2)
        amb = (a1["Area2"] > 0.5 * a1["Area"]) & ((a1["RT2"] - a1["RT"]).abs() > 0.3)      # 候補のピークが 2 本あって決められない組は外す
        a1 = a1[~amb].reset_index()
        col = colrec(cid_, match)
        for exp, s in a1.groupby("EXPID"):
            r0 = s.iloc[0]
            mid = "wz_%s_%02d" % (oid, int(exp))
            buf = "bicarb" if r0["pHs"] in (6.8, 10.5) else "acetate" if r0["pHs"] in (4.9, 5.8) else "formate"
            add_method("wiczling_osf", mid, dict(col), NOTE["wiczling_osf"], length_mm=50, id_mm=3.0, particle_um=2.5, temp_C=float(r0["Temp"]), flow_mL_min=float(r0["F"]),
                       t0_min=float(r0["to"]), pH_weak=float(r0["pHs"]), frac_acn=1.0 if str(r0["Mod"]).strip() == "ACN" else 0.0,
                       add_formic=0, add_acetic=0, add_tfa=0, add_ammonium=int(buf != "bicarb"), add_basic=int(buf == "bicarb"), add_phosphate=0,
                       gradient=json.dumps([[float(r0["td"]), float(r0["fio"])], [float(r0["td"] + r0["tg"]), float(r0["fik"])]]), isocratic_phi=np.nan, fixed_cal=1)
            for metid, rt in zip(s["METID"], s["RT"] - r0["te"]):
                r = f_wz.get(str(int(metid)))
                if r is not None:
                    rows.append((mid, reg(r), float(rt), np.nan))
    # 2e. Zhang ら（CC0）
    for p in zfiles:
        nn = os.path.basename(p)[8:10]
        wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
        S = [[("" if c is None else str(c).replace("\n", " ").strip()) for c in r] for r in wb["LC setups"].iter_rows(values_only=True)]
        kv = {r[0]: r[1] for r in S if r and r[0]}
        dims = [float(x) for x in re.findall(r"[\d.]+", kv.get("Column dimensions", ""))][:3]
        L, ID = max(dims[0], dims[1]), min(dims[0], dims[1])
        A_, B_ = kv.get("Mobile phase A", ""), kv.get("Mobile phase B", "")
        phiA = 0.1 if "methanol=90:10" in A_.replace(" ", "").lower() else 0.0
        gi = next(i for i, r in enumerate(S) if r and r[0].startswith("Gradient"))
        g = []
        for r in S[gi + 1:]:
            try:
                g.append((float(r[1]), float(r[2]), float(r[3])))
            except (ValueError, IndexError):
                pass
        bmax = max(x[2] for x in g); im = max(i for i, x in enumerate(g) if x[2] == bmax)
        prog = [[t, round(phiA + (1 - phiA) * b / 100.0, 4)] for t, _, b in g[:im + 1]]
        al = A_.lower()
        pH = 6.8 if "ammonium acetate" in al and "formic" not in al else 2.7 if ("0.1% formic" in al and "ammonium" not in al) else 3.0 if "0.1% formic" in al else 3.7
        col = colrec(*ZHANG_COLUMN[nn])
        note = "A: %s / B: %s。pH は添加剤からの見積もり。死時間は表の記載" % (A_[:60], B_[:40])
        if col["column_match"] == "approx":
            note += "。カラムの性質は、HSM2 の表で同じ固定相の粒径違いの製品の値（%s）を使った" % col["column_name"]
        add_method("zhang2024_scidb", "zhang_%s" % nn, col, note, length_mm=L, id_mm=ID, particle_um=dims[2],
                   temp_C=float(kv.get(next(k for k in kv if k.startswith("Column temperature")), "nan")), flow_mL_min=g[0][1], t0_min=float(kv["Dead time (min)"]),
                   pH_weak=pH, frac_acn=1.0 if "acetonitrile" in B_.lower() else 0.0, add_formic=int("formic" in al), add_acetic=0, add_tfa=0,
                   add_ammonium=int("ammonium" in al), add_basic=0, add_phosphate=0, gradient=json.dumps(prog), isocratic_phi=np.nan)
        R = list(wb["RT"].iter_rows(values_only=True))
        hd = [("" if c is None else str(c).replace("\n", " ")) for c in R[0]]
        i_smi = next(i for i, h in enumerate(hd) if "SMILES" in h); i_rt = next(i for i, h in enumerate(hd) if h.startswith("RT"))
        for r in R[1:]:
            try:
                rt = float(r[i_rt])
            except (TypeError, ValueError):
                continue
            ft = f_zh.get(str(r[i_smi]).replace("\xa0", "").strip()) if r[i_smi] else None
            if ft is not None:
                rows.append(("zhang_%s" % nn, reg(ft), rt, np.nan))

    # ---- 書き出し（同梱の base_data/ と同じ形）
    os.makedirs(OUT, exist_ok=True)
    order = ["method_id", "source", "column_id", "column_name", "column_match", "phase"] + HSM_COLS + ["logkEB", "length_mm", "id_mm", "particle_um", "temp_C", "flow_mL_min",
             "t0_min", "pH_weak", "frac_acn", "add_formic", "add_acetic", "add_tfa", "add_ammonium", "add_basic", "add_phosphate", "gradient", "isocratic_phi", "fixed_cal", "exclude", "note"]
    pd.DataFrame(methods)[order].to_csv(os.path.join(OUT, "methods.csv"), index=False, encoding="utf-8", lineterminator="\n")
    R = pd.DataFrame(rows, columns=["method_id", "ck", "rt_min", "logk_iso"])
    R.to_csv(os.path.join(OUT, "rows.csv.gz"), index=False, encoding="utf-8", lineterminator="\n", compression=dict(method="gzip", mtime=0))
    used = sorted(set(R["ck"]))
    # イオン化の規則の数は、標準の SMILES から数え直す（公開の重みを作った時と同じ。SD ファイルの構造から数えると、2 化合物で 1 つずれる）
    C.log("イオン化の規則の数を、標準の SMILES から数え直す")
    f_can = featurize("canonical", [(ck, comp[ck]["smiles"]) for ck in used], work)
    for ck in used:
        comp[ck] = dict(comp[ck], ion=f_can[ck]["ion"])
    T = pd.DataFrame([[ck, comp[ck]["smiles"]] + comp[ck]["d"] + comp[ck]["g"] + comp[ck]["ion"] + [1, comp[ck]["fp"]] for ck in used],
                     columns=["ck", "smiles"] + ["d_" + x for x in desc] + ["g_" + x for x in grp] + ["ion_" + x for x in ion] + ["ion_ok", "fp"])
    T.to_csv(os.path.join(OUT, "compounds.csv.gz"), index=False, encoding="utf-8", lineterminator="\n", compression=dict(method="gzip", mtime=0))
    st = [[nm, r["smiles"], r["inchikey"]] + r["d"] + r["g"] + r["ion"] for nm, r in f_man.items() if nm not in ("arginie", "p-nitrophenol")]
    pd.DataFrame(st, columns=["name", "smiles", "inchikey"] + ["d_" + x for x in desc] + ["g_" + x for x in grp] + ["ion_" + x for x in ion]).to_csv(
        os.path.join(OUT, "stock.csv"), index=False, encoding="utf-8", lineterminator="\n")
    # 名前の索引の材料（出所の表の名前と構造のまま）
    nm = [(n, s, "stock") for n, s in zip(man["name"], man["smiles"])] + [(n, s, "weber") for n, s in zip(web["name"], web["smiles"])]
    for d in sorted(glob.glob(rp("wiczling_osf", "*"))):
        fn, fs_ = os.path.join(d, "4-compounds-names.csv"), os.path.join(d, "5-smiles.csv")
        if os.path.exists(fn) and os.path.exists(fs_):
            a, b = pd.read_csv(fn), pd.read_csv(fs_)
            key = [c for c in a.columns if c in b.columns]
            ncol = [c for c in a.columns if "name" in c.lower()]; scol = [c for c in b.columns if "smiles" in c.lower()]
            mm = a.merge(b, on=key[0])
            nm += [(n, s, "wiczling") for n, s in zip(mm[ncol[0]], mm[scol[0]])]
            break
    R0 = list(openpyxl.load_workbook(zfiles[0], read_only=True, data_only=True)["RT"].iter_rows(values_only=True))
    hd = [("" if c is None else str(c).replace("\n", " ")) for c in R0[0]]
    i_smi = next(i for i, h in enumerate(hd) if "SMILES" in h); i_nm = next((i for i, h in enumerate(hd) if "name" in h.lower()), None)
    if i_nm is not None:
        nm += [(r[i_nm], str(r[i_smi]).replace("\xa0", "").strip(), "zhang") for r in R0[1:] if r[i_smi] and r[i_nm]]
    pd.DataFrame(nm, columns=["name", "smiles", "source"]).to_csv(os.path.join(OUT, "names.csv.gz"), index=False, encoding="utf-8", lineterminator="\n",
                                                                compression=dict(method="gzip", mtime=0))
    for f in ["names_manual.tsv", "names_weber.tsv", "names_exclude.tsv", "retention_form.json"]:
        shutil.copyfile(os.path.join(C.BASE, f), os.path.join(OUT, f))
    C.log("できた: 分析法 %d、測定 %d、化合物 %d → %s" % (len(methods), len(R), len(used), OUT))


if __name__ == "__main__":
    main()
