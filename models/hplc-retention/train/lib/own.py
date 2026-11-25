# SPDX-License-Identifier: MIT
"""自社のデータ（own_data/ の methods.csv・measurements.csv）を、学習の表の形に直す。形の説明は OWN_DATA.md。

構造（SMILES）の特徴量は、ページと同じ js/features.js を Node（RDKit.js）で計算する（scripts/featurize.cjs）。
同じ化合物（InChIKey の前半 14 文字が同じ）が基礎データにあれば、基礎データの特徴量を使う（学習の時と同じ決まり）。
データはこの PC の中だけで扱う。どこにも送らない。
"""
import os, re, json, math
import numpy as np, pandas as pd
from . import common as C
from .table import columns_table

PHASE_WORDS = {"c18": "C18", "c8": "C8", "phenyl": "Phenyl", "cn": "CN", "cyano": "CN", "ep": "EP", "f": "F", "pfp": "F", "short": "short"}
HSM = ["H", "A", "B", "C", "V", "D"]
ADDS = ["formic", "acetic", "tfa", "ammonium", "basic", "phosphate"]


class OwnDataError(SystemExit):
    pass


def _num(v):
    if v is None:
        return np.nan
    if isinstance(v, (int, float, np.integer, np.floating)):
        return float(v)
    s = str(v).strip()
    if not s:
        return np.nan
    try:
        return float(s)
    except ValueError:
        return None          # 数でない


def _program(s, where):
    """'0:0.05;10:0.95;12:0.95'（分:有機溶媒の割合 0〜1）または JSON の [[0,0.05],[10,0.95]] → [[t, phi], ...]"""
    s = str(s).strip()
    if s.startswith("["):
        prog = json.loads(s)
    else:
        prog = []
        for part in re.split(r"[;\s]+", s):
            if not part:
                continue
            t, p = part.split(":")
            prog.append([float(t), float(p)])
    if len(prog) < 2:
        raise OwnDataError("%s: 勾配は 2 点以上（例 0:0.05;10:0.95）" % where)
    ts = [p[0] for p in prog]
    if any(b < a for a, b in zip(ts, ts[1:])):
        raise OwnDataError("%s: 勾配の時刻が増えていない" % where)
    if any(not (0 <= p[1] <= 1) for p in prog):
        raise OwnDataError("%s: 勾配の割合は 0〜1 で書く（%% ではない）" % where)
    return [[float(t), float(p)] for t, p in prog]


def _dead_time(L, idm, flow):
    """js/retention_model.js の deadTime と同じ目安（空隙率 0.65）"""
    return math.pi * (idm / 20.0) ** 2 * (L / 10.0) * 0.65 / flow


def _column(r, where, cols):
    """カラムの性質: columns.json の id か一般名、固定相の種類の言葉（その種類の中央値）、または custom（H〜D と kEB か logkEB を書く）"""
    key = str(r.get("column", "")).strip()
    byid = {c["id"]: c for c in cols["columns"]}
    byname = {c["name"].lower(): c for c in cols["columns"]}
    if key in byid or key.lower() in byname:
        c = byid.get(key) or byname[key.lower()]
        return dict(column_id=c["id"], column_name=c["name"], column_match="exact", phase=c["phase"], logkEB=float(c["logkEB"]),
                    **{k: float(c[k]) for k in HSM})
    if key.lower() in PHASE_WORDS:
        ph = PHASE_WORDS[key.lower()]
        sel = [c for c in cols["columns"] if c["phase"] == ph]
        if not sel:
            raise OwnDataError("%s: columns.json に固定相 %s のカラムが無い" % (where, ph))
        med = {k: float(np.median([c[k] for c in sel])) for k in HSM + ["logkEB"]}
        return dict(column_id="", column_name="%s（columns.json の %s %d 本の中央値）" % (ph, ph, len(sel)), column_match="phase_median", phase=ph, **med)
    if key.lower() == "custom":
        v = {k: _num(r.get(k)) for k in HSM}
        lk = _num(r.get("logkEB"))
        if lk is None or not np.isfinite(lk):
            k = _num(r.get("kEB"))
            lk = math.log10(k) if (k is not None and np.isfinite(k) and k > 0) else np.nan
        ph = str(r.get("phase", "")).strip()
        ph = PHASE_WORDS.get(ph.lower(), ph)
        bad = [k for k in HSM if v[k] is None or not np.isfinite(v[k])] + ([] if np.isfinite(lk) else ["kEB か logkEB"])
        if bad:
            raise OwnDataError("%s: column が custom の時は %s を書く" % (where, "・".join(bad)))
        return dict(column_id="", column_name=str(r.get("column_name", "") or "自社のカラム"), column_match="custom", phase=ph or "other", logkEB=lk, **v)
    raise OwnDataError("%s: column '%s' が分からない。weights/columns.json の id（例 hsm2_3）か一般名（例 C18 No.3）、"
                       "固定相の種類（C18・C8・Phenyl・CN・EP・F・short）、custom のどれか" % (where, key))


def read_methods(d):
    p = os.path.join(d, "methods.csv")
    if not os.path.exists(p):
        raise OwnDataError("%s が無い" % p)
    T = pd.read_csv(p, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    need = ["method_id", "column", "organic"]
    miss = [c for c in need if c not in T.columns]
    if miss:
        raise OwnDataError("methods.csv に列 %s が無い" % ", ".join(miss))
    cols = columns_table()
    out, seen = [], set()
    for i, r in T.iterrows():
        r = {k: (v.strip() if isinstance(v, str) else v) for k, v in r.items()}
        where = "methods.csv の %d 行目（%s）" % (i + 2, r["method_id"])
        mid = r["method_id"]
        if not mid or mid in seen:
            raise OwnDataError("%s: method_id が空か重なっている" % where)
        seen.add(mid)
        col = _column(r, where, cols)
        org = r["organic"].lower()
        frac_acn = 1.0 if org in ("acn", "acetonitrile", "アセトニトリル") else 0.0 if org in ("meoh", "methanol", "メタノール") else _num(org)
        if frac_acn is None or not (0 <= frac_acn <= 1):
            raise OwnDataError("%s: organic は acn か meoh（または有機溶媒に占めるアセトニトリルの割合 0〜1）" % where)
        phi, grad = _num(r.get("phi", "")), str(r.get("gradient", "") or "").strip()
        if phi is None:
            raise OwnDataError("%s: phi が数でない" % where)
        if np.isfinite(phi) == bool(grad):
            raise OwnDataError("%s: 等組成なら phi（0〜1）、勾配なら gradient の、どちらか 1 つだけを書く" % where)
        if np.isfinite(phi) and not (0 <= phi <= 1):
            raise OwnDataError("%s: phi は 0〜1（%% ではない）" % where)
        nums = {k: _num(r.get(k, "")) for k in ["t0_min", "dwell_min", "pH", "temp_C", "length_mm", "id_mm", "particle_um", "flow_mL_min"]}
        bad = [k for k, v in nums.items() if v is None]
        if bad:
            raise OwnDataError("%s: 数でない欄 %s" % (where, ", ".join(bad)))
        t0, dwell = nums["t0_min"], nums["dwell_min"]
        t0_note = ""
        if not np.isfinite(t0):
            if all(np.isfinite(nums[k]) for k in ["length_mm", "id_mm", "flow_mL_min"]):
                t0 = _dead_time(nums["length_mm"], nums["id_mm"], nums["flow_mL_min"])
                t0_note = "t0 はカラムの寸法と流量からの目安（空隙率 0.65）"
            elif grad:
                raise OwnDataError("%s: 勾配では t0_min（死時間）か、length_mm・id_mm・flow_mL_min が要る" % where)
        fixed = 0
        gjson = ""
        if grad:
            prog = _program(grad, where)
            if np.isfinite(dwell):
                prog = [[t + dwell, p] for t, p in prog]      # 装置の遅れは勾配の表に入れる（基礎データの Wiczling と同じ）。補正の当てはめはしない
                fixed = 1 if not t0_note else 0
            gjson = json.dumps(prog)
        add = str(r.get("additive", "") or "").lower()
        out.append(dict(method_id="own_" + mid, source=C.OWN_SRC, **col,
                        length_mm=nums["length_mm"], id_mm=nums["id_mm"], particle_um=nums["particle_um"], temp_C=nums["temp_C"],
                        flow_mL_min=nums["flow_mL_min"], t0_min=t0 if np.isfinite(t0) else np.nan, pH_weak=nums["pH"], frac_acn=float(frac_acn),
                        **{"add_" + a: int(a in add) for a in ADDS},
                        gradient=gjson, isocratic_phi=phi if np.isfinite(phi) else np.nan, fixed_cal=fixed, exclude=0,
                        note="; ".join(x for x in [str(r.get("note", "") or ""), t0_note] if x)))
    return pd.DataFrame(out)


def read_measurements(d, md):
    p = os.path.join(d, "measurements.csv")
    if not os.path.exists(p):
        raise OwnDataError("%s が無い" % p)
    T = pd.read_csv(p, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    for c in ["method_id", "smiles"]:
        if c not in T.columns:
            raise OwnDataError("measurements.csv に列 %s が無い" % c)
    mrow = md.set_index("method_id")
    out, skipped = [], []
    for i, r in T.iterrows():
        where = "measurements.csv の %d 行目" % (i + 2)
        mid = "own_" + r["method_id"].strip()
        if mid not in mrow.index:
            raise OwnDataError("%s: method_id '%s' が methods.csv に無い" % (where, r["method_id"]))
        m = mrow.loc[mid]
        smi = r["smiles"].strip()
        if not smi:
            raise OwnDataError("%s: smiles が空" % where)
        rt, lk, k = (_num(r.get(x, "")) for x in ("rt_min", "log_k", "k"))
        if None in (rt, lk, k):
            raise OwnDataError("%s: rt_min・log_k・k に数でない値" % where)
        if np.isfinite(m["isocratic_phi"]):
            if not np.isfinite(lk):
                if np.isfinite(k):
                    lk = math.log10(k) if k > 0 else np.nan
                elif np.isfinite(rt) and np.isfinite(m["t0_min"]):
                    lk = math.log10((rt - m["t0_min"]) / m["t0_min"]) if rt > m["t0_min"] else np.nan
            if not np.isfinite(lk):
                skipped.append((where, "等組成の保持が分からない（log_k・k、または rt_min と t0_min が要る。rt が t0 以下）")); continue
            out.append(dict(method_id=mid, smiles=smi, rt_min=np.nan, logk_iso=float(lk), name=r.get("name", "").strip()))
        else:
            if not np.isfinite(rt):
                skipped.append((where, "勾配には rt_min が要る")); continue
            # rt が t0 以下（保持されていない）の行も残す。lkw を逆算できないので学習には使われないが、評価には入る（基礎データと同じ扱い）
            out.append(dict(method_id=mid, smiles=smi, rt_min=float(rt), logk_iso=np.nan, name=r.get("name", "").strip()))
    return pd.DataFrame(out, columns=["method_id", "smiles", "rt_min", "logk_iso", "name"]), skipped


def featurize(smiles, work):
    """SMILES の一覧 → {smiles: 特徴量}（Node + RDKit.js + js/features.js）。読めない構造は入らない。"""
    os.makedirs(work, exist_ok=True)
    pin, pout = os.path.join(work, "own_structures.tsv"), os.path.join(work, "own_features.jsonl")
    uniq = list(dict.fromkeys(smiles))
    with open(pin, "w", encoding="utf-8", newline="\n") as f:
        for i, s in enumerate(uniq):
            f.write("%d\t%s\n" % (i, s))
    C.node("featurize.cjs", "tsv", pin, pout)
    res = {}
    for line in open(pout, encoding="utf-8"):
        r = json.loads(line)
        if r.get("ok") and r.get("inchikey"):
            res[uniq[int(r["id"])]] = r
    return res


def load(d, work, base_comp):
    """→ (methods, rows, 新しい化合物 {ck: 特徴量}, 名前 [(name, smiles)], 報告)"""
    md = read_methods(d)
    meas, skipped = read_measurements(d, md)
    feats = featurize(list(meas["smiles"]), work)
    new, cks, names, bad = {}, [], [], []
    for s, nm in zip(meas["smiles"], meas["name"]):
        r = feats.get(s)
        if r is None:
            cks.append(None); bad.append(s); continue
        ck = r["inchikey"][:14]
        cks.append(ck)
        if ck not in base_comp and ck not in new:
            new[ck] = dict(smiles=r["smiles"], d=np.asarray(r["d"], float), g=np.asarray(r["g"], float), ion=np.asarray(r["ion"], float), ion_ok=True, fp=r["fp"])
        if nm:
            names.append((nm, s))
    meas["ck"] = cks
    rows = meas[meas["ck"].notna()][["method_id", "ck", "rt_min", "logk_iso"]].reset_index(drop=True)
    used = set(rows["method_id"])
    empty = [m for m in md["method_id"] if m not in used]
    md = md[md["method_id"].isin(used)].reset_index(drop=True)
    rep = dict(n_methods=int(len(md)), n_rows=int(len(rows)), n_compounds=int(rows["ck"].nunique()), n_new_compounds=len(new),
               n_known_compounds=int(rows["ck"].nunique() - len(new)), unreadable_smiles=sorted(set(bad)),
               skipped_rows=[dict(row=w, reason=x) for w, x in skipped], methods_without_rows=empty)
    return md, rows, new, names, rep
