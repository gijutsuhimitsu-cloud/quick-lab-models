"""溶解度のモデル（勾配ブースティングの木）を学習して、JavaScript で読める JSON に書き出す。

入力（train/build_tables.cjs が書く）:
  <work>/train_25C.csv     1 行 = 溶質 × 溶媒（25 °C 付近の実測の平均）。有機溶媒のモデル（tm・notm）の学習に使う
  <work>/train_25C_aq.csv  水の行（BigSolDB の水・AqSolDB・自社のデータの水）。水のモデル（tm_w・notm_w）の学習に使う
出力:
  <out>/solubility_gbt.json   4 つのモデルと、検証での外れ（cv）。weights/solubility_gbt.json と同じ形
  <out>/reference_gbt.json    照合の見本（学習側 = scikit-learn の予測）。tests/check.cjs <out> で JavaScript と比べる
  <out>/gbt_cv.json           検証のまとめ

使い方（models/solubility/ で）:
  python train/train_gbt.py                     既定（公開の版と同じ設定。CPU で数分）
  python train/train_gbt.py --quick             試し用（木を減らす。値は本番と違う）
  python train/train_gbt.py --work DIR --out DIR --jobs N
必要な物: Python 3.9 以上、numpy・pandas・scikit-learn・joblib。

検証は溶質ごとに 5 つの組に分け、4 組で学習して残りの 1 組で外れを測る（同じ溶質が学習と検証の両方に入らない）。
GPU を使いたい時：GradientBoostingRegressor の所を LightGBM・XGBoost の GPU 版に差し替えられる（このファイルは CPU だけ）。
ただし書き出し（export）は scikit-learn の木の形を前提にしているので、差し替えた時は export も書き直すこと。
"""
import argparse, io, json, os, sys

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GroupKFold

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
VERSION = "0.1"

ap = argparse.ArgumentParser()
ap.add_argument("--work", default=os.environ.get("SOL_WORK", os.path.join(HERE, "work")))
ap.add_argument("--out", default=os.environ.get("SOL_OUT", os.path.join(HERE, "out")))
ap.add_argument("--quick", action="store_true", help="木を 60 本にして速く回す（流れの確かめ用）")
ap.add_argument("--jobs", type=int, default=5, help="並べて走らせる数（検証の 5 組）")
ap.add_argument("--version", default=None, help="meta.version に書く版（既定は 0.1、自社のデータを足した時は 0.1+own）")
args = ap.parse_args()

# 公開の版の設定（ページの第 5 版と同じ。溶質を分けた 5 分割の検証で選んだ）
PARAMS = dict(n_estimators=500, max_depth=5, learning_rate=0.04, subsample=0.8, min_samples_leaf=10, max_features=0.5, random_state=0)
if args.quick:
    PARAMS["n_estimators"] = 60
W_BIG = 3.0            # 水のモデル：BigSolDB の水の行の重み（AqSolDB は 1）
MIN_ROWS, MIN_GROUPS = 20, 5   # これより少ない時はそのモデルを作らない（推論は、ある方のモデルで代わりにする）
META = {"group", "has_tm", "solvent", "y", "src", "aqg", "aqsd", "w", "can", "tms", "tmc"}
os.makedirs(args.out, exist_ok=True)


def stats(e):
    e = np.asarray(e, float)
    if len(e) == 0:
        return {"n": 0}
    a = np.sort(np.abs(e))
    n = len(a)
    q = lambda p: float(round(a[int(p * (n - 1))], 3))
    return {"n": int(n), "bias": float(round(np.mean(e), 3)), "p50": q(0.5), "p68": q(0.68), "p90": q(0.9)}


def export(m, cols):
    """左の子はいつも親の次の番号（深さ優先で作る木）なので書かない。r が -1 の所が葉。葉の f・t と、枝の v は 0 にして縮める。
    しきい値は 6 桁（5 桁に丸めると、分かれ目の上の値が反対へ行くことがある）。"""
    trees = []
    for est in m.estimators_[:, 0]:
        t = est.tree_
        leaf = t.children_right < 0
        assert all(t.children_left[i] == i + 1 for i in range(t.node_count) if not leaf[i])
        trees.append({"f": [0 if leaf[i] else int(t.feature[i]) for i in range(t.node_count)],
                      "t": [0 if leaf[i] else round(float(t.threshold[i]), 6) for i in range(t.node_count)],
                      "r": [int(v) for v in t.children_right],
                      "v": [round(float(t.value[i, 0, 0]), 5) if leaf[i] else 0 for i in range(t.node_count)]})
    return {"init": round(float(m.init_.constant_[0][0]), 5), "lr": m.learning_rate, "cols": cols, "trees": trees}


def weights_or_none(w):
    """重みがすべて 1 なら None（scikit-learn に重みを渡さない。公開の版と同じ結果にするため）"""
    return None if w is None or np.all(w == 1.0) else w


def fit(Xa, ya, wa):
    return GradientBoostingRegressor(**PARAMS).fit(Xa, ya, sample_weight=weights_or_none(wa))


OWN_IN_REF = [False]   # 照合の見本に自社の行を使ったか（reference_gbt.json の own_rows）


def sample_rows(rows, src, tms, n=40):
    """照合の見本の行：公開データの行から等間隔に n 行。自社のデータの行と、融点が手元の融点の表から来た行（tms = mp）は入れない（自社のデータだけの時は自社の行）"""
    ok = (src[rows] != "own") & (tms[rows] != "mp")
    if not ok.any():                  # 自社のデータだけで学習した時は、自社の行で照合する（照合の値に自社の構造が入る。run_all.cjs は tests/own/ へ写す）
        ok = tms[rows] != "mp"
        OWN_IN_REF[0] = True
    if not ok.any():
        ok = np.ones(len(rows), bool)
    pub = rows[ok]
    return pub[:: max(1, len(pub) // n)][:n]


def ref_entry(df, X_, k, cols, m):
    e = {"rows": X_[k].tolist(), "pred": m.predict(X_[k][:, cols]).tolist(), "solvent": df["solvent"].values[k].tolist()}
    if "can" in df.columns:
        e["can"] = df["can"].values[k].tolist()
        e["tmc"] = [None if v in ("", None) else float(v) for v in df["tmc"].values[k]]
    return e


def folds_of(X, y, g):
    k = min(5, len(np.unique(g)))
    if k < 2:
        return []
    return list(GroupKFold(n_splits=k).split(X, y, g))


# ---------- 有機溶媒のモデル（BigSolDB のすべての溶媒の行で学習。ページは水以外に使う） ----------
d = pd.read_csv(os.path.join(args.work, "train_25C.csv"), keep_default_na=False)
names = [c for c in d.columns if c not in META]
ix = {n: i for i, n in enumerate(names)}
X, y, g = d[names].values.astype(float), d["y"].values.astype(float), d["group"].values.astype(int)
src = d["src"].values if "src" in d.columns else np.array(["bigsol"] * len(d))
w_all = d["w"].values.astype(float) if "w" in d.columns else np.ones(len(d))
own = src == "own"
tms = d["tms"].values if "tms" in d.columns else np.array([""] * len(d))
water = d["sv_water"].values.astype(float) == 1
SETS = {"tm": {"rows": d["has_tm"].values.astype(int) == 1, "drop": ["formula0"]},
        "notm": {"rows": np.ones(len(d), bool), "drop": ["crystal", "formula"]}}
out, cv, ref = {}, {}, {"note": "学習側（scikit-learn）の予測。tests/check.cjs が JavaScript の計算と比べる。rows は入力の列（names の順）", "names": names, "models": {}}
for key, S in SETS.items():
    cols = [ix[n] for n in names if n not in S["drop"]]
    rows = np.where(S["rows"])[0]
    if len(rows) < MIN_ROWS or len(np.unique(g[rows])) < MIN_GROUPS:
        print(key, "行が少ない（%d 行）。このモデルは作らない" % len(rows))
        continue
    Xs, ys, gs, ws = X[rows][:, cols], y[rows], g[rows], w_all[rows]
    pred = np.full(len(rows), np.nan)
    folds = folds_of(Xs, ys, gs)
    fitted = Parallel(n_jobs=args.jobs)(delayed(lambda a, b: fit(Xs[a], ys[a], ws[a]).predict(Xs[b]))(tr, te) for tr, te in folds)
    for (tr, te), p in zip(folds, fitted):
        pred[te] = p
    e = pred - ys
    base = X[rows][:, ix["formula" if key == "tm" else "formula0"]] - ys          # 比べる相手：式のまま
    w = water[rows]
    o = own[rows]
    # 実測 1 点から換算した時：同じ溶質のほかの溶媒の実測に、モデルの差を足す（有機溶媒どうし）
    pe = []
    for _, sub in pd.DataFrame({"g": gs[~w], "p": pred[~w], "y": ys[~w]}).groupby("g"):
        p, yy = sub["p"].values, sub["y"].values
        if len(p) > 1 and not np.isnan(p).any():
            dm = (yy[:, None] + (p[None, :] - p[:, None])) - yy[None, :]
            pe.append(dm[~np.eye(len(p), dtype=bool)])
    ok = ~np.isnan(e)
    cv[key] = {"organic": stats(e[ok & ~w & ~o]) if o.any() else stats(e[ok & ~w]), "water": stats(e[ok & w & ~o]) if o.any() else stats(e[ok & w]),
               "organic_formula_baseline": stats(base[ok & ~w]), "water_formula_baseline": stats(base[ok & w]),
               "one_point_organic": stats(np.concatenate(pe)) if pe else {"n": 0}}
    if o.any():
        cv[key]["own_organic"], cv[key]["own_water"] = stats(e[ok & ~w & o]), stats(e[ok & w & o])
        if (~o).sum() == 0:          # 自社のデータだけで学習した時は、幅に自社のデータの外れを使う
            cv[key]["organic"], cv[key]["water"] = cv[key]["own_organic"], cv[key]["own_water"]
    print(key, json.dumps({k: v.get("p68") for k, v in cv[key].items()}, ensure_ascii=False), flush=True)
    m = fit(Xs, ys, ws)
    out[key] = export(m, cols)
    ref["models"][key] = ref_entry(d, X, sample_rows(rows, src, tms), cols, m)

# ---------- 水のモデル（BigSolDB の水は重み W_BIG、AqSolDB の水は 1）。有機溶媒のモデルは上のまま ----------
AQ = os.path.join(args.work, "train_25C_aq.csv")
meta_w = None
if os.path.exists(AQ):
    a = pd.read_csv(AQ, keep_default_na=False)
    a = a[a["sv_water"].astype(float) == 1].reset_index(drop=True)
    XA, yA, gA = a[names].values.astype(float), a["y"].values.astype(float), a["group"].values.astype(int)
    srcA = a["src"].values
    bigA, aqA, ownA = srcA == "bigsol", srcA == "aqsol", srcA == "own"
    wA = np.where(bigA, W_BIG, 1.0) * (a["w"].values.astype(float) if "w" in a.columns else 1.0)
    for key, S in SETS.items():
        if key not in cv:
            cv[key] = {}
        cols = [ix[n] for n in names if n not in S["drop"]]
        rows = np.where(a["has_tm"].astype(int).values == 1)[0] if key == "tm" else np.arange(len(a))
        if len(rows) < MIN_ROWS or len(np.unique(gA[rows])) < MIN_GROUPS:
            print(key + "_w", "行が少ない（%d 行）。このモデルは作らない" % len(rows))
            continue
        Xs, ys, gs, ws = XA[rows][:, cols], yA[rows], gA[rows], wA[rows]
        folds = folds_of(Xs, ys, gs)
        res = Parallel(n_jobs=args.jobs)(delayed(lambda p, q: fit(Xs[p], ys[p], ws[p]).predict(Xs[q]))(tr, te) for tr, te in folds)
        pred = np.full(len(rows), np.nan)
        for (tr, te), r in zip(folds, res):
            pred[te] = r
        e, bb, aq, oo = pred - ys, bigA[rows], aqA[rows], ownA[rows]
        ok = ~np.isnan(e)
        cv[key]["water_one_model"] = cv[key].get("water", {"n": 0})        # 上の 1 つのモデルでの水（BigSolDB の水だけで測った値）
        cv[key]["water"] = stats(e[ok & ~oo]) if (~oo).any() else stats(e[ok & oo])   # 水だけのモデル、水の全部（BigSolDB と AqSolDB）で測った値
        cv[key]["water_bigsol"], cv[key]["water_aqsol"] = stats(e[ok & bb]), stats(e[ok & aq])
        if oo.any():
            cv[key]["own_water"] = stats(e[ok & oo])
        print(key + "_w", json.dumps({k: cv[key][k].get("p68") for k in ("water", "water_bigsol", "water_aqsol")}, ensure_ascii=False), flush=True)
        m = fit(Xs, ys, ws)
        out[key + "_w"] = export(m, cols)
        tmsA = a["tms"].values if "tms" in a.columns else np.array([""] * len(a))
        ref["models"][key + "_w"] = ref_entry(a, XA, sample_rows(rows, srcA, tmsA), cols, m)
    meta_w = {"n_rows": int(len(a)), "n_solutes": int(a["group"].nunique()), "n_aqsol": int(aqA.sum()), "weight_bigsol": W_BIG,
              "source": "BigSolDB 2.2 の水 + AqSolDB (Sorkun et al.), Harvard Dataverse, doi:10.7910/DVN/OVHAW8, CC0 1.0"}
    if ownA.any():
        meta_w["n_own"] = int(ownA.sum())

out["names"] = names
out["cv"] = cv
n_own = int(own.sum())
version = args.version or (VERSION + ("+own" if n_own or (meta_w and meta_w.get("n_own")) else ""))
out["meta"] = {"name": "solubility_gbt", "version": version, "made": pd.Timestamp.now().strftime("%Y-%m-%d"), "params": PARAMS,
               "n_rows": int((~own).sum()), "n_solutes": int(pd.Series(g[~own]).nunique()),
               "source": "BigSolDB 2.2 (Krasnov et al.), Zenodo, doi:10.5281/zenodo.22648301, CC BY 4.0", "water": meta_w,
               "license": "CC BY 4.0（出典は ATTRIBUTION.md）", "quick": bool(args.quick)}
if n_own:
    out["meta"]["n_own_rows"] = n_own
json.dump(cv, open(os.path.join(args.out, "gbt_cv.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
ref["own_rows"] = bool(OWN_IN_REF[0])
json.dump(ref, open(os.path.join(args.out, "reference_gbt.json"), "w", encoding="utf-8"), separators=(",", ":"))
dst = os.path.join(args.out, "solubility_gbt.json")
with open(dst, "w", encoding="utf-8") as f:
    json.dump(out, f, separators=(",", ":"), ensure_ascii=False)
print("書いた", dst, os.path.getsize(dst), "B")
