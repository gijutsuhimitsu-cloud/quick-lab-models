"""融点の推算モデル（勾配ブースティングの木）を学習して、JavaScript で読める JSON に書き出す。

入力: <base>/tm_train.csv.gz（train/base_data。融点の公開データ（CC0 1.0）から作った表。作り直しは train/rebuild/bradley_tm.cjs）
      --check-solubility の時だけ: <work>/train_25C.csv・train_solutes.json（train/build_tables.cjs）、weights/formula.json
出力: <out>/tm_model.json（weights/tm_model.json と同じ形）、<out>/reference_tm.json（照合の見本）、<out>/tm_cv.json
使い方（models/solubility/ で）:
  python train/train_tm.py                       融点のモデルだけ（1 分ほど）
  python train/train_tm.py --check-solubility    推算した融点が溶解度の見込みに効くかも測る（数分）
融点の外れは、溶解度の学習の溶質（in_bigsol = 1）を学習から外して、その溶質で測る。書き出すモデルは全部の行で学習し直した物。
"""
import argparse, io, json, os, sys

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GroupKFold

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--work", default=os.environ.get("SOL_WORK", os.path.join(HERE, "work")))
ap.add_argument("--out", default=os.environ.get("SOL_OUT", os.path.join(HERE, "out")))
ap.add_argument("--weights", default=os.environ.get("SOL_WEIGHTS", os.path.join(os.path.dirname(HERE), "weights")))
ap.add_argument("--base", default=os.environ.get("SOL_BASE", os.path.join(HERE, "base_data")))
ap.add_argument("--check-solubility", action="store_true")
ap.add_argument("--quick", action="store_true", help="木を 40 本にする（流れの確かめ用）")
args = ap.parse_args()
P_TM = dict(n_estimators=250, max_depth=5, learning_rate=0.08, subsample=0.7, min_samples_leaf=5, random_state=0)
P_SOL = dict(n_estimators=300, max_depth=4, learning_rate=0.06, subsample=0.8, min_samples_leaf=5, random_state=0)
if args.quick:
    P_TM["n_estimators"] = 40
    P_SOL["n_estimators"] = 40
META = {"group", "has_tm", "solvent", "y", "src", "aqg", "aqsd", "w", "can", "tms", "tmc"}
os.makedirs(args.out, exist_ok=True)


def stats(e):
    a = np.sort(np.abs(e)); n = len(a); q = lambda p: float(round(a[int(p * (n - 1))], 3))
    return {"n": int(n), "bias": float(round(np.mean(e), 3)), "p50": q(0.5), "p68": q(0.68), "p90": q(0.9)}


def export(m, cols):
    """左右の子の番号（l・r。葉は -1）を書く形。入力は 4 桁に丸めた記述子と整数の数なので、しきい値は 5 桁で足りる"""
    trees = []
    for est in m.estimators_[:, 0]:
        t = est.tree_
        trees.append({"f": [int(v) if v >= 0 else 0 for v in t.feature], "t": [round(float(v), 5) for v in t.threshold],
                      "l": [int(v) for v in t.children_left], "r": [int(v) for v in t.children_right], "v": [round(float(v), 3) for v in t.value[:, 0, 0]]})
    return {"init": round(float(m.init_.constant_[0][0]), 3), "lr": m.learning_rate, "cols": cols, "trees": trees}


# ---- 融点のモデル ----
t = pd.read_csv(os.path.join(args.base, "tm_train.csv.gz"))
fx = [c for c in t.columns if c not in ("in_bigsol", "mpC")]
tr, te = t[t.in_bigsol == 0], t[t.in_bigsol == 1]
res = {"n_train": int(len(tr))}
if len(te):
    m0 = GradientBoostingRegressor(**P_TM).fit(tr[fx].values, tr.mpC.values)
    res["tm_holdout_bigsol_solutes_degC"] = stats(m0.predict(te[fx].values) - te.mpC.values)
    res["tm_baseline_mean_degC"] = stats(tr.mpC.mean() - te.mpC.values)
    print("融点", json.dumps(res, ensure_ascii=False), flush=True)

# ---- 推算した融点は溶解度に効くか（融点のある溶質 × 溶媒で）----
if args.check_solubility and len(te):
    K = json.load(open(os.path.join(args.weights, "formula.json"), encoding="utf-8"))["K"]
    sol = json.load(open(os.path.join(args.work, "train_solutes.json"), encoding="utf-8"))
    tm_pred = m0.predict(np.array([s["x"] for s in sol]))
    d = pd.read_csv(os.path.join(args.work, "train_25C.csv"), keep_default_na=False)
    names = [c for c in d.columns if c not in META]
    ix = {n: i for i, n in enumerate(names)}
    d = d[(d.has_tm.astype(int) == 1) & (d.group.astype(int) < len(sol))].reset_index(drop=True)
    X, y, g, w = d[names].values.astype(float), d.y.values.astype(float), d.group.values.astype(int), d.sv_water.values.astype(float) == 1
    crystal = lambda tm: np.where(tm <= 25, 0.0, -56.5 * (tm - 25) / (2.302585 * 8.314462 * 298.15))
    Xp = X.copy()
    c = crystal(tm_pred[g])
    Xp[:, ix["crystal"]] = c
    Xp[:, ix["formula"]] = np.where(w, 0.5 - X[:, ix["logPow"]] + c, K["k0"] - K["kP"] * X[:, ix["logPow"]] + K["kC"] * c + K["kT"] * X[:, ix["transfer"]])
    c_tm = [ix[n] for n in names if n != "formula0"]
    c_no = [ix[n] for n in names if n not in ("crystal", "formula")]
    p = {k: np.zeros(len(y)) for k in ("true", "pred", "none")}
    for a, b in GroupKFold(n_splits=5).split(X, y, g):
        mt = GradientBoostingRegressor(**P_SOL).fit(X[a][:, c_tm], y[a])
        p["true"][b] = mt.predict(X[b][:, c_tm]); p["pred"][b] = mt.predict(Xp[b][:, c_tm])
        p["none"][b] = GradientBoostingRegressor(**P_SOL).fit(X[a][:, c_no], y[a]).predict(X[b][:, c_no])
    res["solubility_organic"] = {"measured_Tm": stats((p["true"] - y)[~w]), "predicted_Tm_into_measured_Tm_model": stats((p["pred"] - y)[~w]), "no_Tm_model": stats((p["none"] - y)[~w])}
    print("solubility_organic", json.dumps(res["solubility_organic"], ensure_ascii=False))
json.dump(res, open(os.path.join(args.out, "tm_cv.json"), "w", encoding="utf-8"), indent=1)

# ---- 書き出す：全部の行で学習し直す ----
m1 = GradientBoostingRegressor(**P_TM).fit(t[fx].values, t.mpC.values)
out = {"model": export(m1, list(range(len(fx)))), "cv": res.get("tm_holdout_bigsol_solutes_degC"), "features": fx,
       "meta": {"name": "tm_model", "version": "0.1", "n": int(len(t)), "params": P_TM, "quick": bool(args.quick), "license": "CC BY 4.0（出典は ATTRIBUTION.md）",
                "what": "融点 / °C の推算。入力は js/features.js の desc 41 個 ＋ 官能基の数 56 個（solute.x）。cv は溶解度の学習の溶質を外して測った外れ（°C）",
                "source": "Jean-Claude Bradley Open Melting Point Dataset (Bradley, Williams, Lang), figshare, doi:10.6084/m9.figshare.1031637, CC0 1.0（出所の欄が非営利の薬のデータベースと試薬会社のカタログの行を除く）"}}
json.dump(out, open(os.path.join(args.out, "tm_model.json"), "w", encoding="utf-8"), separators=(",", ":"), ensure_ascii=False)
# 照合の見本：溶解度の学習の溶質の行（融点そのものは書かない。入力と予測だけ）
k = np.arange(len(t))[:: max(1, len(t) // 40)][:40]
json.dump({"note": "学習側（scikit-learn）の融点の予測。rows は入力の列（features の順）", "rows": t[fx].values[k].tolist(), "pred": m1.predict(t[fx].values[k]).tolist()},
          open(os.path.join(args.out, "reference_tm.json"), "w", encoding="utf-8"), separators=(",", ":"))
print("書いた", os.path.join(args.out, "tm_model.json"), os.path.getsize(os.path.join(args.out, "tm_model.json")), "B")
