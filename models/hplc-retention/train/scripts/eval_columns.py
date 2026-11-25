# SPDX-License-Identifier: MIT
"""学習に使っていないカラムでの評価（任意。木を 5 回学習し直すので、木 900 本だと 5〜10 分）。
  python train/scripts/eval_columns.py [名前（既定 free）] [--strict] [--trees N]
勾配の分析法をカラムで 5 つに分け、1 組ずつ外して学習し、外した分析法の保持時間を予測する（外した分析法は補正なし＝死時間の倍率 1・遅れ 0）。
  --strict: 同じカラム（columns.json の id）の行を、すべての出所（等組成・HSM2 の再構成・SMRT も）から外す。
前提: python train/run.py を同じ名前で実行してある（train/work/<名前>/calibration.json を使う）。
出力: train/out/<名前>/metrics_unseen_columns[_strict].json・eval_unseen_columns[_strict].csv
注意: 分け方はカラムの id（一般名）のハッシュで決める。以前の非公開の作業の表（製品名で分けた）とは分け方が違うので、数字は少し違う。
"""
import os, sys, json, hashlib, argparse
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from lib import common as C
from lib.runs import data_for_run
from lib.data import fit, metrics_rt

ap = argparse.ArgumentParser()
ap.add_argument("name", nargs="?", default="free")
ap.add_argument("--strict", action="store_true")
ap.add_argument("--trees", type=int, default=0)
A = ap.parse_args()
O = os.path.join(C.OUT, A.name)
D, cfg, cal_all = data_for_run(A.name)
trees = A.trees or cfg["trees"]
sources = sorted(set(D.src[D.in_tr]))
md, rows = D.md, D.rows
tm = set(md["method_id"][md["source"].isin(sources)])
gm = md[md["method_id"].isin(tm) & ~md["source"].isin(["smrt_figshare", "hsm2_mendeley", "weber_isocratic_mendeley"]) & md["gradient"].fillna("").str.startswith("[")]
md["colkey"] = [i if isinstance(i, str) and i else (c if isinstance(c, str) else "") for i, c in zip(md["column_id"], md["column_name"])]
gm = md[md["method_id"].isin(set(gm["method_id"]))]
cols = sorted(set(gm["colkey"]))
fold = {c: int(hashlib.md5(c.encode()).hexdigest()[:6], 16) % 5 for c in cols}
res = []
for k in range(5):
    keys = {c for c in cols if fold[c] == k}
    held = set(gm[gm["colkey"].isin(keys)]["method_id"])
    held_all = set(md[md["colkey"].isin(keys)]["method_id"]) if A.strict else held
    is_held = rows["method_id"].isin(held).values
    is_out = rows["method_id"].isin(held_all).values
    cal = {m: v for m, v in cal_all.items() if m not in held}
    y = D.targets(cal)
    tr = D.in_tr & ~is_out & ~D.test_c & np.isfinite(y)
    yp = fit(D.X[tr], y[tr], D.w[tr], D.mono, trees).predict(D.X)
    sel = is_held & D.is_grad
    rtp = D.predict_rt(yp, cal, sel)
    for tag, s in (("all_compounds", sel), ("unseen_compounds", sel & D.test_c)):
        r = metrics_rt(rows[s].reset_index(drop=True), rtp[s], md)
        r["fold"] = k; r["set"] = tag
        res.append(r)
    C.log("分割 %d: 外したカラム %d、分析法 %d" % (k, len(keys), len(held)))
R = pd.concat(res)
SUF = "_strict" if A.strict else ""
R.to_csv(os.path.join(O, "eval_unseen_columns%s.csv" % SUF), index=False)
out = {}
for tag, g in R.groupby("set"):
    out[tag] = dict(n_methods=int(len(g)), n_meas=int(g["n"].sum()), n_columns=int(g["column"].nunique()), median_method_mae_min=float(g["mae_min"].median()),
                    median_rel_mae=float(g["rel_mae"].median()), median_spearman=float(g["spearman"].median()), spearman_q25=float(g["spearman"].quantile(0.25)))
out["trees"] = trees
C.dump_json(out, os.path.join(O, "metrics_unseen_columns%s.json" % SUF), compact=False)
print(json.dumps(out, indent=1, ensure_ascii=False))
