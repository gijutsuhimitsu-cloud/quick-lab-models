# SPDX-License-Identifier: MIT
"""学習し直さずに、重みを別のデータで評価する（自社のデータの形のフォルダ。形は OWN_DATA.md）。
  python train/scripts/eval_external.py <データのフォルダ> [--weights <重みのフォルダ>] [--cal <calibration.json>]
    既定の重みは weights/（公開の重み）。例: python train/scripts/eval_external.py train/examples/zhang_dataset10
出す数字（勾配の分析法ごと）: 保持時間の誤差（分・分析時間に対する割合）・溶出順の相関・分析時間の 5%／10% 以内の割合・偏り。
  - all: 全部の化合物 / unseen: 公開の重みの学習に使っていない化合物（InChIKey のハッシュで 10%）/ seen: それ以外
  - 補正なし（死時間は表の値、遅れ 0。遅れを書いた分析法は、その値）。--cal を渡すと、同じカラムの学習の分析法の補正の中央値を使った値も出す
出力: train/out/eval_<フォルダ名>.json（--out で変えられる）。データはどこにも送らない。
"""
import os, sys, json, argparse
import numpy as np, pandas as pd
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from lib import common as C, table, own as OWN, trees as TR, chromlib2 as c2
from lib.data import Data, is_test_compound, par_of

ap = argparse.ArgumentParser()
ap.add_argument("data")
ap.add_argument("--weights", default=C.WEIGHTS)
ap.add_argument("--cal", default="", help="学習の時の分析法の補正（train/work/<名前>/calibration.json）。同じカラムの補正の中央値を参考に使う")
ap.add_argument("--out", default="")
A = ap.parse_args()
table.load_par()
model = json.load(open(os.path.join(A.weights, "model_free.json"), encoding="utf-8"))
for k in ("acn", "meoh"):
    c2.PAR[k] = dict(model["meta"]["par"][k])          # 重みに書いてある式の定数を使う
bmd, _, comp = table.load_base()
name = os.path.basename(os.path.normpath(A.data))
work = os.path.join(C.WORK, "eval_" + name)
md, rows, new, names, rep = OWN.load(A.data, work, comp)
comp = {k: v for k, v in comp.items() if k in set(rows["ck"])}
comp.update(new)
D = Data(md, rows, comp, sources=[], holdout=0.1)
assert D.X.shape[1] == model["n_features"]
lkw = TR.predict(model["baseline"], model["trees"], D.X)
cal_all = {k: tuple(v) for k, v in json.load(open(A.cal, encoding="utf-8")).items()} if A.cal else {}


def stats(a, p, tg):
    if len(a) < 5:
        return dict(n=int(len(a)))
    e = np.abs(a - p)
    return dict(n=int(len(a)), mae_min=float(e.mean()), medae_min=float(np.median(e)), rel_mae=float(e.mean() / tg), spearman=float(spearmanr(a, p).statistic),
                within_5pct_run=float(np.mean(e <= 0.05 * tg)), within_10pct_run=float(np.mean(e <= 0.10 * tg)), bias_min=float(np.mean(p - a)))


out = dict(data=name, weights=os.path.relpath(A.weights, C.MODEL), model=model["meta"]["name"], n_trees=len(model["trees"]), data_report=rep, methods={})
unseen = np.array([is_test_compound(c) for c in D.rows["ck"]])
for mid, ix in D.groups.items():
    m = D.mrow.loc[mid]
    sel = np.zeros(len(D.rows), bool); sel[ix] = True
    r = dict(column=m["column_name"], n=int(len(ix)))
    if np.isfinite(m["isocratic_phi"]):
        e = c2.logk(lkw[ix], m["isocratic_phi"], par_of(m)) - D.rows["logk_iso"].values[ix]
        r["all"] = dict(n=int(len(ix)), logk_mae=float(np.mean(np.abs(e))), logk_bias=float(np.mean(e)))
    else:
        tg = max(t for t, _ in json.loads(m["gradient"]))
        a = D.rows["rt_min"].values
        p0 = D.predict_rt(lkw, {}, sel)
        r["no_calibration"] = {k: stats(a[s], p0[s], tg) for k, s in (("all", sel), ("unseen", sel & unseen), ("seen", sel & ~unseen))}
        if cal_all and m["column_id"]:
            same = set(bmd["method_id"][bmd["column_id"] == m["column_id"]])
            cc = [v for k, v in cal_all.items() if k in same]
            if cc:
                cs, cd = float(np.median([c[0] for c in cc])), float(np.median([c[1] for c in cc]))
                p1 = D.predict_rt(lkw, {mid: (cs, cd)}, sel)
                r["calibration_from_same_column"] = dict(t0_scale=cs, dwell_min=cd, n_methods=len(cc), **{k: stats(a[s], p1[s], tg) for k, s in (("all", sel), ("unseen", sel & unseen))})
    out["methods"][mid] = r
p = A.out or os.path.join(C.OUT, "eval_%s.json" % name)
C.dump_json(out, p, compact=False)
for mid, r in out["methods"].items():
    if "no_calibration" in r:
        s = r["no_calibration"]
        print("%s（%s）: 全部 %d 化合物 %.2f 分（%.1f%%、順位 %.3f）／学習に使っていない化合物 %d 個 %.1f%%" % (
            mid, r["column"], s["all"]["n"], s["all"]["mae_min"], 100 * s["all"]["rel_mae"], s["all"]["spearman"], s["unseen"]["n"], 100 * s["unseen"].get("rel_mae", np.nan)))
        if "calibration_from_same_column" in r:
            q = r["calibration_from_same_column"]
            print("  同じカラムの補正の中央値（死時間 ×%.2f・遅れ %.1f 分）: 全部 %.1f%%" % (q["t0_scale"], q["dwell_min"], 100 * q["all"]["rel_mae"]))
    else:
        print("%s（%s）: 等組成 %d 点、log k の誤差 %.3f" % (mid, r["column"], r["all"]["n"], r["all"]["logk_mae"]))
print("書いた:", p)
