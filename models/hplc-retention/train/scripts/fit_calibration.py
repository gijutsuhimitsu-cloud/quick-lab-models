# SPDX-License-Identifier: MIT
"""「測定値 1 点で合わせ直す」規則の係数を、木のモデルの交差検証の外れから当てはめ直す（任意。CPU で数分）。
  python train/scripts/fit_calibration.py [名前（既定 free）]
前提: python train/run.py を同じ名前で実行してある（train/work/<名前>/ の oof_tree.npy・fold_c.npy・calibration.json を使う）。

ページ（js/retention_model.js の calibrate）は、利用者が 1 点測った保持時間から lkw の外れ d を求め、同じ化合物の別の条件の予測に w·d を足す:
  w = w0 ×（カラムが別なら col）×（有機溶媒が別なら solv）×（floor + (1 − floor)·exp(−dq/tau)）、dq = |Δ負の電荷| + |Δ正の電荷|
同じ化合物・同じ出所（＝同じ装置）の中で、合わせ直しに使う 1 行と当てる行の組を作り、交差検証の外れ（自分を学習に使っていないモデルの外れ）が
どれだけ持ち越されるかに、上の式を最小二乗で当てはめる（化合物 × 出所ごとに同じ重み）。化合物の 4 分割の外でも確かめる。
出力: train/out/<名前>/calibration_rule.json（係数）と calibration_rule.md（表）。
  係数をページで使うには、uncertainty_free.json の meta.calibration に {w0, col, solv, floor, tau} を入れる（retention_model.js が読む。無ければ既定の値）。
"""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from lib import common as C
from lib.runs import data_for_run
from lib.data import Y_FLOOR, Y_CEIL

name = sys.argv[1] if len(sys.argv) > 1 else "free"
W, O = os.path.join(C.WORK, name), os.path.join(C.OUT, name)
D, cfg, cal = data_for_run(name)
oof = np.load(os.path.join(W, "oof_tree.npy"), allow_pickle=False)
foldc = np.load(os.path.join(W, "fold_c.npy"), allow_pickle=False).astype(int)
assert len(oof) == len(D.rows) and len(foldc) == len(D.cks), "work/ の物と学習の表が合わない（run.py を同じ設定で実行し直す）"
YR = D.targets(cal, clip=False)
p = np.clip(oof, -3, 14)
ing = D.in_tr & np.isfinite(YR) & np.isfinite(p)
unc = ing & (YR < Y_CEIL - 1e-6) & (YR > Y_FLOOR + 1e-6)          # 目的変数の上限・下限に当たった行は、合わせ直しに使わない
res = YR - p
SRC = sorted(set(D.src)); src = np.array([SRC.index(s) for s in D.src])
MN = D.mnames
colsig = np.round(D.M[:, :7], 4)                                  # H, A, B, C, V, D, logkEB
_, colid = np.unique(colsig, axis=0, return_inverse=True)
colid = np.asarray(colid).ravel()[D.mi]
acn = D.M[D.mi, MN.index("frac_acn")]
IF = D.X[:, -10:]
qneg, qpos = IF[:, 0], IF[:, 1]
ci = D.ci
grp = ci.astype(np.int64) * 16 + src


def draw(seed):
    """化合物 × 出所ごとに、合わせ直しに使う行を 1 つ選ぶ（打ち切りでない行から）。各行にとっての合わせ直しの行（無ければ −1）"""
    rng = np.random.default_rng(seed)
    u = rng.random(len(YR)); u[~unc] = -1
    o = np.lexsort((u, grp))
    last = np.r_[grp[o][1:] != grp[o][:-1], True]
    best = o[last]
    ok = u[best] >= 0
    c = pd.Series(best[ok], index=grp[best[ok]]).reindex(grp).values
    c = np.where(np.isnan(c), -1, c).astype(np.int64)
    c[c == np.arange(len(YR))] = -1
    c[~ing] = -1
    return c


rows_ = []
for seed in range(6):
    c = draw(seed); j = np.where((c >= 0) & unc)[0]; cc = c[j]
    rows_.append(pd.DataFrame(dict(j=j, c=cc, rj=res[j], rc=res[cc], src=src[j], ck=ci[j], seed=seed, same_col=colid[j] == colid[cc], same_solv=acn[j] == acn[cc],
                                   dq=np.abs(qneg[j] - qneg[cc]) + np.abs(qpos[j] - qpos[cc]), same_method=D.mi[j] == D.mi[cc])))
T = pd.concat(rows_, ignore_index=True)
T = T[~T["same_method"]].reset_index(drop=True)
T["wt"] = 1.0 / T.groupby(["ck", "src", "seed"])["rj"].transform("size")
C.log("組 %d（化合物 × 出所 %d）" % (len(T), T.groupby(["ck", "src"]).ngroups))


def wrule(g, th):
    w0, col, solv, flo, tau = th
    return w0 * np.where(g["same_col"], 1.0, col) * np.where(g["same_solv"], 1.0, solv) * (flo + (1 - flo) * np.exp(-np.asarray(g["dq"]) / tau))


from scipy.optimize import least_squares


def fitrule(mask):
    g = T[mask]; sw = np.sqrt(g.wt.values)
    return least_squares(lambda th: sw * (g.rj.values - wrule(g, th) * g.rc.values), [0.85, 0.85, 0.7, 0.6, 0.2],
                         bounds=([0, 0, 0, 0, 0.02], [1, 1, 1, 1, 3])).x


mae = lambda g, e: float((g.wt * np.abs(e)).sum() / g.wt.sum())
th = fitrule(np.ones(len(T), bool))
fc = foldc[T["ck"].values]
e_cv = np.zeros(len(T)); per_fold = []
for k in range(4):
    m = fc == k
    if not m.any():
        continue
    thk = fitrule(~m); per_fold.append([float(x) for x in thk])
    e_cv[m] = T.rj.values[m] - wrule(T[m], thk) * T.rc.values[m]
DEFAULT = (0.87, 0.88, 0.72, 0.68, 0.06)          # js/retention_model.js の CAL_DEFAULT
T["e_none"] = T.rj; T["e_full"] = T.rj - T.rc; T["e_rule"] = e_cv; T["e_default"] = T.rj - wrule(T, DEFAULT) * T.rc
L = ["# 測定値 1 点で合わせ直す規則（%s、木のモデルの交差検証の外れから）" % name, "",
     "係数: w0 %.3f・カラムが別 %.3f・溶媒が別 %.3f・下限 %.3f・tau %.3f（ページの既定は 0.87・0.88・0.72・0.68・0.06）" % tuple(th), "",
     "分割ごとに当てはめ直した係数: " + "；".join("(%s)" % ", ".join("%.3f" % x for x in f) for f in per_fold), "",
     "lkw の誤差（化合物 × 出所ごとに同じ重み）", "", "| 出所 | 化合物 | 合わせ直す前 | 外れを全部足す | ページの既定の係数 | 規則（分割の外で当てはめ） |", "|---|---|---|---|---|---|"]
by = {}
for s in sorted(T["src"].unique()):
    g = T[T["src"] == s]
    by[SRC[s]] = dict(n_compounds=int(g.ck.nunique()), none=mae(g, g.e_none), full=mae(g, g.e_full), page_default=mae(g, g.e_default), rule_cv=mae(g, g.e_rule))
    L.append("| %s | %d | %.3f | %.3f | %.3f | %.3f |" % (SRC[s], g.ck.nunique(), mae(g, g.e_none), mae(g, g.e_full), mae(g, g.e_default), mae(g, g.e_rule)))
by["all"] = dict(n_compounds=int(T.ck.nunique()), none=mae(T, T.e_none), full=mae(T, T.e_full), page_default=mae(T, T.e_default), rule_cv=mae(T, T.e_rule))
L.append("| 全部 | %d | %.3f | %.3f | %.3f | %.3f |" % (T.ck.nunique(), mae(T, T.e_none), mae(T, T.e_full), mae(T, T.e_default), mae(T, T.e_rule)))
T["w"] = wrule(T, th); T["wb"] = pd.cut(T["w"], [0, .3, .4, .5, .6, .7, .8, 1.0])
L += ["", "合わせ直したあとの外れの大きさ ≒ 前 × sqrt(1 − w²)（ページの幅の縮め方）の確かめ", "", "| w の範囲 | 組 | 前の RMS | あとの RMS | 比 | sqrt(1−w²) |", "|---|---|---|---|---|---|"]
for b, g in T.groupby("wb", observed=True):
    r0 = np.sqrt((g.wt * g.rj ** 2).sum() / g.wt.sum()); r1 = np.sqrt((g.wt * g.e_rule ** 2).sum() / g.wt.sum())
    L.append("| %s | %d | %.3f | %.3f | %.2f | %.2f |" % (b, len(g), r0, r1, r1 / r0, np.sqrt(1 - ((g.wt * g.w).sum() / g.wt.sum()) ** 2)))
rule = dict(w0=float(th[0]), col=float(th[1]), solv=float(th[2]), floor=float(th[3]), tau=float(th[4]))
C.dump_json(dict(name=name, model="tree", rule=rule, rounded={k: round(v, 2) for k, v in rule.items()}, per_fold=per_fold, by_source=by,
                 n_pairs=int(len(T)), page_default=dict(w0=0.87, col=0.88, solv=0.72, floor=0.68, tau=0.06),
                 note="w = w0 × (カラムが別なら col) × (溶媒が別なら solv) × (floor + (1 − floor)·exp(−dq/tau))。ページで使うには uncertainty_free.json の meta.calibration に入れる"),
             os.path.join(O, "calibration_rule.json"), compact=False)
open(os.path.join(O, "calibration_rule.md"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("\n".join(L))
