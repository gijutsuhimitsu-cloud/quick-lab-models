# SPDX-License-Identifier: MIT
"""木のモデル（scikit-learn の HistGradientBoostingRegressor）を、ページ用の JSON の形に書き出す・読む。

形（js/retention_model.js の evalTrees と同じ）: trees[k] = { f: 特徴量の番号（葉は -1）, t: しきい値, m: 欠測を左へ送るか,
  l・r: 子の番号, v: 葉の値（小数 6 桁） }。x[f] <= t なら左。予測 = baseline + Σ 葉の値。
pickle は使わない（読み込むだけでプログラムが動く形式は配らない、という決まり）。
"""
import numpy as np


def dump(model):
    """学習済みのモデル → (baseline, trees)。しきい値は float64 のまま（丸めると枝が変わる）。欠測だけを分ける枝の ±inf は ±1e300。"""
    trees = []
    for (pred,) in model._predictors:
        nd = pred.nodes
        leaf = nd["is_leaf"].astype(bool)
        trees.append(dict(
            f=[-1 if l else int(f) for l, f in zip(leaf, nd["feature_idx"])],
            t=[0 if l else (1e300 if x == np.inf else -1e300 if x == -np.inf else float(x)) for l, x in zip(leaf, nd["num_threshold"])],
            m=[int(x) for x in nd["missing_go_to_left"]],
            l=[int(x) for x in nd["left"]], r=[int(x) for x in nd["right"]],
            v=[round(float(x), 6) if l else 0 for l, x in zip(leaf, nd["value"])]))
    return float(np.ravel(model._baseline_prediction)[0]), trees


def predict(baseline, trees, X):
    """JSON の木で予測する（numpy でまとめて。Python と JavaScript の照合や、書き出した重みの評価に使う）。"""
    X = np.asarray(X, np.float64)
    n = len(X)
    out = np.full(n, float(baseline))
    rows = np.arange(n)
    for t in trees:
        f = np.asarray(t["f"], np.int64); th = np.asarray(t["t"], np.float64); ml = np.asarray(t["m"], bool)
        lc = np.asarray(t["l"], np.int64); rc = np.asarray(t["r"], np.int64); v = np.asarray(t["v"], np.float64)
        node = np.zeros(n, np.int64)
        act = rows
        while len(act):
            nd = node[act]
            fi = f[nd]
            inner = fi >= 0
            act, nd, fi = act[inner], nd[inner], fi[inner]
            if not len(act):
                break
            x = X[act, fi]
            left = np.where(np.isnan(x), ml[nd], x <= th[nd])
            node[act] = np.where(left, lc[nd], rc[nd])
        out += v[node]
    return out
