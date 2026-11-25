# SPDX-License-Identifier: MIT
"""グラフのモデル用の予測の幅（確かさの目安。この写しでは未実行。GPU の交差検証の結果が要る）。木のモデル用（run.py の予測の幅）と同じ考え方・同じ形。
  python train/gpu/train_uncertainty_gnn.py <交差検証のモデルの名前（run_gpu_cv.sh の gSD）> <名前（free など）> <データのフォルダ> <runs のフォルダ> <重みのフォルダ>
  - 外れは、交差検証（化合物のかたまりの 4 分割）での「自分を学習に使っていないモデル」の予測と、実測から逆算した lkw の差。
  - その大きさを、ページで作れる 153 個の入力から分位点回帰（木）で当てる。幅のモデル自身も同じ 4 分割で作り、分割の外で 68%・90% が入るように倍率を決める。
出力: <重みのフォルダ>/uncertainty_free_gnn.json（retention_model.js の create の 3 つ目に、木の版の uncertainty_free.json の代わりに渡す）、
      <重みのフォルダ>/../uncertainty_gnn_metrics.json
"""
import os, sys, json, glob
import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
KIND, GROUP = sys.argv[1], sys.argv[2]
P, RUNS, WD = (os.path.abspath(x) for x in sys.argv[3:6])
NF = 4
R = np.load(os.path.join(P, "rows.npz"))
meta = json.load(open(os.path.join(P, "meta.json"), encoding="utf-8"))
ci, mi = R["ci"], R["mi"]
X = np.hstack([R["XC"][ci], R["IC"][ci], R["M"][mi], R["IF"]]).astype(np.float32)
YR = R["y_raw_" + GROUP]; fold = R["fold_" + GROUP]; src = R["src"]
oof = np.full(len(YR), np.nan, np.float32)
for k in range(NF):
    ps = [np.load(os.path.join(d, "pred.npy")) for d in sorted(glob.glob(os.path.join(RUNS, "%s_%s_f%d_s*" % (KIND, GROUP, k))))]
    assert ps, (KIND, k)
    oof[fold == k] = np.mean(ps, 0)[fold == k]
oof = np.clip(oof, -3.0, 14.0)
# 打ち切りの行（教師が上限・下限）は、外れの大きさが分からないので使わない
tr = R["in_" + GROUP] & np.isfinite(YR) & np.isfinite(oof) & (YR < meta["y_ceil"] - 1e-6) & (YR > meta["y_floor"] + 1e-6)
res = np.abs(YR - oof)
w = R["w"]


def qmodel(q):
    return HistGradientBoostingRegressor(loss="quantile", quantile=q, max_iter=200, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=50, l2_regularization=1.0,
                                         random_state=0, early_stopping=False)


out = dict(kind=KIND, group=GROUP, n_rows=int(tr.sum()), oof_mae=float(np.mean(res[tr])), oof_q68=float(np.quantile(res[tr], .68)), oof_q90=float(np.quantile(res[tr], .9)), quantiles={})
models = {}
for q in [0.68, 0.90]:
    s_oof = np.full(len(YR), np.nan)
    for k in range(NF):
        a = tr & (fold != k); b = tr & (fold == k)
        s_oof[b] = qmodel(q).fit(X[a], res[a], sample_weight=w[a]).predict(X[b])
    s_oof = np.maximum(s_oof, 0.02)
    scale = float(np.quantile(res[tr] / s_oof[tr], q))
    s = s_oof * scale
    m = qmodel(q).fit(X[tr], res[tr], sample_weight=w[tr])
    models[q] = dict(model=m, scale=scale)
    order = np.argsort(s[tr]); parts = np.array_split(order, 5)
    bins = [dict(mean_predicted=float(s[tr][ix].mean()), mean_abs_error=float(res[tr][ix].mean()), coverage=float(np.mean(res[tr][ix] <= s[tr][ix]))) for ix in parts]
    by = {meta["sources"][k]: float(np.mean(res[tr & (src == k)] <= s[tr & (src == k)])) for k in range(len(meta["sources"])) if (tr & (src == k)).sum() > 50}
    out["quantiles"]["%.2f" % q] = dict(scale=scale, coverage=float(np.mean(res[tr] <= s[tr])), coverage_by_source=by, spearman=float(spearmanr(s[tr], res[tr]).statistic),
                                        bins=bins, median_width=float(np.median(s[tr])), q10_width=float(np.quantile(s[tr], .1)), q90_width=float(np.quantile(s[tr], .9)))
    print("q=%.2f scale %.3f | 分割の外で入る割合 %.3f | 順位の相関 %.3f | 幅の中央値 %.2f | 出所ごと %s | %s" % (
        q, scale, out["quantiles"]["%.2f" % q]["coverage"], out["quantiles"]["%.2f" % q]["spearman"], np.median(s[tr]),
        {k[:6]: round(v, 2) for k, v in by.items()}, " ".join("%.2f→%.2f" % (b["mean_predicted"], b["mean_abs_error"]) for b in bins)), flush=True)
json.dump(out, open(os.path.join(WD, "..", "uncertainty_gnn_metrics.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)


def dump(m):
    trees = []
    for (pred,) in m._predictors:
        nd = pred.nodes; leaf = nd["is_leaf"].astype(bool)
        trees.append(dict(f=[-1 if l else int(f) for l, f in zip(leaf, nd["feature_idx"])],
                          t=[0 if l else (1e300 if x == np.inf else -1e300 if x == -np.inf else float(x)) for l, x in zip(leaf, nd["num_threshold"])],
                          m=[int(x) for x in nd["missing_go_to_left"]], l=[int(x) for x in nd["left"]], r=[int(x) for x in nd["right"]],
                          v=[round(float(x), 6) if l else 0 for l, x in zip(leaf, nd["value"])]))
    return dict(baseline=float(np.ravel(m._baseline_prediction)[0]), trees=trees)


# 学習に入っている化合物（ストックの多く）のための倍率: 本番のモデルは全部の行で学習しているので、交差検証の外れより小さい。
# 木のモデルの版と同じ倍率（約 0.5）を、目安としてそのまま使う。
old = json.load(open(os.path.join(WD, "uncertainty_free.json"), encoding="utf-8"))["meta"]
ex = dict(meta=dict(name=GROUP, feature_version=meta["fv"], n_features=int(X.shape[1]), floor=0.02, **{"for": "gnn"}, cv_model=KIND,
                    what="グラフのモデルの lkw の予測の外れの大きさの見積もり。s68: 外れの 68% がこの値より小さい、s90: 90% がこの値より小さい",
                    check=dict(coverage_68=out["quantiles"]["0.68"]["coverage"], coverage_90=out["quantiles"]["0.90"]["coverage"], spearman=out["quantiles"]["0.68"]["spearman"]),
                    seen_scale=old.get("seen_scale"), seen_stock=old.get("seen_stock")),
          s68=dict(scale=models[0.68]["scale"], **dump(models[0.68]["model"])), s90=dict(scale=models[0.90]["scale"], **dump(models[0.90]["model"])))
p = os.path.join(WD, "uncertainty_free_gnn.json")
json.dump(ex, open(p, "w", encoding="utf-8"), ensure_ascii=False, separators=(",", ":"))
print("written", p, "%.2f MB" % (os.path.getsize(p) / 1e6))
