# SPDX-License-Identifier: MIT
"""学習の表（1 行 = 1 測定）と入力の行列、保持時間 ↔ lkw の変換、分析法ごとの補正。

入力（特徴量の第 4 版、153 個）: 化合物の記述子 41・官能基 56・イオン化の規則ごとの数 30、分析法 16
（カラムの性質 H, A, B, C, V, D, logkEB、固定相の種類 7、アセトニトリルの割合、温度）、pH から作るイオン化の量 10。
生の pH と添加剤は入力に入れない（pH はイオン化の量を通してだけ入る）。イオン化の量には向きの制約（電荷が増えると保持は上がらない）。
目的変数: 見かけの log kw（lkw）。log10 k = lkw − (a + b·lkw)·φ/(1 + B·φ)。勾配の保持時間からは数値積分で逆算する。
"""
import json, hashlib
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from . import chromlib2 as c2, ionfeat

FV = 4                                   # 特徴量の版（重みの meta.feature_version）
HSM_COLS = ["H", "A", "B", "C", "V", "D"]
PHASES = ["C18", "C8", "Phenyl", "CN", "EP", "F", "short"]
COND = ["frac_acn", "pH_weak", "temp_C", "add_formic", "add_acetic", "add_tfa", "add_ammonium", "add_basic", "add_phosphate"]
N_DESC, N_GRP = 41, 56
Y_FLOOR, Y_CEIL = -1.0, 9.0
W_ISO = 3.0                              # 等組成の行の重み
S_GRID = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]      # 分析法ごとの補正: 死時間の倍率
D_GRID = [0.0, 0.2, 0.5, 1.0, 2.0]                # 分析法ごとの補正: 装置の遅れ（分）
SMRT_FIXED = dict(s=1.0, d=3.8)                   # SMRT は補正を固定（遅れ 3.8 分は溶出の始まりからの見積もり）


def is_test_compound(ck, frac=0.1):
    """評価用に学習から外す化合物（InChIKey の前半のハッシュで 10%）。どの分析法の行も外れる。"""
    h = int(hashlib.md5(ck.encode()).hexdigest()[:8], 16)
    return (h % 1000) < frac * 1000


def comp_matrix(comp):
    cks = sorted(comp)
    idx = {c: i for i, c in enumerate(cks)}
    X = np.zeros((len(cks), N_DESC + N_GRP), np.float32)
    for c, i in idx.items():
        X[i, :N_DESC] = comp[c]["d"]
        X[i, N_DESC:] = comp[c]["g"]
    return cks, idx, X


def method_matrix(md):
    M = np.zeros((len(md), len(HSM_COLS) + 1 + len(PHASES) + len(COND)), np.float32)
    for j, c in enumerate(HSM_COLS + ["logkEB"]):
        M[:, j] = md[c].astype(float).values
    for j, p in enumerate(PHASES):
        M[:, len(HSM_COLS) + 1 + j] = (md["phase"] == p).astype(float).values
    for j, c in enumerate(COND):
        M[:, len(HSM_COLS) + 1 + len(PHASES) + j] = md[c].astype(float).values
    return M, HSM_COLS + ["logkEB"] + ["phase_" + p for p in PHASES] + COND


def par_of(m):
    return c2.PAR["acn" if not (m["frac_acn"] < 0.5) else "meoh"]


def metrics_rt(rows, rt_pred, md):
    """分析法ごとの保持時間の誤差と順番"""
    res = []
    mrow = md.set_index("method_id")
    for mid, ix in rows.groupby("method_id").indices.items():
        a, p = rows["rt_min"].values[ix], rt_pred[ix]
        ok = np.isfinite(a) & np.isfinite(p)
        if ok.sum() < 5:
            continue
        a, p = a[ok], p[ok]
        prog = json.loads(mrow.loc[mid, "gradient"])
        tg = max(t for t, _ in prog) if prog else np.nan
        res.append(dict(method_id=mid, source=mrow.loc[mid, "source"], column=mrow.loc[mid, "column_name"], n=int(ok.sum()),
                        mae_min=float(np.mean(np.abs(a - p))), medae_min=float(np.median(np.abs(a - p))),
                        run_min=float(tg), rel_mae=float(np.mean(np.abs(a - p)) / tg) if tg > 0 else np.nan,
                        spearman=float(spearmanr(a, p).statistic) if ok.sum() > 2 else np.nan))
    return pd.DataFrame(res)


class Data:
    """md: 分析法の表、rows: 測定の表、comp: {ck: 特徴量}、sources: 学習に入れる出所、holdout: 評価用に外す化合物の割合、
    own_weight: 自社のデータ（出所 own）の行の重みの倍率"""

    def __init__(self, md, rows, comp, sources, holdout=0.1, own_weight=1.0):
        md = md[~md["exclude"].fillna(False).astype(bool)].reset_index(drop=True)
        self.md, self.comp = md, comp
        self.cks, cidx, self.XC = comp_matrix(comp)
        self.M, self.mnames = method_matrix(md)
        midx = {m: i for i, m in enumerate(md["method_id"])}
        rows = rows[rows["method_id"].isin(midx)].reset_index(drop=True)
        self.rows = rows
        self.ci = rows["ck"].map(cidx).values
        if np.isnan(self.ci.astype(float)).any():
            raise SystemExit("測定の表に、化合物の表に無い ck がある")
        self.ci = self.ci.astype(np.int64)
        self.mi = rows["method_id"].map(midx).values.astype(np.int64)
        self.fv = FV
        # 並び: 化合物（記述子 41・官能基 56）、イオン化の規則ごとの数、分析法 16、イオン化の量 10
        self.IC = np.zeros((len(self.cks), len(ionfeat.rules())), np.float32)
        for i, c in enumerate(self.cks):
            if comp[c].get("ion_ok", True):
                self.IC[i] = comp[c]["ion"]
        pka, typ, perm = ionfeat.groups(self.IC)
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                IF = ionfeat.features(pka[self.ci], typ[self.ci], perm[self.ci], md["pH_weak"].values[self.mi].astype(float), ionfeat.CAP_NO_GROUP)
        keep = [i for i, n in enumerate(self.mnames) if n != "pH_weak" and not n.startswith("add_")]
        self.M = self.M[:, keep]; self.mnames = [self.mnames[i] for i in keep]
        IF = IF[:, ionfeat.keep(FV)]
        self.X = np.hstack([self.XC[self.ci], self.IC[self.ci], self.M[self.mi], IF])
        self.mono = np.zeros(self.X.shape[1], int)
        self.mono[self.X.shape[1] - IF.shape[1]:] = ionfeat.signs(FV)
        self.is_grad = rows["logk_iso"].isna().values
        self.holdout = holdout
        self.test_c = np.array([is_test_compound(c, holdout) for c in rows["ck"]]) if holdout > 0 else np.zeros(len(rows), bool)
        self.src = md["source"].values[self.mi]
        self.in_tr = np.isin(self.src, list(sources))
        self.w = np.where(self.is_grad, 1.0, W_ISO) * np.where(self.src == "own", float(own_weight), 1.0)
        self.mrow = md.set_index("method_id")
        self.groups = rows.groupby("method_id").indices
        self._cache = {}

    def target_method(self, mid, s=1.0, d=0.0):
        key = (mid, round(s, 3), round(d, 3))
        if key in self._cache:
            return self._cache[key]
        m = self.mrow.loc[mid]
        ix = self.groups[mid]
        par = par_of(m)
        if np.isfinite(m["isocratic_phi"]):
            y = c2.lkw_from_isocratic(self.rows["logk_iso"].values[ix], m["isocratic_phi"], par)
        else:
            y = c2.lkw_from_gradient(self.rows["rt_min"].values[ix], m["t0_min"] * s, json.loads(m["gradient"]), par, d, n=800, iters=40)
        self._cache[key] = y
        return y

    def targets(self, cal, clip=True):
        y = np.full(len(self.rows), np.nan)
        for mid, ix in self.groups.items():
            c = cal.get(mid, (1.0, 0.0))
            y[ix] = self.target_method(mid, c[0], c[1])
        return np.clip(y, Y_FLOOR, Y_CEIL) if clip else y

    def predict_rt(self, lkw, cal, sel=None):
        out = np.full(len(self.rows), np.nan)
        for mid, ix in self.groups.items():
            m = self.mrow.loc[mid]
            if np.isfinite(m["isocratic_phi"]):
                continue
            if sel is not None:
                ix = ix[sel[ix]]
                if not len(ix):
                    continue
            c = cal.get(mid, (1.0, 0.0))
            out[ix] = c2.rt_from_lkw(lkw[ix], m["t0_min"] * c[0], json.loads(m["gradient"]), par_of(m), c[1], n=800)
        return out

    def calibrate(self, lkw_model, use, methods):
        """分析法ごとに、死時間の倍率と装置の遅れを、予測の保持時間が実測に合うように格子から選ぶ（学習用の化合物だけで）。"""
        cal = {}
        for mid in methods:
            m = self.mrow.loc[mid]
            if np.isfinite(m["isocratic_phi"]):
                continue
            if mid == "smrt":
                cal[mid] = (SMRT_FIXED["s"], SMRT_FIXED["d"])
                continue
            if int(m.get("fixed_cal", 0) or 0) == 1:       # 遅れと死時間が記録されている分析法は補正しない
                continue
            ix = self.groups[mid]
            ix = ix[use[ix]]
            if len(ix) < 8:
                cal[mid] = (1.0, 0.0)
                continue
            if len(ix) > 300:
                ix = ix[np.random.default_rng(0).choice(len(ix), 300, replace=False)]
            rt = self.rows["rt_min"].values[ix]
            prog, par = json.loads(m["gradient"]), par_of(m)
            best = (1e9, 1.0, 0.0)
            for s in S_GRID:
                for d in D_GRID:
                    p = c2.rt_from_lkw(lkw_model[ix], m["t0_min"] * s, prog, par, d, n=400)
                    e = np.nanmedian(np.abs(p - rt)) + 1e-4 * (abs(s - 1) + d)
                    if e < best[0]:
                        best = (e, s, d)
            cal[mid] = (best[1], best[2])
        return cal


def fit(X, y, w, mono=None, trees=900, seed=0):
    from sklearn.ensemble import HistGradientBoostingRegressor
    m = HistGradientBoostingRegressor(max_iter=int(trees), learning_rate=0.08, max_leaf_nodes=63, min_samples_leaf=20,
                                      l2_regularization=1.0, random_state=seed, early_stopping=False,
                                      monotonic_cst=(mono if mono is not None and len(mono) == X.shape[1] else None))
    return m.fit(X, y, sample_weight=w)


def iso_metrics(D, lkw_pred, sel):
    """等組成の行の log k の誤差（アセトニトリルの定数で。Weber・HSM2 用）"""
    if not sel.any():
        return {}
    r = D.rows[sel]
    phi = D.md["isocratic_phi"].values[D.mi[sel]]
    lk = c2.logk(lkw_pred[sel], phi, c2.PAR["acn"])
    e = lk - r["logk_iso"].values
    df = pd.DataFrame(dict(mid=r["method_id"].values, e=e, obs=r["logk_iso"].values, pred=lk))
    sp = [spearmanr(g["obs"], g["pred"]).statistic for _, g in df.groupby("mid") if len(g) >= 6]
    return dict(n=int(sel.sum()), logk_mae=float(np.mean(np.abs(e))), logk_bias=float(np.mean(e)), logk_rmse=float(np.sqrt(np.mean(e ** 2))),
                within_factor_1p5=float(np.mean(np.abs(e) < np.log10(1.5))), within_factor_2=float(np.mean(np.abs(e) < np.log10(2))),
                within_factor_3=float(np.mean(np.abs(e) < np.log10(3))),
                spearman_median_per_condition=float(np.median(sp)) if sp else None, n_conditions=len(sp))


def summ(r):
    if r is None or not len(r):
        return {}
    w = r["n"].values
    return dict(n_methods=int(len(r)), n_meas=int(w.sum()), mae_min_weighted=float(np.average(r["mae_min"], weights=w)),
                median_method_mae_min=float(r["mae_min"].median()), median_rel_mae=float(r["rel_mae"].median()),
                median_spearman=float(r["spearman"].median()), spearman_q25=float(r["spearman"].quantile(0.25)))
