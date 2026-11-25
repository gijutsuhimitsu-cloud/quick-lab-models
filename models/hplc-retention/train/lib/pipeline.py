# SPDX-License-Identifier: MIT
"""学習から書き出しまで（1 つの Python の中で続けて実行する。途中のモデルを pickle に保存しないため）。

  fit_main        木 900 本の本体。分析法ごとの補正（死時間の倍率・装置の遅れ）を 3 回まわして決める。評価 A〜C
  export          重みの JSON（model_free.json）と照合用の値（reference_free.json）、ストックの化合物
  uncertainty     予測の幅（分位点回帰の木 200 本 × 2）。化合物のかたまりの 4 分割の外れから
  seen_flags      学習データに入っている化合物の幅の倍率と、ストックの印
  name_index      名前 → 構造の索引
"""
import os, json, time, hashlib, collections, re
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from . import common as C, chromlib2 as c2, ionfeat, cluster, trees as TR
from .data import (Data, fit, iso_metrics, summ, metrics_rt, comp_matrix, FV, W_ISO, Y_FLOOR, Y_CEIL, SMRT_FIXED, HSM_COLS)

LIC_FREE = ("学習データ: METLIN SMRT〈Apache 2.0〉、HSM2 database・Weber の等組成データ・MCaRT・Wiczling グループの 5 カラムのデータ〈CC BY 4.0〉、"
            "Zhang ら〈CC0 1.0〉。RepoRT から読んだ保持値・構造は入っていない。出典の表示とライセンスの一覧は ATTRIBUTION.md")


def license_text(opt):
    if not opt.own:
        return LIC_FREE
    if opt.own_only:
        return "自社のデータだけで学習した重み（公開の基礎データは使っていない）。扱いはデータの持ち主が決める"
    return LIC_FREE.replace("学習データ: ", "学習データ: 公開の基礎データ（") .replace("。RepoRT", "）に、自社のデータを足した。RepoRT") + \
        "。配る時は ATTRIBUTION.md を付け、自社のデータを足したことを書き足す"


# ---------------------------------------------------------------- 本体
def fit_main(D, opt):
    md, rows = D.md, D.rows
    in_tr = D.in_tr
    is_web = D.src == "weber_isocratic_mendeley"
    w = D.w
    tr_methods = sorted(set(rows["method_id"][in_tr]))
    cal = {"smrt": (SMRT_FIXED["s"], SMRT_FIXED["d"])}
    log = dict(name=opt.name, par=c2.PAR, trees=opt.trees, holdout=opt.holdout, sources=sorted(set(D.src[in_tr])), rounds=[])
    model = None
    for rnd in range(3):
        t = time.time()
        y = D.targets(cal)
        ok = np.isfinite(y)
        tr = in_tr & ~D.test_c & ok
        model = fit(D.X[tr], y[tr], w[tr], D.mono, opt.trees)
        yp = model.predict(D.X)
        te = in_tr & D.test_c & ok
        e = float(np.mean(np.abs(y[te] - yp[te]))) if te.any() else float("nan")
        C.log("本体 %d 回目: 学習の行 %d、評価用の化合物の lkw の誤差 %.3f（%.0f 秒）" % (rnd + 1, tr.sum(), e, time.time() - t))
        log["rounds"].append(dict(round=rnd, n_train=int(tr.sum()), heldout_lkw_mae=e))
        if rnd < 2:
            cal = D.calibrate(yp, in_tr & ~D.test_c, tr_methods)
    y = D.targets(cal)
    yp = model.predict(D.X)
    out = opt.outdir
    # A. 学習に使っていない化合物（勾配の分析法の保持時間）
    selA = in_tr & D.test_c & D.is_grad
    if selA.any():
        rtp = D.predict_rt(yp, cal, selA)
        rA = metrics_rt(rows[selA].reset_index(drop=True), rtp[selA], md)
        rA.to_csv(os.path.join(out, "eval_heldout_compounds.csv"), index=False)
        log["heldout_compounds"] = summ(rA)
        log["heldout_compounds_by_source"] = {s: summ(g) for s, g in rA.groupby("source")}
    # B. 学習に使っていない分析法（基礎データだけの時は無い）
    selB = (~in_tr) & D.is_grad
    if selB.any():
        rtp = D.predict_rt(yp, {}, selB)
        rB = metrics_rt(rows[selB].reset_index(drop=True), rtp[selB], md)
        rB.to_csv(os.path.join(out, "eval_unseen_methods.csv"), index=False)
        log["unseen_methods"] = summ(rB)
    # C. 等組成（Weber）: 化合物を 5 つに分けて外す交差検証
    if (in_tr & is_web).any() and not opt.skip_weber_cv:
        cks = sorted(set(rows["ck"][is_web]))
        fold = {c: int(hashlib.md5(c.encode()).hexdigest()[:6], 16) % 5 for c in cks}
        pred = np.full(len(rows), np.nan)
        for k in range(5):
            held = rows["ck"].map(lambda c: fold.get(c, -1) == k).values     # その化合物は、どの出所の行も外す
            tr = in_tr & ~held & np.isfinite(y)
            pred[is_web & held] = fit(D.X[tr], y[tr], w[tr], D.mono, opt.trees).predict(D.X[is_web & held])
        log["isocratic_weber_cv_by_compound"] = iso_metrics(D, pred, is_web)
        log["isocratic_weber_cv_by_compound"]["n_compounds"] = len(cks)
        C.log("等組成（Weber）の交差検証: log k の誤差 %.3f" % log["isocratic_weber_cv_by_compound"]["logk_mae"])
    elif is_web.any() and not (in_tr & is_web).any():
        log["isocratic_weber_external"] = iso_metrics(D, yp, is_web)
    iso_h = in_tr & ~D.is_grad & (D.src == "hsm2_mendeley")
    if iso_h.any():
        log["hsm2_fit_logk_mae"] = iso_metrics(D, yp, iso_h)["logk_mae"]
    log["calibration"] = {k: v for k, v in cal.items() if v != (1.0, 0.0)}
    log["features"] = dict(n_comp=int(D.XC.shape[1]), method=D.mnames)
    C.dump_json(log, os.path.join(out, "metrics.json"), compact=False)
    C.dump_json({k: list(v) for k, v in cal.items()}, os.path.join(opt.workdir, "calibration.json"), compact=False)
    return model, cal, log, y, yp


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def _intify(v):
    v = float(v)
    return int(v) if v.is_integer() else v


def stock_compounds():
    """base_data/stock.csv（こちらで作った名前 → 構造の表の 102 種と、その特徴量）"""
    from .table import feature_names, read_csv
    desc, grp, ion = feature_names()
    S = read_csv(os.path.join(C.BASE, "stock.csv"), dtype={"name": str, "smiles": str, "inchikey": str}, keep_default_na=False)
    return [dict(name=r["name"], smiles=r["smiles"], inchikey=r["inchikey"], desc=[_intify(r["d_" + x]) for x in desc],
                 groups=[int(r["g_" + x]) for x in grp], ion=[int(r["ion_" + x]) for x in ion]) for _, r in S.iterrows()]


# ---------------------------------------------------------------- 書き出し
def export(D, model, log, opt):
    W = opt.weights_dir
    base, trees = TR.dump(model)
    cks, cidx, XC = D.cks, {c: i for i, c in enumerate(D.cks)}, D.XC
    _ok = D.in_tr
    mw = XC[[cidx[c] for c in sorted(set(D.rows["ck"][_ok]))], 0]          # 分子量の範囲は、学習に入っている化合物だけから
    _ok2 = _ok & (D.src != "smrt_figshare")                                 # SMRT（分子量 180〜500 に集中）に引きずられないよう、SMRT 以外の範囲と合わせた広い方
    mw2 = XC[[cidx[c] for c in sorted(set(D.rows["ck"][_ok2]))], 0] if _ok2.any() else mw
    out = dict(meta=dict(name=opt.name, version=2, license=license_text(opt), feature_version=FV,
                         target="lkw: log10 k = lkw - (a + b*lkw) * phi/(1 + B*phi)", par=c2.PAR,
                         metrics=dict(heldout_compounds=log.get("heldout_compounds"), isocratic=log.get("isocratic_weber_cv_by_compound"))),
               baseline=base, n_comp=int(XC.shape[1]), n_ion=len(ionfeat.rules()),
               n_features=int(XC.shape[1] + len(D.mnames) + len(ionfeat.rules()) + len(ionfeat.keep(FV))), method_features=D.mnames,
               range=dict(mw=[float(min(np.quantile(mw, 0.005), np.quantile(mw2, 0.005))), float(max(np.quantile(mw, 0.995), np.quantile(mw2, 0.995)))]),
               trees=trees)
    if opt.own:
        out["meta"]["trained_on"] = dict(base=not opt.own_only, own_rows=int((D.src == "own").sum()), own_weight=opt.own_weight)
    out["meta"]["ion_rules"] = [[n, t, p] for n, t, p in ionfeat.rules()]
    out["meta"]["ion_const"] = dict(max_groups=ionfeat.MAX_GROUPS, max_rep=ionfeat.MAX_REP, beta_shift=ionfeat.BETA_SHIFT, ph_default=ionfeat.PH_DEFAULT,
                                    cap_no_group=ionfeat.CAP_NO_GROUP, keep=ionfeat.keep(FV), monotone=ionfeat.signs(FV),
                                    amines=list(ionfeat.AMINES), names=[ionfeat.NAMES[i] for i in ionfeat.keep(FV)])
    out = _clean(out)
    p = os.path.join(W, "model_free.json")
    os.makedirs(W, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    C.log("重み %s（%.2f MB、木 %d 本）" % (os.path.relpath(p, C.TRAIN), os.path.getsize(p) / 1e6, len(trees)))
    # カラムの表は学習しない（HSM2 の値を写しただけ）。同じ物をそろえて置く
    cols = json.load(open(os.path.join(C.WEIGHTS, "columns.json"), encoding="utf-8"))
    with open(os.path.join(W, "columns.json"), "w", encoding="utf-8") as f:
        json.dump(cols, f, ensure_ascii=False, separators=(",", ":"))
    stock = stock_compounds()
    # 照合用: 学習の表からランダムな 300 行の入力と Python の予測、勾配の計算、ページの入口（predictRetention）を通した例
    rng = np.random.default_rng(0)
    X = D.X[rng.integers(0, len(D.X), 300)]
    assert X.shape[1] == out["n_features"]
    yp = model.predict(X)
    progs = [[[0, 0.05], [10, 0.95], [12, 0.95]], [[0, 0.1], [2, 0.1], [20, 0.6], [25, 1.0]], [[0, 0.3]]]
    rt = []
    for i, pr in enumerate(progs):
        lk = np.linspace(0.2, 6.0, 12)
        r = 0.8 * (1 + 10 ** c2.logk(lk, pr[0][1], c2.PAR["acn"])) if len(pr) == 1 else c2.rt_from_lkw(lk, 0.8, pr, c2.PAR["acn"], 0.3 * i, n=20000)
        rt.append(dict(prog=pr, t0=0.8, dwell=0.3 * i if len(pr) > 1 else 0.0, lkw=[float(v) for v in lk], rt=[float(v) for v in r]))
    PH = ["C18", "C8", "Phenyl", "CN", "EP", "F", "short"]
    ADDS = ["formic", "acetic", "tfa", "ammonium", "basic", "phosphate"]
    conds = [dict(organic="acn", phi=0.5, pH=2.8, tempC=35, additive="phosphate"), dict(organic="acn", phi=0.3, pH=3.0, tempC=40, additive="formic"),
             dict(organic="meoh", phi=0.6, pH=6.8, tempC=30, additive="ammonium"), dict(organic="acn", phi=0.4, pH=9.7, tempC=25, additive="ammonium"),
             dict(organic="acn", phi=0.2, additive="formic"),
             dict(organic="acn", gradient=[[0, 0.05], [10, 0.95], [12, 0.95]], pH=2.7, tempC=40, additive="formic", t0=0.6, dwell=0.3),
             dict(organic="meoh", gradient=[[2.1, 0.05], [32.1, 0.8]], pH=4.9, tempC=25, additive="ammonium", t0=0.532)]
    dflt = {d["id"] for d in cols["defaults"]}
    colsel = [c for c in cols["columns"] if c["id"] in dflt]
    e2e, Xe = [], []
    for si in range(0, len(stock), 3):
        s = stock[si]
        for c in colsel:
            for cd in conds:
                mv = {k: c[k] for k in HSM_COLS}; mv["logkEB"] = c["logkEB"]
                mv.update({"phase_" + p: (1.0 if c["phase"] == p else 0.0) for p in PH})
                mv.update(frac_acn=(0.0 if cd["organic"] == "meoh" else 1.0), pH_weak=cd.get("pH", np.nan), temp_C=cd.get("tempC", np.nan))
                mv.update({"add_" + a: (1.0 if cd.get("additive", "formic") == a else 0.0) for a in ADDS})
                v = list(s["desc"]) + list(s["groups"]) + list(s["ion"]) + [mv[n] for n in D.mnames]
                g3 = ionfeat.groups(np.array([s["ion"]], float))
                with np.errstate(all="ignore"):
                    import warnings
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        v += [float(z) for z in ionfeat.features(g3[0], g3[1], g3[2], np.array([cd.get("pH", np.nan)], float), ionfeat.CAP_NO_GROUP)[0][ionfeat.keep(FV)]]
                Xe.append(v); e2e.append(dict(stock=si, column=c["id"], cond=cd))
    lk = model.predict(np.array(Xe, np.float32))
    for e, l in zip(e2e, lk):
        cd = e["cond"]; par = c2.PAR["meoh" if cd["organic"] == "meoh" else "acn"]
        e["lkw"] = float(l)
        if "gradient" in cd:
            e["tR"] = float(c2.rt_from_lkw(np.array([l]), cd["t0"], cd["gradient"], par, cd.get("dwell", 0.0), n=20000)[0])
        else:
            e["k"] = float(10 ** c2.logk(l, cd["phi"], par))
    with open(os.path.join(W, "reference_free.json"), "w") as f:
        json.dump(dict(x=[[None if not np.isfinite(v) else float(v) for v in row] for row in X], y=[float(v) for v in yp], rt=rt, e2e=_clean(e2e)), f)
    return stock


# ---------------------------------------------------------------- 予測の幅
def _qmodel(q, n):
    from sklearn.ensemble import HistGradientBoostingRegressor
    return HistGradientBoostingRegressor(loss="quantile", quantile=q, max_iter=int(n), learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=50,
                                         l2_regularization=1.0, random_state=0, early_stopping=False)


def uncertainty(D, model, cal, y, opt, NF=4):
    ok = np.isfinite(y)
    in_tr, w = D.in_tr, D.w
    tr = in_tr & ~D.test_c & ok
    te = in_tr & D.test_c & ok
    in_c = np.zeros(len(D.cks), bool); in_c[np.unique(D.ci[in_tr])] = True
    non_smrt = np.zeros(len(D.cks), bool); non_smrt[np.unique(D.ci[in_tr & (D.src != "smrt_figshare")])] = True
    fold_c, n_clu, n_used = cluster.folds(D.cks, D.comp, in_c, ~non_smrt, NF)
    C.log("予測の幅: 化合物のかたまり %d 個（%d 化合物から。この版の化合物は %d）" % (n_clu, n_used, in_c.sum()))
    fold = fold_c[D.ci]
    t0 = time.time()
    oof = np.full(len(y), np.nan)
    for k in range(NF):
        a = tr & (fold != k); b = tr & (fold == k)
        if b.any():
            oof[b] = fit(D.X[a], y[a], w[a], D.mono, opt.trees).predict(D.X[b])
        C.log("  分割 %d の外れ: 行 %d、lkw の誤差 %.3f（%.0f 秒）" % (k, b.sum(), np.mean(np.abs(oof[b] - y[b])) if b.any() else np.nan, time.time() - t0))
    np.save(os.path.join(opt.workdir, "oof_tree.npy"), oof, allow_pickle=False)
    np.save(os.path.join(opt.workdir, "fold_c.npy"), fold_c.astype(np.int8), allow_pickle=False)
    res = np.abs(y - oof)
    out = dict(name=opt.name, n_rows=int(tr.sum()), oof_mae=float(np.mean(res[tr])), oof_q68=float(np.quantile(res[tr], .68)), oof_q90=float(np.quantile(res[tr], .9)), quantiles={})
    models = {}
    p_final = model.predict(D.X)
    err_te = np.abs(y - p_final)
    ns = D.src != "smrt_figshare"
    for q in [0.68, 0.90]:
        s_oof = np.full(len(y), np.nan)
        for k in range(NF):
            a = tr & (fold != k); b = tr & (fold == k)
            if b.any():
                s_oof[b] = _qmodel(q, opt.qtrees).fit(D.X[a], res[a], sample_weight=w[a]).predict(D.X[b])
        s_oof = np.maximum(s_oof, 0.02)
        scale = float(np.quantile(res[tr] / s_oof[tr], q))          # 分割の外で q の割合が幅に入る倍率
        m = _qmodel(q, opt.qtrees).fit(D.X[tr], res[tr], sample_weight=w[tr])
        s = np.maximum(m.predict(D.X), 0.02) * scale
        models[q] = dict(model=m, scale=scale)
        qq = out["quantiles"]["%.2f" % q] = dict(scale=scale, oof_coverage_before_scale=float(np.mean(res[tr] <= s_oof[tr])))
        if te.any():
            order = np.argsort(s[te]); parts = np.array_split(order, 5)
            qq.update(test_coverage=float(np.mean(err_te[te] <= s[te])), test_coverage_not_smrt=float(np.mean(err_te[te & ns] <= s[te & ns])) if (te & ns).any() else None,
                      test_spearman=float(spearmanr(s[te], err_te[te]).statistic),
                      test_spearman_not_smrt=float(spearmanr(s[te & ns], err_te[te & ns]).statistic) if (te & ns).sum() > 2 else None,
                      test_bins=[dict(mean_predicted=float(s[te][ix].mean()), mean_abs_error=float(err_te[te][ix].mean()), coverage=float(np.mean(err_te[te][ix] <= s[te][ix])), n=int(len(ix))) for ix in parts],
                      median_width=float(np.median(s[te])), q10_width=float(np.quantile(s[te], .1)), q90_width=float(np.quantile(s[te], .9)))
            if (te & (D.src == "own")).any():
                o = te & (D.src == "own")
                qq["test_coverage_own"] = float(np.mean(err_te[o] <= s[o]))
        um = (~in_tr) & ok
        if um.sum() > 100:
            qq["unseen_methods_coverage"] = float(np.mean(err_te[um] <= s[um]))
            qq["unseen_methods_spearman"] = float(spearmanr(s[um], err_te[um]).statistic)
        C.log("  幅 q=%.2f: 倍率 %.3f、評価用の化合物で幅に入る割合 %s" % (q, scale, ("%.3f" % qq["test_coverage"]) if "test_coverage" in qq else "—"))
    C.dump_json(out, os.path.join(opt.outdir, "uncertainty_metrics.json"), compact=False)

    def dump(m):
        b, t = TR.dump(m)
        return dict(baseline=b, trees=t)
    q68, q90 = out["quantiles"]["0.68"], out["quantiles"]["0.90"]
    ex = dict(meta=dict(name=opt.name, feature_version=FV, n_features=int(D.X.shape[1]), floor=0.02,
                        what="lkw の予測の外れの大きさの見積もり。s68: 外れの 68% がこの値より小さい、s90: 90% がこの値より小さい",
                        check=dict(test_coverage_68=q68.get("test_coverage"), test_coverage_90=q90.get("test_coverage"), test_spearman=q68.get("test_spearman"))),
              s68=dict(scale=models[0.68]["scale"], **dump(models[0.68]["model"])), s90=dict(scale=models[0.90]["scale"], **dump(models[0.90]["model"])))
    rng = np.random.default_rng(0)
    ix = rng.integers(0, len(D.X), 300)
    chk = dict(X=[[None if not np.isfinite(v) else float(v) for v in r] for r in D.X[ix]],
               s68=[float(v) for v in np.maximum(models[0.68]["model"].predict(D.X[ix]), 0.02) * models[0.68]["scale"]],
               s90=[float(v) for v in np.maximum(models[0.90]["model"].predict(D.X[ix]), 0.02) * models[0.90]["scale"]])
    with open(os.path.join(opt.weights_dir, "reference_uncertainty_free.json"), "w", encoding="utf-8") as f:
        json.dump(_clean(chk), f, allow_nan=False)
    return ex, models


def seen_flags(D, model, y, ex, models, stock, opt):
    """学習データに入っている化合物の幅の倍率（meta.seen_scale）と、ストックの印（seen・meta.seen_stock）"""
    ok = np.isfinite(y)
    tr = D.in_tr & ~D.test_c & ok
    err = np.abs(y - model.predict(D.X))
    sc = []
    for q in [0.68, 0.90]:
        s = np.maximum(models[q]["model"].predict(D.X[tr]), 0.02) * models[q]["scale"]
        sc.append(float(np.quantile(err[tr] / s, q)))
    ex["meta"]["seen_scale"] = dict(s68=sc[0], s90=sc[1], note="学習データに入っている化合物（compound.seen）では、幅にこの倍率を掛ける")
    seen = set(D.rows["ck"].values[tr])
    for c in stock:
        c["seen"] = {opt.name: bool((c.get("inchikey") or "")[:14] in seen)}
    ex["meta"]["seen_stock"] = sorted({(c.get("inchikey") or "")[:14] for c in stock if (c.get("inchikey") or "")[:14] in seen})
    with open(os.path.join(opt.weights_dir, "uncertainty_free.json"), "w", encoding="utf-8") as f:
        json.dump(_clean(ex), f, ensure_ascii=False, separators=(",", ":"), allow_nan=False)     # データが少ないと、確かめの数字が nan になることがある
    with open(os.path.join(opt.weights_dir, "stock_compounds.json"), "w", encoding="utf-8") as f:
        json.dump(dict(feature_version=FV, compounds=stock), f, ensure_ascii=False, separators=(",", ":"))
    C.log("学習済みの化合物の幅の倍率 %.3f・%.3f、ストック %d 種のうち学習データにある物 %d" % (sc[0], sc[1], len(stock), sum(c["seen"][opt.name] for c in stock)))


# ---------------------------------------------------------------- 名前の索引
def _norm(n):
    n = str(n).replace("\xa0", " ").strip().lower()
    return re.sub(r"\s+", " ", n)


def _good(n):
    if not (2 <= len(n) <= 70):
        return False
    if re.fullmatch(r"[\d\W_]+", n):
        return False
    if re.search(r"inchi|^cid|^hmdb|^chebi|^cas[ :]|^unk|unknown|^nan$|^na$", n):
        return False
    return True


def name_index(opt, own_names=()):
    """base_data/names.csv.gz（出所の表の名前と構造）＋自社のデータの名前 → name_index_free.json（RDKit.js で読めない行は、あとで Node が外す）"""
    from .table import read_csv
    T = read_csv(os.path.join(C.BASE, "names.csv.gz"), dtype=str, keep_default_na=False)
    EX = read_csv(os.path.join(C.BASE, "names_exclude.tsv"), sep="\t", dtype=str, keep_default_na=False)
    excl = set(EX["name"])
    cand = collections.defaultdict(collections.Counter)
    use_base = not opt.own_only
    rows = list(zip(T["name"], T["smiles"], T["source"])) if use_base else []
    rows += [(n, s, "own") for n, s in own_names]
    for name, smi, src in rows:
        n = _norm(name)
        if n in excl:
            continue
        if _good(n) and isinstance(smi, str) and 1 < len(smi) < 300:
            cand[n][(smi.strip(), src)] += 1
    srcs = ({"stock", "weber", "wiczling", "zhang"} if use_base else set()) | ({"own"} if own_names else set())
    pri = {"stock": 0, "weber": 1, "wiczling": 2, "zhang": 3, "own": 4}
    out = []
    for n in sorted(cand):
        c = collections.Counter({k: v for k, v in cand[n].items() if k[1] in srcs})
        if c:
            (smi, src), _ = sorted(c.items(), key=lambda kv: (pri[kv[0][1]], -kv[1]))[0]
            out.append([n, smi, src])
    lic = "CC BY 4.0（出所: Weber ら・Wiczling グループ〈CC BY 4.0〉、Zhang ら〈CC0〉、こちらで作った表）"
    if own_names:
        lic += "＋自社のデータの名前（出所の印 own）"
    js = dict(meta=dict(name="free" if not opt.own else opt.name, license=lic, n=len(out), sources=sorted(srcs),
                        note="名前は小文字。構造は出所の表の SMILES のまま（正規化していない）。ページでは RDKit.js で読み直して使う"), names=out)
    p = os.path.join(opt.weights_dir, "name_index_free.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(js, f, ensure_ascii=False, separators=(",", ":"))
    C.log("名前の索引: %d 件（RDKit.js で読めない行は次の段で外す）" % len(out))


# ---------------------------------------------------------------- 自社のデータの評価
def own_report(D, model, cal, y, opt):
    """自社のデータの行で、公開の重み（weights/model_free.json、補正なし）と、いま作った重み（補正あり）を比べる。
    評価用の化合物（学習に使っていない 10%）と、全部の行（学習に使った化合物を含む）。"""
    own = D.src == "own"
    if not own.any():
        return None
    shipped = json.load(open(os.path.join(C.WEIGHTS, "model_free.json"), encoding="utf-8"))
    p_base = TR.predict(shipped["baseline"], shipped["trees"], D.X[own])
    lk_base = np.full(len(y), np.nan); lk_base[own] = p_base
    lk_new = model.predict(D.X)
    res = dict(note="before = 公開の重み（weights/model_free.json）・分析法の補正なし。after = この学習の重み・この学習で決めた補正（遅れを書いた分析法は補正なし）",
               methods={})
    for mid in sorted(set(D.rows["method_id"][own])):
        ix = D.groups[mid]
        m = D.mrow.loc[mid]
        r = dict(n_rows=int(len(ix)), n_heldout_rows=int(D.test_c[ix].sum()), calibration=list(cal.get(mid, (1.0, 0.0))))
        sel = np.zeros(len(y), bool); sel[ix] = True
        if np.isfinite(m["isocratic_phi"]):
            par = c2.PAR["acn" if not (m["frac_acn"] < 0.5) else "meoh"]
            obs = D.rows["logk_iso"].values
            for tag, s in (("heldout", sel & D.test_c), ("all", sel)):
                if s.sum() >= 3:
                    eb = np.abs(c2.logk(lk_base[s], m["isocratic_phi"], par) - obs[s]); ea = np.abs(c2.logk(lk_new[s], m["isocratic_phi"], par) - obs[s])
                    r[tag] = dict(n=int(s.sum()), logk_mae_before=float(eb.mean()), logk_mae_after=float(ea.mean()))
        else:
            prog = json.loads(m["gradient"]); tg = max(t for t, _ in prog)
            rb = D.predict_rt(lk_base, {}, sel); ra = D.predict_rt(lk_new, cal, sel)
            obs = D.rows["rt_min"].values
            for tag, s in (("heldout", sel & D.test_c), ("all", sel)):
                s = s & np.isfinite(rb) & np.isfinite(ra)
                if s.sum() >= 3:
                    eb, ea = np.abs(rb[s] - obs[s]), np.abs(ra[s] - obs[s])
                    r[tag] = dict(n=int(s.sum()), run_min=tg, mae_min_before=float(eb.mean()), mae_min_after=float(ea.mean()),
                                  rel_mae_before=float(eb.mean() / tg), rel_mae_after=float(ea.mean() / tg),
                                  spearman_before=float(spearmanr(obs[s], rb[s]).statistic) if s.sum() > 3 else None,
                                  spearman_after=float(spearmanr(obs[s], ra[s]).statistic) if s.sum() > 3 else None)
        res["methods"][mid] = r
    C.dump_json(res, os.path.join(opt.outdir, "own_eval.json"), compact=False)
    return res


def write_report(opt, log, unc_metrics_path, own_res, own_rep):
    L = ["# 学習の結果（%s）" % opt.name, "",
         "- 木の数 %d、評価用に外した化合物の割合 %.2f、学習に入れた出所: %s" % (opt.trees, opt.holdout, ", ".join(log.get("sources", []))),
         "- 重み: `weights/`（models/hplc-retention/ で `node tests/check.cjs train/out/%s/weights` を実行すると、JavaScript との一致を確かめられる）" % opt.name, ""]
    hc = log.get("heldout_compounds")
    if hc:
        L += ["## 学習に使っていない化合物（勾配の分析法の保持時間）", "", "| 出所 | 分析法 | 測定 | 保持時間の誤差（分析時間に対する割合、分析法ごとの中央値） | 順位の相関 |", "|---|---|---|---|---|"]
        for s, v in sorted(log.get("heldout_compounds_by_source", {}).items()):
            L.append("| %s | %d | %d | %.1f%% | %.2f |" % (s, v["n_methods"], v["n_meas"], 100 * v["median_rel_mae"], v["median_spearman"]))
        L.append("| 全体 | %d | %d | %.1f%% | %.2f |" % (hc["n_methods"], hc["n_meas"], 100 * hc["median_rel_mae"], hc["median_spearman"]))
        L.append("")
    iw = log.get("isocratic_weber_cv_by_compound")
    if iw:
        L += ["- 等組成（Weber、化合物の 5 分割）: log k の誤差 %.3f、2 倍以内 %.0f%%" % (iw["logk_mae"], 100 * iw["within_factor_2"]), ""]
    if os.path.exists(unc_metrics_path):
        u = json.load(open(unc_metrics_path, encoding="utf-8"))["quantiles"]
        if "test_coverage" in u["0.68"]:
            L += ["- 予測の幅に入る割合（評価用の化合物）: 68%% の幅に %.1f%%、90%% の幅に %.1f%%" % (100 * u["0.68"]["test_coverage"], 100 * u["0.90"]["test_coverage"]), ""]
    if own_rep:
        L += ["## 自社のデータ", "", "- 分析法 %d、測定 %d、化合物 %d（基礎データにも入っている化合物 %d）" % (own_rep["n_methods"], own_rep["n_rows"], own_rep["n_compounds"], own_rep["n_known_compounds"])]
        if own_rep["unreadable_smiles"]:
            L.append("- 読めなかった構造 %d 個（外した）: %s" % (len(own_rep["unreadable_smiles"]), ", ".join(own_rep["unreadable_smiles"][:10])))
        if own_rep["skipped_rows"]:
            L.append("- 外した行 %d（理由は `own_eval.json` の隣の `own_data_report.json`）" % len(own_rep["skipped_rows"]))
        L.append("")
    if own_res:
        L += ["| 分析法 | 行 | 評価用の化合物の行 | 前（公開の重み・補正なし） | 後（この重み） |", "|---|---|---|---|---|"]
        for mid, r in own_res["methods"].items():
            h = r.get("heldout") or {}
            if "rel_mae_before" in h:
                L.append("| %s | %d | %d | 保持時間の誤差 %.2f 分（%.1f%%） | %.2f 分（%.1f%%） |" % (mid, r["n_rows"], h["n"], h["mae_min_before"], 100 * h["rel_mae_before"], h["mae_min_after"], 100 * h["rel_mae_after"]))
            elif "logk_mae_before" in h:
                L.append("| %s | %d | %d | log k の誤差 %.3f | %.3f |" % (mid, r["n_rows"], h["n"], h["logk_mae_before"], h["logk_mae_after"]))
            else:
                L.append("| %s | %d | %d | （評価用の化合物が 3 行未満） | |" % (mid, r["n_rows"], r["n_heldout_rows"]))
        L += ["", "「評価用の化合物」は、どの学習にも使っていない化合物（InChIKey のハッシュで %.0f%%）。全部の行の値は `own_eval.json`。" % (100 * opt.holdout), ""]
    open(os.path.join(opt.outdir, "REPORT.md"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
