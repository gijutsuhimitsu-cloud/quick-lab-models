# SPDX-License-Identifier: MIT
"""保持の式の定数（a, b, B）を、基礎データから当てはめ直す（任意。1 分以内）。同梱の base_data/retention_form.json と比べるだけで、書き換えない。
  python train/scripts/fit_form.py
  log10 k = lkw − (a + b·lkw)·φ/(1 + B·φ)
  - アセトニトリル: Weber らの等組成の実測（105 化合物 × 割合 × 温度）。化合物・温度ごとに lkw だけが自由で、残差の 2 乗和が最小の定数
  - メタノール: Wiczling グループの、同じ化合物・カラム・pH・温度を勾配時間 30 分と 90・270 分で測った組。30 分の保持時間から lkw を逆算し、
    90・270 分の保持時間を予測した相対誤差の中央値が最小の定数（当てはめは 2 カラム × pH 3 水準、検査は全部）
出力: train/out/form/retention_form.json（使うには base_data/retention_form.json と入れ替えて学習し直す。重みの meta.par も変わる）
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from lib import common as C, table, chromlib2 as c2

md, rows, _ = table.load_base()
ship = json.load(open(os.path.join(C.BASE, "retention_form.json"), encoding="utf-8"))
# ---- アセトニトリル（Weber）
w = rows.merge(md[["method_id", "source", "isocratic_phi", "temp_C"]], on="method_id")
w = w[w["source"] == "weber_isocratic_mendeley"]
groups = [(s["isocratic_phi"].values, s["logk_iso"].values) for _, s in w.groupby(["ck", "temp_C"])]


def sse(p, detail=False):
    a, b, B = p
    tot, n, errs = 0.0, 0, []
    for phi, y in groups:
        u = phi / (1 + B * phi)
        A = 1 - b * u
        lkw = np.sum(A * (y + a * u)) / np.sum(A * A)
        r = lkw * A - a * u - y
        tot += np.sum(r * r); n += len(y); errs.extend(r)
    return (np.sqrt(tot / n), np.array(errs)) if detail else tot


r = minimize(sse, [4.2, 1.8, 1.5], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-6, maxiter=4000))
rmse, errs = sse(r.x, True)
acn = dict(a=float(r.x[0]), b=float(r.x[1]), B=float(r.x[2]))
C.log("アセトニトリル: a %.4f, b %.4f, B %.4f（同梱 %.4f, %.4f, %.4f）、log k の RMSE %.4f（%d 点）" % (
    acn["a"], acn["b"], acn["B"], ship["acn"]["a"], ship["acn"]["b"], ship["acn"]["B"], rmse, len(w)))
# ---- メタノール（Wiczling）
wz = md[md["source"] == "wiczling_osf"].copy()
wz["prog"] = wz["gradient"].map(json.loads)
wz["tg"] = wz["prog"].map(lambda p: round(p[1][0] - p[0][0]))
rr = rows[rows["method_id"].isin(set(wz["method_id"]))]
by = {m: g.set_index("ck")["rt_min"] for m, g in rr.groupby("method_id")}
COLS_FIT = ["hsm2_303", "hsm2_484"]          # 当てはめに使う 2 カラム（columns.json の id）
PHS_FIT = [2.5, 5.8, 9.7]


def pairs(mod_acn, cols, phs):
    out = []
    sub = wz[(wz["frac_acn"] == mod_acn) & wz["column_id"].isin(cols) & wz["pH_weak"].isin(phs)]
    for _, g in sub.groupby(["column_id", "pH_weak", "temp_C"]):
        a = g[g["tg"] == 30]
        if not len(a):
            continue
        a = a.iloc[0]
        for _, b in g[g["tg"] > 30].iterrows():
            s1, s2 = by.get(a["method_id"]), by.get(b["method_id"])
            if s1 is None or s2 is None:
                continue
            s1 = s1[~s1.index.duplicated()]; s2 = s2[~s2.index.duplicated()]
            ck = s1.index.intersection(s2.index)
            if len(ck) < 20:
                continue
            out.append((a["prog"], a["t0_min"], s1.loc[ck].values, b["prog"], b["t0_min"], s2.loc[ck].values))
    return out


def err(par, P, n=400):
    e = []
    for pa, t0a, rta, pb, t0b, rtb in P:
        lk = c2.lkw_from_gradient(rta, t0a, pa, par, 0.0, n=n, iters=32)
        pr = c2.rt_from_lkw(lk, t0b, pb, par, 0.0, n=n)
        ok = np.isfinite(pr) & (rta > 1.3 * t0a)
        e.extend(np.abs(pr[ok] - rtb[ok]) / rtb[ok])
    return float(np.median(e)), len(e)


Pm = pairs(0.0, COLS_FIT, PHS_FIT)
C.log("メタノール: 当てはめの組 %d" % len(Pm))


def obj(x):
    a, b, B = x
    if B < 0 or b <= 0 or b * (1.0 / (1.0 + B)) > 0.95:       # 1 − b·u(1) > 0 を保つ
        return 9.9
    return err(dict(a=a, b=b, B=B), Pm)[0]


r2 = minimize(obj, [acn["a"], acn["b"], acn["B"]], method="Nelder-Mead", options=dict(xatol=0.02, fatol=2e-4, maxiter=90))
meoh = dict(a=float(r2.x[0]), b=float(r2.x[1]), B=float(r2.x[2]))
ALL = sorted(wz["column_id"].unique()); ALLPH = sorted(wz["pH_weak"].unique())
chk = dict(meoh_fitted=err(meoh, pairs(0.0, ALL, ALLPH)), acn=err(acn, pairs(1.0, ALL, ALLPH)))
C.log("メタノール: a %.4f, b %.4f, B %.4f（同梱 %.4f, %.4f, %.4f）、勾配時間を変えた予測の相対誤差の中央値 %.4f（%d 点）" % (
    meoh["a"], meoh["b"], meoh["B"], ship["meoh"]["a"], ship["meoh"]["b"], ship["meoh"]["B"], chk["meoh_fitted"][0], chk["meoh_fitted"][1]))
out = dict(acn=dict(acn, note="Weber らの等組成の実測に当てはめた値"), meoh=dict(meoh, note="Wiczling グループのデータ（メタノール、勾配時間 30 分 → 90・270 分の予測）で決めた"),
           fit=dict(n_points=int(len(w)), rmse_logk=float(rmse)), gradient_time_check={k: dict(median_rel_error=v[0], n=v[1]) for k, v in chk.items()},
           shipped=dict(acn={k: ship["acn"][k] for k in "abB"}, meoh={k: ship["meoh"][k] for k in "abB"}))
C.dump_json(out, os.path.join(C.OUT, "form", "retention_form.json"), compact=False)
print(json.dumps({k: out[k] for k in ("acn", "meoh", "shipped")}, ensure_ascii=False, indent=1))
