# SPDX-License-Identifier: MIT
"""保持のモデル（第 2 版の式）。js/retention_model.js と同じ式。

  u(phi)   = phi / (1 + B*phi)
  S        = a + b * lkw                      （傾きは log kw と連動させる）
  log10 k  = lkw - S * u(phi) = lkw * (1 - b*u) - a*u

a, b, B は有機溶媒ごとの定数（PAR）。化合物ごとの未知数は lkw（見かけの log kw）だけ。
勾配: (1/t0) ∫_0^{tR - t0} dt / k(phi(t)) = 1 を数値積分で解く。
"""
import numpy as np

# Weber らの等組成の実測（C18、アセトニトリル／10 mM リン酸、25〜75 °C、105 化合物）への当てはめ。
# 値は base_data/retention_form.json を読む（scripts/fit_form.py で当てはめ直せる。なければ下の初期値）。
PAR = {"acn": dict(a=4.231, b=1.808, B=1.496), "meoh": dict(a=4.231, b=1.808, B=1.496)}
LN10 = np.log(10.0)


def load_par(path):
    import json, os
    if os.path.exists(path):
        j = json.load(open(path, encoding="utf-8"))
        for k in ("acn", "meoh"):
            if k in j:
                PAR[k] = {x: float(j[k][x]) for x in ("a", "b", "B")}
    return PAR


def u_of(phi, B):
    return phi / (1.0 + B * phi)


def logk(lkw, phi, par):
    u = u_of(np.asarray(phi, float), par["B"])
    return lkw * (1.0 - par["b"] * u) - par["a"] * u


def lkw_from_isocratic(logk_obs, phi, par):
    u = u_of(np.asarray(phi, float), par["B"])
    return (np.asarray(logk_obs, float) + par["a"] * u) / (1.0 - par["b"] * u)


def _grid(prog, t_dwell, t_max, n):
    pts = sorted((float(t) + t_dwell, float(p)) for t, p in prog)
    t = np.linspace(0.0, max(t_max, 1e-6), n)
    phi = np.interp(t, [p[0] for p in pts], [p[1] for p in pts])
    return t, phi, pts[-1][1]


def lkw_from_gradient(tr, t0, prog, par, t_dwell=0.0, n=1200, chunk=4000, lo=-3.0, hi=14.0, iters=44):
    """保持時間（分、配列）→ lkw。tr <= t0 は nan。二分法（F は lkw について単調減少）。"""
    tr = np.asarray(tr, float)
    out = np.full(tr.shape, np.nan)
    ok = np.where(tr > t0)[0]
    if not len(ok):
        return out
    tau_all = tr[ok] - t0
    t, phi, _ = _grid(prog, t_dwell, float(tau_all.max()), n)
    u = u_of(phi, par["B"])
    dt = t[1] - t[0]
    c1 = LN10 * par["a"] * u                 # exp(c1 - LN10*lkw*(1-b u))
    c2 = LN10 * (1.0 - par["b"] * u)
    for s in range(0, len(ok), chunk):
        ix = ok[s:s + chunk]
        tau = tau_all[s:s + chunk]
        pos = tau / dt
        j = np.minimum(pos.astype(int), n - 2)
        fr = pos - j
        a_lo = np.full(len(ix), lo)
        a_hi = np.full(len(ix), hi)
        for _ in range(iters):
            mid = 0.5 * (a_lo + a_hi)
            f = np.exp(c1[None, :] - mid[:, None] * c2[None, :])          # 1/k
            cum = np.cumsum((f[:, 1:] + f[:, :-1]) * 0.5 * dt, axis=1)
            cum = np.hstack([np.zeros((len(ix), 1)), cum])
            r = np.arange(len(ix))
            val = cum[r, j] + fr * (cum[r, j + 1] - cum[r, j])            # ∫_0^tau dt/k
            big = val > t0                                                # まだ速すぎる → lkw を上げる
            a_lo = np.where(big, mid, a_lo)
            a_hi = np.where(big, a_hi, mid)
        out[ix] = 0.5 * (a_lo + a_hi)
    return out


def rt_from_lkw(lkw, t0, prog, par, t_dwell=0.0, n=1200, chunk=4000, t_max=None):
    """lkw（配列）→ 保持時間（分）。"""
    lkw = np.asarray(lkw, float)
    out = np.full(lkw.shape, np.nan)
    tm = t_max if t_max else max(float(max(p[0] for p in prog)) + t_dwell + 5.0, 10.0)
    t, phi, phi_end = _grid(prog, t_dwell, tm, n)
    u = u_of(phi, par["B"])
    dt = t[1] - t[0]
    c1 = LN10 * par["a"] * u
    c2 = LN10 * (1.0 - par["b"] * u)
    ue = u_of(phi_end, par["B"])
    for s in range(0, len(lkw), chunk):
        x = lkw[s:s + chunk]
        f = np.exp(c1[None, :] - x[:, None] * c2[None, :]) / t0
        cum = np.cumsum((f[:, 1:] + f[:, :-1]) * 0.5 * dt, axis=1)
        cum = np.hstack([np.zeros((len(x), 1)), cum])
        reached = cum[:, -1] >= 1.0
        j = np.minimum((cum < 1.0).sum(axis=1) - 1, n - 2)
        j = np.maximum(j, 0)
        r = np.arange(len(x))
        den = np.maximum(cum[r, j + 1] - cum[r, j], 1e-300)
        tau = t[j] + (1.0 - cum[r, j]) / den * dt
        # 表の終わりまでに出ない物は、最後の割合のまま続くとして延長
        rate_end = np.exp(LN10 * par["a"] * ue - x * LN10 * (1.0 - par["b"] * ue)) / t0
        tau = np.where(reached, tau, t[-1] + (1.0 - cum[:, -1]) / np.maximum(rate_end, 1e-300))
        out[s:s + chunk] = t0 + tau
    return out


if __name__ == "__main__":
    par = PAR["acn"]
    prog = [(0, 0.05), (3, 0.05), (5, 0.5), (20, 0.85), (23, 0.85)]
    lk = np.array([0.3, 1.5, 3.0, 4.5, 6.0])
    tr = rt_from_lkw(lk, 1.1, prog, par, 0.5, n=4000)
    back = lkw_from_gradient(tr, 1.1, prog, par, 0.5, n=4000)
    print(tr, back)
    assert np.allclose(back, lk, atol=2e-3), back - lk
    iso = rt_from_lkw(lk, 1.0, [(0, 0.4)], par, 0.0, t_max=4000, n=200000)
    k = 10 ** logk(lk, 0.4, par)
    print(iso, 1.0 * (1 + k))
    assert np.allclose(iso[:4], (1 + k)[:4], rtol=2e-3)
    assert np.allclose(lkw_from_isocratic(logk(lk, 0.3, par), 0.3, par), lk)
    print("ok")
