# SPDX-License-Identifier: MIT
"""イオン化の特徴量（特徴量の第 2 版から。いまの重みは第 4 版）。js/retention_model.js の ionFeatures と同じ計算。
規則（名前・種類・代表的な pKa）は js/features.js の ION の表から読む（写しを持たない）。
  入力: 規則ごとの部位の数（化合物ごと）、pH（測定ごと）
  出力: [負の電荷, 正の電荷, 正味, 合計, log10(中性の割合), 酸の数, 塩基の数, いちばん強い酸の pKa, いちばん強い塩基の pKa,
         pKa(酸) − pH, pH − pKa(塩基), 正 × 負]
"""
import os, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ["ion_qneg", "ion_qpos", "ion_qnet", "ion_qabs", "ion_logneutral", "ion_n_acid", "ion_n_base", "ion_pka_acid_min", "ion_pka_base_max",
         "ion_pka_a_minus_pH", "ion_pH_minus_pka_b", "ion_zwitter"]
AMINES = ("amine_primary", "amine_secondary", "amine_tertiary")
MAX_GROUPS, MAX_REP, BETA_SHIFT, PH_DEFAULT = 8, 3, -1.5, 3.0
_rules = None


def rules():
    """[(名前, 種類, pKa)]。種類は acid / base / mod / perm+"""
    global _rules
    if _rules is None:
        s = open(os.path.join(HERE, "..", "..", "js", "features.js"), encoding="utf-8").read()
        s = s[s.index("var ION = ["):s.index("];", s.index("var ION = ["))]
        _rules = [(m.group(1), m.group(2), float(m.group(3)))
                  for m in re.finditer(r"\['([a-z0-9_+]+)',\s*.*?,\s*'(acid|base|mod|perm\+)',\s*(-?[\d.]+)\]", s)]
    return _rules


def groups(counts):
    """counts: (n, 規則の数) → pKa の表 (n, 8)、種類 (n, 8。酸 −1・塩基 +1)、永久の正電荷の数 (n)"""
    R = rules(); names = [r[0] for r in R]
    n = len(counts)
    pka = np.full((n, MAX_GROUPS), np.nan); typ = np.zeros((n, MAX_GROUPS)); perm = np.zeros(n)
    beta = counts[:, names.index("amine_beta_hetero")] > 0
    fill = np.zeros(n, int)
    for j, (nm, t, pk) in enumerate(R):
        c = counts[:, j].astype(int)
        if t == "perm+":
            perm += c
            continue
        if t not in ("acid", "base"):
            continue
        for rep in range(1, MAX_REP + 1):
            ix = np.where((c >= rep) & (fill < MAX_GROUPS))[0]
            v = np.full(len(ix), pk)
            if nm in AMINES:
                v = v + BETA_SHIFT * beta[ix]
            pka[ix, fill[ix]] = v; typ[ix, fill[ix]] = -1 if t == "acid" else 1
            fill[ix] += 1
    return pka, typ, perm


# 第 4 版: 正味の電荷（2）と両性の積（11）を外し、残りに向きの制約を付ける（イオン化すると逆相の保持は下がる）
KEEP4 = [0, 1, 3, 4, 5, 6, 7, 8, 9, 10]
SIGN4 = {0: -1, 1: -1, 3: -1, 4: 1, 9: 1, 10: 1}      # −1: 増えると lkw は上がらない、+1: 増えると lkw は下がらない


def keep(fv):
    """その版で使うイオン化の特徴量の番号"""
    return list(KEEP4) if fv >= 4 else list(range(len(NAMES)))


def signs(fv):
    """keep(fv) の並びに合わせた向きの制約（0 は制約なし）"""
    return [SIGN4.get(k, 0) if fv >= 4 else 0 for k in keep(fv)]


CAP_NO_GROUP = 10.0      # 第 3 版: 酸（塩基）の基が無い時の「pKa − pH」（「pH − pKa」）は、pH によらない定数にする


def features(pka, typ, perm, pH, cap=None):
    """行ごと（pka・typ・perm は行に合わせて並べた物）。cap を与えると、基の無い側の差を定数にする。"""
    pH = np.where(np.isfinite(pH), pH, PH_DEFAULT)[:, None]
    with np.errstate(over="ignore", invalid="ignore"):
        fa = np.nan_to_num(np.where(typ < 0, 1.0 / (1.0 + 10.0 ** (pka - pH)), 0.0))
        fb = np.nan_to_num(np.where(typ > 0, 1.0 / (1.0 + 10.0 ** (pH - pka)), 0.0))
        qneg = fa.sum(1); qpos = fb.sum(1) + perm
        logneu = np.log10(np.maximum(1 - fa, 1e-6)).sum(1) + np.log10(np.maximum(1 - fb, 1e-6)).sum(1)
        logneu = np.where(perm > 0, -6.0, logneu)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)          # 酸（塩基）の基が無い行は全部 nan（下で 14・0 に置き換える）
            pa = np.nanmin(np.where(typ < 0, pka, np.nan), axis=1); pb = np.nanmax(np.where(typ > 0, pka, np.nan), axis=1)
    pa = np.where(np.isfinite(pa), pa, 14.0); pb = np.where(np.isfinite(pb), pb, 0.0)
    na = (typ < 0).sum(1); nb = (typ > 0).sum(1)
    da = pa - pH[:, 0]; db = pH[:, 0] - pb
    if cap is not None:
        da = np.where(na == 0, cap, da); db = np.where(nb == 0, cap, db)
    return np.column_stack([qneg, qpos, qpos - qneg, qpos + qneg, logneu, na, nb, pa, pb, da, db, qpos * qneg]).astype(np.float32)
