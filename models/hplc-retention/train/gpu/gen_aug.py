# SPDX-License-Identifier: MIT
"""蒸留のための「分子 × 条件」の組を作る（手元の CPU。教師のグラフのモデルに予測させて、小さいモデルの学習に足す行）。
  python train/gpu/gen_aug.py <データのフォルダ> [1 化合物あたりの数（既定 12）] [名前（既定 free）] [all]
  - 化合物: その版の学習に入っている化合物（all を付けると評価用の 10% も入れる。本番用の全部の行での学習はこちら）
  - 条件: 3 分の 2 は学習データにある分析法から（出所ごとに重みを付けて選ぶ。SMRT の 1 条件は実測で足りているので入れない）、
          3 分の 1 は HSM2 の 551 カラム × ランダムな pH（2〜10.5）・溶媒（アセトニトリル 6 割）・温度（25〜50 °C）
  → 実測は増えない。SMRT にしか無い分子を、ほかのカラム・pH で教師がどう予測するかを、小さいモデルに写すための行。
出力: <データのフォルダ>/aug.npz（ci, mi, IF, M_ext＝分析法の表の後ろに合成の条件を足した物, pH）
"""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from lib import ionfeat

P = os.path.abspath(sys.argv[1])
K = int(sys.argv[2]) if len(sys.argv) > 2 else 12
GROUP = sys.argv[3] if len(sys.argv) > 3 else "free"
ALL = len(sys.argv) > 4 and sys.argv[4] == "all"
R = np.load(os.path.join(P, "rows.npz"))
meta = json.load(open(os.path.join(P, "meta.json"), encoding="utf-8"))
rng = np.random.default_rng(0)
src_m = np.array(meta["method_source"])
ci, mi = R["ci"], R["mi"]
use = R["in_" + GROUP] & (True if ALL else ~R["test_c"])
comps = np.unique(ci[use])
mn = meta["mnames"]
M = R["M"]
pH_m = np.array([np.nan if v is None else v for v in meta["method_pH"]], float)
used_m = np.unique(mi[use])
wsrc = {"wiczling_osf": 0.40, "hsm2_mendeley": 0.20, "zhang2024_scidb": 0.15, "mcart_msd_zenodo": 0.10, "weber_isocratic_mendeley": 0.15, "own": 0.30}
pm = np.array([wsrc.get(src_m[m], 0.0) / max((src_m[used_m] == src_m[m]).sum(), 1) for m in used_m])
if pm.sum() <= 0:
    pm = np.ones(len(used_m))
pm /= pm.sum()
k_real = (2 * K) // 3; k_syn = K - k_real
hs = np.where(src_m == "hsm2_mendeley")[0]
if not len(hs):
    raise SystemExit("合成の条件に使う HSM2 のカラムが無い（自社のデータだけの版では、蒸留の組を作らない）")
n_syn = 4000
base = rng.choice(hs, n_syn)
M_syn = M[base].copy()
M_syn[:, mn.index("frac_acn")] = (rng.random(n_syn) < 0.6).astype(np.float32)
M_syn[:, mn.index("temp_C")] = rng.uniform(25, 50, n_syn)
pH_syn = rng.uniform(2.0, 10.5, n_syn)
M_ext = np.vstack([M, M_syn]).astype(np.float32)
pH_ext = np.concatenate([pH_m, pH_syn])
a_ci = np.repeat(comps, K)
a_mi = np.concatenate([rng.choice(used_m, (len(comps), k_real), p=pm), len(M) + rng.integers(0, n_syn, (len(comps), k_syn))], axis=1).ravel()
pka, typ, perm = ionfeat.groups(R["IC"])
with np.errstate(all="ignore"):
    IF = ionfeat.features(pka[a_ci], typ[a_ci], perm[a_ci], pH_ext[a_mi], ionfeat.CAP_NO_GROUP)[:, meta["ion_keep"]].astype(np.float32)
    chk = rng.choice(np.where(use)[0], min(2000, int(use.sum())), replace=False)        # 元の行で同じ計算をして、書き出した IF と一致するか
    IF0 = ionfeat.features(pka[ci[chk]], typ[ci[chk]], perm[ci[chk]], pH_m[mi[chk]], ionfeat.CAP_NO_GROUP)[:, meta["ion_keep"]]
d = np.nanmax(np.abs(np.nan_to_num(IF0) - np.nan_to_num(R["IF"][chk])))
assert d < 1e-4, d
np.savez_compressed(os.path.join(P, "aug.npz"), ci=a_ci.astype(np.int32), mi=a_mi.astype(np.int32), IF=IF, M_ext=M_ext, pH=pH_ext[a_mi].astype(np.float32))
print("化合物 %d、組 %d（学習データの分析法 %d・合成の条件 %d）、照合 %.1e" % (len(comps), len(a_ci), len(comps) * k_real, len(comps) * k_syn, d))
print(pd.Series(src_m[a_mi[a_mi < len(M)]]).value_counts().to_dict())
