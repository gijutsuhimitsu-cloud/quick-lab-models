# SPDX-License-Identifier: MIT
"""GPU で学習するグラフのモデルのための、学習の表の書き出し（手元の CPU。PyTorch・RDKit は要らない）。
  python train/gpu/export_graph_data.py [名前（既定 free）] [出力のフォルダ（既定 train/work/<名前>/gpu_data）]
前提: python train/run.py を同じ名前で実行してある（分析法ごとの補正・予測の幅の 4 分割・木の重みを使う）。
出力（train_teacher.py がそのまま読む形）:
  rows.npz      XC（記述子 41＋官能基 56）・IC（イオン化の規則 30）・M（分析法 16）・ci・mi・IF（イオン化の量 10）・w・test_c・is_grad・src、
                in_<名前>（学習に入る行）・fold_<名前>・foldc_<名前>（化合物のかたまりの 4 分割。予測の幅と同じ）、
                y_raw_<名前>（実測から逆算した lkw。上限・下限で切る前）、tree_<名前>（木の重みの予測。比べる基準）
  compounds.tsv ck と SMILES（graphs_js.cjs がグラフを作る）
  meta.json     出所・分析法の名前・pH など。cal_<名前>.json は分析法ごとの補正
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from lib import common as C, ionfeat, trees as TR
from lib.runs import data_for_run
from lib.data import FV, Y_FLOOR, Y_CEIL

name = sys.argv[1] if len(sys.argv) > 1 else "free"
OUT = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.join(C.WORK, name, "gpu_data")
os.makedirs(OUT, exist_ok=True)
D, cfg, cal = data_for_run(name)
foldc = np.load(os.path.join(C.WORK, name, "fold_c.npy"), allow_pickle=False).astype(np.int8)
assert len(foldc) == len(D.cks)
srcs = sorted(set(D.src))
yr = D.targets(cal, clip=False)
model = json.load(open(os.path.join(C.OUT, name, "weights", "model_free.json"), encoding="utf-8"))
tree = TR.predict(model["baseline"], model["trees"], D.X)
n_ion = len(ionfeat.keep(FV))
np.savez_compressed(os.path.join(OUT, "rows.npz"),
                    XC=D.XC.astype(np.float32), IC=D.IC.astype(np.float32), M=D.M.astype(np.float32),
                    ci=D.ci.astype(np.int32), mi=D.mi.astype(np.int32), IF=D.X[:, -n_ion:].astype(np.float32),
                    w=D.w.astype(np.float32), test_c=D.test_c, is_grad=D.is_grad, src=np.array([srcs.index(s) for s in D.src], np.int8),
                    **{"in_" + name: D.in_tr, "fold_" + name: foldc[D.ci], "foldc_" + name: foldc, "y_raw_" + name: yr.astype(np.float32), "tree_" + name: tree.astype(np.float32)})
with open(os.path.join(OUT, "compounds.tsv"), "w", encoding="utf-8", newline="\n") as f:
    for c in D.cks:
        f.write("%s\t%s\n" % (c, D.comp[c]["smiles"]))
json.dump({k: list(v) for k, v in cal.items()}, open(os.path.join(OUT, "cal_%s.json" % name), "w"))
json.dump(dict(sources=srcs, mnames=D.mnames, ion_keep=ionfeat.keep(FV), ion_signs=[int(x) for x in ionfeat.signs(FV)], y_floor=Y_FLOOR, y_ceil=Y_CEIL,
               fv=FV, n_rows=int(len(D.rows)), n_comp=int(len(D.cks)), n_methods=int(len(D.md)), method_ids=list(D.md["method_id"]),
               method_source=list(D.md["source"]), method_pH=[None if not np.isfinite(v) else float(v) for v in D.md["pH_weak"].astype(float)]),
          open(os.path.join(OUT, "meta.json"), "w", encoding="utf-8"), ensure_ascii=False)
C.log("書き出した: 行 %d、化合物 %d、分析法 %d → %s" % (len(D.rows), len(D.cks), len(D.md), OUT))
