# SPDX-License-Identifier: MIT
"""graphs_js.cjs の出力（ブラウザと同じコードで作ったグラフ）を、学習スクリプトの形（graphs_js.npz）にまとめる。
  python gpu/graphs_js_to_npz.py <データのフォルダ>
結合は向きつきに直す（1 本を両方向に）。"""
import os, sys, json
import numpy as np
P = sys.argv[1]
m = json.load(open(os.path.join(P, "gjs_meta.json")))
rd = lambda n, dt: np.fromfile(os.path.join(P, n), dtype=dt)
atom = rd("gjs_atom.bin", np.uint8).reshape(-1, 8); bond = rd("gjs_bond.bin", np.uint8).reshape(-1, 3)
a, b = rd("gjs_a.bin", np.int16), rd("gjs_b.bin", np.int16)
aptr, bptr = rd("gjs_aptr.bin", np.int32).astype(np.int64), rd("gjs_bptr.bin", np.int32).astype(np.int64)
assert len(atom) == m["n_atoms"] and len(bond) == m["n_bonds"] and aptr[-1] == len(atom) and bptr[-1] == len(bond)
for k, d in enumerate(m["atom_dims"]):
    assert atom[:, k].max() < d, (k, atom[:, k].max(), d)
for k, d in enumerate(m["bond_dims"]):
    assert bond[:, k].max() < d, (k, bond[:, k].max(), d)
esrc = np.column_stack([a, b]).ravel(); edst = np.column_stack([b, a]).ravel()
np.savez_compressed(os.path.join(P, "graphs_js.npz"), atom=atom, bond=np.repeat(bond, 2, axis=0), esrc=esrc.astype(np.int16), edst=edst.astype(np.int16),
                    aptr=aptr, bptr=bptr * 2, adims=np.array(m["atom_dims"]), bdims=np.array(m["bond_dims"]), version=np.array(m["version"]))
print("molecules", m["n_mol"], "unparsed", m["unparsed"], "atoms", len(atom), "directed bonds", 2 * len(bond))
for k, d in enumerate(m["atom_dims"]):
    print("  atom feature", k, np.bincount(atom[:, k], minlength=d).tolist())
for k, d in enumerate(m["bond_dims"]):
    print("  bond feature", k, np.bincount(bond[:, k], minlength=d).tolist())
