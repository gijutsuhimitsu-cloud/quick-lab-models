// SPDX-License-Identifier: MIT
// 学習データの全化合物のグラフを、ブラウザと同じコード（js/graph_features.js ＋ 同梱の RDKit.js）で作る（手元の Node。8 万分子で数分）。
//   node train/gpu/graphs_js.cjs <データのフォルダ（compounds.tsv がある所）>
// 出力: <フォルダ>/gjs_atom.bin（uint8, 原子 × 8）・gjs_bond.bin（uint8, 結合 × 3）・gjs_a.bin・gjs_b.bin（int16, 分子の中の原子の番号）・
//       gjs_aptr.bin・gjs_bptr.bin（int32, 分子ごとの始まり）・gjs_meta.json。あとで graphs_js_to_npz.py が graphs_js.npz にまとめる。
const fs = require('fs'), path = require('path');
const MODEL = path.join(__dirname, '..', '..');
const G = require(path.join(MODEL, 'js', 'graph_features.js'));
const initRDKitModule = require(path.join(MODEL, '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'));
if (!process.argv[2]) { console.error('使い方: node graphs_js.cjs <データのフォルダ>'); process.exit(2); }
const dir = path.resolve(process.argv[2]);
(async () => {
  const R = await initRDKitModule();
  const lines = fs.readFileSync(path.join(dir, 'compounds.tsv'), 'utf8').split('\n').filter(l => l.length);
  const atom = [], bond = [], ea = [], eb = [], aptr = [0], bptr = [0];
  let bad = 0, na = 0, nbd = 0, maxAtoms = 0;
  const t0 = Date.now();
  for (let i = 0; i < lines.length; i++) {
    const smi = lines[i].split('\t')[1];
    const g = G.fromSmiles(R, smi);
    if (g) {
      atom.push(g.atom); bond.push(g.bond); ea.push(Int16Array.from(g.a)); eb.push(Int16Array.from(g.b));
      na += g.n; nbd += g.nb; if (g.n > maxAtoms) maxAtoms = g.n;
    } else bad++;
    aptr.push(na); bptr.push(nbd);
    if (i % 20000 === 0) console.log(i, 'of', lines.length, ((Date.now() - t0) / 1000).toFixed(0) + 's');
  }
  const cat = (arrs, T, n) => { const o = new T(n); let p = 0; for (const a of arrs) { o.set(a, p); p += a.length; } return Buffer.from(o.buffer, o.byteOffset, o.byteLength); };
  fs.writeFileSync(path.join(dir, 'gjs_atom.bin'), cat(atom, Uint8Array, na * 8));
  fs.writeFileSync(path.join(dir, 'gjs_bond.bin'), cat(bond, Uint8Array, nbd * 3));
  fs.writeFileSync(path.join(dir, 'gjs_a.bin'), cat(ea, Int16Array, nbd));
  fs.writeFileSync(path.join(dir, 'gjs_b.bin'), cat(eb, Int16Array, nbd));
  fs.writeFileSync(path.join(dir, 'gjs_aptr.bin'), Buffer.from(Int32Array.from(aptr).buffer));
  fs.writeFileSync(path.join(dir, 'gjs_bptr.bin'), Buffer.from(Int32Array.from(bptr).buffer));
  fs.writeFileSync(path.join(dir, 'gjs_meta.json'), JSON.stringify({ version: G.VERSION, atom_dims: G.ATOM_DIMS, bond_dims: G.BOND_DIMS, n_mol: lines.length, n_atoms: na, n_bonds: nbd, unparsed: bad, max_atoms: maxAtoms }));
  console.log('molecules', lines.length, 'unparsed', bad, 'atoms', na, 'bonds', nbd, 'max atoms', maxAtoms, ((Date.now() - t0) / 1000).toFixed(0) + 's');
})();
