// SPDX-License-Identifier: MIT
// グラフのモデルの JavaScript の推論（js/gnn_model.js ＋ graph_features.js ＋ features.js）を、PyTorch の予測と照合する（この写しでは未実行）。
//   node train/gpu/check_gnn.cjs <train_teacher.py --export 1 の出力のフォルダ（manifest.json・weights.bin・check.json）>
//   1. 記述子は Python 側の値を使い、グラフは SMILES から JavaScript で作る → lkw と z の差
//   2. 記述子も JavaScript（features.js）で作る → 端から端までの差
//   3. 速さ（分子ごとの encode、条件ごとの lkw）
const fs = require('fs'), path = require('path');
const MODEL = path.join(__dirname, '..', '..'), EX = path.join(MODEL, 'js');
const GNN = require(path.join(EX, 'gnn_model.js')), CG = require(path.join(EX, 'graph_features.js')), F = require(path.join(EX, 'features.js'));
const initRDKitModule = require(path.join(MODEL, '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'));
const dir = path.resolve(process.argv[2]);
const now = () => Number(process.hrtime.bigint()) / 1e6;
(async () => {
  const R = await initRDKitModule();
  const man = JSON.parse(fs.readFileSync(path.join(dir, 'manifest.json'), 'utf8'));
  const buf = fs.readFileSync(path.join(dir, 'weights.bin'));
  const W = new Float32Array(buf.buffer, buf.byteOffset, buf.byteLength / 4);
  const chk = JSON.parse(fs.readFileSync(path.join(dir, 'check.json'), 'utf8'));
  const M = GNN.create(man, W);
  let d1 = 0, dz = 0, d2 = 0, dd = 0, n = 0, bad = 0, worst = null;
  const bySrc = {};
  for (const r of chk.rows) {
    const g = CG.fromSmiles(R, r.smiles);
    if (!g) { bad++; continue; }
    const z = M.encode(g, r.desc);
    const p = M.lkw(z, r.M, r.IF);
    const e1 = Math.abs(p - r.lkw);
    if (e1 > d1) { d1 = e1; worst = r.smiles; }
    for (let i = 0; i < 8; i++) dz = Math.max(dz, Math.abs(z[i] - r.z8[i]));
    // 端から端まで: 記述子も JavaScript で
    const m = R.get_mol(r.smiles); m.remove_hs_in_place();
    const f = F.compute(R, m); m.delete();
    const desc = [].concat(Array.from(f.desc), Array.from(f.groups), Array.from(f.ion));
    for (let i = 0; i < desc.length; i++) { const a = desc[i], b = r.desc[i]; if (isFinite(a) && isFinite(b)) dd = Math.max(dd, Math.abs(a - b) / (1 + Math.abs(b))); }
    const p2 = M.lkw(M.encode(g, desc), r.M, r.IF);
    d2 = Math.max(d2, Math.abs(p2 - r.lkw));
    bySrc[r.source] = Math.max(bySrc[r.source] || 0, e1);
    n++;
  }
  console.log('rows', n, 'unparsed', bad, '| max |JS - PyTorch| lkw', d1.toExponential(2), '(z の最初の 8 個', dz.toExponential(2) + ')', '| 端から端まで', d2.toExponential(2), '| 記述子の相対差', dd.toExponential(2));
  console.log('worst', worst, JSON.stringify(Object.fromEntries(Object.entries(bySrc).map(([k, v]) => [k, +v.toExponential(1)]))));
  // 速さ
  const smis = chk.rows.slice(0, 200).map(r => r.smiles), gs = [], ds = [];
  let t = now();
  for (const s of smis) gs.push(CG.fromSmiles(R, s));
  const tGraph = (now() - t) / smis.length;
  chk.rows.slice(0, 200).forEach(r => ds.push(r.desc));
  for (let i = 0; i < 20; i++) M.encode(gs[i], ds[i]);
  t = now(); const zs = [];
  for (let i = 0; i < gs.length; i++) zs.push(M.encode(gs[i], ds[i]));
  const tEnc = (now() - t) / gs.length;
  const r0 = chk.rows[0];
  for (let i = 0; i < 2000; i++) M.lkw(zs[i % zs.length], r0.M, r0.IF);
  t = now();
  for (let i = 0; i < 20000; i++) M.lkw(zs[i % zs.length], r0.M, r0.IF);
  const tHead = (now() - t) / 20000;
  const atoms = gs.reduce((a, g) => a + g.n, 0) / gs.length;
  const out = { dir: path.basename(dir), params: man.n_params, bytes: buf.byteLength, max_abs_diff_lkw: d1, max_abs_diff_end_to_end: d2, rows: n,
                graph_ms_per_molecule: tGraph, encode_ms_per_molecule: tEnc, head_ms_per_prediction: tHead, mean_atoms: atoms, node: process.version };
  fs.writeFileSync(path.join(dir, 'js_check.json'), JSON.stringify(out, null, 1));
  console.log('speed: graph %s ms/mol, encode %s ms/mol (平均 %s 原子), head %s ms/予測'.replace('%s', tGraph.toFixed(2)).replace('%s', tEnc.toFixed(2)).replace('%s', atoms.toFixed(0)).replace('%s', tHead.toFixed(4)));
})();
