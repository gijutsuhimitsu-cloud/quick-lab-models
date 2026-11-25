// SPDX-License-Identifier: MIT
// グラフのモデルを、ページ用の重みのフォルダに置く（この写しでは未実行。GPU で作った重みがまだ無いため）。
//   node train/gpu/make_gnn_export.cjs <train_teacher.py --export 1 の出力のフォルダ> <重みのフォルダ（train/out/<名前>/weights）>
//   1. gnn_free.json（並びと正規化の定数）・gnn_free.bin（重み、float32）・gnn_free_weights.json（同じ重みの base64。.bin を置けない所用）
//   2. ストックの化合物の z（分子ごとの表現）の事前計算 stock_z_free.json（ストックだけなら RDKit を読まずに使える）
//   3. js/retention_model.js を通した確認（木とグラフの予測の例、pH の向き、事前計算の z とその場の計算の一致）
// ファイル名は weights/ のほかのファイルと同じく「_free」で固定。中の group（＝重みの meta.name）で版を見分ける。
'use strict';
const fs = require('fs'), path = require('path');
const MODEL = path.join(__dirname, '..', '..'), JS = path.join(MODEL, 'js');
const CR = require(path.join(JS, 'retention_model.js'));
const GNN = require(path.join(JS, 'gnn_model.js')), CG = require(path.join(JS, 'graph_features.js'));
const initRDKitModule = require(path.join(MODEL, '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'));
if (process.argv.length < 4) { console.error('使い方: node make_gnn_export.cjs <書き出しのフォルダ> <重みのフォルダ>'); process.exit(2); }
const dir = path.resolve(process.argv[2]), W = path.resolve(process.argv[3]);
(async () => {
  const R = await initRDKitModule();
  const man = JSON.parse(fs.readFileSync(path.join(dir, 'manifest.json'), 'utf8'));
  const model = JSON.parse(fs.readFileSync(path.join(W, 'model_free.json'), 'utf8'));
  if (man.group !== model.meta.name) throw new Error('グラフのモデルの group（' + man.group + '）と、木の重みの meta.name（' + model.meta.name + '）が違う');
  man.source_dir = path.basename(dir);
  fs.writeFileSync(path.join(W, 'gnn_free.json'), JSON.stringify(man));
  fs.copyFileSync(path.join(dir, 'weights.bin'), path.join(W, 'gnn_free.bin'));
  const buf = fs.readFileSync(path.join(W, 'gnn_free.bin'));
  fs.writeFileSync(path.join(W, 'gnn_free_weights.json'), JSON.stringify({ note: 'gnn_free.bin を base64 にした物', encoding: 'base64', bytes: buf.byteLength, data: buf.toString('base64') }));
  const Wt = new Float32Array(buf.buffer, buf.byteOffset, buf.byteLength / 4);
  const columns = JSON.parse(fs.readFileSync(path.join(W, 'columns.json'), 'utf8'));
  const stock = JSON.parse(fs.readFileSync(path.join(W, 'stock_compounds.json'), 'utf8'));
  const up = path.join(W, 'uncertainty_free_gnn.json');
  const unc = fs.existsSync(up) ? JSON.parse(fs.readFileSync(up, 'utf8')) : null;
  const G = GNN.create(man, Wt);
  const P = CR.create(model, columns, unc, G);
  const zs = {}; let bad = 0;
  for (const c of stock.compounds) {
    const g = CG.fromSmiles(R, c.smiles);
    if (!g) { bad++; continue; }
    P.embed(c, g);
    const key = String(c.inchikey || c.name);
    if (zs[key]) throw new Error('鍵が重なった: ' + key);
    zs[key] = Array.from(c.z).map(v => +v.toFixed(5));
  }
  fs.writeFileSync(path.join(W, 'stock_z_free.json'), JSON.stringify({ meta: { name: model.meta.name, zdim: G.zdim, gnn_params: man.n_params, key: 'InChIKey の全体（27 文字。立体の違う分子を分けるため）', note: 'compound.z = Float32Array.from(z[compound.inchikey]) として使う' }, z: zs }));
  console.log('ストックの z:', Object.keys(zs).length, '/', stock.compounds.length, '（グラフを作れない', bad + '）| 重み', man.n_params, '個', (buf.byteLength / 1e6).toFixed(2), 'MB');
  const col = columns.defaults[0].id;
  const iso = { organic: 'acn', phi: 0.4, pH: 3, tempC: 35, length_mm: 150, id_mm: 4.6, particle_um: 5, flow: 1.0 };
  let dmax = 0;
  for (const c of stock.compounds) {
    const a = P.predictRetention(c, col, iso).lkw;
    const c2 = Object.assign({}, c, { z: Float32Array.from(zs[String(c.inchikey || c.name)]) });
    dmax = Math.max(dmax, Math.abs(P.predictRetention(c2, col, iso).lkw - a));
  }
  console.log('事前計算の z（小数 5 桁）とその場の計算の lkw の差の最大:', dmax.toExponential(2));
  // pH の向き: pH を 2 → 11 に上げた時に lkw の向きが途中で変わる化合物（酸と塩基の両方を持つ物は変わってよい）。グラフだけで変わる物が 0 なら、向きの決まりは守られている
  let viol = 0, violTree = 0, onlyGnn = 0;
  for (const c of stock.compounds) {
    let prev = null, dirn = 0, ok = true, prevT = null, dirT = 0, okT = true;
    for (let pH = 2; pH <= 11.01; pH += 0.25) {
      const r = P.predictRetention(c, col, Object.assign({}, iso, { pH }));
      if (prev !== null && Math.abs(r.lkw - prev) > 1e-5) { const d = r.lkw > prev ? 1 : -1; if (dirn === 0) dirn = d; else if (d !== dirn) ok = false; }
      if (prevT !== null && Math.abs(r.lkwTree - prevT) > 1e-5) { const d = r.lkwTree > prevT ? 1 : -1; if (dirT === 0) dirT = d; else if (d !== dirT) okT = false; }
      prev = r.lkw; prevT = r.lkwTree;
    }
    if (!ok) viol++; if (!okT) violTree++; if (!ok && okT) onlyGnn++;
  }
  console.log('pH の向きが途中で変わる化合物: グラフ', viol, '・木', violTree, '/', stock.compounds.length, '| グラフだけで変わる:', onlyGnn);
})();
