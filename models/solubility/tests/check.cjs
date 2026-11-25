// SPDX-License-Identifier: MIT
// JavaScript の予測が、学習側（Python・scikit-learn）の予測と合うかを確かめる。
//   node tests/check.cjs              （models/solubility/ で実行。公開の重み weights/ と、照合用の値 tests/reference_*.json）
//   node tests/check.cjs <フォルダ>     （学習し直した重みを確かめる時。重みも照合用の値も、そのフォルダの物を使う。
//                                       train/train_gbt.py・train_tm.py が書いた train/out など）
// 確かめること:
//   1. 木の計算：同じ入力の行で、JS（js/solfeat.js）と scikit-learn の予測の差（4 つの溶解度のモデル・融点のモデル）
//   2. 構造から：照合の行の SMILES から JS で入力を作り直し（RDKit.js）、入力の差と予測の差
//   3. 入口（js/solubility.js）：サイトのページ（第 5 版）の計算と同じ値になるか（tests/reference_api.json。公開の重みの時だけ）。
//      SMILES から足した分子では、融点がページと同じ順（溶解度の照合の表 → 融点の照合の表 tm_lookup.json）で入るかも確かめる
// 合格の目安：1 は 1e-4 未満（融点は 0.05 °C 未満）、2 は入力 1e-4 未満・予測 1e-3 未満、3 は 1e-9 未満。
'use strict';
const fs = require('fs'), path = require('path');
const HERE = path.join(__dirname, '..');
const DIR = process.argv[2] ? path.resolve(process.argv[2]) : null;
// 照合の値：フォルダを渡せばそこ。渡さなければ tests/own/（自社のデータだけで学習して --install した時。git に入らない）か tests/
const OWNREF = path.join(__dirname, 'own');
const W = DIR || path.join(HERE, 'weights'), REF = DIR || (fs.existsSync(path.join(OWNREF, 'reference_gbt.json')) ? OWNREF : __dirname);
const SM = require(path.join(HERE, 'js', 'solubility.js'));
const SF = require(path.join(HERE, 'js', 'solfeat.js'));
const FE = require(path.join(HERE, 'js', 'features.js'));
const AB = require(path.join(HERE, 'js', 'abraham.js'));
const initRDKitModule = require(path.join(HERE, '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'));
const read = p => JSON.parse(fs.readFileSync(p, 'utf8'));
// 重みの見分け（4 つのモデルの最初の木と init から作る短い印）。入口の照合は、この印が同じ重みの時だけ行う
const fingerprint = G => { let h = 2166136261; const s = JSON.stringify(['tm', 'notm', 'tm_w', 'notm_w'].map(k => G[k] ? [G[k].init, G[k].trees.length, G[k].trees[0], G[k].trees[G[k].trees.length - 1]] : null)); for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619) >>> 0; } return h.toString(16); };
const has = p => fs.existsSync(p);
let ok = true;
const report = (name, v, lim, n) => { const pass = v < lim; if (!pass) ok = false; console.log((pass ? '合格  ' : '不合格') + '  ' + name + '  最大の差 ' + v.toExponential(2) + '（目安 < ' + lim + '）' + (n !== undefined ? '  ' + n + ' 件' : '')); };

(async () => {
  // 重みは W から。学習し直した重みのフォルダに無い表（溶媒の表など）は weights/ の物を使う
  const wfile = f => has(path.join(W, f)) ? path.join(W, f) : path.join(HERE, 'weights', f);
  const S = SM.create({ gbt: read(wfile('solubility_gbt.json')), formula: read(wfile('formula.json')), mix: read(wfile('mix_fit.json')), solvents: read(wfile('solvents.json')),
    abraham: read(wfile('abraham_solute_model.json')), tm: has(wfile('tm_model.json')) ? read(wfile('tm_model.json')) : null, examples: read(wfile('example_compounds.json')),
    lookup: has(wfile('lookup_25C.json')) ? read(wfile('lookup_25C.json')) : null, tmLookup: has(wfile('tm_lookup.json')) ? read(wfile('tm_lookup.json')) : null });
  const G = read(wfile('solubility_gbt.json')), formula = read(wfile('formula.json')), ABM = read(wfile('abraham_solute_model.json'));
  console.log('重み ' + W + '（版 ' + (G.meta && G.meta.version) + '）　照合の値 ' + REF);

  // 1. 木の計算
  const rg = path.join(REF, 'reference_gbt.json');
  if (DIR && !has(path.join(DIR, 'solubility_gbt.json'))) { console.log('不合格  ' + DIR + ' に solubility_gbt.json がない（学習し直した重みが無いまま、公開の重みで照合してしまうため止める）'); ok = false; }
  if (!has(rg)) { console.log('reference_gbt.json がない'); ok = false; }
  const ref = has(rg) ? read(rg) : { models: {} };
  if (!Object.keys(ref.models).length) { console.log('不合格  照合の値にモデルが 1 つもない'); ok = false; }
  for (const key of Object.keys(ref.models)) {
    const r = ref.models[key]; let d = 0;
    if (!G[key] || !(r.rows && r.rows.length) || !(r.pred && r.pred.length === r.rows.length)) { console.log('不合格  ' + key + '：重みにこのモデルがない、または照合の行が空・数が合わない'); ok = false; continue; }
    r.rows.forEach((row, i) => { d = Math.max(d, Math.abs(SF.evalModel(G[key], row) - r.pred[i])); });
    report('木の計算 ' + key.padEnd(6), d, 1e-4, r.rows.length);
  }
  const rt = path.join(REF, 'reference_tm.json');
  if (has(rt) && has(wfile('tm_model.json'))) {
    const r = read(rt), TM = read(wfile('tm_model.json')); let d = 0;
    r.rows.forEach((row, i) => { d = Math.max(d, Math.abs(SF.evalModel(TM.model, row) - r.pred[i])); });
    report('融点のモデル（°C）', d, 0.05, r.rows.length);
  }

  // 2. 構造から入力を作り直す
  const R = await initRDKitModule();
  const SV = new Map(S.solvents.map(s => [s.name, s]));
  for (const key of Object.keys(ref.models)) {
    const r = ref.models[key]; if (!r.can || !G[key]) continue;
    let dx = 0, dy = 0, n = 0, nBig = 0;
    r.can.forEach((can, i) => {
      const sv = SV.get(r.solvent[i]); if (!sv) return;
      const m = R.get_mol(can); m.remove_hs_in_place(); const f = FE.compute(R, m), d = JSON.parse(m.get_descriptors()); m.delete();
      const x = f.desc.concat(f.groups), abr = AB.predictSolute(ABM, x); abr.V = AB.mcgowanV(R, can);
      const sol = { M: d.amw, logPow: d.CrippenClogP, TmC: key.startsWith('tm') ? r.tmc[i] : null, abr, x };
      if (key.startsWith('tm') && r.tmc[i] === null) return;
      const row = SF.build(sol, sv, formula.K, formula.tm_default_C, 25);
      const cols = G[key].cols; for (const c of cols) dx = Math.max(dx, Math.abs(row[c] - r.rows[i][c]));
      const e = Math.abs(SF.evalModel(G[key], row) - r.pred[i]); dy = Math.max(dy, e); if (e > 1e-3) nBig++; n++;
    });
    report('構造から ' + key.padEnd(6) + ' 入力', dx, 1e-4, n);
    report('構造から ' + key.padEnd(6) + ' 予測', dy, 1e-3, n);
  }

  // 3. 入口の計算（サイトのページ 第 5 版と同じか）
  const ra = path.join(REF, 'reference_api.json');
  const A = has(ra) ? read(ra) : null;
  if (A && A.weights_fingerprint !== fingerprint(G)) console.log('（入口の照合は飛ばす：reference_api.json は公開の重み（印 ' + A.weights_fingerprint + '）の値。今の重みの印は ' + fingerprint(G) + '）');
  else if (A) {
    let d = 0, n = 0, vbad = 0, pbad = 0, nPrep = 0, nTmLookup = 0;
    const V = { yes: 'stays_dissolved', no: 'precipitates', edge: 'edge' };
    for (const c of A.cases) {
      if (c.kind === 'prepare') {
        const p = S.prepare(R, c.smiles); nPrep++; if (!p.inLookup && p.tmSrc === 'mp') nTmLookup++;
        if (p.smiles !== c.can || p.TmC !== c.TmC || p.tmSrc !== c.tmSrc || p.tmHint !== c.tmHint || p.inLookup !== c.inLookup) { pbad++; if (pbad <= 3) console.log('  違う: ' + c.smiles + ' 融点 ' + p.TmC + '/' + c.TmC + ' 出所 ' + p.tmSrc + '/' + c.tmSrc); }
        for (const o of c.out) { const e = S.estimate(p, o.solvent); d = Math.max(d, Math.abs(e.logS - o.logS), Math.abs(e.band - o.band)); n++; }
        continue;
      }
      const sol = S.withTm(S.fromExample(c.example), c.TmC);
      if (c.kind === 'estimate' || c.kind === 'reference') for (const o of c.out) {
        const e = S.estimate(sol, o.solvent, c.ref ? { reference: c.ref } : {});
        d = Math.max(d, Math.abs(e.logS - o.logS), Math.abs(e.band - o.band)); n++;
      }
      if (c.kind === 'mix') for (const o of c.out) { const e = S.mix(sol, o.a, o.b, o.f); d = Math.max(d, Math.abs(Math.log10(e.mgPerML) - o.y), Math.abs(e.band - o.band)); n++; }
      if (c.kind === 'dilution') for (const o of c.out) {
        const e = S.dilution(sol, { stockSolvent: o.stock, stockMgPerML: o.conc, stockUL: o.take, diluentUL: o.vol, diluent: { water: 1 - o.acn, acn: o.acn } });
        d = Math.max(d, Math.abs(Math.log10(e.solubilityMgPerML) - o.y), Math.abs(e.band - o.band), Math.abs(e.factor - o.factor) / o.factor); if (e.verdict !== V[o.v]) vbad++; n++;
      }
    }
    report('入口（ページ 第 5 版との差）', d, 1e-9, n);
    if (nPrep) { const okP = pbad === 0 && nTmLookup > 0; if (!okP) ok = false; console.log((okP ? '合格  ' : '不合格') + '  SMILES から足した分子 ' + nPrep + ' 件：構造・融点・出所・推算の目安がページと同じ（違う ' + pbad + ' 件）。融点の照合の表から融点が入った分子 ' + nTmLookup + ' 件'); }
    if (vbad) { ok = false; console.log('不合格  希釈の判定が違う ' + vbad + ' 件'); }
  } else console.log('（reference_api.json がないので、入口の照合は飛ばす）');

  // 4. 速さの目安
  const sol = S.fromExample(S.examples[0].id), t0 = Date.now(); let k = 0;
  while (Date.now() - t0 < 300) { S.estimateAll(sol); k++; }
  console.log('速さ  16 溶媒の見込み 1 回あたり ' + (300 / k).toFixed(2) + ' ms（' + k + ' 回）');
  console.log(ok ? 'すべて合格' : '不合格がある');
  process.exit(ok ? 0 : 1);
})().catch(e => { console.error(e); process.exit(1); });
