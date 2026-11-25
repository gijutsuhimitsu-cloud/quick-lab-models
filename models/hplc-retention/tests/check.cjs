// JavaScript の予測が、学習側（Python）の予測と合うかを確かめる。
//   node tests/check.cjs              （models/hplc-retention/ で実行）
//   node tests/check.cjs <重みのフォルダ>   （学習し直した重みを確かめる時。照合用の値はそのフォルダの reference_*.json）
// 合格の目安：lkw の差 1e-4 未満、保持時間の相対差 1e-3 未満、予測の幅の差 1e-4 未満。
'use strict';
const fs = require('fs'), path = require('path');
const HERE = path.join(__dirname, '..');
const W = process.argv[2] ? path.resolve(process.argv[2]) : path.join(HERE, 'weights');
const REF = process.argv[2] ? W : __dirname;
const CR = require(path.join(HERE, 'js', 'retention_model.js'));
const F = require(path.join(HERE, 'js', 'features.js'));
const initRDKitModule = require(path.join(HERE, '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'));
const read = (p) => JSON.parse(fs.readFileSync(p, 'utf8'));

(async () => {
  const model = read(path.join(W, 'model_free.json'));
  const unc = read(path.join(W, 'uncertainty_free.json'));
  const columns = read(path.join(W, 'columns.json'));
  const stock = read(path.join(W, 'stock_compounds.json')).compounds;
  const chk = read(path.join(REF, 'reference_free.json'));
  let ok = true;

  // 1. 木の出力（lkw）
  let dLkw = 0;
  chk.x.forEach((row, i) => { const x = row.map(v => v === null ? NaN : v); dLkw = Math.max(dLkw, Math.abs(CR.evalTrees(model, x) - chk.y[i])); });
  // 2. 勾配の保持時間の計算
  let dRt = 0;
  for (const c of chk.rt) c.lkw.forEach((lk, i) => { dRt = Math.max(dRt, Math.abs(CR.retentionTime(lk, c.t0, c.prog, model.meta.par.acn, c.dwell, 8000) - c.rt[i]) / c.rt[i]); });
  // 3. ページと同じ入口（predictRetention）を通した値
  const P = CR.create(model, columns, unc);
  let dE = 0, dErel = 0, nE = 0;
  for (const e of (chk.e2e || [])) {
    const r = P.predictRetention(stock[e.stock], e.column, e.cond);
    dE = Math.max(dE, Math.abs(r.lkw - e.lkw));
    if (e.tR != null) dErel = Math.max(dErel, Math.abs(r.tR - e.tR) / e.tR);
    if (e.k != null) dErel = Math.max(dErel, Math.abs(r.k - e.k) / e.k);
    nE++;
  }
  console.log('lkw  max |JS - Python| = ' + dLkw.toExponential(2) + '  (' + chk.x.length + ' 行)');
  console.log('保持時間  最大の相対差 = ' + dRt.toExponential(2));
  console.log('predictRetention  lkw ' + dE.toExponential(2) + ' / k・tR の相対差 ' + dErel.toExponential(2) + '  (' + nE + ' 件)');
  if (!(dLkw < 1e-4 && dRt < 1e-3 && dE < 1e-4 && dErel < 1e-3)) ok = false;

  // 4. 予測の幅
  const up = path.join(REF, 'reference_uncertainty_free.json');
  if (fs.existsSync(up)) {
    const cu = read(up); let d68 = 0, d90 = 0;
    cu.X.forEach((row, i) => {
      const x = row.map(v => v === null ? NaN : v);
      d68 = Math.max(d68, Math.abs(Math.max(CR.evalTrees(unc.s68, x), unc.meta.floor) * unc.s68.scale - cu.s68[i]));
      d90 = Math.max(d90, Math.abs(Math.max(CR.evalTrees(unc.s90, x), unc.meta.floor) * unc.s90.scale - cu.s90[i]));
    });
    console.log('予測の幅  max |JS - Python| = s68 ' + d68.toExponential(2) + ', s90 ' + d90.toExponential(2));
    if (!(d68 < 1e-4 && d90 < 1e-4)) ok = false;
  }

  // 5. 構造（SMILES）から特徴量を作る道も、ストックの事前計算と同じ値になるか
  const R = await initRDKitModule();
  let dS = 0, nS = 0;
  for (const c of stock.slice(0, 40)) {
    const m = R.get_mol(c.smiles); if (!m || !m.is_valid()) continue;
    m.remove_hs_in_place(); const f = F.compute(R, m); m.delete();
    const cond = { organic: 'acn', phi: 0.4, pH: 3, tempC: 30 };
    const a = P.predictRetention({ desc: f.desc, groups: f.groups, ion: f.ion }, columns.defaults[0].id, cond).lkw;
    const b = P.predictRetention({ desc: c.desc, groups: c.groups, ion: c.ion }, columns.defaults[0].id, cond).lkw;
    dS = Math.max(dS, Math.abs(a - b)); nS++;
  }
  console.log('SMILES から作った特徴量とストックの事前計算  lkw の差 ' + dS.toExponential(2) + '  (' + nS + ' 分子)');
  if (!(dS < 1e-6)) ok = false;

  console.log(ok ? '合格' : '不合格');
  process.exit(ok ? 0 : 1);
})();
