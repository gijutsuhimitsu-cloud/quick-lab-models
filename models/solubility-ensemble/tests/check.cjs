// SPDX-License-Identifier: MIT
// 呼び出し口（js/solubility_ensemble.js）の値が、サイトのページの計算（tests/reference.json）と一致するかを確かめる。
// 使い方: node tests/check.cjs [重みのフォルダ]   （models/solubility-ensemble/ で。Node 18 以上）
'use strict';
const fs = require('fs'), path = require('path');
const ROOT = path.resolve(__dirname, '..');
const E = require(path.join(ROOT, 'js', 'solubility_ensemble.js'));
const REF = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference.json'), 'utf8'));
const S = E.loadNode(process.argv[2] ? path.resolve(process.argv[2]) : undefined);

let bad = 0, n = 0, maxDiff = 0;
const show = [];
function near(a, b, tol, what) { n++; const d = Math.abs(a - b); if (d > maxDiff && tol < 1e-3) maxDiff = d; if (!(d <= tol)) { bad++; if (show.length < 8) show.push(what + ' ' + a + ' / ' + b); } }
function same(a, b, what) { n++; if (a !== b) { bad++; if (show.length < 8) show.push(what + ' ' + a + ' / ' + b); } }
const byId = id => S.fromExample(id);
const refs = list => list.map(r => ({ solvent: r[0], mgPerML: r[1] }));

for (const c of REF.presets) {
  let sol = byId(c.id); if (c.noTm) sol = Object.assign(S.withTm(sol, null), { smiles: sol.smiles });
  for (const sv in c.expect) {
    const x = c.expect[sv], e = S.estimate(sol, sv, { references: refs(c.refs) }), w = 'preset ' + c.id + ' ' + sv + (c.noTm ? ' 融点なし' : '') + ' ref' + c.refs.length;
    if (x.oo) same(e.outOfScope, x.oo, w);
    else { near(e.logS, x.logS, 1e-9, w + ' logS'); near(e.band, x.band, 1e-9, w + ' band'); }
  }
}
for (const c of REF.mix) {
  const m = S.mix(byId(c.id), c.a, c.b, c.f, { references: refs(c.refs) }), w = 'mix ' + c.id + ' ' + c.a + '/' + c.b + ' ' + c.f;
  if (c.expect.na) same(m.mgPerML, null, w);
  else if (!c.expect.one) { same(m.onePhase, false, w); same(m.mgPerML, null, w); }
  else { near(Math.log10(m.mgPerML), c.expect.y, 1e-9, w + ' y'); near(m.band, c.expect.band, 1e-9, w + ' band'); }
}
for (const c of REF.dilution) {
  const p = c.p, d = S.dilution(byId(c.id), { stockSolvent: p.sv, stockMgPerML: p.conc, stockUL: p.take, diluentUL: p.vol, diluent: { water: 1 - p.acn, acn: p.acn } }), w = 'dilution ' + c.id + ' ' + JSON.stringify(p);
  if (c.expect.na) same(d.solubilityMgPerML, null, w);
  else {
    near(Math.log10(d.solubilityMgPerML), c.expect.y, 1e-9, w + ' y'); near(d.band, c.expect.band, 1e-9, w + ' band'); near(d.concMgPerML, c.expect.conc, 1e-9, w + ' conc');
    same(d.verdict, c.expect.v === 'yes' ? 'stays_dissolved' : c.expect.v === 'no' ? 'precipitates' : 'edge', w + ' verdict');
  }
}

require(path.resolve(ROOT, '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'))().then(R => {
  for (const c of REF.molecules) {
    const x = c.expect, w = 'mol ' + c.smiles;
    let sol; try { sol = S.prepare(R, c.smiles); } catch (e) { same('error ' + e.code, 'ok', w); continue; }
    same(sol.TmC, x.tm, w + ' 融点'); same(sol.M, x.M, w + ' M'); same(sol.outOfScope, x.oo, w + ' 対象外');
    for (const sv in x.tmv) {
      const [y, band, lit] = x.tmv[sv], e = S.estimate(sol, sv);
      if (y === null) same(e.mgPerML, null, w + ' ' + sv);
      else { near(Math.log10(e.mgPerML), y, 0.0051, w + ' ' + sv + ' y'); near(e.band, band, 0.0051, w + ' ' + sv + ' band'); }
      if (lit !== null) near(e.measured ? Math.log10(e.measured.mgPerML) : NaN, lit, 0.0051, w + ' ' + sv + ' 実測');
    }
  }
  for (const c of REF.errors) { let code = 'ok'; try { S.prepare(R, c.smiles); } catch (e) { code = e.code; } same(code, c.code, 'error ' + JSON.stringify(c.smiles.slice(0, 20))); }
  if (n < 1000) { console.log('照合が少なすぎる:', n); bad++; }
  console.log('照合', n, '件、不一致', bad, '件、最大の差（例の固体）', maxDiff.toExponential(2));
  if (bad) { console.log(show.join('\n')); process.exit(1); }
  console.log('OK');
});
