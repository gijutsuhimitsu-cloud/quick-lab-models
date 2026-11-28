// SPDX-License-Identifier: MIT
// 使い方の例（Node）。 node examples/node_example.cjs   （models/solubility-ensemble/ で）
'use strict';
const path = require('path');
const E = require(path.join(__dirname, '..', 'js', 'solubility_ensemble.js'));
const S = E.loadNode();
const f = v => v === null || v === undefined ? '—' : v.toPrecision(3);
const line = e => `  ${e.solvent.padEnd(8)} ${f(e.mgPerML)} mg/mL` + (e.mgPerML === null ? `（${e.reason}）` : `（約 ${e.factor.toFixed(1)} 倍の幅）`) + (e.verdict ? ' ' + e.verdict : '') + (e.measured ? `  文献の実測 ${f(e.measured.mgPerML)}` : '');

const caf = S.fromExample('カフェイン');
console.log('カフェイン、目標 10 mg/mL');
S.estimateAll(caf, { target: 10 }).forEach(e => console.log(line(e)));

console.log('メタノールで 8 mg/mL と測った時の、アセトニトリルの見込み');
console.log(line(S.estimate(caf, 'acn', { references: [{ solvent: 'meoh', mgPerML: 8 }] })));

const m = S.mix(caf, 'acn', 'water', 0.3);
console.log('アセトニトリル 70 % / 水 30 %:', f(m.mgPerML), 'mg/mL');

const d = S.dilution(caf, { stockSolvent: 'dmso', stockMgPerML: 50, stockUL: 10, diluentUL: 90, diluent: { water: 0.7, acn: 0.3 } });
console.log('DMSO の原液 50 mg/mL を 10 倍に希釈:', d.verdict, '（薄めた後', f(d.concMgPerML), 'mg/mL、溶解度の見込み', f(d.solubilityMgPerML), 'mg/mL）');

require(path.join(__dirname, '..', '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'))().then(R => {
  const c = S.prepare(R, 'CC(=O)Oc1ccccc1C(=O)O', { name: 'アスピリン' });
  console.log('SMILES から:', c.name, '融点', c.TmC, '°C');
  S.estimateAll(c, { solvents: ['water', 'etoh', 'acn', 'toluene'] }).forEach(e => console.log(line(e)));
});
