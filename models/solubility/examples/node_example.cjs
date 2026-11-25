// SPDX-License-Identifier: MIT
// いちばん短い使い方（Node）。models/solubility/ で:  node examples/node_example.cjs
// Node 18 以上。ほかに入れる物はない（構造を読む RDKit.js はリポジトリの vendor/rdkit を使う）。
'use strict';
const path = require('path');
const SolubilityModel = require('../js/solubility.js');
const initRDKitModule = require('../../../vendor/rdkit/RDKit_minimal.js');

const S = SolubilityModel.loadNode(path.join(__dirname, '..', 'weights'));
const mg = v => v >= 100 ? v.toFixed(0) : v >= 1 ? v.toFixed(2) : v.toExponential(1);
const JA = { yes: '溶ける見込み', edge: '境目', no: '溶けない見込み', stays_dissolved: '析出しない見込み', precipitates: '析出する見込み' };

// 1. 例の固体（RDKit.js なしで使える）。目標 10 mg/mL で、主な 16 溶媒
const caf = S.fromExample('カフェイン');
console.log('■ ' + caf.name + '（融点 ' + caf.TmC + ' °C）、目標 10 mg/mL');
for (const e of S.estimateAll(caf, { target: 10 })) {
  const sv = S.solvent(e.solvent);
  console.log('  ' + (sv.ja || sv.name).padEnd(14) + ' 見込み ' + mg(e.mgPerML).padStart(7) + ' mg/mL（' + mg(e.range[0]) + '〜' + mg(e.range[1]) + '）' +
    (e.measured ? '  実測 ' + mg(e.measured.mgPerML) + '（doi:' + e.measured.doi + '）' : '') + '  → ' + JA[e.verdict] + '（' + e.basis + '）');
}

// 2. 手元の実測 1 点（メタノールに 25 mg/mL）から、ほかの溶媒を換算する
const r = S.estimate(caf, 'acn', { reference: { solvent: 'meoh', mgPerML: 25 } });
console.log('\n■ メタノールの実測 25 mg/mL から換算したアセトニトリル: ' + mg(r.mgPerML) + ' mg/mL（約 ' + r.factor.toFixed(1) + ' 倍の幅）');

// 3. 混合溶媒（DMSO と水。水の体積分率 30 %）と、HPLC に打つ前の希釈
const m = S.mix(caf, 'dmso', 'water', 0.3, { target: 10 });
console.log('■ DMSO 70 % / 水 30 %: ' + mg(m.mgPerML) + ' mg/mL（' + mg(m.range[0]) + '〜' + mg(m.range[1]) + '）→ ' + JA[m.verdict]);
const d = S.dilution(caf, { stockSolvent: 'dmso', stockMgPerML: 50, stockUL: 20, diluentUL: 180, diluent: { water: 0.8, acn: 0.2 } });
console.log('■ DMSO 50 mg/mL の原液 20 µL を、水 80 % / アセトニトリル 20 % 180 µL で薄める: ' + mg(d.concMgPerML) + ' mg/mL（' + d.factor + ' 倍）、この組成での溶解度 ' + mg(d.solubilityMgPerML) + ' mg/mL → ' + JA[d.verdict]);

// 4. 自分の分子（SMILES）。RDKit.js（wasm）を読む。構造はどこにも送らない
//    融点は 入れた値 ＞ 溶解度の照合の表 ＞ 融点の照合の表（weights/tm_lookup.json）の順で自動で入る。名前は渡した名前か、なければ SMILES
initRDKitModule().then(R => {
  const SRC = { user: '入れた値', sol: '溶解度の表（BigSolDB）', mp: '融点の公開データ' };
  const c = S.prepare(R, 'CC(C)(C)c1ccc(O)cc1', { name: '4-tert-ブチルフェノール' });
  console.log('\n■ ' + c.name + '（' + c.smiles + '）: 溶解度の照合の表に ' + (c.inLookup ? 'ある' : 'ない') + '、融点 ' + (c.TmC === null ? '不明' : c.TmC + ' °C（' + SRC[c.tmSrc] + '）'));
  const tm = S.estimateTm(c);
  console.log('  融点の推算（目安・計算には使わない）: 約 ' + Math.round(tm.TmC) + ' °C（外れの約 7 割が ± ' + Math.round(tm.band) + ' °C）');
  const none = S.withTm(c, null);
  for (const id of ['water', 'etoh', 'toluene', 'hexane']) {
    const e = S.estimate(c, id, { target: 10 }), e0 = S.estimate(none, id, { target: 10 });
    console.log('  ' + id.padEnd(8) + ' 融点あり ' + mg(e.mgPerML).padStart(7) + ' mg/mL（約 ' + e.factor.toFixed(1) + ' 倍の幅）／融点を使わないと ' + mg(e0.mgPerML).padStart(7) + '（約 ' + e0.factor.toFixed(1) + ' 倍）' + (e0.capped ? ' ※純粋な溶質の濃度で切った' : ''));
  }
  const u = S.prepare(R, 'CCOC(=O)c1ccc(NC(=O)C2CC2)cc1');
  console.log('■ 名前を渡さない分子: name = ' + u.name + '、融点 ' + (u.TmC === null ? '不明（推算の目安 約 ' + u.tmHint + ' °C）' : u.TmC + ' °C'));
});
