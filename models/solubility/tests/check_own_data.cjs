// SPDX-License-Identifier: MIT
// 自社のデータの読み取り（train/lib/own_data.cjs）を、例のファイル train/examples/own_data_example.csv で確かめる。Python は要らない。
//   node tests/check_own_data.cjs [自社のデータの CSV]
// 確かめること：読めた行・読まなかった行（理由つき）、単位の換算（mg/mL と mol/L が同じ log10 mol/L になる）、溶媒の名前の対応。
'use strict';
const path = require('path');
const HERE = path.join(__dirname, '..');
const C = require(path.join(HERE, 'train', 'lib', 'common.cjs'));
const OWN = require(path.join(HERE, 'train', 'lib', 'own_data.cjs'));
const file = process.argv[2] ? path.resolve(process.argv[2]) : path.join(HERE, 'train', 'examples', 'own_data_example.csv');
let ok = true; const check = (name, cond, extra) => { if (!cond) ok = false; console.log((cond ? '合格  ' : '不合格') + '  ' + name + (extra ? '  ' + extra : '')); };

C.rdkit().then(R => {
  const solvents = C.readJSON(path.join(HERE, 'weights', 'solvents.json')).solvents, model = C.readJSON(path.join(HERE, 'weights', 'abraham_solute_model.json'));
  // 名前で引ける表は、学習（train/build_tables.cjs）と同じ：例の固体の日本語の名前だけ
  const nameToSmiles = new Map(C.readJSON(path.join(HERE, 'weights', 'example_compounds.json')).compounds.map(e => [e.name.toLowerCase(), e.smiles]));
  const own = OWN.readOwn([file], { R, solvents, model, nameToSmiles });
  const raw = C.readTable(file).rows.filter(r => r.some(c => String(c).trim()));
  console.log('ファイル ' + file + '：' + raw.length + ' 行 → 溶質 × 溶媒 ' + own.rows.length + ' 組（溶質 ' + own.solutes.size + '）、読まなかった行 ' + own.skipped.length);
  own.skipped.forEach(s => console.log('  読まなかった  ' + path.basename(s.line) + '  ' + s.why));
  check('学習に使える行が 1 つ以上ある', own.rows.length > 0, own.rows.length ? '' : '（どの行も読めなかった。上の理由を見て直す）');
  check('どの行も、読めたか理由つきで飛ばしたか', own.rows.reduce((a, r) => a + r.n, 0) + own.skipped.length === raw.length);
  check('値は log10 mol/L で -12〜2 の間', own.rows.every(r => r.logS > -12 && r.logS < 2));
  // 単位の換算：同じ量を mg/mL と mol/L で書いた 2 行は同じ値になる
  const M = 151.165, a = OWN.UNITS['mg/ml'](10, M), b = OWN.UNITS['mol/l'](10 / M), c = OWN.UNITS['mmol/l'](10000 / M), d = OWN.UNITS['ug/ml'](10000, M);
  check('単位の換算（mg/mL・mol/L・mmol/L・µg/mL）', Math.abs(a - b) < 1e-12 && Math.abs(a - c) < 1e-12 && Math.abs(a - d) < 1e-12);
  // 溶媒の指定：id・英語の名前・日本語の名前・SMILES のどれでも同じ溶媒になる
  const sv = OWN.solventResolver(R, solvents), ids = ['acn', 'Acetonitrile', 'アセトニトリル', 'MeCN'].map(k => (sv(k) || {}).id).concat([(sv('', 'CC#N') || {}).id]);
  check('溶媒の指定（id・名前・日本語・略号・SMILES）', ids.every(x => x === 'acn'), ids.join(','));
  if (!process.argv[2]) {
    check('例のファイル：50 行のうち 47 行を使い、40 °C の 3 行を飛ばす', own.rows.reduce((s, r) => s + r.n, 0) === 47 && own.skipped.length === 3 && own.skipped.every(s => /20〜30/.test(s.why)));
    check('例のファイル：溶質 11 種', own.solutes.size === 11);
  }
  console.log(ok ? 'すべて合格' : '不合格がある');
  process.exit(ok ? 0 : 1);
}).catch(e => { console.error(e); process.exit(1); });
