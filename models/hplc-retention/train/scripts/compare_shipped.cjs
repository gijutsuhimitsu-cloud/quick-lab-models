// SPDX-License-Identifier: MIT
// 作り直した重みが、公開の重み（weights/）を再現しているかを確かめる。
//   node compare_shipped.cjs <作り直した重みのフォルダ>
// 1. 公開の照合用の値（tests/reference_free.json・reference_uncertainty_free.json）を、作り直した重みで計算した時の差
//    （tests/check.cjs と同じ計算。合格の目安は lkw の差 1e-4 未満。公開の重みそのものでも 1e-5 程度の差がある＝葉の値を小数 6 桁に丸めているため）
// 2. 木の形（特徴量・しきい値・葉の値）の一致、3. ファイルの中身が 1 バイトも違わないか
'use strict';
const fs = require('fs'), path = require('path');
const MODEL = path.join(__dirname, '..', '..');
const CR = require(path.join(MODEL, 'js', 'retention_model.js'));
const read = (p) => JSON.parse(fs.readFileSync(p, 'utf8'));
const dir = path.resolve(process.argv[2] || '.');
const SW = path.join(MODEL, 'weights'), ST = path.join(MODEL, 'tests');

function treeDiff(a, b) {
  const r = { n: [a.trees.length, b.trees.length], sameShape: 0, maxT: 0, maxV: 0, baseline: Math.abs(a.baseline - b.baseline) };
  for (let i = 0; i < Math.min(a.trees.length, b.trees.length); i++) {
    const x = a.trees[i], y = b.trees[i];
    const same = ['f', 'l', 'r', 'm'].every(k => x[k].length === y[k].length && x[k].every((v, j) => v === y[k][j]));
    if (!same) continue;
    r.sameShape++;
    x.t.forEach((v, j) => { r.maxT = Math.max(r.maxT, Math.abs(v - y.t[j])); });
    x.v.forEach((v, j) => { r.maxV = Math.max(r.maxV, Math.abs(v - y.v[j])); });
  }
  return r;
}

const model = read(path.join(dir, 'model_free.json')), unc = read(path.join(dir, 'uncertainty_free.json'));
const columns = read(path.join(dir, 'columns.json')), stock = read(path.join(dir, 'stock_compounds.json')).compounds;
const ref = read(path.join(ST, 'reference_free.json'));
let dL = 0;
ref.x.forEach((row, i) => { dL = Math.max(dL, Math.abs(CR.evalTrees(model, row.map(v => v === null ? NaN : v)) - ref.y[i])); });
const P = CR.create(model, columns, unc);
let dE = 0, dR = 0;
for (const e of ref.e2e) {
  const r = P.predictRetention(stock[e.stock], e.column, e.cond);
  dE = Math.max(dE, Math.abs(r.lkw - e.lkw));
  if (e.tR != null) dR = Math.max(dR, Math.abs(r.tR - e.tR) / e.tR);
  if (e.k != null) dR = Math.max(dR, Math.abs(r.k - e.k) / e.k);
}
const cu = read(path.join(ST, 'reference_uncertainty_free.json'));
let d68 = 0, d90 = 0;
cu.X.forEach((row, i) => {
  const x = row.map(v => v === null ? NaN : v);
  d68 = Math.max(d68, Math.abs(Math.max(CR.evalTrees(unc.s68, x), unc.meta.floor) * unc.s68.scale - cu.s68[i]));
  d90 = Math.max(d90, Math.abs(Math.max(CR.evalTrees(unc.s90, x), unc.meta.floor) * unc.s90.scale - cu.s90[i]));
});
console.log('公開の照合用の値（tests/）との差: lkw ' + dL.toExponential(2) + '（300 行）, predictRetention の lkw ' + dE.toExponential(2) +
  '・k/tR の相対差 ' + dR.toExponential(2) + '（' + ref.e2e.length + ' 件）, 幅 s68 ' + d68.toExponential(2) + '・s90 ' + d90.toExponential(2));
const sm = read(path.join(SW, 'model_free.json')), su = read(path.join(SW, 'uncertainty_free.json'));
const tm = treeDiff(model, sm), t68 = treeDiff(unc.s68, su.s68), t90 = treeDiff(unc.s90, su.s90);
const show = (n, t) => console.log('木の形（' + n + '）: 本数 ' + t.n.join(' / ') + '、形が同じ木 ' + t.sameShape + '、しきい値の差の最大 ' + t.maxT.toExponential(2) +
  '、葉の値の差の最大 ' + t.maxV.toExponential(2) + '、baseline の差 ' + t.baseline.toExponential(2));
show('本体', tm); show('幅 s68', t68); show('幅 s90', t90);
const same = {};
for (const f of ['model_free.json', 'uncertainty_free.json', 'stock_compounds.json', 'name_index_free.json', 'columns.json']) {
  same[f] = fs.readFileSync(path.join(dir, f)).equals(fs.readFileSync(path.join(SW, f)));
}
console.log('ファイルが 1 バイトも違わない:', Object.entries(same).map(([k, v]) => k + ' ' + (v ? 'はい' : 'いいえ')).join('、'));
const ok = dL < 1e-4 && dE < 1e-4 && dR < 1e-3 && d68 < 1e-4 && d90 < 1e-4;
console.log(ok ? '再現できている（公開の照合用の値との差が目安の中）' : '再現できていない（差が目安を超えた）');
fs.writeFileSync(path.join(dir, '..', 'compare_shipped.json'), JSON.stringify({ lkw: dL, e2e_lkw: dE, e2e_rel: dR, s68: d68, s90: d90, trees: { model: tm, s68: t68, s90: t90 }, identical_files: same, ok }, null, 1));
process.exit(ok ? 0 : 1);
