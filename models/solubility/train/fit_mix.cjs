// SPDX-License-Identifier: MIT
/* 混合溶媒の補正の係数（weights/mix_fit.json と同じ形）を合わせ直す。
     対数線形        log C(φ) = (1 − φ) log C_A + φ log C_B                 φ は溶媒 B の体積分率
     補正つき        log C(φ) = 上の式 + φ(1 − φ) [A0 + A1 (1 − 2φ)]       A0・A1 は溶媒の組ごとの定数（溶質によらない）
   入力: <base>/mix_series_25C.csv.gz（MixtureSolDB、CC BY 4.0 から作った表。作り直しは train/rebuild/mix_from_raw.cjs）
   使い方（models/solubility/ で）: node train/fit_mix.cjs
   出力: <out>/mix_fit.json
   係数は、溶質ごとに 5 つの組に分け、4 組で合わせて残りの 1 組で外れを測る。溶質が 8 種に満たない組は、組の種類（水＋有機、有機＋有機）の係数を使う。 */
'use strict';
const path = require('path');
const C = require('./lib/common.cjs'), P = C.P;
const BASE = path.resolve(process.env.SOL_BASE || P.BASE);
function stats(e) {
  if (!e.length) return { n: 0 };
  const a = e.map(Math.abs).sort((p, q) => p - q), n = a.length;
  return { n, bias: +(e.reduce((s, v) => s + v, 0) / n).toFixed(3), p50: +a[Math.floor(0.5 * (n - 1))].toFixed(3), p68: +a[Math.floor(0.68 * (n - 1))].toFixed(3), p90: +a[Math.floor(0.9 * (n - 1))].toFixed(3) };
}
function fit2(rows) {          // r ≈ A0·b0 + A1·b1 の最小二乗
  let s00 = 0, s01 = 0, s11 = 0, t0 = 0, t1 = 0;
  for (const r of rows) { s00 += r.b0 * r.b0; s01 += r.b0 * r.b1; s11 += r.b1 * r.b1; t0 += r.b0 * r.r; t1 += r.b1 * r.r; }
  const det = s00 * s11 - s01 * s01; if (!(Math.abs(det) > 1e-9)) return [s00 > 0 ? t0 / s00 : 0, 0];
  return [(t0 * s11 - t1 * s01) / det, (t1 * s00 - t0 * s01) / det];
}
const MAIN = C.MAIN_SOLVENTS.map(m => m[0]);
const t = C.readTable(C.need(path.join(BASE, 'mix_series_25C.csv.gz')));
const rows = t.objs().map(o => ({ solute: o.solute, pair: o.pair, main: o.main === '1', cls: o.cls, phi: +o.phi, r: +o.r, b0: +o.b0, b1: +o.b1 }));
const solutes = [...new Set(rows.map(r => r.solute))]; let seed = 4242; const rnd = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648;
const fold = new Map(solutes.map(s => [s, Math.floor(rnd() * 5)]));
const build = train => {
  const cls = { wo: fit2(train.filter(r => r.cls === 'wo')), oo: fit2(train.filter(r => r.cls === 'oo')) }, by = new Map();
  for (const r of train) { if (!by.has(r.pair)) by.set(r.pair, []); by.get(r.pair).push(r); }
  const pairs = new Map(); for (const [k, v] of by) { const ns = new Set(v.map(r => r.solute)).size; if (ns >= 8) pairs.set(k, fit2(v).concat([ns])); }
  return { cls, pairs };
};
const corr = (M, r) => { const A = M.pairs.get(r.pair) || M.cls[r.cls]; return A[0] * r.b0 + A[1] * r.b1; };
const E = { lin: { all: [], wo: [], oo: [], main: [] }, cor: { all: [], wo: [], oo: [], main: [] } }, mid = { lin: [], cor: [] };
for (let k = 0; k < 5; k++) {
  const M = build(rows.filter(r => fold.get(r.solute) !== k));
  for (const r of rows) if (fold.get(r.solute) === k) {
    const el = -r.r, ec = corr(M, r) - r.r;          // 外れ = 見込み − 実測
    E.lin.all.push(el); E.cor.all.push(ec); E.lin[r.cls].push(el); E.cor[r.cls].push(ec); if (r.main) { E.lin.main.push(el); E.cor.main.push(ec); }
    if (r.phi > 0.35 && r.phi < 0.65) { mid.lin.push(el); mid.cor.push(ec); }
  }
}
const S = o => Object.fromEntries(Object.entries(o).map(([k, v]) => [k, stats(v)]));
const cv = { log_linear: S(E.lin), corrected: S(E.cor), middle_35_65: { log_linear: stats(mid.lin), corrected: stats(mid.cor) } };
const M = build(rows), mainPairs = {};
for (const [k, v] of M.pairs) if (k.split('|').every(x => MAIN.includes(x))) mainPairs[k] = [+v[0].toFixed(3), +v[1].toFixed(3), v[2]];
const out = { meta: { name: 'mix_fit', version: '0.1', made: new Date().toISOString().slice(0, 10), license: 'CC BY 4.0（出典は ATTRIBUTION.md）' },
  source: 'MixtureSolDB (Malikov et al.), Zenodo, doi:10.5281/zenodo.18660057, CC BY 4.0。298.15 ± 5 K、純粋な溶媒 2 つの値がそろった系列の途中の点 ' + rows.length + ' 点',
  equation: 'log C = (1-phi) log C_A + phi log C_B + phi(1-phi)[A0 + A1(1-2phi)]。組の名前は id の辞書順で A|B、phi は B の体積分率',
  pairs: mainPairs, cls: { wo: M.cls.wo.map(v => +v.toFixed(3)), oo: M.cls.oo.map(v => +v.toFixed(3)) }, band: { log_linear: cv.log_linear.all.p68, corrected: cv.corrected.all.p68, corrected_wo: cv.corrected.wo.p68, corrected_oo: cv.corrected.oo.p68 }, cv };
C.writeJSON(path.join(P.OUT, 'mix_fit.json'), out, true);
console.log('点', rows.length, '溶質', solutes.length, '幅（約 7 割）対数線形', out.band.log_linear, '→ 補正つき', out.band.corrected, '（水＋有機', out.band.corrected_wo, '、有機どうし', out.band.corrected_oo, '）主な溶媒の組', Object.keys(mainPairs).length);
console.log('書いた: ' + path.join(P.OUT, 'mix_fit.json'));
