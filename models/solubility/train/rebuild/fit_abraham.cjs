// SPDX-License-Identifier: MIT
/* 溶質の Abraham 記述子（E・S・A・B）を構造から推算する一次式（weights/abraham_solute_model.json）を作り直す。
   入力: work/abraham_solutes_raw.csv（train/rebuild/abraham_xls_to_csv.py が、CC BY 4.0 の表
         「Compounds with known Abraham descriptors」Bradley, Acree, Lang、figshare doi:10.6084/m9.figshare.1176994 から作る）
   使い方（models/solubility/ で）: node train/rebuild/fit_abraham.cjs
   出力: <out>/abraham_solute_model.json（係数と、交差検証での外れ）
         <out>/base_data/abraham_solutes_clean.csv.gz（構造ごとに平均した表。溶質の記述子の照合に使う）
   特徴量は js/features.js の desc 41 個 ＋ 官能基の数 56 個。V は McGowan の式で計算し、表の V と照合する。 */
'use strict';
const fs = require('fs'), path = require('path');
const C = require('../lib/common.cjs'), P = C.P;
const F = C.lib('features.js');
const AB = C.lib('abraham.js');
const parseCSV = C.parseCSV;

const num = s => { const v = parseFloat(s); return (s !== '' && isFinite(v)) ? v : null; };

// 連立一次方程式（対称正定値）をコレスキー分解で解く
function solveSPD(A, b) {
  const n = b.length, L = A.map(r => r.slice());
  for (let j = 0; j < n; j++) {
    for (let k = 0; k < j; k++) { const l = L[j][k]; for (let i = j; i < n; i++) L[i][j] -= L[i][k] * l; }
    const d = Math.sqrt(L[j][j]); for (let i = j; i < n; i++) L[i][j] /= d;
  }
  const y = new Array(n), x = new Array(n);
  for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i][k] * y[k]; y[i] = s / L[i][i]; }
  for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k]; x[i] = s / L[i][i]; }
  return x;
}
function fitRidge(X, y, lam) {          // X は標準化済み。切片は y の平均
  const n = X.length, p = X[0].length, b = y.reduce((a, v) => a + v, 0) / n;
  const A = Array.from({ length: p }, () => new Float64Array(p)), r = new Float64Array(p);
  for (let i = 0; i < n; i++) { const xi = X[i], yi = y[i] - b; for (let j = 0; j < p; j++) { const v = xi[j]; if (v === 0) continue; r[j] += v * yi; const Aj = A[j]; for (let k = j; k < p; k++) Aj[k] += v * xi[k]; } }
  for (let j = 0; j < p; j++) { for (let k = 0; k < j; k++) A[j][k] = A[k][j]; A[j][j] += lam; }
  return { b, w: solveSPD(A.map(a => Array.from(a)), Array.from(r)) };
}
const predict = (m, x) => { let s = m.b; for (let j = 0; j < x.length; j++) s += m.w[j] * x[j]; return s; };
function stats(err) {
  const a = err.map(Math.abs).sort((p, q) => p - q), n = a.length;
  return { n, mae: +(a.reduce((s, v) => s + v, 0) / n).toFixed(4), rmse: +Math.sqrt(err.reduce((s, v) => s + v * v, 0) / n).toFixed(4), p68: +a[Math.floor(0.68 * (n - 1))].toFixed(4), p90: +a[Math.floor(0.9 * (n - 1))].toFixed(4) };
}

C.rdkit().then(R => {
  const rows = parseCSV(fs.readFileSync(C.need(path.join(P.WORK, 'abraham_solutes_raw.csv'), '先に python train/rebuild/abraham_xls_to_csv.py'), 'utf8'));
  const H = rows[0], ix = k => H.indexOf(k);
  const by = new Map(); let nBad = 0, nSkip = 0, nUse = 0;
  for (const r of rows.slice(1)) {
    if ((r[ix('donotuse')] || '').trim()) { nSkip++; continue; }
    const smi = (r[ix('smiles')] || '').trim(), v = ['E', 'S', 'A', 'B', 'V'].map(k => num(r[ix(k)]));
    // 表は欠けた値を −123 で埋めている。取りうる範囲の外の値も外す
    const okRange = v.every(x => x !== null) && v[0] > -1.5 && v[0] < 6 && v[1] > -1 && v[1] < 6 && v[2] > -0.05 && v[2] < 4 && v[3] > -0.05 && v[3] < 5.5 && v[4] > 0 && v[4] < 10;
    if (!smi || smi.indexOf('.') >= 0 || !okRange) { nSkip++; continue; }
    let m = null; try { m = R.get_mol(smi); } catch (e) { m = null; }
    if (!m || !m.is_valid()) { if (m) m.delete(); nBad++; continue; }
    let can, f;
    try { m.remove_hs_in_place(); can = m.get_smiles(); f = F.compute(R, m); } catch (e) { m.delete(); nBad++; continue; }
    m.delete();
    const Vc = AB.mcgowanV(R, can);
    if (!by.has(can)) by.set(can, { can, name: r[ix('name')], x: f.desc.concat(f.groups), Vc, vals: [], mp: [] });
    const o = by.get(can); o.vals.push(v); const mp = num(r[ix('mp')]); if (mp !== null) o.mp.push(mp); nUse++;
  }
  const all = [...by.values()].map(o => {
    const mean = k => o.vals.reduce((s, v) => s + v[k], 0) / o.vals.length;
    return Object.assign(o, { E: mean(0), S: mean(1), A: mean(2), B: mean(3), V: mean(4) });
  });
  // 表の V が、SMILES から計算した McGowan の体積と合わない構造は、名前と構造の対応があやしいので外す（V は構造だけで決まる量）
  const data = all.filter(o => o.Vc !== null && Math.abs(o.Vc - o.V) <= 0.011);
  console.log('行', rows.length - 1, '使った行', nUse, '外した行（印つき・欠け・範囲の外）', nSkip, '読めない構造', nBad, '構造の数', all.length, '→ V が合う構造', data.length);
  {
    const d0 = all.filter(o => o.Vc !== null).map(o => o.Vc - o.V);
    console.log('V の照合（外す前）', JSON.stringify(stats(d0)), '計算できない構造', all.filter(o => o.Vc === null).length);
  }

  // V の照合（McGowan の式）
  const dv = data.filter(o => o.Vc !== null).map(o => o.Vc - o.V), big = data.filter(o => o.Vc !== null && Math.abs(o.Vc - o.V) > 0.01);
  console.log('V の照合', JSON.stringify(stats(dv)), '差が 0.01 を超える構造', big.length, '計算できない構造', data.filter(o => o.Vc === null).length);
  console.log('  差の大きい例', big.slice(0, 6).map(o => o.name + ' ' + o.V + ' / ' + o.Vc.toFixed(4)).join(' ; '));

  // 特徴量の標準化
  const p = data[0].x.length, mean = new Array(p).fill(0), sd = new Array(p).fill(0);
  for (const o of data) for (let j = 0; j < p; j++) mean[j] += o.x[j] / data.length;
  for (const o of data) for (let j = 0; j < p; j++) sd[j] += (o.x[j] - mean[j]) ** 2 / data.length;
  for (let j = 0; j < p; j++) sd[j] = Math.sqrt(sd[j]);
  const Z = data.map(o => o.x.map((v, j) => sd[j] > 0 ? (v - mean[j]) / sd[j] : 0));

  // 5 分割の交差検証（構造ごと。並びは決まった擬似乱数）
  let seed = 12345; const rnd = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648;
  const fold = data.map(() => Math.floor(rnd() * 5));
  const LAMS = [1, 3, 10, 30, 100, 300], targets = {}, cv = {};
  for (const t of ['E', 'S', 'A', 'B']) {
    const y = data.map(o => o[t]); let best = null;
    for (const lam of LAMS) {
      const err = [];
      for (let k = 0; k < 5; k++) {
        const tr = [], te = []; fold.forEach((f, i) => (f === k ? te : tr).push(i));
        const m = fitRidge(tr.map(i => Z[i]), tr.map(i => y[i]), lam);
        for (const i of te) { let v = predict(m, Z[i]); if ((t === 'A' || t === 'B') && v < 0) v = 0; err.push(v - y[i]); }
      }
      const s = stats(err); if (!best || s.rmse < best.s.rmse) best = { lam, s };
    }
    const m = fitRidge(Z, y, best.lam);
    targets[t] = { b: +m.b.toFixed(5), w: m.w.map(v => +v.toFixed(5)) };
    const sy = Math.sqrt(y.reduce((s, v) => s + (v - m.b) ** 2, 0) / y.length);
    cv[t] = Object.assign({ lambda: best.lam, sd_of_values: +sy.toFixed(4), r2: +(1 - (best.s.rmse / sy) ** 2).toFixed(3) }, best.s);
    console.log(t, JSON.stringify(cv[t]));
  }
  const out = {
    meta: { name: 'abraham_solute_linear', version: '0.1', license: 'CC BY 4.0（出典は ATTRIBUTION.md）', feature_version: F.FEATURE_VERSION, n_structures: data.length, made: new Date().toISOString().slice(0, 10),
      source: 'Compounds with known Abraham descriptors (Bradley et al.), figshare, doi:10.6084/m9.figshare.1176994, CC BY 4.0',
      note: '一次式（リッジ回帰）。cv は 5 分割の交差検証での外れ（推算 − 表の値）。A・B は 0 で切る。' },
    features: F.names(), mean: mean.map(v => +v.toFixed(5)), sd: sd.map(v => +v.toFixed(5)), targets, cv
  };
  C.writeJSON(path.join(P.OUT, 'abraham_solute_model.json'), out);
  C.writeTable(path.join(P.OUT, 'base_data', 'abraham_solutes_clean.csv.gz'), ['smiles', 'E', 'S', 'A', 'B', 'V', 'V_mcgowan', 'n_rows'],
    data.map(o => [o.can, o.E.toFixed(3), o.S.toFixed(3), o.A.toFixed(3), o.B.toFixed(3), o.V.toFixed(4), o.Vc === null ? '' : o.Vc.toFixed(4), o.vals.length]));
  console.log('書いた: ' + path.join(P.OUT, 'abraham_solute_model.json') + '、' + path.join(P.OUT, 'base_data', 'abraham_solutes_clean.csv.gz'));
}).catch(e => { console.error(e); process.exit(1); });
