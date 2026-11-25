// SPDX-License-Identifier: MIT
/* BigSolDB 2.2（CC BY 4.0）の元の表から、学習の土台の表と、式の係数（formula.json）を作り直す。
   式の見込み（js/solub.js ＋ abraham.js）と実測を比べて、外れの大きさも測る。
   元のデータ（train/raw/bigsoldb/ に置く。https://zenodo.org/records/22648301 、doi:10.5281/zenodo.22648301）:
     BigSolDBv2.2.csv、BigSolDBv2.2_gamma.csv（融点の列を使う）
   ほかに使う物: 溶媒の係数（work/solvents_abraham.json があればそれ、なければ weights/solvents.json）、
     溶質の記述子の表（<base>/abraham_solutes_clean.csv.gz）、記述子の推算式（<weights>/abraham_solute_model.json）、
     融点の照合の表（<weights>/tm_lookup.json。BigSolDB に融点がない溶質に使う。作り直しは train/rebuild/bradley_tm.cjs）
   使い方（models/solubility/ で）: node train/rebuild/bigsol_from_raw.cjs
   出力: <out>/base_data/bigsol_25C.csv.gz（25 °C 付近の実測を溶質 × 溶媒ごとに平均した表）
         <out>/base_data/bigsol_solutes.csv.gz（溶質ごとの M・log Pow・融点（tmSrc：sol = BigSolDB、mp = 融点の公開データ）・推算した記述子）
   溶質の名前は書かない（配布元の表の名前は写さない。構造で引く）。
         <out>/formula.json（式の 4 つの係数、融点が分からない時に置く値、外れ）、work/validation_bigsol.json
   比べる温度は 298.15 ± 5 K。同じ溶質 × 溶媒の点は log S を平均して 1 点にする。
   環境変数 SOL_BASE・SOL_WEIGHTS で、読む土台の表と重みのフォルダを変えられる（既定 train/base_data・weights）。 */
'use strict';
const fs = require('fs'), path = require('path');
const C = require('../lib/common.cjs'), P = C.P;
const F = C.lib('features.js'), AB = C.lib('abraham.js'), SC = C.lib('solub.js');
const BASE = path.resolve(process.env.SOL_BASE || P.BASE), WTS = path.resolve(process.env.SOL_WEIGHTS || P.WEIGHTS);
const median = C.median;
function stats(err) {
  if (!err.length) return { n: 0 };
  const a = err.map(Math.abs).sort((p, q) => p - q), n = a.length, m = err.reduce((s, v) => s + v, 0) / n;
  return { n, bias: +m.toFixed(3), mae: +(a.reduce((s, v) => s + v, 0) / n).toFixed(3), rmse: +Math.sqrt(err.reduce((s, v) => s + v * v, 0) / n).toFixed(3), p50: +a[Math.floor(0.5 * (n - 1))].toFixed(3), p68: +a[Math.floor(0.68 * (n - 1))].toFixed(3), p90: +a[Math.floor(0.9 * (n - 1))].toFixed(3) };
}
// 合わせる前の式：log S = 0.5 − log Pow + crystal + transfer（水は transfer = 0）
const chain = (s, tr, TmC) => Math.min(0.5 - s.logPow + SC.crystal(TmC, 25) + tr, SC.logSmax(s.M));

C.rdkit().then(R => {
  const canon = s => C.canon(R, s);
  // 溶媒の係数（正規化した SMILES で対応づける。同じ構造が 2 回あれば先の物）
  const svFile = fs.existsSync(path.join(P.WORK, 'solvents_abraham.json')) ? path.join(P.WORK, 'solvents_abraham.json') : path.join(WTS, 'solvents.json');
  const SV = C.readJSON(svFile).solvents.map(s => Object.assign({}, s, { name: s.name || s.en }));
  const svBy = new Map();
  for (const s of SV) { if (!s.smiles || /vol|wet|oil|Air/i.test(s.name)) continue; const c = (() => { let m = null; try { m = R.get_mol(s.smiles); } catch (e) { } if (!m || !m.is_valid()) { if (m) m.delete(); return null; } const v = m.get_smiles(); m.delete(); return v; })(); if (c && !svBy.has(c)) svBy.set(c, s); }
  // 表にある溶質の記述子
  const known = new Map(); { const t = C.readTable(C.need(path.join(BASE, 'abraham_solutes_clean.csv.gz'))); for (const r of t.rows) known.set(r[t.ix('smiles')], { E: +r[t.ix('E')], S: +r[t.ix('S')], A: +r[t.ix('A')], B: +r[t.ix('B')] }); }
  const model = C.readJSON(C.need(path.join(WTS, 'abraham_solute_model.json')));
  // 融点（BigSolDB の gamma の表。溶質ごとの中央値）
  const tmBy = new Map(); { const t = C.readTable(C.need(path.join(P.RAW, 'bigsoldb', 'BigSolDBv2.2_gamma.csv'), 'https://zenodo.org/records/22648301 から取って train/raw/bigsoldb/ に置く')); for (const r of t.rows) { const v = parseFloat(r[t.ix('melting_point_K')]); if (isFinite(v)) { const k = r[t.ix('SMILES_Solute')]; if (!tmBy.has(k)) tmBy.set(k, []); tmBy.get(k).push(v); } } }
  // 融点の公開データの照合の表（weights/tm_lookup.json）。BigSolDB に融点がない溶質に使う
  const TL = C.findTmLookup(WTS), TML = TL ? TL.mp : {};
  console.log('融点の照合の表', TL ? '使う（' + Object.keys(TML).length + ' 構造）' : 'なし（BigSolDB の融点だけ）');
  // 実測（25 °C 付近）を 溶質 × 溶媒 で平均
  const B = C.readTable(C.need(path.join(P.RAW, 'bigsoldb', 'BigSolDBv2.2.csv')));
  const iS = B.ix('SMILES_Solute'), iT = B.ix('Temperature_K'), iV = B.ix('SMILES_Solvent'), iN = B.ix('Solvent'), iL = B.ix('LogS(mol/L)'), iD = B.ix('Source'), iNm = B.ix('Compound_Name');
  const pair = new Map();
  for (const r of B.rows) {
    const T = parseFloat(r[iT]), L = parseFloat(r[iL]); if (!isFinite(T) || !isFinite(L)) continue;
    if (Math.abs(T - 298.15) > 5) continue;
    const k = r[iS] + '\t' + r[iV]; if (!pair.has(k)) pair.set(k, { solute: r[iS], solvSmi: r[iV], solvName: r[iN], L: [], doi: new Set() });
    const o = pair.get(k); o.L.push(L); o.doi.add(r[iD]);
  }
  // 溶質ごとの入力
  const solutes = new Map(); let nNoMol = 0;
  for (const o of pair.values()) {
    if (solutes.has(o.solute)) continue;
    let m = null; try { m = R.get_mol(o.solute); } catch (e) { }
    if (!m || !m.is_valid() || o.solute.indexOf('.') >= 0) { if (m) m.delete(); solutes.set(o.solute, null); nNoMol++; continue; }
    m.remove_hs_in_place(); const can = m.get_smiles(), f = F.compute(R, m), d = JSON.parse(m.get_descriptors()); m.delete();
    const V = AB.mcgowanV(R, can); if (V === null) { solutes.set(o.solute, null); nNoMol++; continue; }
    const pred = AB.predictSolute(model, f.desc.concat(f.groups)), kn = known.get(can) || null;
    const tmS = tmBy.has(o.solute), tm = tmS ? median(tmBy.get(o.solute)) - 273.15 : (typeof TML[can] === 'number' ? TML[can] : null);
    solutes.set(o.solute, { can, M: d.amw, logPow: d.CrippenClogP, TmC: tm, tmSrc: tm === null ? null : (tmS ? 'sol' : 'mp'), pred: Object.assign({ V }, pred), known: kn ? Object.assign({ V }, kn) : null });
  }
  const svCanon = new Map(); const svOf = smi => { if (!svCanon.has(smi)) { let m = null; try { m = R.get_mol(smi); } catch (e) { } let c = null; if (m && m.is_valid()) c = m.get_smiles(); if (m) m.delete(); svCanon.set(smi, c); } return svBy.get(svCanon.get(smi)) || null; };
  const tmAll = [...solutes.values()].filter(s => s && s.TmC !== null).map(s => s.TmC);
  const TM_DEFAULT = median(tmAll);
  console.log('25 °C 付近の 溶質 × 溶媒', pair.size, '溶質', solutes.size, '読めない・塩など', nNoMol, '融点つき', tmAll.length, '融点の中央値 °C', TM_DEFAULT.toFixed(0));

  const E = { water_gse: [], org_pred: [], org_noTm: [], water_noTm: [] }, rowsOut = [], perSolute = new Map();
  let nNoSolv = 0;
  for (const o of pair.values()) {
    const s = solutes.get(o.solute); if (!s) continue;
    const L = o.L.reduce((a, v) => a + v, 0) / o.L.length, sv = svOf(o.solvSmi);
    if (!sv) { nNoSolv++; continue; }
    const isW = sv.name.toLowerCase() === 'water';
    const tr = isW ? 0 : AB.logP(sv, s.pred);
    rowsOut.push([s.can, sv.name, L.toFixed(3), o.L.length, [...o.doi].slice(0, 3).join(';')]);
    if (!perSolute.has(o.solute)) perSolute.set(o.solute, []);
    perSolute.get(o.solute).push({ sv: sv.name, L, tr, isW });
    if (s.TmC !== null) { const p = chain(s, tr, s.TmC) - L; (isW ? E.water_gse : E.org_pred).push(p); }
    (isW ? E.water_noTm : E.org_noTm).push(chain(s, tr, TM_DEFAULT) - L);
  }
  const cal = { org_from_org: [], from_water: [] };
  for (const [k, list] of perSolute) {
    const s = solutes.get(k), cap = SC.logSmax(s.M);
    for (const a of list) for (const b of list) {
      if (a === b) continue;
      const e = Math.min(a.L + (b.tr - a.tr), cap) - b.L;
      if (!a.isW && !b.isW) cal.org_from_org.push(e); if (a.isW && !b.isW) cal.from_water.push(e);
    }
  }
  // ---- 式の係数を実測に合わせる（4 個だけ）： log S = k0 + kP·(−log Pow) + kC·crystal + kT·transfer ----
  const solve = (A, b) => { const n = b.length, M = A.map((r, i) => r.concat([b[i]])); for (let i = 0; i < n; i++) { let p = i; for (let r = i + 1; r < n; r++) if (Math.abs(M[r][i]) > Math.abs(M[p][i])) p = r; [M[i], M[p]] = [M[p], M[i]]; for (let r = 0; r < n; r++) if (r !== i) { const f = M[r][i] / M[i][i]; for (let c = i; c <= n; c++) M[r][c] -= f * M[i][c]; } } return M.map((r, i) => r[n] / r[i]); };
  const ols = rows => { const p = rows[0].x.length, A = Array.from({ length: p }, () => new Array(p).fill(0)), b = new Array(p).fill(0); for (const r of rows) for (let i = 0; i < p; i++) { b[i] += r.x[i] * r.y; for (let j = 0; j < p; j++) A[i][j] += r.x[i] * r.x[j]; } return solve(A, b); };
  let sd2 = 777; const rnd2 = () => (sd2 = (sd2 * 1103515245 + 12345) % 2147483648) / 2147483648;
  const foldOf = new Map(); for (const k of perSolute.keys()) foldOf.set(k, Math.floor(rnd2() * 5));
  const fitRows = [], pairRows = [];
  for (const [k, list] of perSolute) {
    const s = solutes.get(k), f = foldOf.get(k), cap = SC.logSmax(s.M);
    if (s.TmC !== null) for (const a of list) fitRows.push({ f, cap, isW: a.isW, x: [1, -s.logPow, SC.crystal(s.TmC, 25), a.tr], x0: [1, -s.logPow, SC.crystal(TM_DEFAULT, 25), a.tr], y: a.L });
    else for (const a of list) fitRows.push({ f, cap, isW: a.isW, noTm: true, x0: [1, -s.logPow, SC.crystal(TM_DEFAULT, 25), a.tr], y: a.L });
    for (const a of list) for (const b of list) if (a !== b && !a.isW && !b.isW) pairRows.push({ f, cap, dx: b.tr - a.tr, la: a.L, lb: b.L });
  }
  const dot = (K, x) => K.reduce((s, v, i) => s + v * x[i], 0);
  const eFit = { withTm: [], noTm: [], cal: [] };
  for (let k = 0; k < 5; k++) {
    const K = ols(fitRows.filter(r => r.f !== k && r.x && !r.isW));
    for (const r of fitRows) if (r.f === k && !r.isW) { if (r.x) eFit.withTm.push(Math.min(dot(K, r.x), r.cap) - r.y); eFit.noTm.push(Math.min(dot(K, r.x0), r.cap) - r.y); }
    let sxy = 0, sxx = 0; for (const r of pairRows) if (r.f !== k) { sxy += r.dx * (r.lb - r.la); sxx += r.dx * r.dx; }
    const beta = sxy / sxx; for (const r of pairRows) if (r.f === k) eFit.cal.push(Math.min(r.la + beta * r.dx, r.cap) - r.lb);
  }
  const Kall = ols(fitRows.filter(r => r.x && !r.isW)); let sxy = 0, sxx = 0; for (const r of pairRows) { sxy += r.dx * (r.lb - r.la); sxx += r.dx * r.dx; }
  const made = new Date().toISOString().slice(0, 10);
  const fit = { meta: { name: 'formula', version: '0.1', made, license: 'CC BY 4.0（出典は ATTRIBUTION.md）' },
    source: 'BigSolDB 2.2 (Krasnov et al.), Zenodo, doi:10.5281/zenodo.22648301, CC BY 4.0. 298.15 ± 5 K の実測を溶質 × 溶媒ごとに平均した ' + fitRows.filter(r => r.x && !r.isW).length + ' 点（有機溶媒・融点つき）。水は合わせ直さず、一般溶解度式のまま',
    equation: 'log10 S = k0 + kP*(-logPow) + kC*crystal + kT*transfer ; 実測 1 点で合わせる時は log S_b = log S_a + beta*(transfer_b - transfer_a)',
    T_C: 25, K: { k0: +Kall[0].toFixed(4), kP: +Kall[1].toFixed(4), kC: +Kall[2].toFixed(4), kT: +Kall[3].toFixed(4) }, beta: 1, beta_least_squares: +(sxy / sxx).toFixed(4), tm_default_C: +TM_DEFAULT.toFixed(0),
    cv: { organic_with_Tm: stats(eFit.withTm), organic_without_Tm: stats(eFit.noTm), water_with_Tm: stats(E.water_gse), water_without_Tm: stats(E.water_noTm), one_point_organic_from_organic: stats(cal.org_from_org), one_point_organic_from_water: stats(cal.from_water), one_point_beta_least_squares: stats(eFit.cal) } };
  fit.band = { organic_with_Tm: fit.cv.organic_with_Tm.p68, organic_without_Tm: fit.cv.organic_without_Tm.p68, water_with_Tm: fit.cv.water_with_Tm.p68, water_without_Tm: fit.cv.water_without_Tm.p68, one_point: fit.cv.one_point_organic_from_organic.p68, one_point_from_water: fit.cv.one_point_organic_from_water.p68 };
  C.writeJSON(path.join(P.OUT, 'formula.json'), fit, true);
  C.writeJSON(path.join(P.WORK, 'validation_bigsol.json'), { made, n_pairs_25C: pair.size, n_no_solvent_coeff: nNoSolv, tm_default_C: +TM_DEFAULT.toFixed(1), fitted: fit }, true);
  C.writeTable(path.join(P.OUT, 'base_data', 'bigsol_25C.csv.gz'), ['smiles', 'solvent', 'logS_molL', 'n_points', 'doi'], rowsOut);
  // 溶質の表（融点の出所は tmSrc：sol = BigSolDB、mp = 融点の公開データ）
  const sol = [...solutes.values()].filter(Boolean);
  C.writeTable(path.join(P.OUT, 'base_data', 'bigsol_solutes.csv.gz'), ['can', 'M', 'logPow', 'TmC', 'tmSrc', 'V', 'E', 'S', 'A', 'B', 'kE', 'kS', 'kA', 'kB'],
    sol.map(s => { const k = s.known || {}; return [s.can, s.M, s.logPow, s.TmC === null ? '' : s.TmC, s.tmSrc || '', s.pred.V, s.pred.E, s.pred.S, s.pred.A, s.pred.B, s.known ? k.E : '', s.known ? k.S : '', s.known ? k.A : '', s.known ? k.B : '']; }));
  console.log('係数', JSON.stringify(fit.K), '融点の既定', fit.tm_default_C, '\n' + Object.entries(fit.cv).map(([k, v]) => '  ' + k.padEnd(34) + ' n=' + v.n + ' p68=' + v.p68).join('\n'));
  console.log('書いた: ' + path.join(P.OUT, 'formula.json') + '、' + path.join(P.OUT, 'base_data') + ' の bigsol_25C.csv.gz・bigsol_solutes.csv.gz（溶質 ' + sol.length + '、点 ' + rowsOut.length + '）');
}).catch(e => { console.error(e); process.exit(1); });
