// SPDX-License-Identifier: MIT
/* MixtureSolDB（2 成分の混合溶媒での溶解度、CC BY 4.0、doi:10.5281/zenodo.18660057）の元の表から、
   混合溶媒の補正を合わせるための表（25 °C 付近の系列の途中の点と、対数線形からのずれ）を作り直す。
   元のデータ（train/raw/ に置く）:
     mixturesoldb/MixtureSolDB.csv（https://zenodo.org/records/18660057）
     bigsoldb/BigSolDBv2.2_densities.csv（純粋な溶媒の密度。体積分率への換算に使う）
   使い方（models/solubility/ で）: node train/rebuild/mix_from_raw.cjs
   出力: <out>/base_data/mix_series_25C.csv.gz（solute, pair, main, cls, phi, r, b0, b1）。続けて node train/fit_mix.cjs
   選び方：同じ論文・同じ溶質・同じ溶媒の組・同じ温度（298.15 ± 5 K）で、純粋な溶媒 2 つの値がそろっている系列だけ。
   濃度は g L⁻¹（溶液の密度は、溶媒の体積が足し算になると置いた混合溶媒の密度で近似）。体積分率は、純粋な溶媒の密度と分子量から換算。 */
'use strict';
const path = require('path');
const C = require('../lib/common.cjs'), P = C.P;
const WTS = path.resolve(process.env.SOL_WEIGHTS || P.WEIGHTS);

C.rdkit().then(R => {
  const info = new Map();          // 溶媒の SMILES → {can, M}
  const sv = smi => { if (!info.has(smi)) { let m = null; try { m = R.get_mol(smi); } catch (e) { } if (!m || !m.is_valid()) { if (m) m.delete(); info.set(smi, null); } else { const d = JSON.parse(m.get_descriptors()); info.set(smi, { can: m.get_smiles(), M: d.amw }); m.delete(); } } return info.get(smi); };
  const dens = new Map();
  for (const r of C.readTable(C.need(path.join(P.RAW, 'bigsoldb', 'BigSolDBv2.2_densities.csv'))).rows) { const k = r[0].toLowerCase(), T = +r[1], v = +r[2]; if (!isFinite(v)) continue; const o = dens.get(k); if (!o || Math.abs(T - 298.15) < Math.abs(o.T - 298.15)) dens.set(k, { T, v }); }
  // 主な溶媒 16 の id（正規化した SMILES で対応づける）
  const fs = require('fs'), svFile = fs.existsSync(path.join(P.WORK, 'solvents_abraham.json')) ? path.join(P.WORK, 'solvents_abraham.json') : path.join(WTS, 'solvents.json');
  const SV = C.readJSON(C.need(svFile)).solvents;
  const idOf = new Map(C.MAIN_SOLVENTS.map(([id, , en]) => [sv(SV.find(s => s.name === en).smiles).can, id]));
  // 商品名で載っている溶媒は、構造で見分けて化学の名前で書く（並びの向きは元の名前で決めるので、合わせた係数は変わらない）
  const ren = new Map(Object.entries(C.SOLVENT_RENAME_SMILES).map(([smi, nm]) => [sv(smi).can, nm]));
  const t = C.readTable(C.need(path.join(P.RAW, 'mixturesoldb', 'MixtureSolDB.csv'), 'https://zenodo.org/records/18660057 から取って train/raw/mixturesoldb/ に置く')), ix = t.ix;
  const iS = ix('SMILES_Solute'), iT = ix('Temperature_K'), iG = ix('Solubility(g/100g)'), i1 = ix('Solvent1'), i2 = ix('Solvent2'), is1 = ix('SMILES_Solvent1'), is2 = ix('SMILES_Solvent2'), if1 = ix('Fraction_Solvent1'), iFT = ix('Fraction_Type'), iD = ix('Source');
  const groups = new Map(); let nNoDens = 0, nT = 0;
  for (const r of t.rows) {
    const T = +r[iT]; if (!(Math.abs(T - 298.15) <= 5)) continue; nT++;
    const a = sv(r[is1]), b = sv(r[is2]), da = dens.get(r[i1].toLowerCase()), db = dens.get(r[i2].toLowerCase()), S = +r[iG], f1 = +r[if1];
    if (!a || !b || !da || !db) { nNoDens++; continue; }
    if (!(S > 0) || !(f1 >= 0 && f1 <= 1)) continue;
    let wa, wb; const ft = r[iFT].toLowerCase();
    if (ft.indexOf('mol') >= 0) { wa = f1 * a.M; wb = (1 - f1) * b.M; } else if (ft.indexOf('vol') >= 0) { wa = f1 * da.v; wb = (1 - f1) * db.v; } else { wa = f1; wb = 1 - f1; }
    const va = wa / da.v, vb = wb / db.v, phiA = va / (va + vb), rho = (wa + wb) / (va + vb);
    const logC = Math.log10(S / (100 + S) * rho * 1000);          // g L⁻¹
    const ka = idOf.get(a.can) || r[i1].toLowerCase(), kb = idOf.get(b.can) || r[i2].toLowerCase(), swap = ka > kb;
    const da_ = ren.get(a.can) || ka, db_ = ren.get(b.can) || kb;
    const pair = swap ? db_ + '|' + da_ : da_ + '|' + db_, phi = swap ? phiA : 1 - phiA;
    const key = [r[iD], r[iS], pair, T.toFixed(2)].join('\t');
    if (!groups.has(key)) groups.set(key, { solute: r[iS], pair, main: idOf.has(a.can) && idOf.has(b.can), water: ka === 'water' || kb === 'water', pts: [] });
    groups.get(key).pts.push({ phi, logC });
  }
  const rows = [];
  for (const g of groups.values()) {
    const e0 = g.pts.filter(p => p.phi < 1e-6), e1 = g.pts.filter(p => p.phi > 1 - 1e-6); if (!e0.length || !e1.length) continue;
    const c0 = e0.reduce((s, p) => s + p.logC, 0) / e0.length, c1 = e1.reduce((s, p) => s + p.logC, 0) / e1.length;
    for (const p of g.pts) if (p.phi > 0.02 && p.phi < 0.98) rows.push([g.solute, g.pair, g.main ? 1 : 0, g.water ? 'wo' : 'oo', p.phi, p.logC - ((1 - p.phi) * c0 + p.phi * c1), p.phi * (1 - p.phi), p.phi * (1 - p.phi) * (1 - 2 * p.phi)]);
  }
  C.writeTable(path.join(P.OUT, 'base_data', 'mix_series_25C.csv.gz'), ['solute', 'pair', 'main', 'cls', 'phi', 'r', 'b0', 'b1'], rows);
  console.log('298 K 付近の行', nT, '密度がない溶媒を含む行', nNoDens, '系列', groups.size, '途中の点', rows.length, '→', path.join(P.OUT, 'base_data', 'mix_series_25C.csv.gz'));
}).catch(e => { console.error(e); process.exit(1); });
