// SPDX-License-Identifier: MIT
/* AqSolDB（水への溶解度、CC0 1.0、doi:10.7910/DVN/OVHAW8）の元の表から、学習の土台の表を作り直す。
   元のデータ: Harvard Dataverse の記録ページから curated-solubility-dataset.csv を手で取って、train/raw/aqsoldb/ に置く
              （このリポジトリは取りに行く道具を持たない）。
   使い方（models/solubility/ で）: node train/rebuild/aqsol_from_raw.cjs [AqSolDB の CSV]
   出力: <out>/base_data/aqsol_rows.csv.gz（can, logS, sd, grp, src, id, occ, inBig。同じ構造は log S を平均して 1 行。id は AqSolDB の番号）
         <out>/base_data/aqsol_solutes.csv.gz（BigSolDB に無い溶質の入力。bigsol_solutes と同じ列。融点は <weights>/tm_lookup.json から）
   溶質の名前は書かない（配布元の表の名前は写さない。構造で引く）。
   外す物：塩と混合物（「.」を含む）、電荷の合計が 0 でない物、読めない構造、McGowan の体積が出ない元素、重い原子 80 個超。
   BigSolDB の溶質の表（<base>/bigsol_solutes.csv.gz）が要る（同じ溶質を同じ番号にそろえるため）。 */
'use strict';
const fs = require('fs'), path = require('path');
const C = require('../lib/common.cjs'), P = C.P;
const F = C.lib('features.js'), AB = C.lib('abraham.js');
const BASE = path.resolve(process.env.SOL_BASE || P.BASE), WTS = path.resolve(process.env.SOL_WEIGHTS || P.WEIGHTS);
const SRC = process.argv[2] || path.join(P.RAW, 'aqsoldb', 'curated-solubility-dataset.csv');
// 角かっこの中の電荷を足す（[N+](=O)[O-] は 0）
const netCharge = smi => { let q = 0; for (const m of smi.matchAll(/\[[^\]]*?([+-]+)(\d*)\]/g)) { const s = m[1][0] === '+' ? 1 : -1; q += s * (m[2] ? +m[2] : m[1].length); } return q; };

C.rdkit().then(R => {
  const t = C.readTable(C.need(SRC, 'Harvard Dataverse（doi:10.7910/DVN/OVHAW8）から curated-solubility-dataset.csv を取って train/raw/aqsoldb/ に置く')), H = t.H, ix = t.ix;
  const iId = ix('ID'), iN = ix('Name'), iS = ix('SMILES'), iL = ix('Solubility'), iSD = ix('SD'), iG = ix('Group'), iO = ix('Ocurrences');
  const model = C.readJSON(C.need(path.join(WTS, 'abraham_solute_model.json')));
  const big = new Set(C.readTable(C.need(path.join(BASE, 'bigsol_solutes.csv.gz'))).rows.map(r => r[0]));
  const TL = C.findTmLookup(WTS), TML = TL ? TL.mp : {};
  const drop = { salt: 0, charge: 0, parse: 0, elem: 0, big: 0, value: 0 };
  const by = new Map(), newSol = new Map();
  for (const r of t.rows) {
    if (r.length < H.length) continue;
    const smi = r[iS], L = parseFloat(r[iL]);
    if (!isFinite(L)) { drop.value++; continue; }
    if (smi.indexOf('.') >= 0) { drop.salt++; continue; }
    if (netCharge(smi) !== 0) { drop.charge++; continue; }
    let m = null; try { m = R.get_mol(smi); } catch (e) { m = null; }
    if (!m || !m.is_valid()) { if (m) m.delete(); drop.parse++; continue; }
    m.remove_hs_in_place(); const can = m.get_smiles(), d = JSON.parse(m.get_descriptors());
    if (d.NumHeavyAtoms > 80) { m.delete(); drop.big++; continue; }
    const V = AB.mcgowanV(R, can); if (V === null) { m.delete(); drop.elem++; continue; }
    if (!big.has(can) && !newSol.has(can)) {
      const f = F.compute(R, m), pred = AB.predictSolute(model, f.desc.concat(f.groups));
      const tm = typeof TML[can] === 'number' ? TML[can] : null;
      newSol.set(can, { can, M: d.amw, logPow: d.CrippenClogP, TmC: tm, tmSrc: tm === null ? null : 'mp', pred: Object.assign({ V }, pred) });
    }
    m.delete();
    if (!by.has(can)) by.set(can, { can, L: [], sd: [], grp: [], src: [], id: [], occ: 0 });
    const o = by.get(can); o.L.push(L); o.sd.push(parseFloat(r[iSD]) || 0); o.grp.push(r[iG]); o.src.push(r[iId][0]); o.id.push(r[iId]); o.occ += parseInt(r[iO], 10) || 1;
  }
  const out = [...by.values()].map(o => [o.can, +(o.L.reduce((a, v) => a + v, 0) / o.L.length).toFixed(4), +Math.max(...o.sd).toFixed(3), o.grp.sort().pop(), o.src.join(''), o.id.join(' '), o.occ, big.has(o.can) ? 1 : 0]);
  C.writeTable(path.join(P.OUT, 'base_data', 'aqsol_rows.csv.gz'), ['can', 'logS', 'sd', 'grp', 'src', 'id', 'occ', 'inBig'], out);
  C.writeTable(path.join(P.OUT, 'base_data', 'aqsol_solutes.csv.gz'), ['can', 'M', 'logPow', 'TmC', 'tmSrc', 'V', 'E', 'S', 'A', 'B'],
    [...newSol.values()].map(s => [s.can, s.M, s.logPow, s.TmC === null ? '' : s.TmC, s.tmSrc || '', s.pred.V, s.pred.E, s.pred.S, s.pred.A, s.pred.B]));
  console.log('行', t.rows.length, '使う構造', out.length, '（BigSolDB にもある', out.filter(o => o[7]).length, '）新しい溶質', newSol.size, '融点つき', [...newSol.values()].filter(s => s.TmC !== null).length, '外した', JSON.stringify(drop));
}).catch(e => { console.error(e); process.exit(1); });
