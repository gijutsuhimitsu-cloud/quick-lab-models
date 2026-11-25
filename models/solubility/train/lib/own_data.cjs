// SPDX-License-Identifier: MIT
// 自社のデータ（CSV）を読んで、学習の行（溶質 × 溶媒の log10 mol/L）にする。形は train/OWN_DATA.md。
// 返す物：{ rows: [{can, name, solvent（溶媒の表の name）, solventId, water, logS, n}], solutes: Map(can → {can, name, M, logPow, TmC, abr, x}), skipped: [{line, why}] }
'use strict';
const C = require('./common.cjs');
const F = C.lib('features.js'), AB = C.lib('abraham.js');

const T_MIN = 20, T_MAX = 30;            // 25 °C のモデルなので、298.15 ± 5 K に当たる行だけを使う（公開データと同じ）
// 単位 → log10(mol/L) への換算。M は分子量 g/mol
const UNITS = {
  'mg/ml': (v, M) => Math.log10(v / M), 'g/l': (v, M) => Math.log10(v / M), 'mg/l': (v, M) => Math.log10(v / 1000 / M),
  'ug/ml': (v, M) => Math.log10(v / 1000 / M), 'µg/ml': (v, M) => Math.log10(v / 1000 / M), 'μg/ml': (v, M) => Math.log10(v / 1000 / M),
  'mol/l': v => Math.log10(v), 'm': v => Math.log10(v), 'mmol/l': v => Math.log10(v / 1000), 'mm': v => Math.log10(v / 1000),
  'log10(mol/l)': v => v, 'logs': v => v, 'log10_mol/l': v => v,
};
const normUnit = u => String(u || '').trim().toLowerCase().replace(/\s+/g, '').replace('l-1', '/l').replace('ml-1', '/ml').replace('mg/ml', 'mg/ml');

// 溶媒の指定（id・英語の名前・日本語の名前・SMILES）→ 溶媒の表の 1 件
function solventResolver(R, solvents) {
  const by = new Map();
  for (const s of solvents) { [s.id, s.name, s.ja, s.short].filter(Boolean).forEach(k => { const key = String(k).toLowerCase(); if (!by.has(key)) by.set(key, s); }); }
  const byCan = new Map(); for (const s of solvents) { const c = C.canon(R, s.smiles); if (c && !byCan.has(c)) byCan.set(c, s); }
  return (name, smiles) => {
    if (name && by.has(String(name).trim().toLowerCase())) return by.get(String(name).trim().toLowerCase());
    if (smiles) { const c = C.canon(R, smiles); if (c && byCan.has(c)) return byCan.get(c); }
    if (name) { const c = C.canon(R, name); if (c && byCan.has(c)) return byCan.get(c); }
    return null;
  };
}

// 構造から、学習と同じ溶質の入力を作る。読めない・扱わない構造は理由を返す
function soluteFromSmiles(R, model, smi) {
  if (!smi) return { why: 'SMILES がない' };
  if (smi.indexOf('.') >= 0) return { why: '塩・混合物（「.」を含む）は扱わない' };
  let m = null; try { m = R.get_mol(smi); } catch (e) { m = null; }
  if (!m || !m.is_valid()) { if (m) m.delete(); return { why: 'SMILES が読めない' }; }
  m.remove_hs_in_place();
  const can = m.get_smiles(), d = JSON.parse(m.get_descriptors()), f = F.compute(R, m); m.delete();
  if (d.NumHeavyAtoms > 80) return { why: '重い原子が 80 個を超える' };
  const V = AB.mcgowanV(R, can); if (V === null) return { why: '体積を出せない元素を含む' };
  const x = f.desc.concat(f.groups);
  return { can, M: d.amw, logPow: d.CrippenClogP, x, abr: Object.assign({ V }, AB.predictSolute(model, x)) };
}

// files: CSV の場所の配列。ctx = {R, solvents, model, nameToSmiles: Map(小文字の名前 → SMILES)}
function readOwn(files, ctx) {
  const { R, solvents, model, nameToSmiles } = ctx, svOf = solventResolver(R, solvents);
  const acc = new Map(), solutes = new Map(), tms = new Map(), skipped = [];
  for (const file of files) {
    const t = C.readTable(file), H = t.H.map(h => h.toLowerCase()), col = (...ks) => { for (const k of ks) { const i = H.indexOf(k.toLowerCase()); if (i >= 0) return i; } return -1; };
    const iSm = col('solute_smiles', 'smiles'), iNm = col('solute_name', 'name'), iSv = col('solvent'), iSvS = col('solvent_smiles'), iT = col('temperature_c', 't_c', 'temperature'),
      iV = col('value', 'solubility'), iU = col('unit'), iL = col('logs_moll', 'logs'), iTm = col('tm_c', 'melting_point_c');
    if ((iSm < 0 && iNm < 0) || (iSv < 0 && iSvS < 0) || iT < 0 || (iV < 0 && iL < 0)) throw new Error(file + ': 列が足りない（solute_smiles か solute_name、solvent か solvent_smiles、temperature_C、value と unit か logS_molL）');
    t.rows.forEach((r, k) => {
      const line = file + ':' + (k + 2), get = i => (i >= 0 && r[i] !== undefined ? String(r[i]).trim() : '');
      if (r.every(c => String(c).trim() === '') || get(0).startsWith('#')) return;
      let smi = get(iSm); const name = get(iNm);
      if (!smi && name) smi = nameToSmiles.get(name.toLowerCase()) || '';
      if (!smi) return skipped.push({ line, why: '名前「' + name + '」では引けない（引けるのは例の固体の日本語の名前だけ）。SMILES を書く' });
      const sv = svOf(get(iSv), get(iSvS)); if (!sv) return skipped.push({ line, why: '溶媒「' + (get(iSv) || get(iSvS)) + '」が溶媒の表（weights/solvents.json）にない' });
      const T = parseFloat(get(iT)); if (!(T >= T_MIN && T <= T_MAX)) return skipped.push({ line, why: '温度 ' + get(iT) + ' °C は 20〜30 °C の外（このモデルは 25 °C だけ）' });
      const key = smi; let s = solutes.get(key);
      if (!s) { s = soluteFromSmiles(R, model, smi); if (s.why) return skipped.push({ line, why: s.why }); solutes.set(key, s); }
      let logS;
      if (iL >= 0 && get(iL) !== '') logS = parseFloat(get(iL));
      else { const v = parseFloat(get(iV)), u = normUnit(get(iU)), fn = UNITS[u]; if (!fn) return skipped.push({ line, why: '単位「' + get(iU) + '」は読めない（mg/mL・g/L・mg/L・µg/mL・mol/L・mmol/L・log10(mol/L)）' }); if (!(v > 0) && !u.startsWith('log')) return skipped.push({ line, why: '値が正の数でない' }); logS = fn(v, s.M); }
      if (!isFinite(logS)) return skipped.push({ line, why: '値が読めない' });
      const tm = parseFloat(get(iTm)); if (isFinite(tm)) { if (!tms.has(s.can)) tms.set(s.can, []); tms.get(s.can).push(tm); }
      const k2 = s.can + '\t' + sv.name; if (!acc.has(k2)) acc.set(k2, { can: s.can, name: name || s.can, solvent: sv.name, solventId: sv.id, water: !!sv.water, L: [] });
      acc.get(k2).L.push(logS);
    });
  }
  const byCan = new Map(); for (const s of solutes.values()) byCan.set(s.can, s);
  for (const [can, v] of tms) byCan.get(can).TmC = C.median(v);
  const rows = [...acc.values()].map(o => ({ can: o.can, name: o.name, solvent: o.solvent, solventId: o.solventId, water: o.water, logS: +(o.L.reduce((a, b) => a + b, 0) / o.L.length).toFixed(4), n: o.L.length }));
  return { rows, solutes: byCan, skipped };
}

module.exports = { readOwn, soluteFromSmiles, solventResolver, UNITS, T_MIN, T_MAX };
