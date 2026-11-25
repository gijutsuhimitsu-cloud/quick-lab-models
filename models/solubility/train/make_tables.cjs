// SPDX-License-Identifier: MIT
/* 推論のプログラムが読む小さな表を作る：溶媒の表、実測の照合の表、例の固体。
   入力: <base>/bigsol_25C.csv.gz・bigsol_solutes.csv.gz・aqsol_rows.csv.gz・aqsol_solutes.csv.gz（学習の土台の表）
         溶媒の係数（work/solvents_abraham.json があればそれ、なければ weights/solvents.json）
         沸点（train/raw/mixturesoldb/solvent_boiling_point.csv があればそれ、なければ weights/solvents.json の値）
         <weights>/abraham_solute_model.json
   使い方（models/solubility/ で）: node train/make_tables.cjs
   出力: <out>/solvents.json、<out>/lookup_25C.json、<out>/example_compounds.json
   融点は土台の表の値（BigSolDB、なければ融点の公開データ）。出所は s・tmSrc（sol = BigSolDB、mp = 融点の公開データ）。
   溶質の名前は配布元の表から写さない（照合の表は構造で引く。例の固体の日本語の名前は、こちらで付けた物）。 */
'use strict';
const fs = require('fs'), path = require('path');
const C = require('./lib/common.cjs'), P = C.P;
const F = C.lib('features.js');
const BASE = path.resolve(process.env.SOL_BASE || P.BASE), WTS = path.resolve(process.env.SOL_WEIGHTS || P.WEIGHTS);
const r3 = v => +(+v).toFixed(3);
// 例の固体：日本語の名前（こちらで付けた）, 構造（= の後ろが SMILES）。融点があり、主な溶媒の実測が 4 つ以上の物だけを使う
const WANT = [
  ['カフェイン', '=Cn1cnc2c1c(=O)n(C)c(=O)n2C'], ['アセトアミノフェン', '=CC(=O)Nc1ccc(O)cc1'], ['テオフィリン', '=Cn1c(=O)c2[nH]cnc2n(C)c1=O'], ['アセトアニリド', '=CC(=O)Nc1ccccc1'],
  ['メチルパラベン', '=COC(=O)c1ccc(O)cc1'], ['エチルパラベン', '=CCOC(=O)c1ccc(O)cc1'], ['プロピルパラベン', '=CCCOC(=O)c1ccc(O)cc1'], ['ブチルパラベン', '=CCCCOC(=O)c1ccc(O)cc1'],
  ['ナフタレン', '=c1ccc2ccccc2c1'], ['ビフェニル', '=c1ccc(-c2ccccc2)cc1'], ['アントラセン', '=c1ccc2cc3ccccc3cc2c1'], ['ベンゾフェノン', '=O=C(c1ccccc1)c1ccccc1'], ['レゾルシノール', '=Oc1cccc(O)c1'],
  ['安息香酸', '=O=C(O)c1ccccc1'], ['サリチル酸', '=O=C(O)c1ccccc1O'], ['イブプロフェン', '=CC(C)Cc1ccc(cc1)C(C)C(=O)O'], ['アスピリン', '=CC(=O)Oc1ccccc1C(=O)O'], ['バニリン', '=COc1cc(C=O)ccc1O'],
  ['スルファニルアミド', '=Nc1ccc(S(N)(=O)=O)cc1'], ['インドメタシン', '=COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1'], ['フルタミド', '=CC(C)C(=O)Nc1ccc([N+](=O)[O-])c(C(F)(F)F)c1'], ['ナプロキセン', '=COc1ccc2cc([C@H](C)C(=O)O)ccc2c1'], ['クルクミン', '=COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O'], ['コウジ酸', '=O=c1cc(CO)occ1O'],
  ['バニリン酸', '=COc1cc(C(=O)O)ccc1O'], ['5-フルオロウラシル', '=O=c1[nH]cc(F)c(=O)[nH]1'], ['ベンズアニリド', '=O=C(Nc1ccccc1)c1ccccc1'], ['プロゲステロン', '=CC(=O)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C'], ['フルオレン', '=c1ccc2c(c1)Cc1ccccc1-2'], ['アルテミシニン', '=C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]32OO4'],
  ['カルベジロール', '=COc1ccccc1OCCNCC(O)COc1cccc2[nH]c3ccccc3c12'], ['グリセオフルビン', '=COC1=CC(=O)C[C@@H](C)[C@]12Oc1c(Cl)c(OC)cc(OC)c1C2=O'], ['ニメスリド', '=CS(=O)(=O)Nc1ccc([N+](=O)[O-])cc1Oc1ccccc1'], ['ベンズイミダゾール', '=c1ccc2[nH]cnc2c1'], ['エテンザミド', '=CCOc1ccccc1C(N)=O'], ['クロザピン', '=CN1CCN(C2=Nc3cc(Cl)ccc3Nc3ccccc32)CC1'],
  ['プラジカンテル', '=O=C(C1CCCCC1)N1CC(=O)N2CCc3ccccc3C2C1'], ['エダラボン', '=CC1=NN(c2ccccc2)C(=O)C1'], ['グリクラジド', '=Cc1ccc(S(=O)(=O)NC(=O)NN2CC3CCCC3C2)cc1'], ['イソニコチンアミド', '=NC(=O)c1ccncc1'], ['1-ニトロナフタレン', '=O=[N+]([O-])c1cccc2ccccc12'],
  ['イトラコナゾール', '=CCC(C)n1ncn(-c2ccc(N3CCN(c4ccc(OCC5COC(Cn6cncn6)(c6ccc(Cl)cc6Cl)O5)cc4)CC3)cc2)c1=O'], ['オクスカルバゼピン', '=NC(=O)N1c2ccccc2CC(=O)c2ccccc21']];
const N_EXAMPLES = 40, AQ_DOI = '10.7910/DVN/OVHAW8';
const slug = s => s.toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_|_$/g, '');
const solRow = o => ({ can: o.can, M: +o.M, logPow: +o.logPow, TmC: C.num(o.TmC), tmSrc: o.tmSrc || null, abr: { V: +o.V, E: +o.E, S: +o.S, A: +o.A, B: +o.B } });

C.rdkit().then(R => {
  const canon = s => C.canon(R, s);
  // ---- 溶媒の表 ----
  const old = fs.existsSync(path.join(WTS, 'solvents.json')) ? C.readJSON(path.join(WTS, 'solvents.json')).solvents : [];
  const svFile = fs.existsSync(path.join(P.WORK, 'solvents_abraham.json')) ? path.join(P.WORK, 'solvents_abraham.json') : path.join(WTS, 'solvents.json');
  const SV = C.readJSON(C.need(svFile)).solvents;
  const bpFile = path.join(P.RAW, 'mixturesoldb', 'solvent_boiling_point.csv');
  const bp = fs.existsSync(bpFile) ? new Map(C.readTable(bpFile).rows.map(r => [r[0].toLowerCase(), parseFloat(r[1])])) : null;
  const B = C.readTable(C.need(path.join(BASE, 'bigsol_25C.csv.gz'))).objs();
  const nRows = new Map(), nSol = new Map();
  for (const r of B) { nRows.set(r.solvent, (nRows.get(r.solvent) || 0) + 1); }
  { const seen = new Map(); for (const r of B) { if (!seen.has(r.solvent)) seen.set(r.solvent, new Set()); seen.get(r.solvent).add(r.smiles); } for (const [k, v] of seen) nSol.set(k, v.size); }
  const byCan = new Map();
  for (const s of SV) { if (!s.smiles || /vol|wet|oil|Air/i.test(s.name)) continue; let m = null; try { m = R.get_mol(s.smiles); } catch (e) { } if (!m || !m.is_valid()) { if (m) m.delete(); continue; } const c = m.get_smiles(); m.delete(); if (!byCan.has(c)) byCan.set(c, s); }
  const used = [...byCan.values()], mainBy = new Map(C.MAIN_SOLVENTS.map(m => [m[2], m])), ids = new Set();
  const rec = s => {
    const mm = mainBy.get(s.name), id = mm ? mm[0] : slug(s.name); if (ids.has(id)) throw new Error('id が重なる: ' + id); ids.add(id);
    let bpC = null;
    if (mm && bp) { const b = bp.get(mm[3].toLowerCase()); bpC = isFinite(b) ? +(b - 273.15).toFixed(0) : null; }
    else if (mm) { const o = old.find(x => x.id === id); bpC = o ? o.bpC : null; }
    return { id, name: s.name, ja: mm ? mm[1] : null, short: mm ? mm[4] : null, main: !!mm, water: s.name.toLowerCase() === 'water', bpC, smiles: s.smiles, inchikey: s.inchikey || null,
      c: s.c, e: s.e, s: s.s, a: s.a, b: s.b, v: s.v, ref: s.ref || null, year: s.year || null, n_rows: nRows.get(s.name) || 0, n_solutes: nSol.get(s.name) || 0 };
  };
  const main = C.MAIN_SOLVENTS.map(m => { const s = used.find(x => x.name === m[2]); if (!s) throw new Error('主な溶媒の係数がない: ' + m[2]); return rec(s); });
  const others = used.filter(s => !mainBy.has(s.name)).map(rec);
  const solvents = main.concat(others), idByName = new Map(solvents.map(s => [s.name, s.id]));
  C.writeJSON(path.join(P.OUT, 'solvents.json'), { meta: { name: 'solvents', version: '0.1', made: new Date().toISOString().slice(0, 10), license: 'CC BY 4.0（出典は ATTRIBUTION.md）',
    source: 'Abraham の溶媒の係数：Abraham Solvent Master Dataset (AbraLlama), Lang & Lee, figshare, doi:10.6084/m9.figshare.25898476, CC0 1.0。沸点（主な溶媒だけ）：MixtureSolDB の solvent_boiling_point.csv, doi:10.5281/zenodo.18660057, CC BY 4.0。n_rows・n_solutes：BigSolDB 2.2 の 25 °C 付近の学習の点の数',
    equation: 'log10 P(water -> solvent) = c + e*E + s*S + a*A + b*B + v*V', note: 'main は図や例で使う 16 溶媒。n_rows が 0 の溶媒は学習に入っていない（見込みは外挿）。同じ構造が表に 2 回ある時は、学習と同じく先の物だけを残した' }, solvents }, true);
  console.log('溶媒', solvents.length, '（主な 16、学習の点がある', solvents.filter(s => s.n_rows > 0).length, '）');

  // ---- 実測の照合の表（BigSolDB の 25 °C 付近の溶媒ごとの平均 ＋ AqSolDB の水） ----
  const meas = new Map();
  for (const r of B) { const id = idByName.get(r.solvent); if (!id) continue; if (!meas.has(r.smiles)) meas.set(r.smiles, {}); meas.get(r.smiles)[id] = [r3(r.logS_molL), +r.n_points, r.doi.split(';')[0]]; }
  const big = C.readTable(C.need(path.join(BASE, 'bigsol_solutes.csv.gz'))).objs().map(solRow);
  const aqRows = C.readTable(C.need(path.join(BASE, 'aqsol_rows.csv.gz'))).objs();
  const aqSol = new Map(C.readTable(C.need(path.join(BASE, 'aqsol_solutes.csv.gz'))).objs().map(o => [o.can, solRow(o)]));
  // 例の固体を選ぶ（BigSolDB の実測の数で選ぶ。AqSolDB の水を足す前に）
  const byC = new Map(big.map(s => [s.can, s]));
  const mainIds = new Set(main.map(s => s.id)), hits = [];
  for (const [jp, key] of WANT) {
    const s = byC.get(canon(key.slice(1)));
    if (!s) { console.log('  なし', jp); continue; }
    const m = Object.fromEntries(Object.entries(meas.get(s.can) || {}).filter(([k]) => mainIds.has(k))), n = Object.keys(m).length;
    if (s.TmC === null || n < 4) { console.log('  外す', jp, '融点', s.TmC, '主な溶媒の実測', n); continue; }
    const mol = R.get_mol(s.can); mol.remove_hs_in_place(); const f = F.compute(R, mol); mol.delete();
    hits.push({ id: 'p' + hits.length, name: jp, smiles: s.can, M: +s.M.toFixed(2), logPow: r3(s.logPow), TmC: +s.TmC.toFixed(0), tmSrc: s.tmSrc,
      abr: { E: r3(s.abr.E), S: r3(s.abr.S), A: r3(s.abr.A), B: r3(s.abr.B), V: +(+s.abr.V).toFixed(4) }, x: f.desc.concat(f.groups), meas: m });
  }
  let nAdd = 0; const aqNew = [];
  for (const o of aqRows) {
    if (!meas.has(o.can)) meas.set(o.can, {});
    const m = meas.get(o.can); if (m.water) continue;
    m.water = [r3(o.logS), +o.occ || 1, AQ_DOI]; nAdd++;
    if (aqSol.has(o.can)) aqNew.push(aqSol.get(o.can));
  }
  const examples = hits.slice(0, N_EXAMPLES).map((h, i) => Object.assign({}, h, { id: 'p' + i }));
  for (const e of examples) { const w = meas.get(e.smiles); if (w && w.water && !e.meas.water) e.meas.water = w.water; }
  C.writeJSON(path.join(P.OUT, 'example_compounds.json'), { meta: { name: 'example_compounds', version: '0.1', made: new Date().toISOString().slice(0, 10), license: 'CC BY 4.0（出典は ATTRIBUTION.md）',
    note: '例の固体（BigSolDB 2.2 に 25 °C 付近の実測が多い物）。x は js/features.js の入力 97 個、abr は推算した Abraham 記述子、meas は溶媒 id → [log10(mol/L), 点の数, 出典の DOI]。name は日本語の名前（こちらで付けた）。TmC は融点 °C、tmSrc はその出所（sol = BigSolDB、mp = 融点の公開データ）。どれも学習に使った溶質なので、見込みは実測に近く出る' }, compounds: examples });
  const look = {};
  for (const s of big) { const m = meas.get(s.can); if (m && Object.keys(m).length) look[s.can] = { t: s.TmC === null ? null : +s.TmC.toFixed(0), s: s.TmC === null ? null : s.tmSrc, m }; }
  for (const s of aqNew) look[s.can] = { t: s.TmC === null ? null : +s.TmC.toFixed(0), s: s.TmC === null ? null : s.tmSrc, m: meas.get(s.can) };
  C.writeJSON(path.join(P.OUT, 'lookup_25C.json'), { meta: { name: 'lookup_25C', version: '0.1', made: new Date().toISOString().slice(0, 10), license: 'CC BY 4.0（出典は ATTRIBUTION.md）' },
    source: 'BigSolDB 2.2, doi:10.5281/zenodo.22648301, CC BY 4.0（25 °C 付近の溶媒ごとの実測）／AqSolDB, doi:10.7910/DVN/OVHAW8, CC0 1.0（水。BigSolDB に水の値がない溶質だけ）',
    T: '298.15 ± 5 K（AqSolDB は室温付近の値）', key: '正規化した SMILES（RDKit.js 2025.3.4、明示の水素を外した形）', value: 't = 融点 °C（なければ null）、s = 融点の出所（sol = BigSolDB、mp = 融点の公開データ doi:10.6084/m9.figshare.1031637、CC0 1.0）、m = 溶媒 id → [log10(mol/L), 点の数（AqSolDB は元の出所の数）, 出典の DOI]', solutes: look });
  console.log('例の固体', examples.length, '（融点つき', examples.filter(e => e.TmC !== null).length, '）照合の表の溶質', Object.keys(look).length, 'AqSolDB で水を補った溶質', nAdd);
  console.log('書いた: ' + ['solvents.json', 'lookup_25C.json', 'example_compounds.json'].map(f => path.join(P.OUT, f)).join('、'));
}).catch(e => { console.error(e); process.exit(1); });
