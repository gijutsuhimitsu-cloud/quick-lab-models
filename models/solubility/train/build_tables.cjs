// SPDX-License-Identifier: MIT
/* 学習用の表を作る（1 行 = 溶質 × 溶媒、25 °C 付近の実測の平均）。入力は推論と同じ js/solfeat.js で作る。
   使い方（models/solubility/ で）:
     node train/build_tables.cjs                                   公開データだけ
     node train/build_tables.cjs --own train/own_data/my.csv       自社のデータを足す（--own は何回でも）
         --own-mode merge     公開データ ＋ 自社のデータ（既定。自社の行の重みは --own-weight、既定 1）
         --own-mode weighted  同じ。自社の行の重みの既定を 5 にする
         --own-mode only      自社のデータだけ（数百行はないと当たらない）
   入力: <base>/bigsol_25C.csv.gz・bigsol_solutes.csv.gz・aqsol_rows.csv.gz・aqsol_solutes.csv.gz（train/base_data）
         <weights>/formula.json・solvents.json・abraham_solute_model.json
         <weights>/tm_lookup.json（土台の表に融点がない溶質に、照合の表の値があれば足す。同梱の土台の表には、もう入っている）
   出力: work/train_25C.csv（有機溶媒のモデル用。BigSolDB のすべての溶媒）、work/train_25C_aq.csv（水のモデル用）、
         work/train_solutes.json（溶質ごとの入力。train_tm.py が使う）、work/own_report.json（自社のデータの読み取りの結果）
   列: group = 溶質の番号（検証の組分けに使う）、has_tm、solvent、y = log10 mol/L、src = bigsol | aqsol | own、w = 重み、
       can = 構造、tms = 融点の出所（sol | mp | own | 空）、tmc = 融点、あとは入力の列（js/solfeat.js の names の順）。 */
'use strict';
const fs = require('fs'), path = require('path');
const C = require('./lib/common.cjs'), P = C.P;
const F = C.lib('features.js'), SF = C.lib('solfeat.js');
const OWN = require('./lib/own_data.cjs');
const BASE = path.resolve(process.env.SOL_BASE || P.BASE), WTS = path.resolve(process.env.SOL_WEIGHTS || P.WEIGHTS);

const argv = process.argv.slice(2), opt = { own: [], mode: 'merge', weight: null };
for (let i = 0; i < argv.length; i++) {
  if (argv[i] === '--own') opt.own.push(path.resolve(argv[++i]));
  else if (argv[i] === '--own-mode') opt.mode = argv[++i];
  else if (argv[i] === '--own-weight') opt.weight = +argv[++i];
  else { console.error('知らない引数: ' + argv[i]); process.exit(2); }
}
if (!['merge', 'weighted', 'only'].includes(opt.mode)) { console.error('--own-mode は merge・weighted・only のどれか'); process.exit(2); }
if (opt.mode === 'only' && !opt.own.length) { console.error('--own-mode only には --own が要る'); process.exit(2); }
const OWN_W = opt.weight !== null ? opt.weight : (opt.mode === 'weighted' ? 5 : 1);
const fmt = x => x.map(v => +(+v).toFixed(5));
const solRow = o => ({ can: o.can, M: +o.M, logPow: +o.logPow, TmC: C.num(o.TmC), tmSrc: o.tmSrc || null, pred: { V: +o.V, E: +o.E, S: +o.S, A: +o.A, B: +o.B } });

C.rdkit().then(R => {
  const fit = C.readJSON(C.need(path.join(WTS, 'formula.json')));
  const solvents = C.readJSON(C.need(path.join(WTS, 'solvents.json'))).solvents;
  const SV = new Map(solvents.map(s => [s.name, Object.assign({}, s, { water: s.name.toLowerCase() === 'water' })]));
  const water = [...SV.values()].find(s => s.water);
  const model = C.readJSON(C.need(path.join(WTS, 'abraham_solute_model.json')));
  const TL = C.findTmLookup(WTS), TML = TL ? TL.mp : null;
  // 溶質：BigSolDB の溶質（番号 0〜）、AqSolDB だけの溶質（その後ろ）
  const solutes = C.readTable(C.need(path.join(BASE, 'bigsol_solutes.csv.gz'))).objs().map(solRow);
  const extra = C.readTable(C.need(path.join(BASE, 'aqsol_solutes.csv.gz'))).objs().map(solRow);
  let nFill = 0;
  const fillTm = s => { if (s.TmC === null && TML && typeof TML[s.can] === 'number') { s.TmC = TML[s.can]; s.tmSrc = 'mp'; nFill++; } };
  solutes.forEach(fillTm); extra.forEach(fillTm);
  console.log('融点つきの溶質 ' + solutes.concat(extra).filter(s => s.TmC !== null).length + '（うち照合の表で足した ' + nFill + '）');
  const all = solutes.concat(extra), idx = new Map();
  const feat = s => { const m = R.get_mol(s.can); m.remove_hs_in_place(); const f = F.compute(R, m); m.delete(); s.x = f.desc.concat(f.groups); s.abr = s.pred; };
  all.forEach((s, i) => { feat(s); idx.set(s.can, i); });

  // 自社のデータ
  let own = { rows: [], solutes: new Map(), skipped: [] };
  if (opt.own.length) {
    const nameToSmiles = new Map();
    if (fs.existsSync(path.join(WTS, 'example_compounds.json'))) for (const e of C.readJSON(path.join(WTS, 'example_compounds.json')).compounds) nameToSmiles.set(e.name.toLowerCase(), e.smiles);       // 名前で引けるのは、例の固体の日本語の名前だけ
    own = OWN.readOwn(opt.own, { R, solvents, model, nameToSmiles });
    for (const [can, s] of own.solutes) {
      if (idx.has(can)) { const t = all[idx.get(can)]; if (typeof s.TmC === 'number') { t.TmC = s.TmC; t.tmSrc = 'own'; } continue; }
      const tmO = typeof s.TmC === 'number' ? s.TmC : (TML && typeof TML[can] === 'number' ? TML[can] : null);      // 自社の融点 ＞ 融点の照合の表
      idx.set(can, all.length); all.push({ can, M: s.M, logPow: s.logPow, TmC: tmO, tmSrc: tmO === null ? null : (typeof s.TmC === 'number' ? 'own' : 'mp'), pred: s.abr, abr: s.abr, x: s.x });
    }
    console.log('自社のデータ', own.rows.length, '組（溶質 × 溶媒）、溶質', own.solutes.size, '（公開データにもある', [...own.solutes.keys()].filter(c => idx.get(c) < solutes.length + extra.length).length, '）、読まなかった行', own.skipped.length);
    own.skipped.slice(0, 12).forEach(s => console.log('  読まなかった', s.line, s.why));
  }

  const names = SF.names(F.names()), head = ['group', 'has_tm', 'solvent', 'y', 'src', 'w', 'can', 'tms', 'tmc'], headAq = ['group', 'has_tm', 'solvent', 'y', 'src', 'aqg', 'aqsd', 'w', 'can', 'tms', 'tmc'];
  const meta = (i, s, sv, y, src, w) => [i, s.TmC === null ? 0 : 1, sv.name, y, src, w, s.can, s.tmSrc || '', s.TmC === null ? '' : s.TmC];
  const rows = [], rowsAq = [];
  const usePublic = opt.mode !== 'only';
  if (usePublic) {
    for (const r of C.readTable(C.need(path.join(BASE, 'bigsol_25C.csv.gz'))).objs()) {
      const i = idx.get(r.smiles), sv = SV.get(r.solvent); if (i === undefined || !sv || i >= solutes.length) continue;
      const s = all[i], x = fmt(SF.build(s, sv, fit.K, fit.tm_default_C, 25));
      rows.push(meta(i, s, sv, r.logS_molL, 'bigsol', 1).concat(x));
      const m = meta(i, s, sv, r.logS_molL, 'bigsol', 1); rowsAq.push(m.slice(0, 5).concat(['', 0]).concat(m.slice(5)).concat(x));
    }
    for (const o of C.readTable(C.need(path.join(BASE, 'aqsol_rows.csv.gz'))).objs()) {
      const i = idx.get(o.can); if (i === undefined) continue;
      const s = all[i], x = fmt(SF.build(s, water, fit.K, fit.tm_default_C, 25));
      const m = meta(i, s, water, +o.logS, 'aqsol', 1); rowsAq.push(m.slice(0, 5).concat([o.grp, +o.sd]).concat(m.slice(5)).concat(x));
    }
  }
  for (const o of own.rows) {
    const i = idx.get(o.can), s = all[i], sv = SV.get(o.solvent), x = fmt(SF.build(s, sv, fit.K, fit.tm_default_C, 25));
    rows.push(meta(i, s, sv, o.logS, 'own', OWN_W).concat(x));
    if (sv.water) { const m = meta(i, s, sv, o.logS, 'own', OWN_W); rowsAq.push(m.slice(0, 5).concat(['', 0]).concat(m.slice(5)).concat(x)); }
  }
  fs.mkdirSync(P.WORK, { recursive: true });
  C.writeTable(path.join(P.WORK, 'train_25C.csv'), head.concat(names), rows);
  C.writeTable(path.join(P.WORK, 'train_25C_aq.csv'), headAq.concat(names), rowsAq);
  C.writeJSON(path.join(P.WORK, 'train_solutes.json'), all.slice(0, solutes.length).map(s => ({ can: s.can, x: s.x })));
  C.writeJSON(path.join(P.WORK, 'own_report.json'), { files: opt.own, mode: opt.mode, weight: OWN_W, rows: own.rows, skipped: own.skipped }, true);
  console.log('train_25C.csv', rows.length, '行（自社', rows.filter(r => r[4] === 'own').length, '）、train_25C_aq.csv', rowsAq.length, '行、入力の列', names.length, '、溶質', all.length);
}).catch(e => { console.error(e); process.exit(1); });
