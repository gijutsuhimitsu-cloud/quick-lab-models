// SPDX-License-Identifier: MIT
/* 融点の公開データ（Bradley、CC0 1.0、doi:10.6084/m9.figshare.1031637）から、
   (1) 融点の照合の表 <out>/tm_lookup.json（weights/tm_lookup.json と同じ形。正規化した SMILES → 融点 °C。同じ構造に複数の値がある時は中央値）
   (2) 融点の推算モデルの学習用の表 <out>/base_data/tm_train.csv.gz（train/base_data/ と同じ形。in_bigsol、mpC、入力 97 列）
   を作る。名前の列は写さない。
   先に: python train/rebuild/bradley_to_csv.py
   使い方（models/solubility/ で）: node train/rebuild/bradley_tm.cjs
   溶解度の学習の溶質の表（<base>/bigsol_solutes.csv.gz）も読む（in_bigsol の印を付けるため。SOL_BASE で変えられる）。
   外す行：train/rebuild/bradley_excluded_keys.txt の key の行（source の列が DrugBank・Aldrich Catalog-Handbook・Sigma-Aldrich の 1,707 行）、donotuse の印つき、「.」を含む構造、読めない構造、
           同じ構造の値の開きが 20 °C を超える物。 */
'use strict';
const fs = require('fs'), path = require('path');
const C = require('../lib/common.cjs'), P = C.P;
const F = C.lib('features.js');
const BASE = path.resolve(process.env.SOL_BASE || P.BASE);
const median = a => { const b = a.slice().sort((x, y) => x - y); return (b[Math.floor((b.length - 1) / 2)] + b[Math.ceil((b.length - 1) / 2)]) / 2; };

C.rdkit().then(R => {
  const ex = new Set(fs.readFileSync(path.join(__dirname, 'bradley_excluded_keys.txt'), 'utf8').split('\n').map(s => s.replace(/#.*/, '').trim()).filter(Boolean));
  const big = new Set(C.readTable(C.need(path.join(BASE, 'bigsol_solutes.csv.gz'))).rows.map(r => r[0]));
  const t = C.readTable(C.need(path.join(P.WORK, 'bradley_mp_raw.csv'), '先に python train/rebuild/bradley_to_csv.py')), ix = t.ix;
  const by = new Map(), exBy = new Map(); let bad = 0, skip = 0, nEx = 0;
  for (const r of t.rows) {
    if (ex.has(String(r[ix('key')]).trim())) { nEx++; const s0 = (r[ix('source')] || '').trim(); exBy.set(s0, (exBy.get(s0) || 0) + 1); continue; }
    const smi = (r[ix('smiles')] || '').trim(), mp = parseFloat(r[ix('mpC')]);
    if ((r[ix('donotuse')] || '').trim() || !smi || smi.indexOf('.') >= 0 || !isFinite(mp) || mp < -200 || mp > 500) { skip++; continue; }
    let m = null; try { m = R.get_mol(smi); } catch (e) { m = null; }
    if (!m || !m.is_valid()) { if (m) m.delete(); bad++; continue; }
    let can, x = null;
    try { m.remove_hs_in_place(); can = m.get_smiles(); if (!by.has(can)) { const f = F.compute(R, m); x = f.desc.concat(f.groups); } } catch (e) { m.delete(); bad++; continue; }
    m.delete();
    if (!by.has(can)) by.set(can, { x, mp: [] });
    by.get(can).mp.push(mp);
  }
  const look = {}, lines = [['in_bigsol', 'mpC'].concat(F.names()).join(',')]; let wide = 0, nBig = 0;
  for (const [can, o] of by) {
    if (Math.max(...o.mp) - Math.min(...o.mp) > 20) { wide++; continue; }
    const mp = median(o.mp); look[can] = +mp.toFixed(1); if (big.has(can)) nBig++;
    lines.push([big.has(can) ? 1 : 0, mp.toFixed(1)].concat(o.x).join(','));
  }
  // 外した行の内訳（元の表の source の列の名前を短くして書く）
  const comma = n => String(n).replace(/\B(?=(\d{3})+(?!\d))/g, ',');
  const short = s => /drugbank/i.test(s) ? 'DrugBank' : /catalog-handbook/i.test(s) ? 'Aldrich Catalog-Handbook' : s;
  const exText = [...exBy].sort((a, b) => b[1] - a[1]).map(([s0, n]) => short(s0) + '（' + comma(n) + ' 行）').join('、');
  C.writeJSON(path.join(P.OUT, 'tm_lookup.json'), { meta: { name: 'tm_lookup', version: '0.1', made: new Date().toISOString().slice(0, 10), license: 'CC BY 4.0（元のデータは CC0 1.0。出典は ATTRIBUTION.md）', n: Object.keys(look).length },
    source: 'Jean-Claude Bradley Open Melting Point Dataset v2 (Jean-Claude Bradley, Antony Williams, Andrew Lang), figshare, doi:10.6084/m9.figshare.1031637, CC0 1.0',
    value: '融点 / °C（同じ構造の値の中央値）', key: '正規化した SMILES（RDKit.js 2025.3.4、明示の水素を外した形）',
    excluded: '元の表の source の列が ' + exText + 'である計 ' + comma(nEx) + ' 行（train/rebuild/bradley_excluded_keys.txt）。ほかの出所（PHYSPROP、Alfa Aesar などほかの試薬会社のカタログ、ハンドブック、ウェブの資料など）の行は使った。使わない印のある行、「.」を含む構造、読めない構造、同じ構造の値の開きが 20 °C を超える物も外した。行ごとの出所（key・source・link）は元のデータにある',
    note: 'CC0 は、元のデータを登録した人が手放せる権利についての宣言で、元の資料（PHYSPROP、Alfa Aesar、ハンドブック、ウェブの資料など）の権利者が個別に許したことや、値を保証することを表すものではない。使い直す時は、必要に応じて元の表の key・source・link を見ること。出所や値の誤りは、その構造と元の表の key を添えて、このリポジトリの管理者に知らせてほしい',
    mp: look });
  fs.mkdirSync(path.join(P.OUT, 'base_data'), { recursive: true });
  fs.writeFileSync(path.join(P.OUT, 'base_data', 'tm_train.csv.gz'), require('zlib').gzipSync(Buffer.from(lines.join('\n') + '\n', 'utf8'), { level: 9 }));
  console.log('行', t.rows.length, '外した行（key の表）', nEx, '（印つき・欠けなど）', skip, '読めない構造', bad, '構造', by.size, '値の開きが 20 °C 超', wide, '残り', lines.length - 1, 'うち溶解度の学習の溶質', nBig);
  console.log('書いた: ' + path.join(P.OUT, 'tm_lookup.json') + ', ' + path.join(P.OUT, 'base_data', 'tm_train.csv.gz'));
}).catch(e => { console.error(e); process.exit(1); });
