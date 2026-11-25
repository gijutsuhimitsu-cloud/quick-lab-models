// SPDX-License-Identifier: MIT
// 名前の索引（<重みのフォルダ>/name_index_free.json）から、RDKit.js で読めない構造の行を外す（何度実行してもよい）。
//   node filter_name_index.cjs <重みのフォルダ>
// 外した行は <重みのフォルダ>/../name_index_dropped.json に書く。
'use strict';
const fs = require('fs'), path = require('path');
const MODEL = path.join(__dirname, '..', '..');
(async () => {
  const dir = path.resolve(process.argv[2] || '.');
  const R = await require(path.join(MODEL, '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'))();
  const p = path.join(dir, 'name_index_free.json');
  const d = JSON.parse(fs.readFileSync(p, 'utf8')), keep = [], out = [];
  for (const row of d.names) {
    let m = null, ok = false;
    try { m = R.get_mol(row[1]); ok = !!m && (!m.is_valid || m.is_valid()); } catch (e) { ok = false; }
    if (m) { try { m.delete(); } catch (e) { /* 何もしない */ } }
    (ok ? keep : out).push(row);
  }
  d.names = keep; d.meta.n = keep.length;
  d.meta.checked = 'RDKit.js で読める構造だけを載せている。出所の表で構造が壊れていた物は直し、名前と構造の食い違いが分かっている物は外してある';
  fs.writeFileSync(p, JSON.stringify(d));
  if (out.length) fs.writeFileSync(path.join(dir, '..', 'name_index_dropped.json'), JSON.stringify(out, null, 1));
  console.log('名前の索引:', keep.length, '件（読めない構造で外した', out.length + '）');
})();
