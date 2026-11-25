// SPDX-License-Identifier: MIT
// 構造 → 標準の SMILES・InChIKey・特徴量。ページと同じ js/features.js を、同梱の RDKit.js（Node）で動かす。
// 学習の特徴量とページの特徴量が食い違わないように、特徴量の計算はこの 1 本だけにしてある。
//   node featurize.cjs tsv <in.tsv> <out.jsonl>            1 列目 id、2 列目 SMILES（見出し行なし）
//   node featurize.cjs sdf <in.sdf> <out.jsonl> <タグ名>     id は SD ファイルのそのタグの値（SMRT は PUBCHEM_COMPOUND_CID）
// 出力の 1 行: { id, ok, smiles, inchikey, d: 記述子 41（小数 4 桁）, g: 官能基 56, fp: 512 ビットの指紋（16 進）, ion: イオン化の規則 30 }
'use strict';
const fs = require('fs');
const path = require('path');
const readline = require('readline');
const MODEL = path.join(__dirname, '..', '..');
const initRDKitModule = require(path.join(MODEL, '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'));
const F = require(path.join(MODEL, 'js', 'features.js'));

function bitsToHex(bits) {
  let out = '';
  for (let i = 0; i < bits.length; i += 4) out += parseInt(bits.substr(i, 4), 2).toString(16);
  return out;
}

function one(R, id, input) {
  let m = null;
  try {
    m = R.get_mol(input);
    if (!m || !m.is_valid()) { if (m) m.delete(); return { id, ok: false }; }
    m.remove_hs_in_place();            // 明示の水素を外してから計算する（ページ側も同じ）
    const smiles = m.get_smiles();
    const inchi = m.get_inchi();
    const ik = inchi ? R.get_inchikey_for_inchi(inchi) : '';
    const f = F.compute(R, m);
    m.delete();
    return { id, ok: true, smiles, inchikey: ik, d: f.desc.map(v => +v.toFixed(4)), g: f.groups, fp: bitsToHex(f.fp), ion: f.ion };
  } catch (e) {
    if (m) try { m.delete(); } catch (_) { /* 何もしない */ }
    return { id, ok: false, err: String(e).slice(0, 80) };
  }
}

async function main() {
  const [mode, inp, outp, tag] = process.argv.slice(2);
  if (!mode || !inp || !outp) { console.error('使い方: node featurize.cjs tsv|sdf <入力> <出力.jsonl> [タグ名]'); process.exit(2); }
  const R = await initRDKitModule();
  if (R.disable_logging) R.disable_logging();
  const out = fs.createWriteStream(outp, { encoding: 'utf8' });
  const rl = readline.createInterface({ input: fs.createReadStream(inp, { encoding: 'utf8' }), crlfDelay: Infinity });
  let n = 0, bad = 0;
  const t0 = Date.now();
  const tick = () => { if (++n % 10000 === 0) console.log(n, '件', '読めない', bad, ((Date.now() - t0) / 1000).toFixed(0) + ' 秒'); };
  if (mode === 'tsv') {
    for await (const line of rl) {
      if (!line.trim()) continue;
      const p = line.split('\t');
      const r = one(R, p[0], p[1]);
      if (!r.ok) bad++;
      out.write(JSON.stringify(r) + '\n');
      tick();
    }
  } else if (mode === 'sdf') {
    let buf = [], tags = {}, cur = null;
    for await (const line of rl) {
      if (line.startsWith('$$$$')) {
        const end = buf.findIndex(l => l.startsWith('M  END'));
        const r = one(R, tags[tag] || String(n), buf.slice(0, end + 1).join('\n'));
        if (!r.ok) bad++;
        out.write(JSON.stringify(r) + '\n');
        buf = []; tags = {}; cur = null;
        tick();
        continue;
      }
      const mt = line.match(/^>\s*<([^>]+)>/);
      if (mt) { cur = mt[1]; continue; }
      if (cur !== null) { if (line.trim()) tags[cur] = line.trim(); else cur = null; }
      buf.push(line);
    }
  } else { console.error('mode は tsv か sdf'); process.exit(2); }
  await new Promise(res => out.end(res));
  console.log('構造の特徴量:', n, '件（読めない', bad + '）', ((Date.now() - t0) / 1000).toFixed(0) + ' 秒');
}
main();
