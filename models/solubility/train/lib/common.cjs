// SPDX-License-Identifier: MIT
// 学習の道具で共通に使う部品：場所、CSV の読み書き（.gz も）、RDKit.js の読み込み。
// 場所はすべてこのファイルからの相対で決める（個人の環境のパスを書かない）。変えたい時は環境変数で渡す：
//   SOL_RAW   元のデータを置いたフォルダ（既定 train/raw）
//   SOL_WORK  途中の物（既定 train/work。.gitignore でリポジトリに入らない）
//   SOL_OUT   できた重み（既定 train/out。同上）
'use strict';
const fs = require('fs'), path = require('path'), zlib = require('zlib');

const TRAIN = path.resolve(__dirname, '..');
const MODEL = path.resolve(TRAIN, '..');
const P = {
  TRAIN, MODEL,
  JS: path.join(MODEL, 'js'),
  WEIGHTS: path.join(MODEL, 'weights'),
  BASE: path.join(TRAIN, 'base_data'),
  RAW: path.resolve(process.env.SOL_RAW || path.join(TRAIN, 'raw')),
  WORK: path.resolve(process.env.SOL_WORK || path.join(TRAIN, 'work')),
  OUT: path.resolve(process.env.SOL_OUT || path.join(TRAIN, 'out')),
  RDKIT: path.resolve(MODEL, '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'),
};

function parseCSV(text) {
  const rows = []; let row = [], cell = '', q = false;
  if (text.charCodeAt(0) === 0xfeff) text = text.slice(1);
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (q) { if (ch === '"') { if (text[i + 1] === '"') { cell += '"'; i++; } else q = false; } else cell += ch; }
    else if (ch === '"') q = true; else if (ch === ',') { row.push(cell); cell = ''; }
    else if (ch === '\n') { row.push(cell); rows.push(row); row = []; cell = ''; } else if (ch !== '\r') cell += ch;
  }
  if (cell !== '' || row.length) { row.push(cell); rows.push(row); }
  return rows;
}
const readText = f => { const b = fs.readFileSync(f); return (f.endsWith('.gz') ? zlib.gunzipSync(b) : b).toString('utf8'); };
// 表を読む。{H: 見出し, rows: [[...]], ix: 名前 → 列の番号, objs(): 1 行ずつの辞書}
function readTable(f) {
  const r = parseCSV(readText(f)).filter(x => !(x.length === 1 && x[0] === ''));
  const H = r[0].map(s => s.trim()), rows = r.slice(1), ix = k => H.indexOf(k);
  return { H, rows, ix, objs: () => rows.map(x => Object.fromEntries(H.map((h, i) => [h, x[i] === undefined ? '' : x[i]]))) };
}
const esc = v => { const s = v === null || v === undefined ? '' : String(v); return /[",\n\r]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; };
// 表を書く。.gz で終わる名前なら gzip で縮める（中身が同じなら毎回同じバイト列になるよう、時刻を入れない）
function writeTable(f, header, rows) {
  const text = [header.map(esc).join(',')].concat(rows.map(r => r.map(esc).join(','))).join('\n') + '\n';
  fs.mkdirSync(path.dirname(f), { recursive: true });
  fs.writeFileSync(f, f.endsWith('.gz') ? zlib.gzipSync(Buffer.from(text, 'utf8'), { level: 9 }) : text);
}
const num = s => { if (s === null || s === undefined || s === '') return null; const v = Number(s); return isFinite(v) ? v : null; };
const readJSON = f => JSON.parse(readText(f));
const writeJSON = (f, o, pretty) => { fs.mkdirSync(path.dirname(f), { recursive: true }); fs.writeFileSync(f, pretty ? JSON.stringify(o, null, 1) : JSON.stringify(o)); };
const need = (f, how) => { if (!fs.existsSync(f)) { console.error('見つからない: ' + f + (how ? '\n  ' + how : '')); process.exit(2); } return f; };
const median = a => { const b = a.slice().sort((x, y) => x - y); return b.length ? b[Math.floor((b.length - 1) / 2)] : null; };

// RDKit.js（リポジトリの vendor/rdkit）を読む
function rdkit() { return require(P.RDKIT)(); }
// 正規化した SMILES（明示の水素を外す）。読めなければ null
function canon(R, smi) {
  let m = null; try { m = R.get_mol(smi); } catch (e) { m = null; }
  if (!m || !m.is_valid()) { if (m) m.delete(); return null; }
  m.remove_hs_in_place(); const c = m.get_smiles(); m.delete(); return c;
}
const lib = name => require(path.join(P.JS, name));
// 融点の照合の表（正規化した SMILES → 融点 °C）。work/tm_lookup.json があればそれ、なければ <重みのフォルダ>/tm_lookup.json。どちらもなければ null
function findTmLookup(weightsDir) {
  for (const f of [path.join(P.WORK, 'tm_lookup.json'), path.join(weightsDir || P.WEIGHTS, 'tm_lookup.json')]) if (fs.existsSync(f)) return { file: f, mp: readJSON(f).mp };
  return null;
}

// 商品名で載っている溶媒を、化学の名前に置き換える（InChIKey で見分ける。元の表の名前は使わない）
const SOLVENT_RENAME = { 'XXJWXESWEXIICW-UHFFFAOYSA-N': 'Diethylene glycol monoethyl ether' };
// 同じ溶媒を構造で見分ける形（MixtureSolDB の溶媒の名前の置き換えに使う。構造は正規化してから比べる）
const SOLVENT_RENAME_SMILES = { 'CCOCCOCCO': 'diethylene glycol monoethyl ether' };

// 主な溶媒 16（ページと同じ id）。[id, 日本語の名前, Abraham の表での名前, 沸点の表（MixtureSolDB）での名前, 略号]
const MAIN_SOLVENTS = [
  ['water', '水', 'water', 'water', null], ['meoh', 'メタノール', 'Methanol', 'methanol', 'MeOH'], ['etoh', 'エタノール', 'Ethanol', 'ethanol', 'EtOH'],
  ['ipa', '2-プロパノール', 'Propan-2-ol', 'isopropanol', 'IPA'], ['acn', 'アセトニトリル', 'Acetonitrile', 'acetonitrile', 'MeCN'], ['acetone', 'アセトン', 'Propanone', 'acetone', null],
  ['etoac', '酢酸エチル', 'Ethyl acetate', 'ethyl acetate', 'EtOAc'], ['thf', 'テトラヒドロフラン', 'Tetrahydrofuran', 'THF', 'THF'], ['dioxane', '1,4-ジオキサン', 'Dioxane', '1,4-dioxane', null],
  ['dmso', 'ジメチルスルホキシド', 'Dimethyl sulfoxide', 'DMSO', 'DMSO'], ['dmf', 'N,N-ジメチルホルムアミド', 'N,N-Dimethylformamide', 'DMF', 'DMF'],
  ['dcm', 'ジクロロメタン', 'Dichloromethane', 'dichloromethane', 'DCM'], ['chcl3', 'クロロホルム', 'Trichloromethane', 'chloroform', 'CHCl₃'],
  ['toluene', 'トルエン', 'Toluene', 'toluene', null], ['hexane', 'ヘキサン', 'Hexane', 'n-hexane', null], ['heptane', 'ヘプタン', 'Heptane', 'n-heptane', null]];

module.exports = { P, parseCSV, readText, readTable, writeTable, num, readJSON, writeJSON, need, median, rdkit, canon, lib, esc, findTmLookup, SOLVENT_RENAME, SOLVENT_RENAME_SMILES, MAIN_SOLVENTS };
