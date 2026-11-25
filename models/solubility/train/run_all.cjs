// SPDX-License-Identifier: MIT
/* 学習し直しの流れを 1 回で走らせる（Windows・macOS・Linux 共通。シェルの書き方に寄らない）。
   使い方（models/solubility/ で）:
     node train/run_all.cjs                                         公開データだけで学習し直す（公開の重みと同じ物ができる。融点の表を手元で作った時）
     node train/run_all.cjs --own train/own_data/my.csv             自社のデータを足す（--own-mode・--own-weight は train/build_tables.cjs と同じ）
     node train/run_all.cjs --quick                                 木を減らして流れだけ確かめる（数十秒）
     node train/run_all.cjs --install                               最後に train/out の重みを weights/ へ写す（照合に合格した時だけ）
   Python と Node の場所は環境変数 PYTHON・NODE（既定は python・node）。
   流れ: build_tables.cjs → train_gbt.py → train_tm.py → fit_mix.cjs → tests/check.cjs train/out */
'use strict';
const fs = require('fs'), path = require('path'), { spawnSync } = require('child_process');
const C = require('./lib/common.cjs'), P = C.P;
const PY = process.env.PYTHON || 'python', NODE = process.env.NODE || process.execPath || 'node';
const argv = process.argv.slice(2), pass = [], flags = new Set();
for (let i = 0; i < argv.length; i++) {
  if (['--own', '--own-mode', '--own-weight'].includes(argv[i])) { pass.push(argv[i], argv[++i]); }
  else if (['--quick', '--install', '--skip-tm'].includes(argv[i])) flags.add(argv[i]);
  else { console.error('知らない引数: ' + argv[i]); process.exit(2); }
}
const run = (cmd, args) => {
  console.log('\n$ ' + [path.basename(cmd)].concat(args.map(a => path.isAbsolute(a) ? path.relative(P.MODEL, a) : a)).join(' '));
  const r = spawnSync(cmd, args, { stdio: 'inherit', cwd: P.MODEL, env: Object.assign({}, process.env, { PYTHONIOENCODING: 'utf-8' }) });
  if (r.error) { console.error('起動できなかった: ' + cmd + '（環境変数 PYTHON・NODE で場所を渡す）\n' + r.error.message); process.exit(1); }
  if (r.status !== 0) { console.error('止まった（終了コード ' + r.status + '）'); process.exit(r.status || 1); }
};
const t0 = Date.now();
run(NODE, [path.join(P.TRAIN, 'build_tables.cjs')].concat(pass));
run(PY, [path.join(P.TRAIN, 'train_gbt.py'), '--work', P.WORK, '--out', P.OUT].concat(flags.has('--quick') ? ['--quick'] : []));
const BASE = path.resolve(process.env.SOL_BASE || P.BASE);
if (!flags.has('--skip-tm') && fs.existsSync(path.join(BASE, 'tm_train.csv.gz'))) run(PY, [path.join(P.TRAIN, 'train_tm.py'), '--work', P.WORK, '--out', P.OUT, '--base', BASE].concat(flags.has('--quick') ? ['--quick'] : []));
else console.log('\n（融点のモデルは学習し直さない）');
run(NODE, [path.join(P.TRAIN, 'fit_mix.cjs')]);
run(NODE, [path.join(P.MODEL, 'tests', 'check.cjs'), P.OUT]);
if (flags.has('--install')) {
  // 今回の実行で書いた物だけを写す（前の実行の残り、たとえば --skip-tm の時の古い融点のモデルは写さない）
  const fresh = f => { const p = path.join(P.OUT, f); return fs.existsSync(p) && fs.statSync(p).mtimeMs >= t0 - 2000; };
  const files = ['solubility_gbt.json', 'tm_model.json', 'mix_fit.json'].filter(fresh);
  for (const f of files) fs.copyFileSync(path.join(P.OUT, f), path.join(P.WEIGHTS, f));
  const refs = ['reference_gbt.json', 'reference_tm.json'].filter(fresh);
  // 自社のデータだけで学習した時は、照合の値に自社の構造が入る。git に入らない tests/own/ に写す（tests/check.cjs は、あればそちらを使う）
  const ownRef = refs.includes('reference_gbt.json') && JSON.parse(fs.readFileSync(path.join(P.OUT, 'reference_gbt.json'), 'utf8')).own_rows === true;
  const refDir = path.join(P.MODEL, 'tests', ownRef ? 'own' : '');
  if (ownRef) fs.mkdirSync(refDir, { recursive: true });
  for (const f of refs) fs.copyFileSync(path.join(P.OUT, f), path.join(refDir, f));
  console.log('\nweights/ へ写した: ' + files.join('、') + '。' + path.relative(P.MODEL, refDir) + '/ へ写した: ' + refs.join('、') + '（tests/reference_api.json は公開の重みの値なので、入口の照合は飛ばされる）');
  if (ownRef) console.log('照合の値に自社の構造が入っているので、tests/own/（.gitignore でリポジトリに入らない）に置いた。weights/ の重みも自社のデータを含むので、公開のリポジトリには入れない。');
}
console.log('\n終わり（' + ((Date.now() - t0) / 1000).toFixed(0) + ' 秒）。できた重み: ' + P.OUT);
