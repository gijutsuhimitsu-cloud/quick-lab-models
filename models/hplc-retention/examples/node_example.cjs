// いちばん短い使い方（Node）。models/hplc-retention/ で実行：
//   node examples/node_example.cjs
//   node examples/node_example.cjs "CC(C)Cc1ccc(cc1)C(C)C(=O)O"    （SMILES を渡すと、その分子で計算）
'use strict';
const fs = require('fs'), path = require('path');
const HERE = path.join(__dirname, '..');
const ChromRetention = require(path.join(HERE, 'js', 'retention_model.js'));
const ChromFeatures = require(path.join(HERE, 'js', 'features.js'));
const initRDKitModule = require(path.join(HERE, '..', '..', 'vendor', 'rdkit', 'RDKit_minimal.js'));
const read = (f) => JSON.parse(fs.readFileSync(path.join(HERE, 'weights', f), 'utf8'));

(async () => {
  const P = ChromRetention.create(read('model_free.json'), read('columns.json'), read('uncertainty_free.json'));
  const columns = read('columns.json');
  const col = columns.defaults[0];                       // 最初の 5 本のうち 1 本目（C18）

  // 分子：SMILES から特徴量を作る（RDKit.js）
  const RDKit = await initRDKitModule();
  const smiles = process.argv[2] || 'CC(C)Cc1ccc(cc1)C(C)C(=O)O';   // 既定はイブプロフェン
  const mol = RDKit.get_mol(smiles);
  if (!mol || !mol.is_valid()) { console.error('SMILES を読めません: ' + smiles); process.exit(1); }
  mol.remove_hs_in_place();
  const f = ChromFeatures.compute(RDKit, mol); mol.delete();
  const compound = { desc: f.desc, groups: f.groups, ion: f.ion };

  // 条件：カラムの寸法と流量は、書かなければ 100 × 2.1 mm・2.7 µm・0.3 mL/min
  const base = { organic: 'acn', pH: 2.8, tempC: 35, length_mm: 150, id_mm: 4.6, particle_um: 5, flow: 1.0 };
  console.log('分子: ' + smiles + '\nカラム: ' + col.name + '\n');

  // 等組成：有機溶媒の割合を変えて
  for (const phi of [0.3, 0.5, 0.7]) {
    const r = P.predictRetention(compound, col.id, Object.assign({ phi }, base));
    console.log('アセトニトリル ' + (phi * 100) + '%  k = ' + r.k.toFixed(2) + '（7 割の幅 ' + r.uncertainty.k68.map(v => v.toFixed(2)).join('〜') + '）  tR = ' + r.tR.toFixed(2) + ' 分');
  }
  // 勾配：0 分 5% → 20 分 95%
  const g = P.predictRetention(compound, col.id, Object.assign({ gradient: [[0, 0.05], [20, 0.95], [25, 0.95]] }, base));
  console.log('\n勾配 5→95%（20 分）  tR = ' + g.tR.toFixed(2) + ' 分（7 割の幅 ' + g.uncertainty.tR68.map(v => v.toFixed(2)).join('〜') + ' 分）');

  // pH：酸（イブプロフェンなど）は pH を上げると早く出る
  const p7 = P.predictRetention(compound, col.id, Object.assign({}, base, { phi: 0.5, pH: 7 }));
  console.log('アセトニトリル 50%・pH 7  k = ' + p7.k.toFixed(2));

  // 自分の装置で 1 回測った値で合わせ直す（同じ分子のほかの条件の予測が寄る）
  P.calibrate(compound, col.id, Object.assign({ phi: 0.5 }, base), { tR: 6.0 });
  const c = P.predictRetention(compound, col.id, Object.assign({ phi: 0.3 }, base));
  console.log('\n50% で tR = 6.0 分と測れた場合 → 30% の予測 k = ' + c.k.toFixed(2) + '（7 割の幅 ' + c.uncertainty.k68.map(v => v.toFixed(2)).join('〜') + '）');
})();
