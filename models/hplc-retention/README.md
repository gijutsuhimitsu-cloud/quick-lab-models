# HPLC の保持のモデル（逆相）

*English: predicts reversed-phase HPLC retention (retention factor k and retention time) from a molecule's 2D structure, the column (551 columns described by Hydrophobic Subtraction Model parameters, generic names only), organic modifier (acetonitrile / methanol), isocratic fraction or gradient program, pH and temperature. Runs in the browser or Node in ~0.1 ms per prediction, with a 68 % / 90 % uncertainty band, optional one-point calibration with your own measurement, and a mixture chromatogram helper. Weights CC BY 4.0, code MIT. See MODEL_CARD.md for accuracy and limits.*

分子の構造と、カラム・移動相の条件から、逆相 HPLC の保持（保持係数 k・保持時間 tR）を見積もります。

- 動かせる条件：カラム（551 本。名前は「C18 No.3」のような一般名）、有機溶媒（アセトニトリル・メタノール）、割合（等組成）か勾配のプログラム、pH、温度、カラムの寸法と流量。
- 予測ごとに「このくらい外れそう」という幅（外れの約 7 割・約 9 割が入る幅）が付きます。
- 自分の装置で 1 回測った値を入れると、同じ分子のほかの条件の予測を寄せられます。
- 何個かの分子を混ぜた時の、ピークの並び・重なり・順番の確かさも出せます（`js/chromatogram.js`）。
- 精度の目安：学習に使っていない化合物で、勾配の保持時間の誤差が分析時間の 7.4%（中央値）、溶出の順番の相関 0.86。くわしくは [MODEL_CARD.md](MODEL_CARD.md)。

## ファイル

| 場所 | 中身 | ライセンス |
|---|---|---|
| `js/retention_model.js` | 予測の入口（`ChromRetention.create(...).predictRetention(...)`） | MIT |
| `js/features.js` | 構造 → 特徴量（記述子 41・官能基 56・酸塩基の基 30）。学習側も同じファイルを使う | MIT |
| `js/chromatogram.js` | 混合物のクロマトグラム（`ChromMixture.run`・`trace`） | MIT |
| `js/gnn_model.js`・`js/graph_features.js` | グラフのモデルの計算（重みは未同梱。`train/gpu/` で作れる） | MIT |
| `weights/model_free.json` | 決定木 900 本の係数（2.7 MB、gzip で 0.79 MB） | CC BY 4.0 |
| `weights/uncertainty_free.json` | 予測の幅の係数 | CC BY 4.0 |
| `weights/columns.json` | カラム 551 本のパラメータ。`defaults` に C18・C8・フェニル・シアノ・極性基内包型の 5 本 | CC BY 4.0 |
| `weights/stock_compounds.json` | ストックの化合物 102 種（名前・SMILES・特徴量の事前計算） | CC BY 4.0 |
| `weights/name_index_free.json` | 化合物の名前 → SMILES の索引（759 件。別名には医薬品の商品名も含む。構造を探す手がかりで、製品との同一性は示さない） | CC BY 4.0 |
| `examples/` | Node の例・ブラウザの最小のページ | MIT |
| `tests/check.cjs` | JavaScript と学習側（Python）の予測の照合 | MIT |
| `train/` | 学習し直す道具・自社データの足し方 | MIT |
| `train/base_data/` | 学習の表（公開データを整えた物、約 20 万行） | CC BY 4.0（元のライセンス：METLIN SMRT は Apache 2.0、Zhang らは CC0。`ATTRIBUTION.md`） |

構造の計算には、リポジトリの直下の `vendor/rdkit/`（RDKit.js、BSD 3-Clause）を使います。ストックの化合物だけなら RDKit.js は要りません。

## すぐ試す

Node（v18 以上）で、このフォルダから：

```
node examples/node_example.cjs
node examples/node_example.cjs "CCc1ccccc1"
node tests/check.cjs
```

ブラウザ：リポジトリの直下で `python -m http.server 8000 --bind 127.0.0.1` を動かし、http://127.0.0.1:8000/models/hplc-retention/examples/ を開きます。

## 使い方

```html
<script src="vendor/rdkit/RDKit_minimal.js"></script>   <!-- SMILES から計算する時だけ -->
<script src="models/hplc-retention/js/features.js"></script>
<script src="models/hplc-retention/js/retention_model.js"></script>
<script>
  // model・columns・uncertainty は weights/ の JSON を fetch で読んだ物
  const P = ChromRetention.create(model, columns, uncertainty);   // 3 つ目は省ける（幅が出なくなる）

  // 分子：ストックの化合物はそのまま渡せる。SMILES からは次のとおり
  const RDKit = await initRDKitModule({ locateFile: () => 'vendor/rdkit/RDKit_minimal.wasm' });
  const mol = RDKit.get_mol('CC(C)Cc1ccc(cc1)C(C)C(=O)O'); mol.remove_hs_in_place();
  const f = ChromFeatures.compute(RDKit, mol); mol.delete();
  const compound = { desc: f.desc, groups: f.groups, ion: f.ion };

  const r = P.predictRetention(compound, columns.defaults[0].id, {
    organic: 'acn', phi: 0.5, pH: 2.8, tempC: 35,
    length_mm: 150, id_mm: 4.6, particle_um: 5, flow: 1.0
  });
  // r.k（等組成の時）、r.tR（分）、r.t0、r.uncertainty.tR68 = [早い側, 遅い側] など
</script>
```

Node では `require('./js/retention_model.js')` のように読めます（`examples/node_example.cjs`）。

### 条件（`cond`）

| 項目 | 意味 | 省いた時 |
|---|---|---|
| `organic` | `'acn'`（アセトニトリル）か `'meoh'`（メタノール） | アセトニトリル |
| `phi` | 有機溶媒の割合（0〜1、等組成） | 0.5 |
| `gradient` | `[[分, 割合], ...]` の折れ線。あれば `phi` より優先 | 等組成 |
| `pH` | 水側の pH。酸・塩基の基を持つ分子の保持だけが動く | 3.0 |
| `tempC` | カラムの温度（°C） | 欠測として扱う |
| `length_mm`・`id_mm`・`particle_um`・`flow` | カラムの長さ・内径（mm）・粒径（µm）・流量（mL/min） | 100・2.1・2.7・0.3 |
| `t0` | 死時間（分）。書けば寸法からの計算より優先 | 寸法と流量から（空隙率 0.65） |
| `dwell` | 装置の遅れ時間（分） | 0 |
| `additive` | 受け取るが、結果は変わらない（添加剤は pH を決める物として扱う） | — |

`organic` が `'acn'`・`'meoh'` 以外の時、割合が 0〜1 の外の時（`50` のように % で書いた時も）、`gradient` の時間が戻る時は、計算せずにエラーを投げます。

### 戻り値

`k`（保持係数。等組成の時）・`tR`（保持時間、分）・`t0`・`lkw`（水 100% に外挿した log k）・`S`・`N`（段数の目安）・`sigma`（ピークの標準偏差、分）・`uncertainty`（`s68`・`s90`・`tR68`・`tR90`・`k68`・`k90`・`seen`）。

### 自分で測った値で合わせ直す

```js
P.calibrate(compound, columnId, condMeasured, { tR: 7.42 });   // 等組成なら { k: 3.1 } でもよい
const r = P.predictRetention(compound, otherColumnId, otherCond); // 自動で寄せられる（r.calibration）
P.clearCalibration(compound);                                     // 消す
```

戻り値の `ok` が false の時は `reason`（`not_retained`・`too_early`・`too_late`）を見てください。`ratio90` が 1.5 を大きく超える時は、ピークの取り違えや条件の入れ違いを疑ってください。

### 混合物のクロマトグラム

```js
const m = ChromMixture.run(P, compounds, columnId, cond);  // m.peaks（早い順）・m.pairs（隣どうしの分離度と、順番の確かさ）
const y = ChromMixture.trace(m, { tMax: m.tMax, n: 600 }); // ピークを足し合わせた曲線
```

予測の外れはピークの幅より大きいので、並びは目安です。7 割の幅（`tR68`）を横棒で描く、順番があやしい隣どうしに印を付ける、などの見せ方をおすすめします。

## 学習し直す・自分のデータを足す

[train/README.md](train/README.md) と [train/OWN_DATA.md](train/OWN_DATA.md)。学習し直したら `node tests/check.cjs <新しい重みのフォルダ>` で、ブラウザの計算と学習側の計算が合うかを確かめてから入れ替えてください。

`js/features.js`（構造 → 特徴量の計算）を変えたら、必ず学習し直してください。学習側も同じファイルで特徴量を作るので、変えたまま古い重みを使うと予測がずれます。重みの `meta.feature_version`（いまは 4）は学習の特徴量の組み方の版で、`features.js` の中の `FEATURE_VERSION`（記述子の計算の版、いまは 2）とは別の番号です。`retention_model.js` はこの 2 つを自動では照合しないので、`tests/check.cjs` で確かめてください。

## 出典とライセンス

重み・表は CC BY 4.0、プログラムは MIT です。配る時は [ATTRIBUTION.md](ATTRIBUTION.md) を残してください。
