# 溶けるか・混ざるか（組み合わせ版 0.2）

*English summary: a more accurate version of the solubility model in [../solubility](../solubility/), for 14 common solvents at 25 °C. Organic solvents use an average of residual-learning trees and a solute × solvent product-kernel regression, blended with a small distilled neural network; water uses a tree model with a similar-molecule correction. It runs in the browser or Node with no build step and no server (weights about 11 MB). Typical error on held-out solutes (70 % of cases): organic solvents within about 2.9× with a melting point and 4.6× without; water about 3.6× / 4.5×; with one to three of your own measurements about 2.3× to 1.9×. Same calculation and weights as the page at https://gijutsu-himitsu.com/lab/solubility/detail/ . Retraining scripts for this ensemble are **not** included yet; to retrain on your own data, use the tree model in ../solubility. Code MIT, weights CC BY 4.0.*

[../solubility](../solubility/)（木のモデル）より当たりやすい版です。サイトの[溶けるか・混ざるか](https://gijutsu-himitsu.com/lab/solubility/check/)と[くわしい版](https://gijutsu-himitsu.com/lab/solubility/detail/)（2026-11-27 の版）と、同じ重み・同じ計算です。

## どちらを使うか

| | このフォルダ（組み合わせ版 0.2） | [../solubility](../solubility/)（木のモデル 0.1） |
|---|---|---|
| 溶媒 | 主な 14 種 | 122 種（ヘキサン・ヘプタンを含む） |
| 幅（約 7 割、有機溶媒） | 約 2.9 倍（融点あり）／約 4.6 倍（なし） | 約 3.3 倍／約 5.4 倍 |
| 手元の実測からの換算 | 3 溶媒まで | 1 溶媒 |
| 対象外 | 分子量 500 以上と両性の分子は、有機溶媒の見込みを出さない | なし |
| 重みの大きさ | 約 11 MB | 約 7 MB |
| 手元のデータでの学び直し | **道具は未同梱** | 同梱（`train/`） |

2 つの幅は、外した物（対象外の分子、誤りの疑いのある実測）が違うので、差のすべてがモデルの良さではありません（[MODEL_CARD.md](MODEL_CARD.md)）。

## 溶媒

`water`・`meoh`・`etoh`・`ipa`・`acn`・`acetone`・`etoac`・`thf`・`dioxane`・`dmso`・`dmf`・`dcm`・`chcl3`・`toluene`

## 使い方

### Node（18 以上。ほかに入れる物はありません）

```js
const E = require('./js/solubility_ensemble.js');
const S = E.loadNode();                                    // weights/ を読む

const caf = S.fromExample('カフェイン');                    // 例の固体 37 種（RDKit.js なしで使える）
S.estimate(caf, 'meoh', { target: 10 });                    // 目標 10 mg/mL。mgPerML・range・factor・verdict・measured
S.estimate(caf, 'acn', { references: [{ solvent: 'meoh', mgPerML: 8 }, { solvent: 'water', mgPerML: 20 }] });  // 手元の実測から換算
S.mix(caf, 'acn', 'water', 0.3);                            // アセトニトリル 70 % / 水 30 %（体積分率）
S.dilution(caf, { stockSolvent: 'dmso', stockMgPerML: 50, stockUL: 10, diluentUL: 90, diluent: { water: 0.7, acn: 0.3 } });

// 自分の分子（SMILES）。構造を読むのに RDKit.js（リポジトリの vendor/rdkit）を使う
const R = await require('../../vendor/rdkit/RDKit_minimal.js')();
const c = S.prepare(R, 'CC(=O)Nc1ccc(O)cc1');               // 融点は 入れた値（{ TmC: 169 }）＞ 照合の表 の順で入る
S.estimateAll(c, { target: 10 });                           // 14 溶媒
```

全部を通した例：`node examples/node_example.cjs`。

### ブラウザ

```html
<script src="js/features.js"></script>
<script src="js/abraham.js"></script>
<script src="js/solub.js"></script>
<script src="js/solfeat.js"></script>
<script src="js/knn.js"></script>
<script src="js/orgmodel.js"></script>
<script src="js/solgnn.js"></script>
<script src="js/solubility_ensemble.js"></script>
<script>
  SolubilityEnsemble.load('weights/').then(S => {
    const r = S.estimate(S.fromExample('ナフタレン'), 'etoh', { target: 50 });
  });
</script>
```

- 自分の分子を足す時だけ、RDKit.js（`vendor/rdkit/RDKit_minimal.js`・`.wasm`、6.9 MB）を読んでください（例：`examples/browser.html`）。
- `fetch` で読むので、簡単なサーバーから開きます（リポジトリの直下で `python -m http.server 8000 --bind 127.0.0.1` を動かし、http://127.0.0.1:8000/models/solubility-ensemble/examples/browser.html を開く）。
- 入力した構造や数値は、どこにも送りません。
- 速さの目安：1 分子の 14 溶媒で、初回は 0.1 秒弱、2 回目からは 10 ミリ秒ほどです（開発に使った PC のブラウザでの値）。重みの読み込みは別にかかります。

## 返す値の要点

| 名前 | 意味 |
|---|---|
| `mgPerML`・`logS`・`range` | 見込み（mg/mL と log10 mol/L、溶液 1 L あたり）と、幅の下の端・上の端 |
| `band`・`factor` | 幅（log10 の片側）と、その倍率。「約 3 倍」は、見込みの 1/3〜3 倍に誤差の約 7 割が入るという意味 |
| `verdict`・`basis` | `target` を渡した時の判定（`yes`・`edge`・`no`）と、その根拠（`reference` 手元の実測・`measured` 文献の実測・`model` 見込み） |
| `outOfScope`・`reason` | 有機溶媒の対象外（`mw` 分子量 500 以上、`amph` 両性の分子）。この時 `mgPerML` は null |
| `measured` | 文献の実測（ある時だけ）。`mgPerML` はいつもモデルの見込みで、実測はここに別に入る |
| `mix` の `onePhase` | false なら 2 層に分かれる割合（見込みは出さない） |
| `dilution` の `warnings` | `stock_not_dissolved`（原液の段階で溶け残る見込み）・`may_separate_with_water` |

## 確かめ

```
node tests/check.cjs
```

呼び出し口の値が、サイトのページの計算と一致することを確かめます（例の固体は 10⁻⁹ まで一致、SMILES から足す分子は小数 2 桁まで一致）。

## 学び直し

この組み合わせを手元のデータで学び直す道具は、まだ同梱していません（ネットワークの学習に GPU を使う段があります）。手元のデータを足して学び直したい時は、[../solubility/train/](../solubility/train/) の木のモデルを使ってください。

## ライセンス

プログラムは MIT、重み・表・照合の値は CC BY 4.0 です。配る時は [ATTRIBUTION.md](ATTRIBUTION.md) を一緒に付けてください。
