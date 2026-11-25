# 溶けるか・混ざるか（溶解度の見込み）

*English summary: a tiny solubility model that runs in the browser or Node with no build step. Given a solid (SMILES, or one of 37 example compounds) it estimates the solubility at 25 °C in 122 solvents (16 common ones are highlighted), with an uncertainty band (the range that contains about 70 % of errors on held-out compounds). It also handles binary solvent mixtures (log-linear mixing plus a fitted pair correction) and checks whether a stock solution will precipitate when diluted before HPLC injection. One measured value (or a melting point) narrows the band; melting points of about 19,000 structures are looked up automatically (open CC0 dataset). Typical accuracy on held-out compounds: organic solvents within ~3.3× (with melting point) / ~5.4× (without) for 70 % of cases; water ~3.8× / ~4.9×; mixed solvents ~1.9×. Code MIT, weights CC BY 4.0, trained only on commercially usable public data (BigSolDB 2.2, MixtureSolDB, AqSolDB, Abraham tables, an open melting-point dataset). Retraining kit with your own data included (train/). Docs in Japanese.*

この固体は、この溶媒に、目標の濃度まで溶けるか。混合溶媒ならどうか。DMSO の原液を水で薄めたら析出しないか。——その「だいたい」を、実験の前にすぐ見るためのモデルです。

- 25 °C での溶解度の見込み（mol/L と mg/mL）を、**幅つき**で返します。幅は「学習に使っていない固体で、外れの約 7 割が入る範囲」です。
- 溶媒は 122 種（主な 16 種：水・メタノール・エタノール・2-プロパノール・アセトニトリル・アセトン・酢酸エチル・THF・1,4-ジオキサン・DMSO・DMF・ジクロロメタン・クロロホルム・トルエン・ヘキサン・ヘプタン）。
- 融点が分かれば入れてください。幅が狭くなります（有機溶媒で約 5.4 倍 → 約 3.3 倍）。入れなくても、公開の表にある約 1.9 万構造は自動で引きます。
- 手元で 1 つの溶媒で測った値があれば、ほかの溶媒への換算になり、幅がいちばん狭くなります（約 2.5 倍）。
- 混合溶媒（2 成分）と、HPLC に打つ前の希釈での析出の確かめもできます。
- 計算はブラウザ（または Node）の中だけ。1 回の見込み（16 溶媒）は 0.2 ms ほど。構造はどこにも送りません。

| 版 | 中身 | 元にしたもの |
|---|---|---|
| 0.1 | 溶解度の木のモデル 4 つ（融点あり／なし × 有機溶媒／水）、式の係数、混合溶媒の補正、融点の推算、溶解度と融点の照合の表、例の固体 37 種 | サイトの体験ページ「溶けるか・混ざるか」の第 5 版（2026-11-25）と同じ計算・同じ重み |

## 使い方

### Node（18 以上。ほかに入れる物はありません）

```js
const SolubilityModel = require('./js/solubility.js');
const S = SolubilityModel.loadNode();               // weights/ を読む

const caf = S.fromExample('カフェイン');            // 例の固体（RDKit.js なしで使える）
const e = S.estimate(caf, 'meoh', { target: 10 });   // 目標 10 mg/mL
// e.mgPerML 見込み、e.range [下, 上]（mg/mL）、e.factor 幅の倍率、e.verdict 'yes' | 'edge' | 'no'、e.measured 文献の実測（あれば）

S.estimate(caf, 'acn', { reference: { solvent: 'meoh', mgPerML: 25 } });   // 手元の実測 1 点から換算
S.mix(caf, 'dmso', 'water', 0.3);                                          // DMSO 70 % / 水 30 %（体積分率）
S.dilution(caf, { stockSolvent: 'dmso', stockMgPerML: 50, stockUL: 20, diluentUL: 180, diluent: { water: 0.8, acn: 0.2 } });
// → { concMgPerML, solubilityMgPerML, range, verdict: 'stays_dissolved' | 'edge' | 'precipitates', ... }

// 自分の分子（SMILES）。構造を読むのに RDKit.js（リポジトリの vendor/rdkit）を使う
const R = await require('../../vendor/rdkit/RDKit_minimal.js')();
const c = S.prepare(R, 'CC(=O)Nc1ccc(O)cc1');                  // 融点は 入れた値（{ TmC: 169 }）＞ 照合の表 の順で入る
// c.TmC・c.tmSrc（'user' | 'sol' | 'mp' | null）、c.tmHint（融点が分からない時の推算。計算には使わない）、c.inLookup（溶解度の照合の表にある）
S.estimateAll(c, { target: 10 });                               // 主な 16 溶媒
S.estimateTm(c);                                                // 融点の推算（目安。計算には使わない）
```

全部を通した例：`node examples/node_example.cjs`。

### ブラウザ

```html
<script src="js/features.js"></script>
<script src="js/abraham.js"></script>
<script src="js/solub.js"></script>
<script src="js/solfeat.js"></script>
<script src="js/solubility.js"></script>
<script>
  SolubilityModel.load('weights/').then(S => {
    const r = S.estimate(S.fromExample('ナフタレン'), 'etoh', { target: 50 });
  });
</script>
```

- 自分の分子を足す時だけ、RDKit.js（`vendor/rdkit/RDKit_minimal.js`・`.wasm`、6.9 MB）を読んでください（例：`examples/browser.html`）。
- 重みの JSON は合わせて 3.8 MB（gzip で 0.9 MB）。照合の表 `lookup_25C.json`（1.5 MB）・`tm_lookup.json`（0.7 MB）と融点のモデル `tm_model.json` は、`load('weights/', { lookup: false, tmLookup: false, tm: false })` で読まないこともできます（その時は融点を引かない）。
- `fetch` で読むので、ファイルを直に開くのではなく、簡単なサーバーから開きます（`python -m http.server` など）。

### 返す値の要点

| 名前 | 意味 |
|---|---|
| `logS` | 見込みの log10(mol/L) |
| `mgPerML`・`range` | 見込みと、幅の下の端・上の端（mg/mL） |
| `band`・`factor` | 幅（log10 の片側）と、その倍率（10^band）。「約 3.3 倍」は、見込みの 1/3.3〜3.3 倍に外れの約 7 割が入るという意味 |
| `verdict` | `target` を渡した時の 3 段の判定。`yes` 幅の下の端でも目標より上、`no` 上の端でも目標より下、`edge` 目標が幅の中。文献の実測がある溶媒は、実測で判定する（`basis: 'measured'`） |
| `model` | 使ったモデル（`tm`・`notm`＝有機溶媒、`tm_w`・`notm_w`＝水） |
| `capped` | 見込みが純粋な溶質の濃度（密度 1.2 g/cm³ と置く）を超えたので、そこで切った |
| `extrapolated` | 学習データにその溶媒の点がない（係数の表にだけある溶媒）。幅は当てにならない |
| `parts` | 式の内訳（固体の項・水になじむ度合い・溶媒への移りやすさ・定数・実測に合わせた補正）。log S への寄与 |

## 精度（学習に使っていない固体で測った値）

外れの約 7 割が入る幅を、倍率で書きます（溶質ごとに 5 つの組に分けた検証。詳しくは [MODEL_CARD.md](MODEL_CARD.md)）。

| 場合 | 幅（約 7 割） | 参考：式だけ |
|---|---|---|
| 有機溶媒・融点あり | **約 3.3 倍** | 約 6 倍 |
| 有機溶媒・融点なし | **約 5.4 倍** | 約 10 倍 |
| 水・融点あり | **約 3.8 倍** | 約 14 倍 |
| 水・融点なし | **約 4.9 倍** | 約 26 倍 |
| 実測 1 点から換算（有機溶媒どうし） | **約 2.5 倍** | 約 3.2 倍 |
| 混合溶媒（2 成分、補正つき） | **約 1.9 倍**（両端の溶媒の値が分かっている時の、途中の外れ） | 対数線形だけで約 2.5 倍 |

- 約 1 割は、これより大きく外れます（有機溶媒・融点ありで 10 倍を超える物が 1 割）。
- 融点は、分かれば入れてください。**構造から推算した融点を入れても、幅は狭くなりませんでした**（だから推算値は計算に使いません。`estimateTm` は目安として返すだけ）。
- 例の固体 37 種は学習に使った固体なので、見込みが実測の近くに出ます。新しい分子での当たり方は、上の表の幅で考えてください。

## 向かない物・分かっている弱点

- **温度は 25 °C だけ**です。ほかの温度は計算しません。
- **塩・イオン液体・混合物（「.」を含む構造）、高分子、金属錯体は扱いません**。重い原子が 80 個を超える分子も断ります。
- **pH は考えていません**。酸・塩基は中性の形の見込みです。緩衝液や塩を含む水には使えません。
- **結晶形（多形・水和物・溶媒和物）は区別しません**。同じ物質でも結晶形で溶解度は数倍変わります。
- **3 成分の混合（原液の溶媒・水・アセトニトリルの希釈など）は、2 成分の組の和での近似で、3 成分の実測との照合はしていません**。補正の分の幅を 2 倍に取っています。
- 混合溶媒の補正は溶質によりません。共溶媒で極大を持つ系（ある割合でいちばん溶ける）は、形は出ても大きさは外れることがあります。
- 学習データにない溶媒（`extrapolated: true`）、とても溶けにくい物（10⁻⁶ mol/L 未満）、とても溶けやすい物（純粋な溶質に近い濃度）は外れやすいです。
- 出所どうしでも値は食い違います（BigSolDB と AqSolDB の両方にある水の値で、差の約 7 割が 2 倍以内、平均で AqSolDB が 0.22 log 低い）。

## 融点について

- 融点は、溶解度の見込みにいちばん効く入力です。`prepare` はページと同じ順で融点を入れます：**入れた値 ＞ 溶解度の照合の表（BigSolDB の融点、なければ融点の公開データ）＞ 融点の照合の表（`weights/tm_lookup.json`）**。出所は `tmSrc`（`user`・`sol` = BigSolDB・`mp` = 融点の公開データ）。
- `tm_lookup.json` は、融点の公開データ（Jean-Claude Bradley Open Melting Point Dataset v2、CC0 1.0）から作った 19,287 構造の融点（同じ構造は中央値）です。**元の表の source の列が DrugBank（875 行）、Aldrich Catalog-Handbook（742 行）、Sigma-Aldrich（90 行）である計 1,707 行は外しています**。ほかの出所（PHYSPROP、Alfa Aesar などほかの試薬会社のカタログ、ハンドブック、ウェブの資料など）の行はそのまま使いました。行ごとの出所は元の表の key・source・link の列にあります。
- CC0 は、元のデータセットを登録した人が手放せる権利についての宣言で、元の資料（PHYSPROP、Alfa Aesar、ハンドブック、ウェブの資料など）の権利者が個別に許したことや、値を保証することを表すものではありません。出所や値の誤りに気づいたら、その構造と元の表の key を添えて知らせてください（[ATTRIBUTION.md](ATTRIBUTION.md) の 6）。
- どちらの表にもない時は `TmC: null` で、融点なしのモデルを使います。`tmHint` に構造からの推算（5 °C 刻み、外れの約 7 割が ± 36 °C）を返しますが、計算には使いません（推算した融点を入れても幅は狭くならなかった）。
- 学習の土台の表（`train/base_data/`）にも同じ融点が入っているので、学習し直す時に何も取りに行かなくても、公開の重みと同じ物ができます。

## 名前について

- **溶質の名前は、配布元の表（BigSolDB・AqSolDB など）から写していません**。照合の表は構造（正規化した SMILES）で引きます。`prepare` で足した分子の `name` は、`{ name: '…' }` で渡した名前か、なければ SMILES です。照合の表にある分子かどうかは `inLookup` で分かります（例のページでは、名前の代わりに SMILES を出す）。
- 名前で引けるのは、例の固体 37 種の日本語の名前（こちらで付けた物）と、溶媒の名前（id・英語・日本語・略号）だけです。

## 自分のデータで学習し直す

[train/README.md](train/README.md) と [train/OWN_DATA.md](train/OWN_DATA.md)。CSV を `train/own_data/` に置いて `node train/run_all.cjs --own train/own_data/<ファイル>` で、公開データと合わせた重みができます（CPU で 3 分ほど）。自社のデータはリポジトリに入りません。

## ファイル

| 場所 | 中身 | ライセンス |
|---|---|---|
| `js/solubility.js` | 入口（読み込み・見込み・換算・混合・希釈・判定） | MIT |
| `js/solfeat.js`・`solub.js`・`abraham.js`・`features.js` | 木の計算と入力、式、Abraham の式、構造の記述子（HPLC のモデルと同じ features.js） | MIT |
| `weights/` | 重みと小さな表（JSON。溶解度・融点の照合の表を含む） | こちらの寄与は CC BY 4.0、元のデータの条件も残る（[ATTRIBUTION.md](ATTRIBUTION.md)、直下の LICENSE-DATA.md） |
| `examples/` | Node の例、ブラウザの 1 枚 | MIT |
| `tests/` | JS と学習側（scikit-learn）の照合、サイトのページの計算との照合、自社のデータの読み取り | プログラムは MIT、照合用の値（`reference_*.json`）は重みと同じ扱い |
| `train/` | 学習し直す道具・土台の表・例の CSV | 道具は MIT。土台の表・例の CSV（MixtureSolDB と BigSolDB の CC BY 4.0）・除外の key の一覧は MIT ではなく、データの条件（[train/base_data/README.md](train/base_data/README.md)） |

確かめ方：`node tests/check.cjs`（すべて合格なら、JS の計算が学習側と合っている）、`node tests/check_own_data.cjs`。

## 引用

使ったデータの出典は [ATTRIBUTION.md](ATTRIBUTION.md)。重みや表を配る時は、その出典の表示を残してください。
