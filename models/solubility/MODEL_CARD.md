# モデルカード：溶けるか・混ざるか（溶解度） 0.1

*English summary: gradient-boosted tree ensembles (scikit-learn, exported to JSON) that predict log10 solubility (mol/L) of neutral organic solids at 25 °C. Inputs: 41 RDKit descriptors, 56 functional-group counts, estimated Abraham solute descriptors, Crippen log P, an ideal-solubility (melting point) term, the solvent's Abraham coefficients, and the value of a fitted physical equation. Separate models for organic solvents and water, each with/without melting point. Validation: 5-fold cross-validation grouped by solute. Mixed solvents: log-linear mixing plus a per-pair correction fitted to MixtureSolDB.*

## 何を予測するか

- 25 °C（298.15 ± 5 K のデータで学習）での、中性の有機の固体の溶解度 log10(mol/L)。
- 溶媒は Abraham の溶媒の係数がある 122 種。そのうち 64 種に学習の点があり、16 種を主な溶媒として扱う（`weights/solvents.json` の `n_rows`・`main`）。
- 予測には幅を付ける：溶質を分けた検証での外れの絶対値の 68 パーセンタイル（log10）。

## 仕組み

1. **溶質の入力**（`js/features.js`・`abraham.js`）：RDKit.js で記述子 41 個と官能基の数 56 個を数える。Abraham 記述子 E・S・A・B は一次式（リッジ回帰）で推算、V は McGowan の式で計算。
2. **式**（`js/solub.js`）：固体の項（理想溶解度。ΔS_f = 56.5 J mol⁻¹ K⁻¹ の Walden の規則）、水は一般溶解度式（Yalkowsky）、有機溶媒は `log S = k0 + kP(−log P_ow) + kC·crystal + kT·log P(水→溶媒)` の 4 つの係数を BigSolDB に合わせた物。融点が分からない時は 150 °C（実測のある固体の中央値）と置く。
3. **木のモデル**（`js/solfeat.js`）：上の全部（114 列）を入力に、勾配ブースティングの木 500 本（深さ 5、学習率 0.04、行の 8 割・列の半分を使う、葉に 10 点以上）。
   - `tm`・`notm`：有機溶媒用。BigSolDB 2.2 の 25 °C 付近の 11,488 組（溶質 1,400、溶媒 64）。`tm` は融点の分かる溶質の行だけ（有機溶媒 8,416 組と水 433 組）、`notm` は全部の行で、融点の列を使わない。ページと同じく、これらは有機溶媒の見込みにだけ使う。
   - `tm_w`・`notm_w`：水用。BigSolDB の水（重み 3）と AqSolDB（重み 1）、9,355 行・溶質 9,100。
4. **実測 1 点からの換算**：手元の実測に、モデルが見込む溶媒どうしの差を足す（固体の項が消える）。
5. **混合溶媒**（`weights/mix_fit.json`）：`log C(φ) = (1−φ) log C_A + φ log C_B + φ(1−φ)[A0 + A1(1−2φ)]`。A0・A1 は溶媒の組ごと（溶質 8 種以上の 22 組）、それ以外は「水＋有機」「有機＋有機」の平均。幅は両端の幅と補正の外れを二乗和で合わせる。
6. **融点**：分かっている融点は、入れた値 ＞ 溶解度の照合の表 ＞ 融点の照合の表（`weights/tm_lookup.json`、19,287 構造）の順で使う。融点の推算（`weights/tm_model.json`、木 250 本）は目安として返すだけで、溶解度の計算には使わない（下の「分かったこと」）。

## 評価の仕方

- 溶質ごとに 5 つの組に分け（scikit-learn の GroupKFold）、4 組で学習して残りの 1 組で外れ（見込み − 実測、log10）を測る。同じ溶質が学習と検証の両方に入らない。
- 「実測 1 点から換算」は、同じ溶質のほかの溶媒の実測から見込んだ時の外れ（検証の組の中で）。
- 混合溶媒は、MixtureSolDB の 25 °C 付近の系列（同じ論文・溶質・溶媒の組・温度で、純粋な溶媒 2 つの値がそろった物）の途中の点 29,609 点。溶質ごとに 5 分割。

## 結果

| 場合 | 点の数 | 外れの 68 % 幅（log10） | 倍率 | 90 % 幅 | 式だけ（68 %） |
|---|---|---|---|---|---|
| 有機溶媒・融点あり（`tm`） | 8,416 | 0.523 | 3.3 倍 | 1.01 | 0.787 |
| 有機溶媒・融点なし（`notm`） | 10,851 | 0.733 | 5.4 倍 | 1.31 | 1.009 |
| 水・融点あり（`tm_w`、BigSolDB と AqSolDB） | 4,142 | 0.575 | 3.8 倍 | 1.22 | 1.135 |
| 水・融点なし（`notm_w`） | 9,355 | 0.692 | 4.9 倍 | 1.40 | 1.421 |
| 　水・融点あり、BigSolDB の水だけで測ると | 433 | 0.672 | 4.7 倍 | 1.58 | — |
| 　水・融点なし、BigSolDB の水だけで測ると | 637 | 0.842 | 7.0 倍 | 1.69 | — |
| 実測 1 点から換算（有機溶媒どうし、`tm`） | 87,502 組 | 0.400 | 2.5 倍 | 0.83 | 0.508 |
| 混合溶媒・補正つき（全部） | 29,609 | 0.277 | 1.9 倍 | 0.55 | 対数線形 0.401 |
| 　水＋有機溶媒 | 16,331 | 0.328 | 2.1 倍 | 0.60 | — |
| 　有機溶媒どうし | 13,278 | 0.227 | 1.7 倍 | 0.42 | — |
| 融点の推算（溶解度の学習の溶質 675 種を外して） | 675 | ± 36 °C | — | ± 64 °C | 平均で置くと ± 77 °C |

- 水の幅は `tm_w`・`notm_w` の検証の値（BigSolDB と AqSolDB の水の全部）。BigSolDB の水だけで測ると広い（上の表の字下げの行）。
- 対数線形の混合は、実測より平均で約 2 倍（−0.33 log）低く出た（混ぜると単純な内挿より溶けやすい）。補正でずれの平均は −0.01 に。

## 分かったこと

- 実測の融点は効く（有機溶媒で 0.733 → 0.523）。**構造から推算した融点は効かなかった**（推算値を入れた時の幅は、融点なしのモデルと同じ）。構造から読める分は、融点なしのモデルがすでに使っている。
- 有機溶媒と水を 1 つのモデルで学習し、そこに AqSolDB を足すと、AqSolDB の水は良くなるが有機溶媒が悪くなった。水だけのモデルを別に作ると両方で良くなった。
- 試して採らなかった物：溶質ごとの重み、深さ 6、列を 3 割、構造の指紋のビットを入力に足す（どれも差が 0.01 以内か悪化）。似た分子の外れで直す補正は改善が小さく、検証の甘さもあるので 0.1 では入れていない。

## データ

| データ | 使い方 | ライセンス |
|---|---|---|
| BigSolDB 2.2 | 溶解度の学習（有機溶媒・水）、式の係数、融点（gamma の表）、溶媒の密度 | CC BY 4.0 |
| AqSolDB | 水のモデルの学習、照合の表の水 | CC0 1.0 |
| MixtureSolDB | 混合溶媒の補正の係数、主な溶媒の沸点 | CC BY 4.0 |
| Abraham Solvent Master Dataset | 溶媒の係数 c・e・s・a・b・v | CC0 1.0 |
| Compounds with known Abraham descriptors | 溶質の記述子の推算式の学習 | CC BY 4.0 |
| Jean-Claude Bradley Open Melting Point Dataset | BigSolDB に融点がない溶質の融点（学習と、融点の照合の表 `tm_lookup.json`）、融点の推算モデル | CC0 1.0（v2。元の表の source の列が DrugBank（875 行）、Aldrich Catalog-Handbook（742 行）、Sigma-Aldrich（90 行）である計 1,707 行を除いた。PHYSPROP・Alfa Aesar などほかの出所の行は使った） |

作者・DOI・こちらで変えたことは [ATTRIBUTION.md](ATTRIBUTION.md)。

## 弱点と使い方の注意

- 25 °C だけ。pH・塩・結晶形（多形・水和物）は考えない。塩・混合物・高分子・重い原子 80 個超は断る。
- 学習の溶質は医薬品・有機合成の中間体・天然物に偏る。無機物、界面活性剤、とても大きい分子、フッ素の多い物は外れやすい。
- 学習データに点のない溶媒（`extrapolated: true`）は、溶媒の係数から外挿するだけで、幅は当てにならない。
- 3 成分の混合（希釈の確かめ）は 2 成分の組の和で近似、実測との照合なし。
- データの出所どうしで値が食い違う（BigSolDB と AqSolDB の水で、差の約 7 割が 2 倍以内）。幅の中には、この測定のばらつきも入っている。
- 例の固体（`weights/example_compounds.json`）と照合の表の固体は学習に使った物。そこでの当たり方は、新しい分子より良く見える。

## 版

| 版 | 日付 | 変えたこと |
|---|---|---|
| 0.1 | 2026-11-25 | 最初の公開。サイトの体験ページの第 5 版と同じ重み・同じ融点の照合の表（学習し直して木が完全に一致することを確かめた）。溶質の名前は配布元の表から写していない |
