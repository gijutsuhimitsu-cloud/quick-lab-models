# 出典（溶けるか・混ざるか 0.1）

*English: data sources used to train and build the weights and tables in `weights/` and `train/base_data/`, with authors, DOIs, licenses and what we changed. When you redistribute the weights or tables (or anything derived from them), keep this attribution. All sources allow commercial use without share-alike terms (CC BY 4.0 or CC0 1.0).*

重み（`weights/`）と土台の表（`train/base_data/`）は CC BY 4.0 で配ります（リポジトリの [LICENSE-DATA.md](../../LICENSE-DATA.md)）。配る時は、このページの出典の表示を一緒に付けてください。ライセンスの欄は、どれも配布元の登録情報（DataCite・figshare・Zenodo・Harvard Dataverse のライセンスの欄）で確かめました（2026-11-25）。

**溶質の名前は、配布元の表から写していません**（どの表も、溶質は正規化した SMILES で持つ）。名前があるのは、例の固体 37 種の日本語の名前（こちらで付けた物）と、溶媒の名前（Abraham の溶媒の表の英語の名前と、こちらで付けた日本語の名前・略号）だけです。

## 使ったデータ

### 1. BigSolDB 2.2（溶解度）

- 作者：Krasnov, L. ほか
- 配布元：Zenodo、doi:[10.5281/zenodo.22648301](https://doi.org/10.5281/zenodo.22648301)（版をまとめる DOI は 10.5281/zenodo.15094978）
- 論文：*Scientific Data* (2025)、doi:[10.1038/s41597-025-05559-8](https://doi.org/10.1038/s41597-025-05559-8)（配布元が引用を求めている）
- ライセンス：CC BY 4.0
- 使ったファイル：`BigSolDBv2.2.csv`（溶解度）、`BigSolDBv2.2_gamma.csv`（融点の列）、`BigSolDBv2.2_densities.csv`（溶媒の密度）
- こちらで変えたこと：298.15 ± 5 K の点だけを使い、同じ溶質 × 溶媒の log S を平均して 1 点にした（11,488 組）。構造を RDKit.js で正規化し、塩・混合物と体積を出せない元素を含む物を外した。溶媒は Abraham の溶媒の表と構造で対応づけた（表にない溶媒の点は使っていない）。融点は溶質ごとの中央値を °C にした。各点の元の論文の DOI（Source の列）を、照合の表と土台の表に残した。溶質の名前（Compound_Name の列）は写していない。
- 入っている所：`weights/solubility_gbt.json`（学習）、`formula.json`（式の係数）、`lookup_25C.json`（実測の照合）、`example_compounds.json`、`solvents.json`（学習の点の数）、`train/base_data/bigsol_25C.csv.gz`・`bigsol_solutes.csv.gz`

### 2. MixtureSolDB（混合溶媒での溶解度）

- 作者：Malikov, D. ほか
- 配布元：Zenodo、doi:[10.5281/zenodo.18660057](https://doi.org/10.5281/zenodo.18660057)（版をまとめる DOI は 10.5281/zenodo.17846306）
- 論文：*Scientific Data* (2026)、doi:[10.1038/s41597-026-07047-z](https://doi.org/10.1038/s41597-026-07047-z)（配布元が引用を求めている）
- ライセンス：CC BY 4.0
- 使ったファイル：`MixtureSolDB.csv`、`solvent_boiling_point.csv`（この表の値の出所の欄は、公的な化学物質のデータベース）
- こちらで変えたこと：298.15 ± 5 K で、同じ論文・溶質・溶媒の組・温度の系列のうち、純粋な溶媒 2 つの値がそろった物だけを使った。濃度を g L⁻¹、割合を体積分率に換算し（BigSolDB の密度の表を使用。密度のない溶媒の 4,956 行は使っていない）、対数線形からのずれ（29,609 点）に溶媒の組ごとの係数を当てはめた。実測の値そのものは重みに入れていない。沸点は主な 16 溶媒だけ °C にした。
- 入っている所：`weights/mix_fit.json`、`solvents.json`（沸点）、`train/base_data/mix_series_25C.csv.gz`
- 学習し直しの例のデータ（`train/examples/own_data_example.csv`）も、この表の純粋な溶媒の点（公開の学習データに無い溶質 11 種）から作った。名前は「例-01」などに置き換え、濃度は mg/mL か mol/L に換算した。各行の `note` に元の論文の DOI を残した。

### 3. AqSolDB（水への溶解度）

- 作者：Sorkun, M. C.; Khetan, A.; Er, S.
- 配布元：Harvard Dataverse、doi:[10.7910/DVN/OVHAW8](https://doi.org/10.7910/DVN/OVHAW8)（v1.0、`curated-solubility-dataset.csv`）
- 論文：AqSolDB, a curated reference set of aqueous solubility and 2D descriptors for a diverse set of compounds. *Scientific Data* 6, 143 (2019)、doi:[10.1038/s41597-019-0151-1](https://doi.org/10.1038/s41597-019-0151-1)
- ライセンス：CC0 1.0（表示の義務はないが、ここに記す）
- こちらで変えたこと：塩・混合物（1,098）、電荷の合計が 0 でない物（86）、読めない構造（2）、体積を出せない元素（70）、重い原子 80 個超（7）を外し、同じ構造は平均した（8,718 構造）。溶質の名前は写していない（AqSolDB の番号 ID は土台の表に残した）。照合の表には、BigSolDB に水の値がない溶質だけ足した。AqSolDB は 9 つのデータをまとめた物で、元のデータごとの条件は確かめていない（配布元は全体を CC0 で配っている）。
- 入っている所：`weights/solubility_gbt.json`（水のモデル）、`lookup_25C.json`（水）、`train/base_data/aqsol_rows.csv.gz`・`aqsol_solutes.csv.gz`

### 4. Abraham の溶媒の係数（AbraLlama のデータ）

- 作者：Lang, A.; Lee, Y.
- 配布元：figshare、Dataset: AbraLlama: Predicting Abraham Model Solute Descriptors and Modified Solvent Parameters Using Llama（2024）、doi:[10.6084/m9.figshare.25898476](https://doi.org/10.6084/m9.figshare.25898476)（ファイル `20240524 Abraham Solvent Master Dataset.xlsx`）
- ライセンス：CC0 1.0
- こちらで変えたこと：水との混合物・水で飽和した溶媒・油・気相の行を外し、同じ構造の 2 つ目以降を外した（122 溶媒）。商品名で載っていた 1 溶媒の名前を、化学の名前（Diethylene glycol monoethyl ether）に置き換えた（係数は変えていない）。各溶媒の係数の元の文献は、`solvents.json` の `ref`・`year` の欄に表のまま残した。
- 入っている所：`weights/solvents.json`

### 5. Compounds with known Abraham descriptors（溶質の Abraham 記述子）

- 作者：Bradley, J.-C.; Acree, W.; Lang, A.
- 配布元：figshare（2014）、doi:[10.6084/m9.figshare.1176994](https://doi.org/10.6084/m9.figshare.1176994)（ファイル `20140922AbrahamDescriptorsBackUp.xls`）
- ライセンス：CC BY 4.0
- こちらで変えたこと：使わない印のある行、欠けた値（−123 で埋めてある）、取りうる範囲の外の値を外した。表の V が構造から計算した McGowan の体積と合わない 354 構造（名前と構造の対応があやしい）を外し、構造ごとに平均した（1,835 構造。名前の列は写していない）。これに一次式（リッジ回帰）を当てはめ、構造から E・S・A・B を推算する係数を作った。
- 入っている所：`weights/abraham_solute_model.json`（係数）、`train/base_data/abraham_solutes_clean.csv.gz`（平均した表）、`bigsol_solutes.csv.gz` の kE〜kB の列

### 6. Jean-Claude Bradley Open Melting Point Dataset（融点）

- 作者：Bradley, J.-C.; Williams, A.; Lang, A.
- 配布元：figshare（2014）、doi:[10.6084/m9.figshare.1031637](https://doi.org/10.6084/m9.figshare.1031637)（ファイル `BradleyMeltingPointDataset.xlsx`）
- ライセンス：CC0 1.0（表示の義務はないが、ここに記す）
- 中身：いろいろな出所から集めた有機化合物の融点の表（28,645 行）。行ごとの出所は、元のデータの source の列にある（公的なデータベース、試薬会社のカタログ、ハンドブック、論文、ウェブサイトなど）。
- こちらで変えたこと：出所の欄が、利用条件が非営利の薬のデータベース、または試薬会社のカタログである 1,707 行を読んでいない（行の番号は `train/rebuild/bradley_excluded_keys.txt`）。ほかの出所の行はそのまま使った。使わない印のある行、「.」を含む構造、読めない構造、同じ構造の値の開きが 20 °C を超える物も外し、構造を RDKit.js で正規化して、同じ構造は中央値にした（19,287 構造）。名前・出所の列は写していない（融点と構造だけ）。
- 入っている所：`weights/tm_lookup.json`（融点の照合の表）、`weights/tm_model.json`（融点の推算モデル）、`lookup_25C.json`・`example_compounds.json`・`train/base_data/bigsol_solutes.csv.gz`・`aqsol_solutes.csv.gz` の融点（出所の印 `mp`）、`train/base_data/tm_train.csv.gz`（融点のモデルの学習用の表：融点と構造の記述子）、`weights/solubility_gbt.json`（学習に使った融点として）

## 同梱した第三者のプログラム

- RDKit.js（@rdkit/rdkit 2025.3.4-1.0.0）：BSD 3-Clause。リポジトリの `vendor/rdkit/`（ライセンスの原文は `vendor/rdkit/LICENSE`）。構造の読み込み、記述子、正規化した SMILES に使う。

## 式の出典（式そのものを自分で実装した）

- 理想溶解度と融解エントロピーの近似（Walden の規則）。
- 一般溶解度式：Jain, N.; Yalkowsky, S. H. *J. Pharm. Sci.* 90, 234 (2001)。
- Abraham の溶媒和の式（線形自由エネルギー関係）と McGowan の特性体積：Abraham, M. H.; McGowan, J. C. *Chromatographia* 23, 243 (1987)。
- Crippen の log P（RDKit.js の CrippenClogP）：Wildman, S. A.; Crippen, G. M. *J. Chem. Inf. Comput. Sci.* 39, 868 (1999)。
- 混合溶媒の対数線形の式（共溶媒の対数線形モデル）と、溶媒の組ごとの補正（こちらで当てはめた）。

## 使っていない物

- **BigMixSolDB**（CC BY-NC-SA 4.0、非営利のみ・継承つき）。名前が MixtureSolDB・BigSolDB と似ているので、取り違えないように書いておく。ここから作った物は何も入っていない。
- **HSP（ハンセン溶解度パラメータ）**の表と推算法：商用で使える出所が見つからなかった（市販のソフトの表の写し、非営利・改変禁止の論文、購読の論文）。このモデルは HSP を使っていない。
- 非営利限定（NC）・改変禁止（ND）・継承（SA）のデータ、出所が分からないデータ。
