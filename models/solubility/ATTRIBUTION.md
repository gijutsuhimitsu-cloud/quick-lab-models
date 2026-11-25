# 出典（溶けるか・混ざるか 0.1）

*English: data sources used to train and build the weights and tables in `weights/`, `train/base_data/` and `train/examples/`, with authors, versions, DOIs, licenses and what we changed. Our own contributions (trained weights, the selection/processing/added parts of the tables, reference values) are offered under CC BY 4.0; upstream material keeps its original terms (CC BY 4.0 or CC0 1.0). When you redistribute the weights or tables (or anything derived from them), keep this attribution. All sources allow commercial use without share-alike terms.*

こちらが許諾できる物（学習済みの重み、表の選び方・加工・追加の部分、照合用の値）は CC BY 4.0 で提供します。元のデータには、それぞれの配布元の条件（下の各節。CC BY 4.0 か CC0 1.0）が引き続きかかります。CC0 の部分や、権利の対象にならない事実に、こちらが新しい条件を付けることはありません。どこに何の条件がかかるかは、リポジトリの直下の [LICENSE](../../LICENSE)（プログラム）と [LICENSE-DATA.md](../../LICENSE-DATA.md)（重み・表・照合の値）に書いてあります。

重みや表を配る時は、このページ（出典・変更の説明）を一緒に付けてください。ライセンス・版・作者は、配布元の登録情報（DataCite・figshare・Zenodo・Harvard Dataverse）で確かめました（2026-11-25）。

**溶質の名前は、配布元の表から写していません**（どの表も、溶質は正規化した SMILES で持つ）。名前があるのは、例の固体 37 種の日本語の名前（こちらで付けた物）と、溶媒の名前（Abraham の溶媒の表の英語の名前と、こちらで付けた日本語の名前・略号）だけです。

## 使ったデータ

### 1. BigSolDB 2.2（溶解度）

- 作者：Lev Krasnov、Dmitry Malikov、Marina Kiseleva、Sergei Tatarin、Sergey Sosnin、Stanislav Bezzubov
- 版：データの版 2.2（Zenodo の版の表示は v3）
- 配布元：Zenodo、doi:[10.5281/zenodo.22648301](https://doi.org/10.5281/zenodo.22648301)（版をまとめる DOI は 10.5281/zenodo.15094978）
- 論文：*Scientific Data* (2025)、doi:[10.1038/s41597-025-05559-8](https://doi.org/10.1038/s41597-025-05559-8)（配布元が引用を求めている）
- ライセンス：CC BY 4.0（[法文](https://creativecommons.org/licenses/by/4.0/legalcode.ja)）
- 使ったファイル：`BigSolDBv2.2.csv`（溶解度）、`BigSolDBv2.2_gamma.csv`（融点の列）、`BigSolDBv2.2_densities.csv`（溶媒の密度）
- こちらで変えたこと：298.15 ± 5 K の点だけを使い、同じ溶質 × 溶媒の log S を平均して 1 点にした（11,488 組）。構造を RDKit.js で正規化し、塩・混合物と体積を出せない元素を含む物を外した。溶媒は Abraham の溶媒の表と構造で対応づけた（表にない溶媒の点は使っていない）。融点は溶質ごとの中央値を °C にした。溶質の名前（Compound_Name の列）は写していない。
- 元の論文の DOI（Source の列）：土台の表 `bigsol_25C.csv.gz` には、平均した点の元の論文の DOI を **1 組につき最大 3 件**残した（2 件の組が 780、3 件の組が 99）。照合の表 `lookup_25C.json` と例の固体 `example_compounds.json` には、**代表の DOI を 1 件**（先頭の 1 件）だけ書いた。1 組に 4 件以上の論文が寄与していることもあるので、元の論文を全部たどる時は、BigSolDB の Source の列を見てください。
- 入っている所：`weights/solubility_gbt.json`（学習）、`formula.json`（式の係数）、`lookup_25C.json`（実測の照合）、`example_compounds.json`、`solvents.json`（学習の点の数）、`train/base_data/bigsol_25C.csv.gz`・`bigsol_solutes.csv.gz`

### 2. MixtureSolDB（混合溶媒での溶解度）

- 作者：Dmitry Malikov、Lev Krasnov、Marina Kiseleva、Elizaveta Meshcheriakova、Fedor Kuznetsov、Vladimir Elistratov、Matvei Vasiyarov、Sergei Tatarin、Stanislav Bezzubov
- 版：Zenodo の v3
- 配布元：Zenodo、doi:[10.5281/zenodo.18660057](https://doi.org/10.5281/zenodo.18660057)（版をまとめる DOI は 10.5281/zenodo.17846306）
- 論文：*Scientific Data* (2026)、doi:[10.1038/s41597-026-07047-z](https://doi.org/10.1038/s41597-026-07047-z)（配布元が引用を求めている）
- ライセンス：CC BY 4.0
- 使ったファイル：`MixtureSolDB.csv`、`solvent_boiling_point.csv`
- こちらで変えたこと：298.15 ± 5 K で、同じ論文・溶質・溶媒の組・温度の系列のうち、純粋な溶媒 2 つの値がそろった物だけを使った。濃度を g L⁻¹、割合を体積分率に換算し（BigSolDB の密度の表を使用。密度のない溶媒の 4,956 行は使っていない）、対数線形からのずれ（29,609 点）に溶媒の組ごとの係数を当てはめた。実測の値そのものは重みに入れていない。商品名で載っている溶媒 1 つは、構造で見分けて化学の名前にした。
- 沸点：主な 16 溶媒の沸点を °C にして `solvents.json` に入れた。元の表の Source の列は、16 件のうち 14 件が Hazardous Substances Data Bank（米国 National Library of Medicine）、1 件（水）が GESTIS Substance Database、1 件（ヘプタン）が International Chemical Safety Cards（米国 NIOSH 経由）。
- 入っている所：`weights/mix_fit.json`、`solvents.json`（沸点）、`train/base_data/mix_series_25C.csv.gz`
- **学習し直しの例のデータ `train/examples/own_data_example.csv`** も、この表の純粋な溶媒の点（公開の学習データに無い溶質 11 種）と、BigSolDB 2.2 の溶媒の密度から作った。名前は「例-01」などに置き換え、濃度は mg/mL か mol/L に換算した。各行の `note` に元の論文の DOI を残した。このファイルには元のデータの条件（MixtureSolDB と BigSolDB の CC BY 4.0）がかかる（プログラムの MIT ではない）。

### 3. AqSolDB（水への溶解度）

- 作者：Murat Cihan Sorkun、Abhishek Khetan、Süleyman Er
- 版：1.0
- 配布元：Harvard Dataverse、doi:[10.7910/DVN/OVHAW8](https://doi.org/10.7910/DVN/OVHAW8)（`curated-solubility-dataset.csv`）
- 論文：AqSolDB, a curated reference set of aqueous solubility and 2D descriptors for a diverse set of compounds. *Scientific Data* 6, 143 (2019)、doi:[10.1038/s41597-019-0151-1](https://doi.org/10.1038/s41597-019-0151-1)
- ライセンス：CC0 1.0（表示の義務はないが、ここに記す）
- こちらで変えたこと：塩・混合物（1,098）、電荷の合計が 0 でない物（86）、読めない構造（2）、体積を出せない元素（70）、重い原子 80 個超（7）を外し、同じ構造は平均した（8,718 構造）。溶質の名前は写していない（AqSolDB の番号 ID は土台の表に残した）。照合の表には、BigSolDB に水の値がない溶質だけ足した。
- 注意：AqSolDB は 9 つの元のデータをまとめた物で、登録は全体で CC0 1.0。元のデータごとの第三者の権利までは確かめていない。

### 4. Abraham の溶媒の係数（AbraLlama のデータ）

- 作者：Andrew Lang、Youngmin Lee
- 版：figshare の v1
- 配布元：figshare、Dataset: AbraLlama: Predicting Abraham Model Solute Descriptors and Modified Solvent Parameters Using Llama（2024）、doi:[10.6084/m9.figshare.25898476.v1](https://doi.org/10.6084/m9.figshare.25898476.v1)（ファイル `20240524 Abraham Solvent Master Dataset.xlsx`）
- ライセンス：CC0 1.0
- こちらで変えたこと：水との混合物・水で飽和した溶媒・油・気相の行を外し、同じ構造の 2 つ目以降を外した（122 溶媒）。商品名で載っていた 1 溶媒の名前を、化学の名前（Diethylene glycol monoethyl ether）に置き換えた（係数は変えていない）。各溶媒の係数の元の文献は、`solvents.json` の `ref`・`year` の欄に表のまま残した。
- 入っている所：`weights/solvents.json`

### 5. Compounds with known Abraham descriptors（溶質の Abraham 記述子）

- 作者：Jean-Claude Bradley、William Acree、Andrew Lang
- 版：figshare の v1（2014）
- 配布元：figshare、doi:[10.6084/m9.figshare.1176994.v1](https://doi.org/10.6084/m9.figshare.1176994.v1)（ファイル `20140922AbrahamDescriptorsBackUp.xls`）
- ライセンス：CC BY 4.0。登録の正式なライセンスの欄と DOI の登録は CC BY 4.0 で、配布元の説明文には「CC0」とも書いてある。この配布物では、出典・変更の説明と CC BY 4.0 の表示を残す（厳しい方に合わせる）。
- こちらで変えたこと：使わない印のある行、欠けた値（−123 で埋めてある）、取りうる範囲の外の値を外した。表の V が構造から計算した McGowan の体積と合わない 354 構造（名前と構造の対応があやしい）を外し、構造ごとに平均した（1,835 構造。名前の列は写していない）。これに一次式（リッジ回帰）を当てはめ、構造から E・S・A・B を推算する係数を作った。
- 入っている所：`weights/abraham_solute_model.json`（係数）、`train/base_data/abraham_solutes_clean.csv.gz`（平均した表）、`bigsol_solutes.csv.gz` の kE〜kB の列

### 6. Jean-Claude Bradley Open Melting Point Dataset（融点）

- 作者：Jean-Claude Bradley、Antony Williams、Andrew Lang
- 版：figshare の v2（DataCite の登録の版。2026-11-25 に配布元の API から取った時の最新の版で、取ったファイルの番号は 1503990）
- 配布元：figshare、doi:[10.6084/m9.figshare.1031637.v2](https://doi.org/10.6084/m9.figshare.1031637.v2)（ファイル `BradleyMeltingPointDataset.xlsx`）
- ライセンス：CC0 1.0（[宣言文](https://creativecommons.org/publicdomain/zero/1.0/legalcode.ja)。表示の義務はないが、ここに記す）
- 中身：いろいろな出所から集めた有機化合物の融点の表（28,645 行）。行ごとの出所は、元の表の key・source・link の列にある。
- **外した行**：元の表の source の列が **DrugBank（875 行）、Aldrich Catalog-Handbook（742 行）、Sigma-Aldrich（90 行）** である計 1,707 行を読んでいない（行の key は `train/rebuild/bradley_excluded_keys.txt`。この一覧はこちらで選んだ物）。**ほかの出所（PHYSPROP、Alfa Aesar などほかの試薬会社のカタログ、ハンドブック、ウェブの資料など）の行は使った**。
- こちらで変えたこと：使わない印のある行、「.」を含む構造、読めない構造、同じ構造の値の開きが 20 °C を超える物も外し、構造を RDKit.js で正規化して、同じ構造の採用した値の中央値を °C で書いた（19,287 構造の構造と融点）。名前・出所の列は写していない。
- 注意：採用した行には、PHYSPROP、Alfa Aesar、ハンドブック、ウェブの資料などから来た測定値が入っています。CC0 は、元のデータセットを登録した人が手放せる権利についての宣言で、それぞれの元の資料の権利者が個別に許したことや、値を保証することを表すものではありません。使い直す時は、必要に応じて元の表の key・source・link を見てください。**出所や値の誤りに気づいたら、その構造と元の表の key を添えて、このリポジトリの管理者に知らせてください。**
- 入っている所：`weights/tm_lookup.json`（融点の照合の表）、`weights/tm_model.json`（融点の推算モデル）、`lookup_25C.json`・`example_compounds.json`・`train/base_data/bigsol_solutes.csv.gz`・`aqsol_solutes.csv.gz` の融点（出所の印 `mp`）、`train/base_data/tm_train.csv.gz`（融点のモデルの学習用の表：融点と構造の記述子）、`weights/solubility_gbt.json`（学習に使った融点として）

## 同梱した第三者のプログラム

- RDKit.js（@rdkit/rdkit 2025.3.4-1.0.0）は、このフォルダには写さず、リポジトリの `vendor/rdkit/` の物を共有して使う（構造の読み込み、記述子、正規化した SMILES）。本体の BSD 3-Clause と、同梱されている第三者の部品・フォントの表示は `vendor/rdkit/` にある。

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

## 免責

元のデータと、こちらが提供する物には、それぞれの条件に書かれた免責がかかります。この出典の表示は、元の作者・権利者がこのモデルを推奨していることを示すものではありません。
