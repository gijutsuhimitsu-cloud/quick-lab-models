# 出典とライセンス（HPLC の保持のモデル）

*English: data sources used to train the weights in `weights/` and to build `train/base_data/`, with their licenses, copyright notices and the changes we made. Keep this file when you redistribute the weights, the tables, or anything derived from them.*

`weights/` の重み・表と `train/base_data/` の表のうち、こちらで作った部分は CC BY 4.0 で配ります（リポジトリの直下の `LICENSE-DATA.md`）。元のデータには、それぞれの元のライセンスがそのまま及びます（下の各項目）。

ライセンスは 2026-11-24 に DOI の登録機関（DataCite）の記録で確かめ、同じ日に利用者が配布元の画面の License 欄でも確かめました。2026-11-25 に DataCite の記録を確かめ直しています。これは事実の整理で、法律の意見ではありません。

**学習の行（保持値・構造・名前）には、継承つき（CC BY-SA）のデータと非営利限定のデータを使っていません。** 例外は 1 つで、METLIN SMRT の測定条件の記載を、継承つきの RepoRT の記録を経由して読んでいます（下の SMRT の項目）。

## 学習に使ったデータ

### カラムのパラメータ（HSM2 database）

- Rutan, Sarah (2020), “Data for: Improvements in the Predictive Accuracy of the Hydrophobic Subtraction Model of Reversed-Phase Selectivity”, Mendeley Data, V1. https://doi.org/10.17632/hmd87ky8gh.1 。関連論文：Stoll ほか, J. Chromatogr. A 1636 (2021) 461682。
- ライセンス：CC BY 4.0（https://creativecommons.org/licenses/by/4.0/legalcode ）。
- こちらで変えたこと：カラムのパラメータは元の数値のまま転記し、kEB を log10(kEB) に変換、名前の前後の空白と固定相の区分を整理。**製品名とメーカー名を外し、「固定相の種類＋HSM2 データベースの番号」の一般名に替えました**（`weights/columns.json`）。溶質とカラムのパラメータから再構成した保持係数は、独立に測った保持ではありません。

### 保持時間と構造（METLIN SMRT dataset）

- データの登録者：Domingo-Almenara, Xavier（DataCite の表記は Domingo, Xavier）, “The METLIN small molecule dataset for machine learning-based retention time prediction”, figshare, v1. https://doi.org/10.6084/m9.figshare.8038913.v1 。原論文：Domingo-Almenara ほか, Nature Communications 10, 5811 (2019). https://doi.org/10.1038/s41467-019-13680-7 。
- ライセンス：Apache License 2.0（全文は `LICENSES/Apache-2.0.txt`）。配布物に NOTICE のファイルはありません。
- 使ったファイル：`SMRT_dataset.csv`・`SMRT_dataset.sdf`。
- こちらで変えたこと：秒から分への変換、有効な構造を持つ行の選別、構造からの特徴量の計算、学習用の値（水 100% に外挿した保持の強さ）への変換。
- **測定条件（カラムの種類・寸法・勾配・流量）の経路**：論文に書かれた条件を、RepoRT（Kretschmer ほか、CC BY-SA 4.0、https://github.com/michaelwitting/RepoRT 、ライセンス https://github.com/michaelwitting/RepoRT/blob/master/LICENSE ）の記録 0186 の記載を経由して読み、こちらの形で書いた個別の事実の値です。論文の本文（Methods）との 1 項目ずつの照合は、まだしていません。RepoRT の保持値・構造は使っていません。死時間・装置の遅れ・pH は、こちらの見積もりです。

### 等組成の保持係数と、保持の式の定数（Weber ら）

- Weber, Stephen (2020), “Retention data”, Mendeley Data, V3. https://doi.org/10.17632/x292z4gd25.3 。
- ライセンス：CC BY 4.0。
- こちらで変えたこと：正の保持係数だけを選び、同じ条件の繰り返しの log10(k) を平均。化合物名に対応する構造はこちらで書き、特徴量の計算と式の定数の当てはめに使用。pH には見積もりを含みます。

### 保持時間と測定条件（MCaRT）

- Beck, Armen G.; Aggarwal, Pankaj ほか (2025), “Code and data for Uncertainty-Aware Learning of Multiple Conditions as a Framework for Streamlined Retention Time Prediction to Accelerate Method Development”, Zenodo, v1. https://doi.org/10.5281/zenodo.17298529 （記録 https://zenodo.org/records/17298529 ）。
- **著作権の表示：Copyright (c) 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.**
- ライセンス：記録の登録には、CC BY 4.0 と上の著作権の表示の両方が載っています。データは CC BY 4.0 の条件で使い、著作権の表示をここに残します。同じ記録のプログラムには別に MIT の表示があります（こちらはプログラムを使っていません）。
- 使ったファイル：`enamine_480.csv`。
- こちらで変えたこと：行の選別、構造からの特徴量の計算、印の付いた保持値の除外、緩衝液からの pH の見積もり、条件の整形。

### 保持時間・構造・測定条件と、メタノールの式の定数（Wiczling グループ）

“LC-TOF-MS data collected for 300 small molecules” の 5 記録（OSF、2022）。5 記録とも CC BY 4.0。作者は登録の順です。

| 記録 | 作者 |
|---|---|
| https://doi.org/10.17605/OSF.IO/ZQTJ7 | Łukasz Kubik, Julia Jacyna, Wiktoria Struck-Lewicka, Michał Markuszewski, Paweł Wiczling |
| https://doi.org/10.17605/OSF.IO/QBV7J | Wiktoria Struck-Lewicka, Paweł Wiczling, Julia Jacyna, Łukasz Kubik, Michał Markuszewski |
| https://doi.org/10.17605/OSF.IO/2MCNW | Wiktoria Struck-Lewicka, Paweł Wiczling, Julia Jacyna, Łukasz Kubik, Michał Markuszewski |
| https://doi.org/10.17605/OSF.IO/Y6S8P | Paweł Wiczling, Wiktoria Struck-Lewicka, Łukasz Kubik, Michał Markuszewski, Julia Jacyna |
| https://doi.org/10.17605/OSF.IO/EVUJ9 | Paweł Wiczling, Michał Markuszewski, Łukasz Kubik, Julia Jacyna, Wiktoria Struck-Lewicka |

こちらで変えたこと：ピーク面積による代表ピークの選択、あいまいなピークの除外、カラム外の時間の補正、構造からの特徴量の計算。付属の pKa・logP の計算値は使っていません。

### 保持時間・構造・測定条件（Zhang ら）

- Zhang, Yan; Gao, Yan; Liu, Fei; Li, Xiu Qin; Li, Cong Kang; Zhang, Qing He, “Retention time dataset for heterogeneous molecules in reversed-phase liquid chromatography”, Science Data Bank, V3. https://doi.org/10.57760/sciencedb.15823 。
- CC0 1.0（https://creativecommons.org/publicdomain/zero/1.0/legalcode ）。この表示は再現のためで、CC0 の義務ではありません。
- こちらで変えたこと：ファイルの選択、構造からの特徴量の計算、勾配の整理、添加剤からの pH の見積もり。Dataset 10 は学習に使わず、評価と、自社データの形の例（`train/examples/zhang_dataset10/`、元の CC0 のまま）にだけ使いました。

## 表

- `weights/name_index_free.json`（化合物の名前 → 構造の索引、759 件）：上の Weber ら・Wiczling グループ（CC BY 4.0）・Zhang ら（CC0）のデータに付いている名前と、こちらで作った表。**配布元が使った別名には、医薬品の商品名なども含まれます。** 名前は構造を探す手がかりで、製剤・塩の形との同一性や、メーカーの推奨・提携を示すものではありません。
- `weights/stock_compounds.json`（ストックの化合物 102 種の名前・構造・特徴量）：こちらで作った表。

## 測定値 1 点で合わせ直す係数

`js/retention_model.js` の既定の 5 つの係数は、公開の重みと同じ学習データ（上の 6 つ）での交差検証から決めた値です（`train/scripts/fit_calibration.py`）。

## プログラムの部品

- **RDKit**（RDKit.js、@rdkit/rdkit 2025.3.4-1.0.0）：構造の読み込みと記述子の計算。本体は BSD 3-Clause。同梱の WebAssembly には第三者のライブラリとフォントが入っています（表示とライセンス文は `../../vendor/rdkit/`）。

## このファイルを残してください

重み・表（`weights/` と `train/base_data/` の物、またはそれを元に学習し直した物）を配る時は、このファイルを一緒に付けてください。変えた時は、前の変更の表示を残し、新しい変更の要約を書き足してください（例：「自社の測定値を足して学習し直した」）。自社のデータそのものや秘密の内容を明かす必要はありません。元の作者・権利者が、この配布物や派生物を推奨していることを示すものではありません。
