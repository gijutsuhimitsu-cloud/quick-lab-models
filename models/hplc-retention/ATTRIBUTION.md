# 出典とライセンス（HPLC の保持のモデル）

*English: data sources used to train the weights in `weights/`, with their licenses and the changes we made. Keep this file when you redistribute the weights or weights derived from them.*

`weights/` の重みと表は CC BY 4.0 で配ります（リポジトリの直下の `LICENSE-DATA.md`）。中身は、下の公開データから作りました。どのデータも、配布元が商用の利用と再配布を許すライセンス（Apache 2.0・CC BY 4.0・CC0 1.0）を示しています。ライセンスは 2026-11-24 に、DOI の登録機関（DataCite）の記録と配布元の画面の License 欄で確かめました。これは事実の整理で、法律の意見ではありません。

**継承つき（CC BY-SA）のデータ、非営利限定のデータは、学習にも表にも入っていません。**

## 学習に使ったデータ

- **カラムのパラメータと、そこから再構成した保持係数**：Rutan, S. (2020), “Data for: Improvements in the Predictive Accuracy of the Hydrophobic Subtraction Model of Reversed-Phase Selectivity”, Mendeley Data, V1. https://doi.org/10.17632/hmd87ky8gh.1 。関連論文：Stoll ほか, J. Chromatogr. A 1636 (2021) 461682。ライセンス：CC BY 4.0。変更：カラムのパラメータは元の数値のまま転記し、kEB を log10(kEB) に変換、名前の前後の空白と固定相の区分を整理。**製品名とメーカー名を外し、「固定相の種類＋HSM2 データベースの番号」の一般名に替えました**（`weights/columns.json`）。溶質とカラムのパラメータから再構成した保持係数は、独立に測った保持ではありません。
- **保持時間と構造**：Domingo-Almenara, X. ほか, METLIN SMRT dataset, figshare, v1. https://doi.org/10.6084/m9.figshare.8038913.v1 。関連論文：https://doi.org/10.1038/s41467-019-13680-7 。ライセンス：Apache License 2.0（全文は `LICENSES/Apache-2.0.txt`。配布物に NOTICE のファイルはありません）。変更：秒から分への変換、構造からの特徴量の計算、行の選別、学習用の値（水 100% に外挿した保持の強さ）への変換。測定条件（カラムの種類・勾配のプログラム）は、論文に書かれた条件を RepoRT のメタデータ（データセット 0186）の記載を経由して読みました。RepoRT から使ったのはこの条件の記載だけで、RepoRT の保持値・構造は入っていません。学習の道具（`train/`）では、この条件を事実として書き込んであります。論文の本文とはまだ照合していません。死時間と装置の遅れは、こちらの見積もりです。
- **等組成の保持係数と、保持の式の定数**：Weber, S. (2020), “Retention data”, Mendeley Data, V3. https://doi.org/10.17632/x292z4gd25.3 。ライセンス：CC BY 4.0。変更：正の保持係数だけを選び、同じ条件の繰り返しの log10(k) を平均。化合物名に対応する構造はこちらで書き、特徴量の計算と式の定数の当てはめに使用。pH には見積もりを含みます。
- **保持時間と測定条件**：Beck, A. G.; Aggarwal, P. ほか (2025), “Code and data for Uncertainty-Aware Learning of Multiple Conditions as a Framework for Streamlined Retention Time Prediction to Accelerate Method Development”（MCaRT）, Zenodo. https://doi.org/10.5281/zenodo.17298529 。データのライセンス：CC BY 4.0。変更：行の選別、構造からの特徴量の計算、印の付いた保持値の除外、緩衝液からの pH の見積もり。同じ記録のプログラム（MIT）は使っていません。
- **保持時間・構造・測定条件と、メタノールの式の定数**：Kubik, Ł.; Jacyna, J.; Struck-Lewicka, W.; Markuszewski, M.; Wiczling, P. (2022), “LC-TOF-MS data collected for 300 small molecules”（5 カラムの 5 記録）, OSF. https://doi.org/10.17605/OSF.IO/ZQTJ7 、https://doi.org/10.17605/OSF.IO/QBV7J 、https://doi.org/10.17605/OSF.IO/2MCNW 、https://doi.org/10.17605/OSF.IO/Y6S8P 、https://doi.org/10.17605/OSF.IO/EVUJ9 。ライセンス：5 記録とも CC BY 4.0。変更：ピーク面積による代表ピークの選択、あいまいなピークの除外、カラム外の時間の補正、構造からの特徴量の計算。付属の pKa・logP の計算値は使っていません。
- **保持時間・構造・測定条件**：Zhang, Y. ほか, “Retention time dataset for heterogeneous molecules in reversed-phase liquid chromatography”, Science Data Bank, V3. https://doi.org/10.57760/sciencedb.15823 。CC0 1.0。変更：ファイルの選択、構造からの特徴量の計算、勾配の整理、添加剤からの pH の見積もり。Dataset 10 は学習に使わず、評価だけに使用。この表示は再現のためで、CC0 の義務ではありません。

## 表

- `weights/name_index_free.json`（化合物の名前 → 構造の索引）：上の Weber ら・Wiczling グループ（CC BY 4.0）・Zhang ら（CC0）のデータに付いている名前と、こちらで作った表。
- `weights/stock_compounds.json`（ストックの化合物 102 種の名前・構造・特徴量）：こちらで作った表。

## プログラムの部品

- **RDKit**（RDKit.js、@rdkit/rdkit 2025.3.4-1.0.0）：構造の読み込みと記述子の計算。BSD 3-Clause（`../../vendor/rdkit/LICENSE`）。

## このファイルを残してください

重み（`weights/` の物、またはそれを元に学習し直した物）を配る時は、このファイルを一緒に付けてください。自社のデータを足した時は「どの公開データに、何を足したか」を書き足してください。
