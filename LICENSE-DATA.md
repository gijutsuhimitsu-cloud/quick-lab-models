# 重み・表のライセンス（CC BY 4.0）

*English: we offer our licensable contributions to the trained weights, the data tables (`models/*/weights/`, `models/*/train/base_data/`) and the reference values in `models/*/tests/` under the Creative Commons Attribution 4.0 International License (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/legalcode . Upstream material keeps its original terms (Apache-2.0, CC BY 4.0, CC0-1.0); each model's `ATTRIBUTION.md` lists the sources, the changes and the notices. When you share the weights or tables (or anything derived from them), keep that attribution.*

## 何に、どのライセンスか

- **こちらが許諾できる部分**（学習済みの重み、表の組み方・特徴量・こちらで作った表、照合用の値）は、「技術の秘密」から **クリエイティブ・コモンズ 表示 4.0 国際（CC BY 4.0）** で提供します。商用の利用、改変、製品への組み込み、再配布ができます。
  - 要約：https://creativecommons.org/licenses/by/4.0/deed.ja
  - 全文：https://creativecommons.org/licenses/by/4.0/legalcode.ja
- **元のデータ**には、配布元の条件がそのまま及びます。HPLC の METLIN SMRT に由来する部分は Apache License 2.0（各モデルの `LICENSES/`）、CC BY の物は CC BY 4.0、CC0 の物は CC0 1.0 です。元の CC0 の部分や、権利の対象にならない事実に、こちらが新しい条件を付けることはありません。
- **自社データの形の例**も、元のデータの条件のままです。HPLC の `train/examples/zhang_dataset10/` は Zhang らの CC0 のデータを整えた物、溶けるかの `train/examples/own_data_example.csv` は MixtureSolDB（CC BY 4.0）の測定と BigSolDB（CC BY 4.0）の密度から作った物です。
- 同梱の RDKit.js（`vendor/rdkit/`）は、それぞれの部品のライセンスです（`vendor/rdkit/` の中の表示）。

## 配る時にしてほしいこと

- 作者・指定された表示、提供された著作権・ライセンス・免責の表示、元の資料へのリンク、変更したことを、配り方に合った方法で伝えてください。各モデルの `ATTRIBUTION.md` を一緒に付けるのがいちばん簡単です。
- 変えた時は、前の変更の表示を残し、新しい変更の要約を書いてください（例：「自社の測定値を足して学習し直した」）。**自社のデータそのものや、秘密の内容を明かす必要はありません。** 社内だけで使うなら、何も公開しなくて大丈夫です。
- 元の作者・権利者や「技術の秘密」が、あなたの配布物を推奨していることを示すものではありません。

## なぜ CC BY なのか

会社で使う人、自社のデータで学習し直したモデルを製品に入れたい人も、気兼ねなく使えるようにするためです。そのため、学習の行（保持値・溶解度などの値・構造・名前）には、継承の義務や非営利の条件があるデータを採用していません（リポジトリの [README.md](README.md) の「ライセンスと、そう決めた理由」）。

## 保証しないこと

重み・表は「そのまま」の状態で提供します。予測の正しさ、特定の目的への適合は保証しません（CC BY 4.0 の第 5 条。元の資料にも、それぞれの条件の免責があります）。

これは事実の整理で、法律の意見ではありません。
