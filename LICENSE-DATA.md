# 重みと表のライセンス（CC BY 4.0）

*English: all trained weights and data tables in `models/*/weights/` and `models/*/train/base_data/` (plus the reference values in `models/*/tests/`) are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/ . When you share them (or weights you derived from them), keep the attribution in that model's `ATTRIBUTION.md`.*

各モデルの `weights/` に入っている学習済みの重みと小さな表、`train/base_data/` に入っている学習の表（公開データを整えた物）、`tests/` の照合用の値は、**クリエイティブ・コモンズ 表示 4.0 国際（CC BY 4.0）** で配ります。

学習の表のうち、元が Apache 2.0（HPLC の METLIN SMRT）や CC0（Zhang ら・AqSolDB・Bradley の融点など）のデータから来た部分は、元のライセンスでも使えます。どの行がどのデータから来たかは、表の出所の列やファイルの名前と、各モデルの `ATTRIBUTION.md` で分かります。

- ライセンスの要約：https://creativecommons.org/licenses/by/4.0/deed.ja
- ライセンスの全文：https://creativecommons.org/licenses/by/4.0/legalcode.ja

## できること

- 使う・写す・配る・改造する・製品に組み込む・売る。商用でも構いません。
- 自分のデータを足して学習し直した重みを作る。その重みを社内だけで使うなら、公開の義務はありません。

## 守ってほしいこと

- 重みや表（そこから作った物を含む）を配る時は、そのモデルの `ATTRIBUTION.md` の出典の表示を一緒に付けてください。
- 中身を変えた時は、変えたことが分かるようにしてください（例：「自社のデータを足して学習し直した」）。

## なぜ CC BY なのか

学習には、継承の義務がない公開データ（CC BY 4.0・CC0 1.0・Apache 2.0）だけを使いました。そのため、重みにも継承の義務を付けずに、CC BY 4.0 で配れます。会社で使う人、自社のデータで学習し直したモデルを製品に入れたい人も、気兼ねなく使えるようにするためです（リポジトリの [README.md](README.md) の「ライセンスと、そう決めた理由」）。

元のデータのうち、Apache 2.0 のデータ（HPLC のモデルの METLIN SMRT）は、その全文を各モデルの `LICENSES/` に同梱しています。

## 保証しないこと

重みと表は「そのまま」の状態で提供します。予測の正しさ、特定の目的への適合は保証しません（CC BY 4.0 の第 5 条）。
