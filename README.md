# quick-lab-models

逆相 HPLC の保持時間や化合物の溶解度など、実験の前に大まかな見当をつけるための軽量な予測モデル集です。

推論は 1 回あたり 1 ミリ秒未満で完了し、ブラウザ上（クライアントサイド）または Node.js だけで動作します。入力した分子構造や測定条件が外部サーバーへ送信されることはありません。

予測には誤差の幅を併記していますが、あくまで実験前の目安です。最終的な条件決定や精密な解析は専用ソフトや実測で行ってください。コードは MIT、学習済みの重みは CC BY 4.0 としており、社内利用・商用利用・製品への組み込みも可能です。

*English summary: Lightweight prediction models for everyday lab work (reversed-phase HPLC retention, solubility) that run locally in the browser or Node.js (< 1 ms per prediction, no server required). Code is MIT, weights are CC BY 4.0, and all models are trained only on commercially usable, non-share-alike public datasets. Retraining scripts for your own data are included. Documentation is in Japanese.*

## 収録モデル

| モデル | 予測対象 | Version | 内容 |
|---|---|---|---|
| [models/hplc-retention](models/hplc-retention/) | 逆相 HPLC の保持（保持係数・保持時間）。カラム、有機溶媒、組成／グラジエント、pH、温度を指定可能。実測値 1 点での補正にも対応 | 0.1 | 推論コード・重み・検証テスト・実行例・再学習スクリプト |
| [models/solubility](models/solubility/) | 固体の溶解度（混合溶媒にも対応） | 0.1 | 推論コード・重み・検証テスト・実行例・再学習スクリプト |

今後もモデルを追加していく予定です。フォルダ構成や呼び出し方は [docs/MODEL_LAYOUT.md](docs/MODEL_LAYOUT.md) に揃えてあるので、1 つ使えれば他も同じ要領で組み込めます。

## 試しかた・ローカル環境での利用

### ブラウザでそのまま試す

[技術の秘密](https://gijutsu-himitsu.com) 内に、ブラウザから操作できる体験ページを用意しています。計算はすべてブラウザ内で完結し、入力した構造や条件が外部に送られることはありません。

- 逆相 HPLC の保持予測：https://gijutsu-himitsu.com/lab/hplc/measure/
- 溶解度予測：公開準備中

※ 体験ページの HPLC の保持予測は、このリポジトリで配布している重みとは別のモデルです。体験ページでは、継承条件付き（CC BY-SA）の公開データも学習に含めた版（GNN を含む）を使っているため、予測値はこのリポジトリの重みと少し異なります。体験ページのモデルは、このリポジトリには含まれていません。

社内規定などで外部サイトを使いにくい場合は、このリポジトリの推論コードと重みからローカル用の HTML 画面を自前で用意することもできます。画面を自作するのが手間で、ローカル環境で動かせる HTML 一式がそのまま欲しい場合は、ご連絡いただければ ZIP にまとめてお渡しすることも可能です。[お問い合わせフォーム](https://gijutsu-himitsu.com/#contact)またはメール（contact@gijutsu-himitsu.com）からご相談ください。

### プログラムから組み込む

各モデルのフォルダ内に JavaScript の推論コードと最小限の実行例（ブラウザ／Node.js）が入っています。具体的な呼び出し方は各モデルの `README.md` を参照してください。

## 手元のデータで再学習する

各モデルの `train/` ディレクトリに、手元の実測データを加えて再学習するためのスクリプトと手順書を置いています。

1. 自分のデータを所定形式の CSV にまとめ、`train/own_data/` に配置します（フォーマットは各モデルの `train/OWN_DATA.md` を参照）。
2. 手順どおりに学習を実行すると、公開データと手元のデータを合わせた重みファイルが生成されます。
3. 生成された重みを `weights/` に入れ替えれば、既存の推論コードのまま動作します。`tests/` のスクリプトで、学習側（Python）と推論側（JavaScript）の計算結果が一致するかも確認できます。

- `train/own_data/` は `.gitignore` の対象なので、手元のデータがそのままリポジトリに含まれることはありません。
- 整形済みの公開学習データ（`train/base_data/`）を同梱しており、追加のダウンロードなしで公開重みを再現・再学習できます。
- 決定木モデルは一般的な PC の CPU で学習できます（目安：HPLC は約 10 分、溶解度は約 3 分）。GPU 環境があれば、HPLC では GNN（グラフニューラルネットワーク）も学習でき、未知化合物に対する予測精度が上がります。
- 手元のデータで再学習した重みを社内だけで使う場合、公開の義務はありません。外部へ配布する場合は、公開データ分の出典表記を残してください。

## ライセンスとデータ選定の方針

| 対象 | ライセンス |
|---|---|
| プログラム（`js/`・`train/`・`tests/`・`examples/`） | MIT（[LICENSE](LICENSE)） |
| 学習済みの重み・参照テーブル（`weights/`）、学習用ベースデータ（`train/base_data/`） | CC BY 4.0（[LICENSE-DATA.md](LICENSE-DATA.md)） |
| 同梱の RDKit.js（`vendor/rdkit/`） | 本体は BSD 3-Clause。同梱コンポーネント（フォント・InChI 等）は `vendor/rdkit/THIRD_PARTY_NOTICES.md` を参照 |
| 元の公開データ | 各データセットのライセンス（Apache 2.0 / CC BY 4.0 / CC0 1.0。詳細は各モデルの `ATTRIBUTION.md`） |

### 継承条件なし（非 SA）・商用可のデータに限定した理由

公開データセットの中には CC BY-SA（継承）が付いたものもあります。これらを学習に含めると精度はわずかに上がりますが（たとえば HPLC の保持時間誤差は 7.4% から 6.8% に改善）、CC BY-SA の継承義務が学習済みの重みにどこまで及ぶかは法的に明確ではありません。自社データを足して再学習したモデルを製品に組み込んだり社外へ渡したりする場面で扱いにくくなるため、本リポジトリでは学習レコード（保持値・溶解度・構造・化合物名）を **CC BY 4.0 / CC0 / Apache 2.0 のみ** に限定し、非営利限定（NC）や継承付き（SA）のデータは採用していません。

HPLC モデルの METLIN SMRT の測定条件（カラム・グラジエント・流量）は、論文に記載された事実情報を RepoRT 経由で参照しています（詳細は `models/hplc-retention/ATTRIBUTION.md`）。

### 利用・再配布時の注意点

- 社内利用のみの場合、モデルや重みを公開する必要はありません。
- 重みやデータを再配布する場合は、各モデルの `ATTRIBUTION.md` と `vendor/rdkit/` の権利表記を残してください。
- 本リポジトリのプログラム作成には AI（Claude・Codex）を使用しています。「技術の秘密」が許諾可能な範囲の権利を MIT ライセンスで提供します。第三者のコード・フォント・データの権利表記は各配布元に従います。
- 予測結果の正確性は保証しません（MIT および CC BY 4.0 の免責事項のとおりです）。重要な判断の前には必ず実測で確認してください。

## モデル追加時のルール

新しいモデルを追加する際の構成やルールは [docs/MODEL_LAYOUT.md](docs/MODEL_LAYOUT.md) にまとめています。

- 学習データは商用利用可能かつ継承条件のないもの（CC0 / CC BY / Apache 2.0 等）に限定し、配布元のライセンスを確認して出典を明記する
- 特定メーカーの製品名は出さず、実機 UI にも似せない
- ブラウザ上で軽量に推論できる設計にし、予測値とあわせて誤差の幅を出力する
- 重みと再学習用コードをセットで公開する。読み込み時に任意コードが実行され得る形式（pickle 等）は使わない

## 制作・連絡先

技術の秘密（https://gijutsu-himitsu.com）

不具合の報告や追加してほしいモデルの要望は GitHub Issues へどうぞ。ローカル用 HTML 一式の相談などは、[お問い合わせフォーム](https://gijutsu-himitsu.com/#contact)またはメール（contact@gijutsu-himitsu.com）からご連絡ください。
