# quick-lab-models：実験の前に試すモデル集

*English summary: tiny, fast prediction models for everyday lab questions ("how long will this compound stay on my HPLC column?", "will it dissolve in this solvent?"). They run entirely in the browser (or Node) in well under a millisecond per prediction. Free to use, including commercially: code is MIT, trained weights are CC BY 4.0, and every model was trained only on public data that allows commercial use without share-alike terms. You can add your own data and retrain. Docs are in Japanese.*

実験の前に「だいたいこうなりそう」をすぐ見るための、軽くて速い予測モデルの集まりです。

- **速い**：1 回の予測は 1 ミリ秒もかかりません。条件のつまみを動かしたら、すぐ答えが変わります。
- **ブラウザの中だけで動く**：入れた構造や条件を、どこにも送りません。サーバーも要りません。
- **当たりを付けるための物**：精度はだいたい 7〜9 割の当たりです。詳しい解析は、専用のソフトや実際の測定でどうぞ。予測には「このくらい外れそう」という幅も付けています。
- **自由に使ってください**：会社でも、製品に組み込んでも、改造して配っても大丈夫です。お金も届け出も要りません。

## モデルの一覧

| フォルダ | 何を予測するか | 版 | 中身 |
|---|---|---|---|
| [models/hplc-retention](models/hplc-retention/) | 逆相 HPLC の保持（保持係数・保持時間）。カラム・有機溶媒・割合・勾配・pH・温度を動かせる。自分で 1 回測った値で合わせ直せる | 0.1 | 推論・重み・照合・例・学習の道具 |
| [models/solubility](models/solubility/) | 固体が溶媒に溶けるか（溶解度）。混合溶媒も | 0.1 | 推論・重み・照合・例・学習の道具 |

これからも増やしていきます。どのモデルも同じ形（[docs/MODEL_LAYOUT.md](docs/MODEL_LAYOUT.md)）なので、1 つ使えたら全部同じように使えます。

## ライセンスと、そう決めた理由

| 何 | ライセンス |
|---|---|
| プログラム（`js/`・`train/`・`tests/`・`examples/`） | MIT（[LICENSE](LICENSE)） |
| 学習済みの重み・小さな表（`weights/`）と、学習の表（`train/base_data/`） | CC BY 4.0（[LICENSE-DATA.md](LICENSE-DATA.md)） |
| 同梱の RDKit.js（`vendor/rdkit/`） | 本体は BSD 3-Clause。中に入っている部品（フォント・InChI など）はそれぞれのライセンス（`vendor/rdkit/THIRD_PARTY_NOTICES.md`） |
| 元のデータ | それぞれの元のライセンス（Apache 2.0・CC BY 4.0・CC0 1.0。各モデルの `ATTRIBUTION.md`） |

**会社で使いたい人もいると思うので、CC BY にしました。**

公開されているデータの中には、CC BY-SA（継承つき）のデータもあります。これを学習に入れると精度は少し上がります（HPLC のモデルなら、保持時間の誤差が 7.4% から 6.8% に）。ただ、そうして作ったモデルを配る人に「同じ CC BY-SA で公開する」という義務がかかるおそれがあります（学習済みの重みにどこまで及ぶかは、はっきり決まっていません）。自社のデータを足して学習し直したモデルを製品に入れたい、社外に渡したい、という時に、これは困ります。

なので、このリポジトリのモデルの学習の行（保持値・溶解度などの値・構造・名前）は、**継承の義務がないデータ（CC BY 4.0・CC0・Apache 2.0）だけ**にしました。非営利限定（NC）のデータも使っていません。少しの精度より、だれでも気兼ねなく使えることを選びました。1 つだけ、HPLC の METLIN SMRT の測定条件（カラム・勾配・流量）は、論文に書かれた条件を継承つきの RepoRT の記載を経由して読んだ事実の値です（`models/hplc-retention/ATTRIBUTION.md`）。

守ってほしいのは 1 つだけです：**重みや表を配る時は、出典の表示（各モデルの `ATTRIBUTION.md`）と、同梱の部品の表示（`vendor/rdkit/`）を残してください。** 元のデータを作った研究者の方々への表示です。社内だけで使うなら、何も公開しなくて大丈夫です。

このリポジトリのプログラムは、AI（Claude・Codex）を使って作りました。「技術の秘密」は、自分が許諾できるプログラムの権利を MIT ライセンスの条件で提供します。第三者のコード・フォント・データの権利の表示は、それぞれの配布元の物です。

予測は保証しません（MIT と CC BY 4.0 の免責のとおりです）。大事な判断の前には、実際に測ってください。

## 自分のデータで学習し直す

各モデルの `train/` に、学習し直す道具と手順があります。

1. 自分のデータを、決まった形の CSV にして `train/own_data/` に置きます（形は各モデルの `train/OWN_DATA.md`）。
2. 手順どおりに学習すると、公開データと自分のデータを合わせた重みができます。
3. できた重みを `weights/` に入れ替えれば、同じプログラムでそのまま動きます。`tests/` で、ブラウザの計算と学習側の計算が合うかを確かめられます。

- `train/own_data/` の中身は `.gitignore` でリポジトリに入らないようにしてあります。自分のデータは手元から出ません。
- 決定木のモデルは CPU で学習できます（ふつうのパソコンで、HPLC は 10 分ほど、溶けるかは 3 分ほど）。公開データを整えた学習の表を同梱しているので、何もダウンロードしなくても、公開の重みと同じ物を作り直せます。GPU があれば、HPLC はグラフのモデル（分子を原子と結合のつながりとして読む小さなニューラルネットワーク）も学習できて、未知の分子への当たりが良くなります。
- 自分のデータで作った重みは、あなたの物です。社内だけで使うなら公開の義務はありません。配る時は、公開データの分の出典の表示を残してください。

## 新しいモデルを足す時の決まり

[docs/MODEL_LAYOUT.md](docs/MODEL_LAYOUT.md) にまとめています。要点：

- データは、商用で使えて継承のない物だけ（CC0・CC BY・Apache 2.0 など）。配布元のライセンスの欄で確かめて、出典を書く。
- 特定の会社の製品名は出さない・装置に似せない。
- 推論はブラウザで軽く。予測には幅を付ける。
- 重みと学習の道具を一緒に出す。読み込むだけでプログラムが動く形式（pickle など）では配らない。

## 作っている人

技術の秘密（https://gijutsu-himitsu.com）。直してほしい所・足してほしいモデルがあれば、Issue でお知らせください。
