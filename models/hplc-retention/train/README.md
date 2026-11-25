# 学習し直す道具（HPLC の保持のモデル）

*English summary: this folder retrains the reversed-phase HPLC retention model from the bundled base data (processed rows of six open datasets: Apache 2.0 / CC BY 4.0 / CC0), reproduces the published weights in `../weights/`, and lets you add your own measurements (`OWN_DATA.md`) and retrain on CPU. Run `python train/run.py` from `models/hplc-retention/` (about 10 minutes on a 24-thread PC). Optional: a GPU graph model (`gpu/`), a full rebuild of the base data from the original downloads, re-fitting the one-point calibration rule. Code: MIT. Base data: licenses of the sources (see `base_data/README.md`).*

## できること

| やりたいこと | コマンド（`models/hplc-retention/` で実行） | 時間の目安 |
|---|---|---|
| 公開の重み（`weights/`）を、同梱の基礎データから作り直す | `python train/run.py` | 約 10 分 |
| 自社のデータを足して学習し直す | `python train/run.py --own train/own_data` | 約 10 分 |
| 自社のデータだけで学習する／自社の行を重くする | `--own-only`／`--own-weight 5` を足す | 数分〜10 分 |
| 速い試し（木を減らす。精度は落ちる） | `python train/run.py --trees 60 --qtrees 30` | 約 5〜6 分（勾配の逆算の時間は変わらない） |
| 学習し直さずに、重みを自社のデータで評価する | `python train/scripts/eval_external.py <フォルダ>` | 1 分 |
| 測定値 1 点で合わせ直す規則を当てはめ直す（任意） | `python train/scripts/fit_calibration.py` | 数分 |
| 学習に使っていないカラムでの評価（任意） | `python train/scripts/eval_columns.py [--strict]` | 約 10 分 |
| 保持の式の定数（a, b, B）を基礎データから当てはめ直す（任意） | `python train/scripts/fit_form.py` | 1 分以内 |
| 基礎データを、元のデータから作り直す（任意） | `python train/scripts/build_from_raw.py --raw <フォルダ>` | 約 10 分 |
| グラフのモデル（GPU、任意） | [gpu/README.md](gpu/README.md) | GPU で 1〜2 時間 |

時間は、この道具を作った PC（24 スレッドの CPU、GPU なし）での実測です。

## 必要な物

- **Python 3**（この道具は 3.14.6 で確かめた）と、numpy・pandas・scipy・scikit-learn（[requirements.txt](requirements.txt)。`pip install -r train/requirements.txt`）。元のデータから作り直す時だけ openpyxl も。
- **Node.js**（この道具は v24 で確かめた。構造の特徴量の計算と、JavaScript との照合に使う）。RDKit.js はリポジトリに同梱（`vendor/rdkit/`）なので、入れる物はありません。
- メモリ 4 GB 以上。GPU は要りません（木のモデル）。
- Python と Node の場所は、環境変数で渡します（個人の環境のパスを道具に書かないため）：
  - `run.sh` は `PYTHON`（既定 `python`）を使います。`run.py` を直接動かす時は、使いたい Python で `python train/run.py` と打ちます。
  - Node は `NODE`（既定 `node`）。

```bash
# Mac・Linux（bash）
cd models/hplc-retention
NODE=node bash train/run.sh

# Windows（PowerShell）
cd models\hplc-retention
$env:PYTHONIOENCODING = "utf-8"      # 日本語の表示が化ける時
$env:NODE = "node"                   # node が PATH に無ければ、node.exe の場所
python train\run.py
```

## 手順の中身

`run.py` は、次の段を 1 つの Python の中で続けて実行します（途中のモデルを pickle に保存しないため。決まりは `docs/MODEL_LAYOUT.md`）。

| 段 | 中身 | 時間 |
|---|---|---|
| 1. 表を作る | 基礎データ（`base_data/`）と自社のデータを読み、1 行 = 1 測定の表と入力の行列（153 個）を作る。自社のデータの構造は Node で特徴量にする | 数秒〜1 分 |
| 2. 本体を学習 | 木 900 本（scikit-learn の HistGradientBoostingRegressor、向きの制約つき）。目的変数は実測から逆算した lkw。分析法ごとの補正（死時間の倍率・装置の遅れ）を、学習用の化合物だけで 3 回まわして決める | 約 3.5 分（勾配の保持時間 → lkw の逆算を含む） |
| 3. 評価 | 学習に使っていない化合物（InChIKey のハッシュで 10%）の保持時間の誤差、等組成（Weber）の化合物の 5 分割（木を 5 回学習） | 約 2.5 分 |
| 4. 書き出し | `model_free.json`（重み）、`reference_free.json`（JavaScript との照合用の値）、ストックの化合物 | 数秒 |
| 5. 予測の幅 | 構造の似た化合物のかたまりの 4 分割で外れを集め（木 900 本 × 4）、外れの大きさを分位点回帰（木 200 本 × 10）で当てる | 約 3.5 分 |
| 6. 印と索引 | 学習データに入っている化合物の幅の倍率、ストックの `seen`、名前 → 構造の索引 | 数秒 |
| 7. 照合（Node） | 名前の索引から RDKit.js で読めない構造を外す → `tests/check.cjs` で JavaScript と Python の一致 → 公開の重みとの比べ（既定の設定で基礎データだけの時） | 数秒 |

## できる物

`train/out/<名前>/`（名前は既定で `free`、自社のデータがあれば `own`。`--name` で変えられる）。`train/out/`・`train/work/` は `.gitignore` でリポジトリに入りません。

| ファイル | 中身 |
|---|---|
| `weights/model_free.json` | 本体の重み（`weights/` と同じ形・同じファイル名。中の `meta.name` が版の名前） |
| `weights/uncertainty_free.json` | 予測の幅の重み |
| `weights/stock_compounds.json`・`name_index_free.json`・`columns.json` | ストックの化合物（`seen` の印つき）、名前の索引、カラムの表（`weights/` の写し） |
| `weights/reference_free.json`・`reference_uncertainty_free.json` | JavaScript との照合用の値（この学習の表から選んだ 300 行など） |
| `REPORT.md` | 結果の要約（出所ごとの誤差、予測の幅に入る割合、自社のデータの前後の比べ） |
| `metrics.json`・`uncertainty_metrics.json`・`eval_heldout_compounds.csv` | 数字の全部 |
| `own_eval.json`・`own_data_report.json` | 自社のデータの評価と、読めなかった構造・外した行 |
| `compare_shipped.json` | 公開の重みとの比べ（基礎データだけの時） |

## 重みを入れ替える

1. 照合が合格していることを確かめる（`run.py` の最後に出る。あとからなら `node tests/check.cjs train/out/<名前>/weights`）。
2. `train/out/<名前>/weights/` の 5 つのファイル（`model_free.json`・`uncertainty_free.json`・`stock_compounds.json`・`name_index_free.json`・`columns.json`）を、ページが読む場所（このリポジトリなら `weights/`）に写す。照合用の `reference_*.json` は `tests/` 用。
   - **自社のデータを足した時は注意**：`name_index_free.json` には自社の化合物の名前と構造（出所の印 `own`）が、`reference_*.json` には自社の行の特徴量が入ることがあります。社内のフォルダに置くのはかまいませんが、公開のリポジトリ（GitHub のフォークなど）や社外に渡す物には入れないでください。名前の索引は、公開の物（このリポジトリの `weights/name_index_free.json`）のままでも動きます。
3. ページの側の呼び方は変わりません（`ChromRetention.create(model, columns, uncertainty)`）。`meta.name` が違っても、重みと予測の幅とストックの印の名前がそろっていれば動きます。

**公開の重みとして配るのは、基礎データだけで作った物です。** 自社のデータを足した重みは、下の「ライセンス」を見てください。

## 公開の重みの再現

基礎データだけで `python train/run.py` を実行すると、最後に `scripts/compare_shipped.cjs` が `weights/` と比べます。この道具を作った PC（numpy 2.4.6・scikit-learn 1.9.0）では、5 つのファイルとも **1 バイトも違わない物** ができました。別の版の scikit-learn では、木のしきい値の選び方がわずかに変わることがあります。その時も、照合用の値との差（lkw）が 1e-4 未満なら、実用上は同じです（公開の重みそのものでも、葉の値を小数 6 桁に丸めているので、Python との差は 2e-5 程度あります）。

## 評価の道具（任意）

- `scripts/eval_external.py <フォルダ> [--weights <重みのフォルダ>] [--cal train/work/free/calibration.json]`：学習し直さずに、自社のデータの形のフォルダで評価する。例：`python train/scripts/eval_external.py train/examples/zhang_dataset10`（公開の重みで、全部の化合物 4.5%、学習に使っていない化合物 7.7%）。
- `scripts/eval_columns.py [名前] [--strict] [--trees N]`：勾配の分析法をカラムで 5 つに分け、1 組ずつ外して学習し直す。`--strict` は、同じカラムの行をすべての出所（等組成・HSM2 の再構成・SMRT も）から外す。分け方は一般名のカラムの番号で決めるので、以前の非公開の作業の数字とは少し違います。
- `scripts/fit_calibration.py [名前]`：ページの「測定値 1 点で合わせ直す」（`retention_model.js` の `calibrate`）の規則の係数を、木のモデルの交差検証の外れから当てはめ直す。結果は `train/out/<名前>/calibration_rule.json`・`.md`。ページで使うには、`uncertainty_free.json` の `meta.calibration` に `{w0, col, solv, floor, tau}` を入れます（`retention_model.js` が読みます。無ければ既定の値）。

- `scripts/fit_form.py`：保持の式の定数（アセトニトリルは Weber らの等組成の実測、メタノールは Wiczling グループの勾配時間を変えた測定）を、基礎データから当てはめ直して、同梱の `base_data/retention_form.json` と比べる（書き換えない。使う時は入れ替えて学習し直す）。

## 元のデータから作り直す（任意）

同梱の基礎データ（`base_data/`）は、下の 6 つの公開データから作った表です。自分で取得した元のファイルから作り直せます（`scripts/build_from_raw.py`）。この道具を作った PC で、作り直した表と同梱の表は、分析法・測定・化合物・特徴量・名前の表まで同じになりました（`scripts/compare_base_data.py` で確かめられます）。

**取得は手で**（ブラウザで配布元のページを開いて、ライセンスを確かめてから）。自動で取りに行くスクリプトは付けていません（配布元によっては、自動の取得を断っているため）。置き方は `build_from_raw.py` の先頭の説明のとおり。

| 出所（フォルダ名） | 配布元・DOI | ライセンス | 取得するファイル |
|---|---|---|---|
| `smrt_figshare` | METLIN SMRT, figshare doi:10.6084/m9.figshare.8038913.v1（論文 doi:10.1038/s41467-019-13680-7） | Apache 2.0 | `SMRT_dataset.csv`・`SMRT_dataset.sdf` |
| `hsm2_mendeley` | HSM2 database, Mendeley Data doi:10.17632/hmd87ky8gh.1 | CC BY 4.0 | 「Download All」の zip の中の `HSM2_Database.xlsx` |
| `weber_isocratic_mendeley` | Weber, Retention data, Mendeley Data doi:10.17632/x292z4gd25（版 3） | CC BY 4.0 | 「Download All」の zip の中の `RetentionData` |
| `mcart_msd_zenodo` | MCaRT, Zenodo doi:10.5281/zenodo.17298529 | CC BY 4.0 | `MC-Retention.zip`（約 200 MB）の中の `enamine_480.csv` |
| `wiczling_osf` | Kubik ほか, OSF doi:10.17605/OSF.IO/ZQTJ7・QBV7J・2MCNW・Y6S8P・EVUJ9（5 記録） | CC BY 4.0 | 各記録の zip を、記録の番号（`zqtj7` など）のフォルダに展開 |
| `zhang2024_scidb` | Zhang ほか, Science Data Bank doi:10.57760/sciencedb.15823 | CC0 1.0 | `Dataset 01.xlsx` ほか 26 個（番号は `build_from_raw.py` の `ZHANG_COLUMN`）。Dataset 10 は評価専用（`examples/zhang_dataset10/`） |

継承つき（CC BY-SA）のデータ（RepoRT など）の保持値・構造、非営利限定のデータ、出所の分からないデータは、基礎データにも作り直しの手順にも入っていません。例外は 1 つで、SMRT の測定条件（カラム・勾配・流量）は、論文に書かれた条件を RepoRT の記録 0186 の記載を経由して読んだ事実の値です（`base_data/README.md`）。

## ライセンス

- この道具（`train/` のプログラム）：MIT（リポジトリの直下の `LICENSE`）。
- 基礎データ（`base_data/`）と例（`examples/`）：元のデータのライセンスのまま（[base_data/README.md](base_data/README.md)）。配る時は、出典の表示（`../ATTRIBUTION.md`）と、SMRT の Apache 2.0 の全文（`../LICENSES/Apache-2.0.txt`）を付けます。
- **基礎データ＋自社のデータで作った重み**：
  - 社内だけで使うなら、出典の表示の義務はありません（ライセンスの義務は、配る時に生じます）。
  - 外に配るなら、`../ATTRIBUTION.md` を付けて、「どの公開データに、何を足したか」を書き足してください（CC BY 4.0・Apache 2.0 の表示の条件）。自社のデータの部分をどう扱うかは、データの持ち主が決めます。
- **自社のデータだけで作った重み**（`--own-only`）：公開データは特徴量の参照にしか使っていません。扱いは持ち主が決めます。

これは事実の整理で、法律の意見ではありません。

## 分かっている癖（再現のために、そのままにしてある物）

- HSM2 の表で同じ名前の行が 2 つあるカラムが 1 本あり（`hsm2_739` と `hsm2_681`）、学習ではカラムの性質の値（H〜D・logkEB）に、先の行（`hsm2_681`）の値が入っています（`base_data/methods.csv` の `column_match` が `exact_same_name_first`）。
- Wiczling グループの元の表で、名前と構造が合わない 4 化合物（METID 194・289・301・87）の行は、学習の表にそのまま入っています（名前の索引からは外してあります）。
- SMRT の標準の SMILES には、元の SD ファイル由来の原子の番号（`[CH3:0]` のような `:0`）が付いた物があります。特徴量の計算には影響しません。
- 同梱の特徴量は、学習の時と同じ値です。標準の SMILES から `featurize.cjs` で計算し直すと、80,988 化合物のうち 42 化合物で記述子の 1 つが 0.0001 だけ違います（元は SD ファイルの構造から計算したため。小数 4 桁の丸めの境目）。

## ファイル

```
train/
  README.md        この説明
  OWN_DATA.md      自社のデータの形と足し方
  run.py・run.sh   入口
  requirements.txt Python のパッケージ
  lib/             部品（表・特徴量・保持の式・学習・書き出し）
  scripts/         Node の道具（featurize.cjs・filter_name_index.cjs・compare_shipped.cjs）と、任意の評価・作り直しの道具
  base_data/       基礎データ（6 つの公開データから作った表）
  examples/        自社のデータの形の例（Zhang らの Dataset 10 を書き直した物）
  own_data/        自社のデータの置き場（README.md 以外は .gitignore でリポジトリに入らない）
  gpu/             グラフのモデル（GPU、任意）
```
