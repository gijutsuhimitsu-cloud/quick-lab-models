# 学習し直す道具（溶けるか・混ざるか）

*English summary: retraining kit. `node train/run_all.cjs` rebuilds the training tables from the bundled processed public data (train/base_data, CC BY 4.0), trains the four tree models with scikit-learn on CPU (about 1.5 min), refits the mixed-solvent correction and checks that the JavaScript predictions match Python. Add `--own <csv>` to include your own data (see OWN_DATA.md); it stays on your machine. Everything needed is bundled (including the melting-point table), so the shipped weights are reproduced exactly without downloading anything. A full rebuild from the original datasets is documented below (download links and DOIs); AqSolDB must be downloaded by hand. Node and Python are found through the PYTHON / NODE environment variables (defaults: python / node).*

## 必要な物

| 物 | 版の目安 | 使う所 |
|---|---|---|
| Node.js | 18 以上 | 表を作る（構造を読むのはリポジトリの `vendor/rdkit` の RDKit.js。入れる物はない） |
| Python | 3.9 以上（確かめたのは 3.14） | 学習 |
| numpy・pandas・scikit-learn・joblib | 確かめたのは scikit-learn 1.9.0・numpy 2.4・pandas 3.0 | 学習 |
| openpyxl | 3.1 | 元のデータ（xlsx）から作り直す時だけ |

- Python と Node の場所は環境変数 `PYTHON`・`NODE` で渡します（既定は `python`・`node`）。例：`PYTHON=/opt/conda/bin/python node train/run_all.cjs`（Windows の PowerShell では `$env:PYTHON = "C:\…\python.exe"`）。
- コマンドは、どれも `models/solubility/` で走らせます。
- 途中の物は `train/work/`、できた重みは `train/out/` に出ます（どちらも `.gitignore` でリポジトリに入りません。場所は `SOL_WORK`・`SOL_OUT` で変えられます）。

## いちばん短い手順

```
node train/run_all.cjs --quick                                   # 流れだけ確かめる（30 秒ほど。木を減らすので値は本番と違う）
node train/run_all.cjs                                           # 公開データだけで学習し直す（3 分ほど）
node train/run_all.cjs --own train/own_data/my_data.csv          # 自社のデータを足す（形は OWN_DATA.md）
node train/run_all.cjs --own train/own_data/my_data.csv --install   # できた重みを weights/ に入れる
```

`run_all.cjs` がする事（1 つずつ走らせてもよい）：

| 順 | コマンド | 中身 | 時間の目安 |
|---|---|---|---|
| 1 | `node train/build_tables.cjs [--own …]` | `base_data/` と `weights/` の係数から、学習用の表 `work/train_25C.csv`（有機溶媒用）・`train_25C_aq.csv`（水用）を作る | 15 秒 |
| 2 | `python train/train_gbt.py` | 4 つのモデル（融点あり／なし × 有機溶媒／水）を、溶質を分けた 5 分割の検証つきで学習し、`out/solubility_gbt.json` と照合用の値 `out/reference_gbt.json` を書く | 1.5 分（CPU 5 並列） |
| 3 | `python train/train_tm.py` | 融点の推算モデル（`base_data/tm_train.csv.gz` から。`--skip-tm` で飛ばせる） | 30 秒 |
| 4 | `node train/fit_mix.cjs` | 混合溶媒の補正の係数（`out/mix_fit.json`） | 1 秒 |
| 5 | `node tests/check.cjs train/out` | できた重みで、JS の計算が scikit-learn と合うかを確かめる | 2 秒 |

時間は 24 コアの Windows の計算機で測った値です。`train_gbt.py --jobs N` で並列の数を変えられます。

- `--install` は、照合に合格した重み（`solubility_gbt.json`・`tm_model.json`・`mix_fit.json`）を `weights/` へ、照合用の値を `tests/` へ写します。元に戻すときは、git で `weights/` を戻してください。
- 学習し直した重みの `meta.version` は、自社のデータを足すと `0.1+own` になります。

## 融点について

公開の重み（版 0.1）は、BigSolDB に融点がない固体に、融点の公開データ（Jean-Claude Bradley Open Melting Point Dataset、CC0 1.0）の値を足して学習しました。その値は同梱しています：

- `base_data/bigsol_solutes.csv.gz`・`aqsol_solutes.csv.gz` の `TmC`（`tmSrc` が `sol` = BigSolDB、`mp` = 融点の公開データ）
- `base_data/tm_train.csv.gz`（融点のモデルの学習用の表：融点と構造の記述子 97 列、19,287 構造）
- `weights/tm_lookup.json`（融点の照合の表。`build_tables.cjs` は、土台の表に融点がない溶質に、この表の値があれば足す）

出所の欄が、利用条件が非営利の薬のデータベース、または試薬会社のカタログである行（1,707 行）は、元の表から読んでいません（行の番号は `rebuild/bradley_excluded_keys.txt`。外したい行が増えたら、ここに足して下の「元のデータから全部を作り直す」）。

`node train/run_all.cjs` で、公開の重みと**木が 1 本残らず同じ**モデルができることを確かめました（2026-11-25、scikit-learn 1.9.0。版が違うと少し変わることがある）。

## 元のデータから全部を作り直す

ふだんは要りません（`base_data/` の表が、下の手順の結果です）。元のデータを取り直して確かめたい時、新しい版のデータに差し替えたい時に使います。元のデータは `train/raw/` に置きます（`.gitignore` で入らない。置き方は [raw/README.md](raw/README.md)）。

| データ | 取る所 | ライセンス | 置く所 |
|---|---|---|---|
| BigSolDB 2.2 | https://zenodo.org/records/22648301 （doi:10.5281/zenodo.22648301）の `BigSolDBv2.2.csv`・`_gamma.csv`・`_densities.csv` | CC BY 4.0 | `raw/bigsoldb/` |
| MixtureSolDB | https://zenodo.org/records/18660057 （doi:10.5281/zenodo.18660057）の `MixtureSolDB.csv`・`solvent_boiling_point.csv` | CC BY 4.0 | `raw/mixturesoldb/` |
| AqSolDB | Harvard Dataverse の記録ページ（doi:10.7910/DVN/OVHAW8）から `curated-solubility-dataset.csv` を**ブラウザで手で取る** | CC0 1.0 | `raw/aqsoldb/` |
| 融点の公開データ | figshare（doi:10.6084/m9.figshare.1031637）の `BradleyMeltingPointDataset.xlsx` | CC0 1.0 | `raw/bradley_mp/` |
| Abraham 記述子の表 | figshare（doi:10.6084/m9.figshare.1176994）の `20140922AbrahamDescriptorsBackUp.xls` | CC BY 4.0 | `raw/abraham_descriptors/` |
| Abraham の溶媒の係数 | figshare（doi:10.6084/m9.figshare.25898476）の `20240524 Abraham Solvent Master Dataset.xlsx` | CC0 1.0 | `raw/abraham_solvents/` |

このリポジトリにデータを取りに行く道具はありません（配布元の利用条件と取り方の決まりに従って、使う人が取ってください）。

作り直しの順（結果は `train/out/` と `train/out/base_data/` に出る。今ある `weights/`・`base_data/` は書き換えない）：

```
python train/rebuild/abraham_xls_to_csv.py          # Abraham 記述子の表 → work/abraham_solutes_raw.csv
node   train/rebuild/fit_abraham.cjs                # → out/abraham_solute_model.json、out/base_data/abraham_solutes_clean.csv.gz
python train/rebuild/abraham_solvents.py            # 溶媒の係数 → work/solvents_abraham.json
python train/rebuild/bradley_to_csv.py              # 融点 → work/bradley_mp_raw.csv
node   train/rebuild/bradley_tm.cjs                 #      → out/tm_lookup.json、out/base_data/tm_train.csv.gz
# ここから先は、作り直した表と係数を読むように場所を渡す（bash の書き方。PowerShell は $env:SOL_BASE = "train/out/base_data" など）
export SOL_BASE=train/out/base_data SOL_WEIGHTS=train/out
node   train/rebuild/bigsol_from_raw.cjs            # → out/formula.json、out/base_data/bigsol_25C.csv.gz・bigsol_solutes.csv.gz
node   train/rebuild/aqsol_from_raw.cjs             # → out/base_data/aqsol_rows.csv.gz・aqsol_solutes.csv.gz
node   train/make_tables.cjs                        # → out/solvents.json・lookup_25C.json・example_compounds.json
node   train/rebuild/mix_from_raw.cjs               # → out/base_data/mix_series_25C.csv.gz
node   train/fit_mix.cjs                            # → out/mix_fit.json
node   train/build_tables.cjs && python train/train_gbt.py && python train/train_tm.py
```

- `bradley_tm.cjs` は、溶解度の学習の溶質の表（`base_data/bigsol_solutes.csv.gz`）で印を付けるので、`SOL_BASE` を渡す前（同梱の表）に走らせます。
- 2026-11-25 にこの順で作り直し、同梱の `base_data/` の 7 つの表と、`weights/` の表（formula・mix_fit・solvents・lookup_25C・tm_lookup・example_compounds・abraham_solute_model。作った日付を除く）が、同梱の物と同じになることを確かめました。
- 元のデータを取らずに、一部だけ作り直すこともできます：`bigsol_from_raw.cjs`・`aqsol_from_raw.cjs` は、`work/tm_lookup.json` がなければ、`SOL_WEIGHTS` のフォルダ（既定は同梱の `weights/`）の `tm_lookup.json` の融点を使います。
- どの表にも、溶質の名前は写しません（構造だけ）。
- 作り直した物を使う時は、`train/out/` の JSON を `weights/` へ、`train/out/base_data/` を `train/base_data/` へ写します。

## ファイル

| 場所 | 中身 |
|---|---|
| `run_all.cjs` | 学習し直しの流れを 1 回で |
| `build_tables.cjs` | 学習用の表（自社のデータの足し込みも） |
| `train_gbt.py`・`train_tm.py` | 学習と書き出し（JSON。pickle は使わない） |
| `fit_mix.cjs`・`make_tables.cjs` | 混合溶媒の補正、溶媒の表・照合の表・例の固体 |
| `lib/common.cjs`・`lib/own_data.cjs` | 共通の部品、自社のデータの読み取り |
| `rebuild/` | 元のデータから作り直す道具 |
| `base_data/` | 学習の土台の表（CC BY 4.0、[base_data/README.md](base_data/README.md)） |
| `examples/` | 自社のデータの形の例（[examples/README.md](examples/README.md)） |
| `own_data/` | 自社のデータの置き場（リポジトリに入らない） |

## GPU について

このモデルは決定木（勾配ブースティング）なので、CPU で数分で学習できます。GPU は要りません。データが大きく増えて時間がかかるようになったら、`train_gbt.py` の `GradientBoostingRegressor` を LightGBM や XGBoost の GPU 版に差し替えられます。ただし、書き出し（`export`）と JS の計算（`js/solfeat.js` の `evalModel`）は scikit-learn の木の形（左の子が親の次の番号）を前提にしているので、差し替えた時はその 2 つも合わせて書き直し、`node tests/check.cjs train/out` で照合してください。このリポジトリでは差し替えていません。

## 自社のデータを足した重みの扱い

- 自社のデータで作った重みはあなたの物です。社内だけで使うなら公開の義務はありません。
- 配る時は、公開データの分の出典の表示（[../ATTRIBUTION.md](../ATTRIBUTION.md)）を残し、学習し直したことが分かるようにしてください（CC BY 4.0）。
- 木のモデルの葉の値には、学習に使った実測の情報が混ざります。自社のデータの秘密が大事な時は、重みを社外に出す前によく考えてください。
