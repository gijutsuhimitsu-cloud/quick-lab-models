# 自社のデータの形と、足し方

*English summary: put a CSV in train/own_data/ with one measurement per row: solute (SMILES, or a name found in the public tables), solvent (an id/name from weights/solvents.json, or its SMILES), temperature in °C (20–30 °C rows are used; the model is for 25 °C), the solubility value and unit (mg/mL, g/L, mg/L, µg/mL, mol/L, mmol/L or log10(mol/L)), and optionally the melting point. Then run `node train/run_all.cjs --own train/own_data/your.csv` (modes: merge / weighted / only). Your data never leaves your machine.*

## CSV の形

1 行に 1 つの測定。文字コードは UTF-8、1 行目は列の名前。列の順番は自由です。例：[examples/own_data_example.csv](examples/own_data_example.csv)。

| 列 | 要るか | 中身 |
|---|---|---|
| `solute_smiles` | どちらか 1 つ | 溶質の構造（SMILES）。塩・混合物（「.」を含む）と、重い原子 80 個を超える物は読みません |
| `solute_name` | 〃 | 名前だけで引けるのは、例の固体 37 種の日本語の名前（`weights/example_compounds.json`）だけです（公開データの表は溶質の名前を持たない）。社内の名前なら `solute_smiles` も書いてください（名前は記録に残るだけ） |
| `solvent` | どちらか 1 つ | 溶媒。`weights/solvents.json` の `id`（`water`・`meoh`・`acn`・`dmso`・`butan_1_ol` など）、英語の名前（`Methanol`）、日本語の名前（`メタノール`）、略号（`MeCN`）のどれか |
| `solvent_smiles` | 〃 | 溶媒の構造。表の 122 溶媒のどれかと同じ構造なら使えます（表にない溶媒は、係数がないので使えません） |
| `temperature_C` | 要る | 測った温度 °C。**20〜30 °C の行だけを使います**（このモデルは 25 °C。外の行は理由を出して飛ばします） |
| `value` と `unit` | どちらか | 溶解度の値と単位：`mg/mL`・`g/L`・`mg/L`・`µg/mL`（`ug/mL` でも可）・`mol/L`（`M`）・`mmol/L`（`mM`）・`log10(mol/L)` |
| `logS_molL` | 〃 | log10(mol/L) で書く時はこの列だけでよい |
| `Tm_C` | 任意 | 融点 °C。書くと融点ありのモデルに入り、当たりが良くなります。公開データの融点（BigSolDB・融点の照合の表）より優先します。書かなければ、照合の表にある構造は自動で融点が入ります |
| `note` など | 任意 | ほかの列は読みません（メモに使ってください） |

- 濃度は「飽和溶液 1 L あたり」で書いてください。g/100 g 溶媒やモル分率のデータは、密度や分子量で mg/mL に直してから入れます。
- 同じ溶質 × 溶媒の行が何行あっても、log10(mol/L) で平均して 1 点にします（公開データと同じ）。
- 読めなかった行は、行番号と理由を表示し、`train/work/own_report.json` にも書きます。

## 足し方

```
node train/run_all.cjs --own train/own_data/my_data.csv                            # 足す（既定：merge、重み 1）
node train/run_all.cjs --own train/own_data/my_data.csv --own-mode weighted        # 自社の行を重く（既定 5 倍）
node train/run_all.cjs --own train/own_data/my_data.csv --own-mode weighted --own-weight 10
node train/run_all.cjs --own train/own_data/my_data.csv --own-mode only            # 自社のデータだけで学習
node train/run_all.cjs --own a.csv --own b.csv                                     # 何ファイルでも
```

| モード | 学習に使う行 | 向く時 |
|---|---|---|
| `merge`（既定） | 公開データ ＋ 自社のデータ（自社の行の重み 1） | まず試す時。公開データの溶質と同じ物があっても、両方の行を残します |
| `weighted` | 同じ。自社の行の重みを `--own-weight`（既定 5） | 自社の化合物の系列（似た骨格）に寄せたい時 |
| `only` | 自社のデータだけ | 数百行以上あり、公開データと性質がかなり違う時。少ないと当たりません（行が 20 未満か溶質が 5 種未満のモデル（例：融点の分かる溶質が少ない時の融点ありのモデル、水の行が少ない時の水のモデル）は作らず、推論ではある方のモデルで代わりにします） |

- 水の行は水のモデル（`tm_w`・`notm_w`）に、ほかの溶媒の行は有機溶媒のモデル（`tm`・`notm`）に入ります（水の行は有機溶媒のモデルにも入る。公開の重みと同じ作り）。
- 自社の溶質も、検証では溶質ごとに組に分けます。`train/out/gbt_cv.json` の `own_organic`・`own_water` が、**学習に使っていない自社の溶質での外れ**です。幅（`organic`・`water`）は公開データの行で測った値のままです（`only` の時は自社の行で測った値）。
- できた重みは `train/out/`。`node tests/check.cjs train/out` が自動で走り、JS と Python の計算が合うかを確かめます。良ければ `--install` を付けて走らせ直すか、`train/out/` の JSON を `weights/` に写します。

## データはどこにも出ません

- `train/own_data/`・`train/work/`・`train/out/` は `.gitignore` でリポジトリに入りません。
- 学習の道具は、手元のファイルを読んで手元に書くだけです。外へ送る処理はありません。
- ただし、できた重み（`train/out/solubility_gbt.json`）の葉の値には、自社のデータの情報が混ざります。重みを社外に渡す時は、それでよいかを確かめてください。

## 例で試す

```
node tests/check_own_data.cjs                                                    # 例のファイルの読み取りだけ（Python 不要、数秒）
node train/run_all.cjs --own train/examples/own_data_example.csv --quick         # 流れを確かめる（30 秒ほど）
```

例のファイルは、公開データ（MixtureSolDB の純粋な溶媒の点、CC BY 4.0）のうち、学習の土台の表に入っていない溶質 11 種・50 行を、自社のデータに見立てて名前を「例-01」などに置き換えた物です（[examples/README.md](examples/README.md)）。40 °C の 3 行は、温度の外で飛ばされる例として残しています。
