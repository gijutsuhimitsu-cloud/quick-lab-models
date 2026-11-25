# 自社のデータを足して学習し直す

*English summary: put two CSV files (`methods.csv`, `measurements.csv`) in `train/own_data/`, run `python train/run.py --own train/own_data`, and read `train/out/own/REPORT.md`, which compares the published weights with your retrained weights on your held-out compounds. Your data stays on your computer; nothing is uploaded. Weights you make are yours to keep internal.*

## 大事なこと

- **データはこの PC の中だけで扱います。** どこにも送りません（道具はネットにつなぎません）。`train/own_data/`・`train/work/`・`train/out/` は `.gitignore` でリポジトリに入りません。
- できた重みは、社内だけで使ってかまいません。外に配る時だけ、出典の表示が要ります（`README.md` の「ライセンス」）。
- できた名前の索引（`name_index_free.json`）と照合用の値（`reference_*.json`）には、自社の化合物の名前・構造・特徴量が入ることがあります。`weights/`・`tests/` に写したあとで、公開のリポジトリや社外に出さないでください（`README.md` の「重みを入れ替える」）。
- 精度が上がるのは、主に**自社の分析法（カラム・溶媒・勾配・装置）**と、**自社でよく測る化合物に似た分子**です。測っていない条件・似ていない分子の精度は、あまり変わりません。

## ファイルの形

`train/own_data/` に 2 つの CSV（UTF-8、1 行目は見出し）を置きます。例は [examples/zhang_dataset10/](examples/zhang_dataset10/)。

### methods.csv（1 行 = 1 つの分析法）

| 列 | 要る？ | 中身 |
|---|---|---|
| `method_id` | 要る | 分析法の名前（英数字。重ならないように） |
| `column` | 要る | カラム。次のどれか：`weights/columns.json` の `id`（例 `hsm2_3`）か一般名（例 `C18 No.3`）／固定相の種類（`C18`・`C8`・`Phenyl`・`CN`・`EP`・`F`・`short`。その種類のカラムの性質の中央値を使う）／`custom`（下の H〜D を自分で書く） |
| `organic` | 要る | 有機溶媒：`acn`（アセトニトリル）か `meoh`（メタノール） |
| `phi` | どちらか 1 つ | 等組成の時の有機溶媒の割合（0〜1。50% なら 0.5） |
| `gradient` | どちらか 1 つ | 勾配の時の表。`時刻（分）:有機溶媒の割合（0〜1）` を `;` でつなぐ。例 `0:0.05;10:0.95;12:0.95`。注入の時刻を 0 に、装置の遅れは入れない。移動相 A にも有機溶媒が入っている時は、全体の割合で書く |
| `t0_min` | 勾配では要る | 死時間（分。ウラシルなど保持されない物の溶出時間）。空なら `length_mm`・`id_mm`・`flow_mL_min` からの目安（空隙率 0.65） |
| `dwell_min` | なくてよい | 装置の遅れ（分）。書くと、その値で固定する。空なら、学習の中で死時間の倍率（0.6〜1.2）と遅れ（0〜2 分）を、学習用の化合物の保持時間に合うように選ぶ |
| `pH` | なくてよい | 水側の緩衝液の pH。空なら 3.0 として扱う（酸・塩基の化合物の予測が変わる） |
| `temp_C` | なくてよい | カラムの温度（°C） |
| `length_mm`・`id_mm`・`particle_um`・`flow_mL_min` | なくてよい | カラムの寸法と流量（死時間の目安にだけ使う） |
| `additive` | なくてよい | 添加剤の言葉（記録だけ。いまのモデルは使わない） |
| `note` | なくてよい | メモ |
| `H`・`A`・`B`・`C`・`V`・`D`・`kEB`（か `logkEB`）・`phase` | `custom` の時 | 自社のカラムの HSM（疎水性差し引きモデル）のパラメータと固定相の種類 |

### measurements.csv（1 行 = 1 つの測定）

| 列 | 要る？ | 中身 |
|---|---|---|
| `method_id` | 要る | methods.csv の `method_id` |
| `smiles` | 要る | 構造（SMILES）。塩は、測っている形（ふつうは中性の親化合物）で |
| `rt_min` | 勾配では要る | 保持時間（分）。等組成でも、`t0_min` があれば使える |
| `log_k`（か `k`） | 等組成の時 | 保持係数の log10（か保持係数そのもの）。書けば `rt_min` より優先 |
| `name` | なくてよい | 化合物の名前。書くと、できた重みの名前の索引に入る |

- 同じ分析法・同じ化合物の繰り返しは、そのまま何行でも書いてかまいません。
- 保持時間が死時間以下の行（保持されていない物）は、学習には使われませんが、評価には入ります。
- RDKit.js で読めない構造の行は外し、`train/out/<名前>/own_data_report.json` に書きます。

## 学習の仕方

`models/hplc-retention/` で：

```bash
python train/run.py --own train/own_data                     # 基礎データ＋自社のデータ（ふつうはこれ）
python train/run.py --own train/own_data --own-weight 5      # 自社のデータの行を 5 倍の重みで（自社の分析法に寄せたい時）
python train/run.py --own train/own_data --own-only          # 自社のデータだけ（行が数千以上ある時。少ないと、ほかの条件の予測が崩れる）
python train/run.py --own train/own_data --name lab2025      # 出力の名前を変える（train/out/lab2025/）
```

木を減らした速い試し（`--trees 60 --qtrees 30`）で形を確かめてから、本番（既定の 900 本）を回すと楽です。

## 結果の確かめ方

`train/out/<名前>/REPORT.md` の「自社のデータ」の表を見ます。

- **前**：公開の重み（`weights/model_free.json`）で、補正なし（死時間は表の値、遅れ 0）の予測。
- **後**：いま作った重みで、学習の中で決めた補正つきの予測。
- 比べるのは「評価用の化合物」＝どの学習にも使っていない化合物（InChIKey のハッシュで 10%。公開の重みの学習でも外してある化合物）。学習に使った化合物の値（`own_eval.json` の `all`）は、良く見えすぎます。
- 評価用の化合物が少ない時（数十個以下）は、数字のぶれが大きいので、目安にとどめます。
- 照合（`tests/check.cjs`）が「合格」なら、ページの JavaScript でも同じ予測になります。

学習し直さずに、公開の重みが自社のデータでどのくらい当たるかだけを見たい時は：

```bash
python train/scripts/eval_external.py train/own_data
```

## 例（公開データで作った、自社のデータの形の例）

[examples/zhang_dataset10/](examples/zhang_dataset10/) は、Zhang ほかの公開データ（CC0）の Dataset 10（1 つの勾配の分析法、335 測定）を、上の形に書き直した物です。**本当の自社のデータではありません。** 公開の重みの学習には入っていない（評価専用の）ファイルなので、「自社の分析法を 1 つ足す」状況の例になります。

```bash
python train/run.py --own train/examples/zhang_dataset10 --name example
```

この道具を作った PC での結果は、`examples/zhang_dataset10/README.md` に書きました。

## 重みを使う・配る

- 入れ替え方は `README.md` の「重みを入れ替える」。ファイル名は `weights/` と同じです。
- 社内だけで使うなら、出典の表示の義務はありません。外に配る時は `../ATTRIBUTION.md` を付けて、「どの公開データに、何を足したか」を書き足してください。
