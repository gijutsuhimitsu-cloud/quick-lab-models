# グラフのモデル（GPU、任意）

*English summary: optional path to train the small graph neural network (GINE, width 128 × 4 layers, 2.1 MB float32) that `js/gnn_model.js` runs in the browser, distilled from two larger teacher GNNs. Needs a CUDA GPU and PyTorch (no Python RDKit: graphs are built on your PC with the bundled RDKit.js, exactly as in the browser). Our runs on one NVIDIA L4 took about 30–120 minutes per step. **This copy of the GPU scripts has not been run** (only syntax-checked); the CPU data-export steps were run and checked.*

**この写しでは未実行**：GPU の上の学習（`train_teacher.py`・`run_gpu.sh`・`run_gpu_cv.sh`）と、その結果を使う道具（`make_gnn_export.cjs`・`check_gnn.cjs`・`train_uncertainty_gnn.py`）は、この道具を作った PC に PyTorch が無いため動かしていません（Python の構文の確かめだけ）。元の非公開の作業で、同じ作りの物を GPU（NVIDIA L4 1 枚）で動かしたことはあります。手元の CPU の段（下の 1〜3）は、この写しで動かして確かめました。

## 何が良くなるか

木のモデル（`model_free.json`）に、分子ごとの表現 `z` を作る小さなグラフのモデルを足すと、学習に使っていない分子の精度が上がります。元の作業の交差検証（化合物のかたまりの 4 分割、継承つきのデータを含む版）では、保持時間の誤差が SMRT 4.7 → 3.1%、MCaRT 7.3 → 6.0%、Zhang 6.5 → 5.4%（Wiczling 10.0 → 9.8% は差なし）。**基礎データだけの版での数字は、まだ測っていません**（このフォルダの手順で GPU を回すと、`run_gpu_cv.sh` の結果で測れます）。

## 必要な物

- CUDA の GPU（元の作業では NVIDIA L4、24 GB）と PyTorch 1.13 以上（CUDA 版。途中の `best.pt` を `weights_only=True` で、重みだけとして読むため）。numpy。
- Python の RDKit は要りません。グラフは手元の PC で、ブラウザと同じコード（`js/graph_features.js` ＋ 同梱の RDKit.js）で作ります。
- 先に、同じ名前で `python train/run.py` を済ませておく（分析法ごとの補正・予測の幅の 4 分割・木の重みを使います）。

## 手順

`models/hplc-retention/` で（名前は `free` の例）：

| 段 | どこで | コマンド | 中身 | 時間の目安 |
|---|---|---|---|---|
| 1 | 手元 | `python train/gpu/export_graph_data.py free` | 学習の表を `train/work/free/gpu_data/` に書き出す（`rows.npz`・`meta.json`・`compounds.tsv`） | 1〜2 分 |
| 2 | 手元 | `node train/gpu/graphs_js.cjs train/work/free/gpu_data` と `python train/gpu/graphs_js_to_npz.py train/work/free/gpu_data` | 全部の化合物の分子グラフ（原子 8 個・結合 3 個の番号の特徴量）→ `graphs_js.npz` | 数分 |
| 3 | 手元 | `python train/gpu/gen_aug.py train/work/free/gpu_data 12 free all` | 蒸留の「分子 × 条件」の組（化合物あたり 12。学習データの分析法 2/3 ＋ HSM2 のカラム × ランダムな pH・溶媒・温度 1/3）→ `aug.npz` | 1 分 |
| 4 | GPU | `train/gpu/` の中身と `gpu_data/` を GPU の計算機に写し、`GROUP=free DATA=gpu_data RUNS=runs bash run_gpu.sh` | 教師（幅 256・5 層）2 本 → 組の予測 → 小さいモデル（幅 128・4 層）の蒸留と書き出し（`runs/final/`） | 元の作業では L4 で約 60〜70 分 |
| 5 | GPU（任意） | `GROUP=free DATA=gpu_data RUNS=runs bash run_gpu_cv.sh` | 化合物の 4 分割の交差検証（精度の数字と、グラフ用の予測の幅に要る） | 元の作業では L4 で約 2 時間 |
| 6 | 手元 | `node train/gpu/check_gnn.cjs <runs/final>` | JavaScript の推論と PyTorch の予測の照合（元の作業では lkw の差 2e-6） | 数秒 |
| 7 | 手元 | `node train/gpu/make_gnn_export.cjs <runs/final> train/out/free/weights` | `gnn_free.json`・`gnn_free.bin`・`gnn_free_weights.json`・`stock_z_free.json` を重みのフォルダに置き、pH の向きなどを確かめる | 数秒 |
| 8 | 手元（任意） | `python train/gpu/train_uncertainty_gnn.py gSD free <gpu_data> <runs> train/out/free/weights` | グラフ用の予測の幅 `uncertainty_free_gnn.json`（段 5 の結果から） | 数分 |

元の作業の GPU の時間（L4 1 枚）：教師 1 本 7〜30 分、本番用の手順（段 4）約 60〜70 分、交差検証（段 5）約 2 時間。1 回の段はおおむね 30〜120 分でした。お金のかかる計算機を借りる時は、先に費用の見込みを確かめてください。

## ページでの使い方

`js/gnn_model.js`（推論）と `js/graph_features.js`（RDKit.js の分子 → グラフ）が読みます。

```js
var gnn = ChromGNN.create(gnnManifest /* gnn_free.json */, gnnWeightsArrayBuffer /* gnn_free.bin */);
var P = ChromRetention.create(model, columns, uncertaintyGnn /* uncertainty_free_gnn.json か uncertainty_free.json */, gnn);
// ストックの化合物: 事前計算の z を付ける（鍵は InChIKey の全体）
stock.compounds.forEach(function (c) { var z = stockZ.z[String(c.inchikey || c.name)]; if (z) c.z = Float32Array.from(z); });
// 利用者が足した分子: P.embed(compound, ChromGraph.fromSmiles(RDKit, smiles)) で z を 1 回だけ計算する
```

- 配るのは `gnn_free.json`・`gnn_free.bin`（数値だけ）と `stock_z_free.json` です。GPU の上で途中に保存される `best.pt` などの PyTorch の点検用のファイルは、中身が pickle の形式なので配りません（`docs/MODEL_LAYOUT.md` の決まり。手元で読み直す時も、信頼できる自分のファイルだけを読む）。
- `z` の無い化合物は、自動で木のモデルになります。グラフのモデルの `group`（`--group` の値）と、木の重みの `meta.name` が同じでないと、`create` が止まります。
- 出口の形は木のモデルと同じ物理の式で、酸・塩基の pH の向きは必ず守られます（イオン化の量の係数が正になる作り）。

## ファイル

| ファイル | 中身 | どこで |
|---|---|---|
| `export_graph_data.py` | 学習の表の書き出し | 手元の CPU（この写しで実行した） |
| `graphs_js.cjs`・`graphs_js_to_npz.py` | 分子グラフ（ブラウザと同じコード） | 手元（実行した） |
| `gen_aug.py` | 蒸留の組 | 手元（実行した） |
| `train_teacher.py` | 教師・小さいモデルの学習、組の予測、ブラウザ用の書き出し（`--export 1`） | GPU（未実行） |
| `run_gpu.sh`・`run_gpu_cv.sh` | GPU の上の手順 | GPU（未実行） |
| `check_gnn.cjs`・`make_gnn_export.cjs` | JavaScript との照合、重みのフォルダへの書き出し | 手元（未実行。GPU の結果が要る） |
| `train_uncertainty_gnn.py` | グラフ用の予測の幅 | 手元（未実行。段 5 の結果が要る） |
