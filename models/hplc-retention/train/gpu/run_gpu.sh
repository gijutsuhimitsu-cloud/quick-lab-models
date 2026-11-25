#!/bin/bash
# SPDX-License-Identifier: MIT
# GPU の計算機の中で回す: 本番用の「蒸留した小さいグラフのモデル」（この写しでは未実行。元の作業では L4 1 枚で約 60〜70 分）。
#   1. 教師（幅 256・5 層）を、その版の全部の行で 2 本学習
#   2. 教師 2 本の平均で「分子 × 条件」の組（gen_aug.py の aug.npz）に lkw を付ける
#   3. 小さいモデル（幅 128・4 層、重み 53 万個＝2.1 MB）を、実測の行＋組 25 万行（重み 0.5）で学習して、ブラウザ用に書き出す
# 使い方（DATA のフォルダに rows.npz・meta.json・compounds.tsv・graphs_js.npz・aug.npz を置いてから）:
#   GROUP=free DATA=gpu_data RUNS=runs bash run_gpu.sh
# できる物: $RUNS/final/（weights.bin・manifest.json・check.json）。手元に持ち帰って make_gnn_export.cjs に渡す。
set -e
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
G="${GROUP:-free}"; D="${DATA:-gpu_data}"; O="${RUNS:-runs}"
C=(--data "$D" --group "$G" --model gnn --graphs graphs_js.npz --select last)   # 配列にして、空白を含むパスでも分かれないように
mkdir -p "$O/teach"
for s in 0 1; do
  [ -f "$O/teach_s$s/result.json" ] || "$PY" train_teacher.py "${C[@]}" --all 1 --seed $s --epochs 40 --out "$O/teach_s$s"
done
"$PY" train_teacher.py "${C[@]}" --all 1 --aug "$D/aug.npz" --predict "$O/teach_s0,$O/teach_s1" --out "$O/teach"
"$PY" train_teacher.py "${C[@]}" --all 1 --seed 0 --hid 128 --layers 4 --head 256 --epochs 30 --aug "$D/aug.npz" --aug_y "$O/teach/aug_pred.npy" \
  --aug_n 250000 --aug_w 0.5 --export 1 --out "$O/final"
date > "$O/DONE_final"
