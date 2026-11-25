#!/bin/bash
# SPDX-License-Identifier: MIT
# GPU の計算機の中で回す（任意）: 化合物のかたまりの 4 分割での交差検証。グラフのモデルの精度の数字と、予測の幅（train_uncertainty_gnn.py）に要る。
# この写しでは未実行。元の作業では L4 1 枚で約 2 時間。aug.npz は gen_aug.py に all を付けて作った物（分割ごとに、外した化合物の組は使わない）。
#   分割 k ごとに: 教師 2 本（分割 k の化合物を学習に使わない）→ 組を予測 → 小さいモデル（蒸留つき）。止める時期は選ばず、最後の周の重み。
#   GROUP=free DATA=gpu_data RUNS=runs bash run_gpu_cv.sh
# できる物: $RUNS/gSD_<GROUP>_f<k>_s0/pred.npy（その分割の外の予測）。手元で train_uncertainty_gnn.py gSD <GROUP> <DATA> <RUNS> <重みのフォルダ>。
set -e
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
G="${GROUP:-free}"; D="${DATA:-gpu_data}"; O="${RUNS:-runs}"
C="--data $D --group $G --model gnn --graphs graphs_js.npz --select last"
for k in ${FOLDS:-0 1 2 3}; do
  for s in 0 1; do
    [ -f "$O/teach_${G}_f${k}_s$s/result.json" ] || $PY train_teacher.py $C --fold $k --seed $s --epochs 40 --out "$O/teach_${G}_f${k}_s$s"
  done
  mkdir -p "$O/teach_${G}_f$k"
  $PY train_teacher.py $C --fold $k --aug "$D/aug.npz" --predict "$O/teach_${G}_f${k}_s0,$O/teach_${G}_f${k}_s1" --out "$O/teach_${G}_f$k"
  [ -f "$O/gSD_${G}_f${k}_s0/result.json" ] || $PY train_teacher.py $C --fold $k --seed 0 --hid 128 --layers 4 --head 256 --epochs 30 --aug "$D/aug.npz" \
    --aug_y "$O/teach_${G}_f$k/aug_pred.npy" --aug_n 250000 --aug_w 0.5 --out "$O/gSD_${G}_f${k}_s0"
  date > "$O/DONE_f$k"
done
