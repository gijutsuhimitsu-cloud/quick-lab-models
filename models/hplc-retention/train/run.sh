#!/bin/bash
# SPDX-License-Identifier: MIT
# 学習し直す（run.py の入口。引数はそのまま run.py に渡す）。
#   bash train/run.sh                         公開の重みを基礎データから作り直す
#   bash train/run.sh --own train/own_data    自社のデータを足して学習する
# Python と Node の場所は環境変数 PYTHON・NODE（既定 python・node）。
set -e
export PYTHONIOENCODING=utf-8
exec "${PYTHON:-python}" "$(dirname "$0")/run.py" "$@"
