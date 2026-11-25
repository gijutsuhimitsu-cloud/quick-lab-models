"""融点の公開データ（xlsx）を CSV にする。openpyxl が要る。

元のデータ: Bradley, J.-C.; Williams, A.; Lang, A. (2014) Jean-Claude Bradley Open Melting Point Dataset. figshare.
            doi:10.6084/m9.figshare.1031637（CC0 1.0）。ファイル BradleyMeltingPointDataset.xlsx
使い方（models/solubility/ で）:
  python train/rebuild/bradley_to_csv.py [入力の .xlsx] [出力の .csv]
  既定の入力: train/raw/bradley_mp/BradleyMeltingPointDataset.xlsx
  既定の出力: train/work/bradley_mp_raw.csv
変換した表から、train/rebuild/bradley_tm.cjs が融点の照合の表（weights/tm_lookup.json）と学習用の表（train/base_data/tm_train.csv.gz）を作る。
元の表（行ごとの名前・出所の列つき）はリポジトリに入れていない（同梱の表は構造と融点だけ）。
"""
import csv, io, os, sys

import openpyxl

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
TRAIN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("SOL_RAW", os.path.join(TRAIN, "raw"))
WORK = os.environ.get("SOL_WORK", os.path.join(TRAIN, "work"))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RAW, "bradley_mp", "BradleyMeltingPointDataset.xlsx")
DST = sys.argv[2] if len(sys.argv) > 2 else os.path.join(WORK, "bradley_mp_raw.csv")

wb = openpyxl.load_workbook(SRC, read_only=True, data_only=True)
ws = wb.worksheets[0]
rows = list(ws.iter_rows(values_only=True))
os.makedirs(os.path.dirname(os.path.abspath(DST)), exist_ok=True)
with open(DST, "w", newline="", encoding="utf-8") as f:
    csv.writer(f).writerows([["" if c is None else c for c in r] for r in rows])
print("シート", ws.title, "行", len(rows) - 1, "→", DST)
