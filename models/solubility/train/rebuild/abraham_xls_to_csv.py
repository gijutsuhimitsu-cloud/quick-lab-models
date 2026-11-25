"""Abraham 記述子の表（拡張子は .xls だが中身は HTML の表）を CSV にする。標準ライブラリだけで動く。

元のデータ: Bradley, J.-C.; Acree, W.; Lang, A. (2014) Compounds with known Abraham descriptors. figshare.
            doi:10.6084/m9.figshare.1176994（CC BY 4.0）。ファイル 20140922AbrahamDescriptorsBackUp.xls
使い方（models/solubility/ で）:
  python train/rebuild/abraham_xls_to_csv.py [入力の .xls] [出力の .csv]
  既定の入力: train/raw/abraham_descriptors/20140922AbrahamDescriptorsBackUp.xls
  既定の出力: train/work/abraham_solutes_raw.csv
"""
import csv, io, os, sys
from html.parser import HTMLParser

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
TRAIN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("SOL_RAW", os.path.join(TRAIN, "raw"))
WORK = os.environ.get("SOL_WORK", os.path.join(TRAIN, "work"))


class T(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.rows, self.row, self.cell, self.inc = [], None, [], False

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.row = []
        elif tag in ("td", "th") and self.row is not None:
            self.inc, self.cell = True, []

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self.inc:
            self.row.append(" ".join("".join(self.cell).split()))
            self.inc = False
        elif tag == "tr" and self.row is not None:
            if self.row:
                self.rows.append(self.row)
            self.row = None

    def handle_data(self, data):
        if self.inc:
            self.cell.append(data)


def main(src, dst):
    raw = open(src, "rb").read()
    for enc in ("utf-8", "cp1252", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    p = T()
    p.feed(text)
    os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
    with open(dst, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(p.rows)
    print("行", len(p.rows), "列", len(p.rows[0]) if p.rows else 0, "→", dst)


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RAW, "abraham_descriptors", "20140922AbrahamDescriptorsBackUp.xls")
    dst = sys.argv[2] if len(sys.argv) > 2 else os.path.join(WORK, "abraham_solutes_raw.csv")
    main(src, dst)
