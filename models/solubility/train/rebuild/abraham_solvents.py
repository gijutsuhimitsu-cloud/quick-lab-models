"""Abraham の溶媒の係数の表（水 → 溶媒の分配。c・e・s・a・b・v）を JSON にする。openpyxl が要る。

元のデータ: Lang, A.; Lee, Y. (2024) Dataset: AbraLlama: Predicting Abraham Model Solute Descriptors and Modified
            Solvent Parameters Using Llama. figshare. doi:10.6084/m9.figshare.25898476（CC0 1.0）。
            ファイル「20240524 Abraham Solvent Master Dataset.xlsx」のシート All Solvent Data
使い方（models/solubility/ で）:
  python train/rebuild/abraham_solvents.py [入力の .xlsx] [出力の .json]
  既定の入力: train/raw/abraham_solvents/20240524 Abraham Solvent Master Dataset.xlsx
  既定の出力: train/work/solvents_abraham.json
変えたこと: 商品名で載っている溶媒 1 つの名前を、化学の名前に置き換える（InChIKey で見分ける）。係数はそのまま。
"""
import io, json, os, sys

import openpyxl

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
TRAIN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("SOL_RAW", os.path.join(TRAIN, "raw"))
WORK = os.environ.get("SOL_WORK", os.path.join(TRAIN, "work"))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RAW, "abraham_solvents", "20240524 Abraham Solvent Master Dataset.xlsx")
DST = sys.argv[2] if len(sys.argv) > 2 else os.path.join(WORK, "solvents_abraham.json")
RENAME = {"XXJWXESWEXIICW-UHFFFAOYSA-N": "Diethylene glycol monoethyl ether"}     # train/lib/common.cjs の SOLVENT_RENAME と同じ

wb = openpyxl.load_workbook(SRC, read_only=True, data_only=True)
ws = wb["All Solvent Data"]
rows = list(ws.iter_rows(values_only=True))
head = [str(c).strip() if c is not None else "" for c in rows[0]]
ix = {k: head.index(k) for k in ["Solvent", "C.SMILES", "c", "e", "s", "a", "b", "v", "REF", "YEAR", "INCHIKEY"]}
out, bad = [], []
for r in rows[1:]:
    if not r or r[ix["Solvent"]] is None:
        continue
    name = str(r[ix["Solvent"]]).strip()
    try:
        co = {k: float(r[ix[k]]) for k in "cesabv"}
    except (TypeError, ValueError):
        bad.append(name)
        continue
    ik = str(r[ix["INCHIKEY"]] or "").strip()
    rec = {"name": RENAME.get(ik, name), "smiles": str(r[ix["C.SMILES"]] or "").strip(), "inchikey": ik,
           "ref": str(r[ix["REF"]] or "").strip(), "year": str(r[ix["YEAR"]] or "").strip()}
    rec.update({k: round(v, 4) for k, v in co.items()})
    out.append(rec)
meta = {"source": "Abraham Solvent Master Dataset (AbraLlama), figshare, doi:10.6084/m9.figshare.25898476, CC0 1.0",
        "equation": "log10 P(water -> solvent) = c + e*E + s*S + a*A + b*B + v*V", "n": len(out)}
os.makedirs(os.path.dirname(os.path.abspath(DST)), exist_ok=True)
with open(DST, "w", encoding="utf-8") as f:
    json.dump({"meta": meta, "solvents": out}, f, ensure_ascii=False)
print("溶媒", len(out), "係数が読めない行", bad, "→", DST)
