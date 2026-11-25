# SPDX-License-Identifier: MIT
"""2 つの基礎データのフォルダを比べる（build_from_raw.py で作り直した物と、同梱の base_data/）。
  python train/scripts/compare_base_data.py <作り直したフォルダ> [<比べる先（既定 train/base_data）>]
行の並び・化合物・特徴量・分析法の条件が同じかを数える。同じなら、学習の結果も同じになる。
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from lib import common as C
from lib.table import read_csv

A = os.path.abspath(sys.argv[1]); B = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else C.BASE
ok = True


def num_diff(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    both_nan = np.isnan(x) & np.isnan(y)
    if ((np.isnan(x) ^ np.isnan(y))).any():
        return np.inf
    d = np.abs(np.where(both_nan, 0, x - y))
    return float(d.max()) if len(d) else 0.0


ma, mb = read_csv(os.path.join(A, "methods.csv"), dtype=str, keep_default_na=False), read_csv(os.path.join(B, "methods.csv"), dtype=str, keep_default_na=False)
print("分析法: %d / %d、並びが同じ: %s" % (len(ma), len(mb), list(ma["method_id"]) == list(mb["method_id"])))
if list(ma["method_id"]) != list(mb["method_id"]):
    ok = False
else:
    for c in mb.columns:
        if c in ("method_id",):
            continue
        try:
            d = num_diff(pd.to_numeric(ma[c].replace("", np.nan)), pd.to_numeric(mb[c].replace("", np.nan)))
            if d > 0:
                print("  数の列 %-14s 差の最大 %.3g" % (c, d))
                ok &= d < 1e-9
        except (ValueError, TypeError):
            n = int((ma[c] != mb[c]).sum())
            if n:
                print("  文字の列 %-12s 違う行 %d（例: %r / %r）" % (c, n, ma[c][ma[c] != mb[c]].iloc[0][:80], mb[c][ma[c] != mb[c]].iloc[0][:80]))
                ok &= c == "note"
ra, rb = read_csv(os.path.join(A, "rows.csv.gz"), dtype={"method_id": str, "ck": str}), read_csv(os.path.join(B, "rows.csv.gz"), dtype={"method_id": str, "ck": str})
same_rows = len(ra) == len(rb) and (ra["method_id"].values == rb["method_id"].values).all() and (ra["ck"].values == rb["ck"].values).all()
print("測定: %d / %d、並び・分析法・化合物が同じ: %s" % (len(ra), len(rb), same_rows))
if same_rows:
    print("  rt_min の差の最大 %.3g、logk_iso の差の最大 %.3g" % (num_diff(ra["rt_min"], rb["rt_min"]), num_diff(ra["logk_iso"], rb["logk_iso"])))
    ok &= num_diff(ra["rt_min"], rb["rt_min"]) < 1e-9 and num_diff(ra["logk_iso"], rb["logk_iso"]) < 1e-9
else:
    ok = False
    print("  分析法ごとの行の数の違い:", (ra.groupby("method_id").size() - rb.groupby("method_id").size()).replace(0, np.nan).dropna().head(10).to_dict())
ca = read_csv(os.path.join(A, "compounds.csv.gz"), dtype={"ck": str, "smiles": str, "fp": str}, keep_default_na=False).set_index("ck")
cb = read_csv(os.path.join(B, "compounds.csv.gz"), dtype={"ck": str, "smiles": str, "fp": str}, keep_default_na=False).set_index("ck")
common = ca.index.intersection(cb.index)
print("化合物: %d / %d、共通 %d" % (len(ca), len(cb), len(common)))
ok &= len(common) == len(ca) == len(cb)
ca, cb = ca.loc[common], cb.loc[common]
print("  標準の SMILES が違う: %d、指紋が違う: %d" % (int((ca["smiles"] != cb["smiles"]).sum()), int((ca["fp"] != cb["fp"]).sum())))
num = [c for c in cb.columns if c.startswith(("d_", "g_", "ion_"))]
dd = np.abs(ca[num].to_numpy(float) - cb[num].to_numpy(float))
bad = (dd > 0).any(axis=1)
print("  特徴量が 1 つでも違う化合物: %d（差の最大 %.3g）" % (int(bad.sum()), float(dd.max()) if dd.size else 0.0))
if bad.any():
    j = np.where(dd.max(axis=0) > 0)[0]
    print("  違いのある列:", ", ".join("%s(%d)" % (num[i], int((dd[:, i] > 0).sum())) for i in j[:12]))
ok &= (ca["fp"] == cb["fp"]).all() and float(dd.max() if dd.size else 0) <= 1e-4 + 1e-12
for f in ["stock.csv", "names.csv.gz"]:
    a, b = read_csv(os.path.join(A, f), dtype=str, keep_default_na=False), read_csv(os.path.join(B, f), dtype=str, keep_default_na=False)
    same = a.shape == b.shape and (a.values == b.values).all()
    print("%s: %s" % (f, "同じ" if same else "違う（%s / %s）" % (a.shape, b.shape)))
    ok &= same or f == "stock.csv"
print("判定:", "同じ（学習に使う値は同じ）" if ok else "違いがある（上を見る）")
