# SPDX-License-Identifier: MIT
"""同梱の基礎データ（base_data/）を読む。

  methods.csv       1 行 = 1 つの分析法（カラムの性質・溶媒・勾配か等組成の割合・死時間・pH・温度）
  rows.csv.gz       1 行 = 1 つの測定（method_id, ck, rt_min または logk_iso）。並びは学習の時の並びのまま（並びが変わると木がわずかに変わる）
  compounds.csv.gz  1 行 = 1 つの化合物（ck = InChIKey の前半 14 文字）。js/features.js で計算した特徴量（記述子 41・官能基 56・イオン化の規則 30・指紋）
小数は「round_trip」で読む（書いた値と 1 ビットも違わない値に戻すため）。
"""
import os, re, json
import numpy as np, pandas as pd
from . import common as C


def feature_names():
    """js/features.js の DESC・GROUPS・ION の名前（並びは特徴量の並び）"""
    js = open(os.path.join(C.JS, "features.js"), encoding="utf-8").read()
    seg = lambda key: js[js.index(key):js.index("];", js.index(key))]
    desc = re.findall(r"'([A-Za-z0-9]+)'", seg("var DESC = ["))
    grp = re.findall(r"\['([a-z0-9_]+)',\s*'", seg("var GROUPS = ["))
    ion = re.findall(r"\['([a-z0-9_+]+)',\s*.*?,\s*'(?:acid|base|mod|perm\+)'", seg("var ION = ["))
    return desc, grp, ion


def read_csv(path, **kw):
    return pd.read_csv(path, float_precision="round_trip", encoding="utf-8", **kw)


def load_methods(path):
    md = read_csv(path, dtype={"method_id": str, "source": str, "column_id": str, "column_name": str, "gradient": str, "note": str})
    md["column_id"] = md["column_id"].fillna("")
    return md


def compounds_from_frame(T, desc, grp, ion):
    """compounds.csv の形の表 → {ck: dict(smiles, d, g, ion, ion_ok, fp)}"""
    D = T[["d_" + x for x in desc]].to_numpy(np.float64)
    G = T[["g_" + x for x in grp]].to_numpy(np.float64)
    I = T[["ion_" + x for x in ion]].to_numpy(np.float64)
    ok = T["ion_ok"].to_numpy(int)
    comp = {}
    for i, (ck, smi, fp) in enumerate(zip(T["ck"], T["smiles"], T["fp"])):
        comp[ck] = dict(smiles=smi, d=D[i], g=G[i], ion=I[i], ion_ok=bool(ok[i]), fp=fp)
    return comp


def load_base(base=None):
    """→ (methods, rows, comp)"""
    base = base or C.BASE
    desc, grp, ion = feature_names()
    md = load_methods(os.path.join(base, "methods.csv"))
    rows = read_csv(os.path.join(base, "rows.csv.gz"), dtype={"method_id": str, "ck": str})
    T = read_csv(os.path.join(base, "compounds.csv.gz"), dtype={"ck": str, "smiles": str, "fp": str}, keep_default_na=False)
    want = ["ck", "smiles"] + ["d_" + x for x in desc] + ["g_" + x for x in grp] + ["ion_" + x for x in ion] + ["ion_ok", "fp"]
    if list(T.columns) != want:
        raise SystemExit("base_data/compounds.csv.gz の列が js/features.js の並びと合わない（特徴量の並びを変えたら、build_from_raw.py で作り直す）")
    return md, rows, compounds_from_frame(T, desc, grp, ion)


def columns_table():
    return json.load(open(os.path.join(C.WEIGHTS, "columns.json"), encoding="utf-8"))


def load_par():
    from . import chromlib2 as c2
    return c2.load_par(os.path.join(C.BASE, "retention_form.json"))
