# SPDX-License-Identifier: MIT
"""一度 run.py で学習した設定（train/out/<名前>/run_config.json）から、同じ学習の表を作り直す（評価・較正・GPU 用の書き出しで使う）。"""
import os, json
import pandas as pd
from . import common as C, table, own as OWN
from .data import Data


def data_for_run(name):
    """→ (Data, 設定, 分析法ごとの補正 {method_id: (死時間の倍率, 遅れ)})"""
    O, W = os.path.join(C.OUT, name), os.path.join(C.WORK, name)
    p = os.path.join(O, "run_config.json")
    if not os.path.exists(p):
        raise SystemExit("%s が無い。先に python train/run.py --name %s を実行する" % (p, name))
    cfg = json.load(open(p, encoding="utf-8"))["args"]
    table.load_par()
    md, rows, comp = table.load_base()
    sources = list(C.FREE_SRC)
    if cfg.get("own"):
        omd, orows, onew, _, _ = OWN.load(cfg["own"], W, comp)
        comp.update(onew)
        if cfg.get("own_only"):
            used = set(orows["ck"])
            md, rows, sources, comp = omd, orows, [C.OWN_SRC], {k: v for k, v in comp.items() if k in used}
        else:
            md, rows = pd.concat([md, omd], ignore_index=True), pd.concat([rows, orows], ignore_index=True)
            sources.append(C.OWN_SRC)
    D = Data(md, rows, comp, sources, holdout=cfg["holdout"], own_weight=cfg["own_weight"])
    cal = {k: tuple(v) for k, v in json.load(open(os.path.join(W, "calibration.json"), encoding="utf-8")).items()}
    return D, cfg, cal
