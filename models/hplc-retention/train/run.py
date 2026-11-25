# SPDX-License-Identifier: MIT
"""HPLC の保持のモデルを学習し直す（入口。Windows でも Mac・Linux でも同じ）。

  python train/run.py                                  公開の重み（free）を基礎データから作り直す（約 10 分）
  python train/run.py --own train/own_data             基礎データ＋自社のデータで学習する
  python train/run.py --own train/own_data --own-weight 5      自社のデータの行を 5 倍の重みで
  python train/run.py --own train/own_data --own-only          自社のデータだけで学習する
  python train/run.py --trees 60 --qtrees 30           木を減らした速い試し（数分。精度は落ちる）

出力: train/out/<名前>/weights/（weights/ と同じファイル名。そのまま入れ替えられる）、REPORT.md・metrics.json など。
途中の物: train/work/<名前>/。どちらも .gitignore でリポジトリに入らない。自社のデータは、どこにも送らない。
Node の場所は環境変数 NODE（既定 node）。特徴量の計算と、JavaScript との照合に使う。
"""
import os, sys, json, time, argparse, platform
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib import common as C

ap = argparse.ArgumentParser(description="HPLC の保持のモデルを学習し直す")
ap.add_argument("--own", default="", help="自社のデータのフォルダ（methods.csv・measurements.csv。形は OWN_DATA.md）")
ap.add_argument("--own-only", action="store_true", help="自社のデータだけで学習する（基礎データは特徴量の参照にだけ使う）")
ap.add_argument("--own-weight", type=float, default=1.0, help="自社のデータの行の重みの倍率（既定 1）")
ap.add_argument("--name", default="", help="出力の名前（既定: 基礎データだけなら free、自社のデータがあれば own）。重みの meta.name になる")
ap.add_argument("--trees", type=int, default=C.env_int("CHROM_TREES", 900), help="本体の木の数（既定 900。環境変数 CHROM_TREES でも）")
ap.add_argument("--qtrees", type=int, default=C.env_int("CHROM_QTREES", 200), help="予測の幅の木の数（既定 200。環境変数 CHROM_QTREES でも）")
ap.add_argument("--holdout", type=float, default=0.1, help="評価用に学習から外す化合物の割合（既定 0.1。0 なら全部で学習し、評価の数字は出ない）")
ap.add_argument("--skip-weber-cv", action="store_true", help="等組成の交差検証（木 5 回分）を省く")
ap.add_argument("--out", default="", help="出力のフォルダ（既定 train/out/<名前>）")
ap.add_argument("--no-node", action="store_true", help="Node の段（名前の索引の確認・JavaScript との照合）を省く")
opt = ap.parse_args()
if opt.own_only and not opt.own:
    ap.error("--own-only には --own が要る")
opt.own = os.path.abspath(opt.own) if opt.own else ""
opt.name = opt.name or ("own" if opt.own else "free")
opt.outdir = os.path.abspath(opt.out) if opt.out else os.path.join(C.OUT, opt.name)
opt.weights_dir = os.path.join(opt.outdir, "weights")
opt.workdir = os.path.join(C.WORK, opt.name)
for d in (opt.outdir, opt.weights_dir, opt.workdir):
    os.makedirs(d, exist_ok=True)

import numpy as np, pandas as pd, sklearn, scipy
from lib import table, own as OWN, pipeline as PL
from lib.data import Data

C.log("出力: %s" % opt.outdir)
C.log("Python %s・numpy %s・pandas %s・scikit-learn %s・scipy %s（公開の重みは numpy 2.4.6・scikit-learn 1.9.0 で作った）" % (
    platform.python_version(), np.__version__, pd.__version__, sklearn.__version__, scipy.__version__))
table.load_par()
md, rows, comp = table.load_base()
C.log("基礎データ: 分析法 %d、測定 %d、化合物 %d" % (len(md), len(rows), len(comp)))
sources = list(C.FREE_SRC)
own_rep, own_names = None, []
if opt.own:
    omd, orows, onew, own_names, own_rep = OWN.load(opt.own, opt.workdir, comp)
    C.dump_json(own_rep, os.path.join(opt.outdir, "own_data_report.json"), compact=False)
    C.log("自社のデータ: 分析法 %d、測定 %d、化合物 %d（新しい化合物 %d）、読めない構造 %d、外した行 %d" % (
        own_rep["n_methods"], own_rep["n_rows"], own_rep["n_compounds"], own_rep["n_new_compounds"], len(own_rep["unreadable_smiles"]), len(own_rep["skipped_rows"])))
    if not len(orows):
        raise SystemExit("自社のデータに使える行が無い（own_data_report.json を見る）")
    comp.update(onew)
    if opt.own_only:
        md, rows, sources = omd, orows, [C.OWN_SRC]
        used = set(orows["ck"])
        comp = {k: v for k, v in comp.items() if k in used}
    else:
        md = pd.concat([md, omd], ignore_index=True)
        rows = pd.concat([rows, orows], ignore_index=True)
        sources.append(C.OWN_SRC)
json.dump(dict(args={k: v for k, v in vars(opt).items()}, sources=sources, versions=dict(python=platform.python_version(), numpy=np.__version__,
          pandas=pd.__version__, sklearn=sklearn.__version__, scipy=scipy.__version__), started=time.strftime("%Y-%m-%d %H:%M:%S")),
          open(os.path.join(opt.outdir, "run_config.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

D = Data(md, rows, comp, sources, holdout=opt.holdout, own_weight=opt.own_weight)
C.log("学習の表: 行 %d（学習に入る出所の行 %d）、入力 %d 個" % (len(D.rows), int(D.in_tr.sum()), D.X.shape[1]))
model, cal, log, y, yp = PL.fit_main(D, opt)
stock = PL.export(D, model, log, opt)
ex, qmodels = PL.uncertainty(D, model, cal, y, opt)
PL.seen_flags(D, model, y, ex, qmodels, stock, opt)
PL.name_index(opt, own_names)
own_res = PL.own_report(D, model, cal, y, opt)
PL.write_report(opt, log, os.path.join(opt.outdir, "uncertainty_metrics.json"), own_res, own_rep)
if not opt.no_node:
    C.log("名前の索引から、RDKit.js で読めない構造を外す")
    C.node("filter_name_index.cjs", opt.weights_dir)
    C.log("JavaScript の計算と Python の計算の照合（tests/check.cjs）")
    C.node(os.path.join(C.TESTS, "check.cjs"), opt.weights_dir)
    if not opt.own and opt.trees == 900 and opt.qtrees == 200 and opt.holdout == 0.1 and not opt.skip_weber_cv:
        C.log("公開の重み（weights/）との比べ（既定の設定で基礎データだけの時）")
        import subprocess
        r = subprocess.run([C.NODE, os.path.join(C.SCRIPTS, "compare_shipped.cjs"), opt.weights_dir])
        if r.returncode != 0:
            C.log("注意: 公開の重みと差がある（scikit-learn・numpy の版が違うと起きうる）。上の tests/check.cjs が合格なら、重みとしては使える")
    else:
        C.log("設定が既定と違うので、公開の重みとの比べは省いた（比べるのは、基礎データだけ・木 900 本・幅の木 200 本・評価用 10% の時）")
C.log("できた。結果の要約: %s" % os.path.join(opt.outdir, "REPORT.md"))
