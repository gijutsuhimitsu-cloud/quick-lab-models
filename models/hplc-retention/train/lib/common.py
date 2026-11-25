# SPDX-License-Identifier: MIT
"""学習の道具で共通に使う場所と設定。

場所はこのファイルからの相対で決める（個人の環境のパスを書かない）。
Node の場所は環境変数 NODE（既定 `node`）。Python は、この道具を動かしている Python をそのまま使う
（`run.sh` は環境変数 PYTHON を見る）。
"""
import os, sys, json, time, subprocess

TRAIN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # models/hplc-retention/train
MODEL = os.path.dirname(TRAIN)                                          # models/hplc-retention
REPO = os.path.dirname(os.path.dirname(MODEL))
JS = os.path.join(MODEL, "js")
WEIGHTS = os.path.join(MODEL, "weights")
TESTS = os.path.join(MODEL, "tests")
VENDOR_RDKIT = os.path.join(REPO, "vendor", "rdkit")
BASE = os.path.join(TRAIN, "base_data")
SCRIPTS = os.path.join(TRAIN, "scripts")
WORK = os.path.join(TRAIN, "work")         # 途中の物（.gitignore で除外）
OUT = os.path.join(TRAIN, "out")           # できた重み（.gitignore で除外）

NODE = os.environ.get("NODE", "node")

# 版ごとの学習の出所。公開の重み（free）は、この 6 つだけで学習する（継承つき〈CC BY-SA〉のデータは入れない）
FREE_SRC = ["smrt_figshare", "hsm2_mendeley", "mcart_msd_zenodo", "wiczling_osf", "zhang2024_scidb", "weber_isocratic_mendeley"]
OWN_SRC = "own"

_t0 = time.time()


def log(*a):
    print("[%6.0fs]" % (time.time() - _t0), *a, flush=True)


def node(script, *args, cwd=None):
    """Node のスクリプトを実行する（失敗したら止める）。"""
    cmd = [NODE, os.path.join(SCRIPTS, script) if not os.path.isabs(script) else script] + [str(a) for a in args]
    r = subprocess.run(cmd, cwd=cwd)
    if r.returncode != 0:
        raise SystemExit("Node のスクリプトが失敗した: %s（NODE=%s）" % (" ".join(cmd[1:]), NODE))


def clean(o):
    """JSON に書けない値（nan・inf）を null にする。タプルはリストに。"""
    import math
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


def dump_json(obj, path, compact=True, ascii=False):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        if compact:
            json.dump(clean(obj), f, ensure_ascii=ascii, separators=(",", ":"), allow_nan=False)
        else:
            json.dump(clean(obj), f, ensure_ascii=ascii, indent=1, allow_nan=False)


def env_int(name, default):
    v = os.environ.get(name, "")
    return int(v) if v.strip() else int(default)
