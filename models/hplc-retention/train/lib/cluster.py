# SPDX-License-Identifier: MIT
"""指紋（features.js の 512 ビットの Morgan 指紋）の類似度で、化合物をかたまりに分ける。予測の幅の 4 分割に使う。"""
import hashlib
import numpy as np


def leader_cluster(bits, thr=0.5, seed=0):
    """bits: (n, 64) uint8。順に見て、既存の代表との Tanimoto が thr 以上ならその組、なければ新しい代表。"""
    pc = np.array([bin(i).count("1") for i in range(256)], np.uint8)
    n = len(bits)
    order = np.random.default_rng(seed).permutation(n)
    cnt = pc[bits].sum(axis=1).astype(np.int32)
    lead_bits = np.zeros((0, bits.shape[1]), np.uint8); lead_cnt = np.zeros(0, np.int32); n_lead = 0
    lab = np.full(n, -1, np.int64)
    for i in order:
        if n_lead:
            inter = pc[lead_bits & bits[i]].sum(axis=1).astype(np.int32)
            tani = inter / np.maximum(lead_cnt + cnt[i] - inter, 1)
            j = int(tani.argmax())
            if tani[j] >= thr:
                lab[i] = j
                continue
        lab[i] = n_lead; n_lead += 1
        lead_bits = np.vstack([lead_bits, bits[i][None, :]]); lead_cnt = np.append(lead_cnt, cnt[i])
    return lab


def folds(cks, comp, in_set, src_of_comp_is_smrt_only, nfold=4, n_smrt=15000, thr=0.5):
    """化合物ごとの分割の番号。in_set の化合物だけでかたまりを作る（その版の学習に入っていない化合物は使わない）。
    SMRT だけにある化合物は n_smrt 個までをかたまりに入れ、残りはハッシュで分ける。"""
    cks = list(cks)
    rng = np.random.default_rng(0)
    cand = np.where(in_set)[0]
    sm = cand[src_of_comp_is_smrt_only[cand]]
    drop = set(rng.choice(sm, max(len(sm) - n_smrt, 0), replace=False).tolist())
    used = np.array([i for i in cand if i not in drop])
    fp = np.zeros((len(used), 64), np.uint8)
    for k, i in enumerate(used):
        fp[k] = np.frombuffer(bytes.fromhex(comp[cks[i]]["fp"]), np.uint8)
    lab = leader_cluster(fp, thr)
    fold_of_cluster = np.array([int(hashlib.md5(str(c).encode()).hexdigest()[:6], 16) % nfold for c in range(lab.max() + 1)])
    out = np.array([int(hashlib.md5(ck.encode()).hexdigest()[8:14], 16) % nfold for ck in cks])
    out[used] = fold_of_cluster[lab]
    return out, int(lab.max() + 1), int(len(used))
