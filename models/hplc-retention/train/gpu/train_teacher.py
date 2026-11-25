# SPDX-License-Identifier: MIT
"""グラフのモデル（教師と、蒸留した小さいモデル）の学習（GPU。PyTorch だけ。RDKit は要らない＝グラフは手元の Node で graphs_js.cjs が作る）。
  python train_teacher.py --data gpu_data --graphs graphs_js.npz --model gnn --group free --seed 0 --epochs 40 --out runs/teach_s0
  手順は README.md と run_gpu.sh。この写しは、元の作業の物から説明と既定の値だけを直した物（この写しでは未実行。PyTorch の無い PC で構文だけ確かめた）
出口は今のモデルと同じ: 機械学習が出すのは lkw 1 つ。保持の式と勾配の積分は、あとで手元の同じ道具で計算する。
  lkw = f0(z, c, u) + Σ_j s_j · softplus(g_j(z, c)) · q_j      z: 分子の表現、c: カラムと条件、q_j: 向きの決まったイオン化の量（s_j = ±1）、u: 向きの決まっていない量
  → 酸・塩基の pH の向きは、木のモデルの向きの制約と同じく、必ず守られる。
目的変数が上限（9）・下限（−1）に当たっている行は打ち切りとして、片側だけの損失にする。
分け方: test_c（10% の化合物）は学習に使わない。学習の化合物の 5% を、止める時期を選ぶための確認用に取る。
"""
import os, sys, json, time, argparse, hashlib
import numpy as np
import torch, torch.nn as nn, torch.nn.functional as F

ap = argparse.ArgumentParser()
ap.add_argument("--data", default="gpu_data"); ap.add_argument("--model", default="gnn"); ap.add_argument("--group", default="free")
ap.add_argument("--seed", type=int, default=0); ap.add_argument("--epochs", type=int, default=40); ap.add_argument("--bs", type=int, default=1024)
ap.add_argument("--lr", type=float, default=1e-3); ap.add_argument("--hid", type=int, default=256); ap.add_argument("--layers", type=int, default=5)
ap.add_argument("--drop", type=float, default=0.1); ap.add_argument("--out", default="runs/tmp"); ap.add_argument("--max_rows", type=int, default=0)
ap.add_argument("--save_z", type=int, default=0)
ap.add_argument("--aug", default="")          # 蒸留用の「分子 × 条件」の組（gen_aug.py の aug.npz）
ap.add_argument("--aug_y", default="")        # その組への教師の予測（.npy）。これがあれば、学習に足す（生徒の学習）
ap.add_argument("--aug_w", type=float, default=0.25)
ap.add_argument("--fold", type=int, default=-1)       # 交差検証: この番号の分割の化合物を学習から外す（rows.npz の fold_<group>）
ap.add_argument("--select", default="val")           # val: 確認用の誤差がいちばん小さい周 / last: 最後の周
ap.add_argument("--graphs", default="graphs.npz")     # グラフのファイル（graphs.npz＝Python の RDKit、graphs_js.npz＝ブラウザと同じコードで作った物）
ap.add_argument("--head", type=int, default=512)      # 出口の全結合の幅
ap.add_argument("--all", type=int, default=0)         # 1: 評価用の化合物を残さず、その版の全部の行で学習（本番用。精度は交差検証で見る）
ap.add_argument("--aug_n", type=int, default=0)       # 学習に足す組の数の上限（0 は全部）
ap.add_argument("--rt_w", type=float, default=0.0)    # > 0: 勾配の行は、実測の保持時間との差（分析時間に対する割合）を主な損失にする（rt_table_<group>.npz が要る）
ap.add_argument("--lkw_w", type=float, default=0.2)   # --rt_w の時、逆算した lkw との差に残す重み（勾配の行）
ap.add_argument("--export", type=int, default=0)      # 1: ブラウザ用に重みを書き出す（weights.bin・manifest.json・check.json）
ap.add_argument("--predict", default="")      # 学習済みのフォルダ（カンマ区切り）。--aug の組を予測して --out/aug_pred.npy に平均を書く（教師の予測）
A = ap.parse_args()
torch.manual_seed(A.seed); np.random.seed(A.seed)
dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
os.makedirs(A.out, exist_ok=True)
R = np.load(os.path.join(A.data, "rows.npz"))
G = np.load(os.path.join(A.data, A.graphs))
meta = json.load(open(os.path.join(A.data, "meta.json"), encoding="utf-8"))
cks = [l.split("\t")[0] for l in open(os.path.join(A.data, "compounds.tsv"), encoding="utf-8")]
smis = [l.rstrip("\n").split("\t")[1] for l in open(os.path.join(A.data, "compounds.tsv"), encoding="utf-8")]
YF, YC = meta["y_floor"], meta["y_ceil"]
signs = np.array(meta["ion_signs"])
jc, ju = np.where(signs != 0)[0], np.where(signs == 0)[0]

ci, mi = R["ci"].astype(np.int64), R["mi"].astype(np.int64)
y_raw = R["y_raw_" + A.group].astype(np.float32)
in_tr = R["in_" + A.group] if ("in_" + A.group) in R.files else np.ones(len(ci), bool)
test_c = R["test_c"]
val_c = np.array([int(hashlib.md5(("val" + c).encode()).hexdigest()[:8], 16) % 1000 < 50 for c in cks])[ci]
ok = np.isfinite(y_raw)
tr = in_tr & ~test_c & ~val_c & ok
va = in_tr & ~test_c & val_c & ok
if A.fold >= 0:
    # 交差検証: 外した分割の化合物は学習に使わない。経過を見るための行は、外した分割から 8,000 行（止める時期の選択には --select last を使う）
    fold = R["fold_" + A.group]
    tr = in_tr & ok & (fold != A.fold)
    held = np.where(in_tr & ok & (fold == A.fold))[0]
    va = np.zeros(len(ci), bool); va[np.random.default_rng(1).choice(held, min(8000, len(held)), replace=False)] = True
    test_c = fold == A.fold
if A.all:
    tr = in_tr & ok
    va = np.zeros(len(ci), bool); va[np.random.default_rng(1).choice(np.where(tr)[0], 8000, replace=False)] = True     # 経過を見るだけ（学習に入っている行）
    test_c = np.zeros(len(ci), bool)
if A.max_rows:
    ix = np.where(tr)[0]; keep = np.random.default_rng(0).choice(ix, min(A.max_rows, len(ix)), replace=False); tr = np.zeros(len(ci), bool); tr[keep] = True
cens = np.where(y_raw >= YC - 1e-6, 1, np.where(y_raw <= YF + 1e-6, -1, 0)).astype(np.int8)
y = np.clip(np.nan_to_num(y_raw, nan=0.0), YF, YC)


def robust(x, m):
    """裾の重い量は符号つきの log1p にしてから、学習の行の平均と幅でそろえる。"""
    x = np.nan_to_num(x.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    x = np.sign(x) * np.log1p(np.abs(x))
    mu, sd = x[m].mean(0), x[m].std(0) + 1e-6
    return ((x - mu) / sd).astype(np.float32), mu, sd


trc = np.zeros(len(cks), bool); trc[np.unique(ci[tr])] = True
Xd, d_mu, d_sd = robust(np.hstack([R["XC"], R["IC"]]), trc)
trm = np.zeros(R["M"].shape[0], bool); trm[np.unique(mi[tr])] = True
Mm = np.nan_to_num(R["M"].astype(np.float32), nan=0.0)
m_mu, m_sd = Mm[trm].mean(0), Mm[trm].std(0) + 1e-6
Mm = (Mm - m_mu) / m_sd
IFr = np.nan_to_num(R["IF"].astype(np.float32), nan=0.0)
n0 = len(ci)
w_all = R["w"].astype(np.float32)
if A.aug:
    # 平均と幅は実測の学習の行から決めたまま、組の行を後ろにつなぐ
    AU = np.load(A.aug)
    Mm = (np.nan_to_num(AU["M_ext"].astype(np.float32), nan=0.0) - m_mu) / m_sd
    IFr = np.vstack([IFr, np.nan_to_num(AU["IF"].astype(np.float32), nan=0.0)])
    ci = np.concatenate([ci, AU["ci"].astype(np.int64)]); mi = np.concatenate([mi, AU["mi"].astype(np.int64)])
    na = len(AU["ci"])
    ya = np.load(A.aug_y).astype(np.float32) if A.aug_y else np.zeros(na, np.float32)
    y = np.concatenate([y, ya]); cens = np.concatenate([cens, np.zeros(na, np.int8)])
    w_all = np.concatenate([w_all, np.full(na, A.aug_w, np.float32)])
    aug_ok = np.ones(na, bool)
    if A.fold >= 0:
        aug_ok = R["foldc_" + A.group][AU["ci"]] != A.fold          # 外した分割の化合物の組は、生徒の学習に使わない
    if A.aug_n and aug_ok.sum() > A.aug_n:
        keep_a = np.random.default_rng(2).choice(np.where(aug_ok)[0], A.aug_n, replace=False); aug_ok = np.zeros(na, bool); aug_ok[keep_a] = True
else:
    na = 0
tr0 = np.where(tr)[0]
q_sd = IFr[tr0][:, jc].std(0) + 1e-6                     # 向きの決まった量は、原点を動かさず幅だけそろえる
Q = IFr[:, jc] / q_sd
u_mu, u_sd = IFr[tr0][:, ju].mean(0), IFr[tr0][:, ju].std(0) + 1e-6
Uu = (IFr[:, ju] - u_mu) / u_sd
T = lambda a, dt=torch.float32: torch.as_tensor(a, dtype=dt, device=dev)
Xd_t, Mm_t, Q_t, U_t = T(Xd), T(Mm), T(Q), T(Uu)
y_t, w_t, cens_t = T(y), T(w_all), T(cens, torch.int8)
ci_t, mi_t = T(ci, torch.long), T(mi, torch.long)
sg_t = T(signs[jc])
atom_t, bond_t = T(G["atom"], torch.long), T(G["bond"], torch.long)
es_t, ed_t = T(G["esrc"], torch.long), T(G["edst"], torch.long)
aptr_t, bptr_t = T(G["aptr"], torch.long), T(G["bptr"], torch.long)
fp_t = torch.as_tensor(G["fp"], device=dev) if A.model == "fpmlp" else None
A_DIMS, B_DIMS = ([int(x) for x in G["adims"]], [int(x) for x in G["bdims"]]) if "adims" in G.files else ([17, 7, 5, 5, 6, 2, 2, 3], [5, 2, 2, 3])


def ranges(ptr, u):
    st = ptr[u]; cnt = ptr[u + 1] - st
    seg = torch.repeat_interleave(torch.arange(len(u), device=dev), cnt)
    new0 = torch.cumsum(cnt, 0) - cnt
    idx = torch.repeat_interleave(st, cnt) + torch.arange(int(cnt.sum()), device=dev) - torch.repeat_interleave(new0, cnt)
    return idx, seg, new0, cnt


class GNN(nn.Module):
    def __init__(s, d, L, drop):
        super().__init__()
        s.ae = nn.ModuleList([nn.Embedding(n, d) for n in A_DIMS])
        s.be = nn.ModuleList([nn.ModuleList([nn.Embedding(n, d) for n in B_DIMS]) for _ in range(L)])
        s.mlp = nn.ModuleList([nn.Sequential(nn.Linear(d, 2 * d), nn.SiLU(), nn.Linear(2 * d, d)) for _ in range(L)])
        s.nrm = nn.ModuleList([nn.LayerNorm(d) for _ in range(L)])
        s.eps = nn.Parameter(torch.zeros(L)); s.drop = drop
        s.out = nn.Sequential(nn.Linear(3 * d, d), nn.SiLU())

    def forward(s, u):
        ai, aseg, a0, acnt = ranges(aptr_t, u)
        bi, bseg, _, _ = ranges(bptr_t, u)
        af, bf = atom_t[ai], bond_t[bi]
        src, dst = es_t[bi] + a0[bseg], ed_t[bi] + a0[bseg]
        h = sum(e(af[:, k]) for k, e in enumerate(s.ae))
        for l in range(len(s.mlp)):
            e = sum(emb(bf[:, k]) for k, emb in enumerate(s.be[l]))
            m = torch.zeros_like(h).index_add_(0, dst, F.silu(h[src] + e))
            h = h + F.dropout(s.mlp[l](s.nrm[l]((1 + s.eps[l]) * h + m)), s.drop, s.training)
        n = len(u)
        sm = torch.zeros(n, h.shape[1], device=dev).index_add_(0, aseg, h)
        mean = sm / acnt.clamp(min=1).unsqueeze(1)
        mx = torch.full((n, h.shape[1]), -1e9, device=dev).scatter_reduce(0, aseg.unsqueeze(1).expand_as(h), h, "amax", include_self=True)
        mx = torch.where(mx < -1e8, torch.zeros_like(mx), mx)
        return s.out(torch.cat([sm / 8.0, mean, mx], 1))


class Net(nn.Module):
    def __init__(s, kind, d, L, drop):
        super().__init__()
        s.kind = kind
        nd = Xd.shape[1]
        s.desc = nn.Sequential(nn.Linear(nd, 256), nn.SiLU(), nn.Dropout(drop), nn.Linear(256, 128), nn.SiLU())
        if kind == "gnn":
            s.enc = GNN(d, L, drop); zd = d + 128
        elif kind == "fpmlp":
            s.enc = nn.Sequential(nn.Linear(2048, 1024), nn.SiLU(), nn.Dropout(drop), nn.Linear(1024, d), nn.SiLU()); zd = d + 128
        else:
            s.enc = None; zd = 128
        s.zd = zd
        nc, nu, nq = Mm.shape[1], Uu.shape[1], Q.shape[1]
        H = A.head
        s.f0 = nn.Sequential(nn.Linear(zd + nc + nu, H), nn.SiLU(), nn.Dropout(drop), nn.Linear(H, H // 2), nn.SiLU(), nn.Linear(H // 2, 1))
        s.g = nn.Sequential(nn.Linear(zd + nc, H // 2), nn.SiLU(), nn.Linear(H // 2, nq))

    def encode(s, u):
        z = s.desc(Xd_t[u])
        if s.kind == "gnn":
            z = torch.cat([s.enc(u), z], 1)
        elif s.kind == "fpmlp":
            z = torch.cat([s.enc(torch.log1p(fp_t[u].float())), z], 1)
        return z

    def head(s, z, c, uu, q):
        base = s.f0(torch.cat([z, c, uu], 1)).squeeze(1)
        coef = F.softplus(s.g(torch.cat([z, c], 1)) - 2.0)
        return base + (coef * sg_t * q).sum(1)

    def forward(s, rows):
        u, inv = torch.unique(ci_t[rows], return_inverse=True)
        z = s.encode(u)[inv]
        return s.head(z, Mm_t[mi_t[rows]], U_t[rows], Q_t[rows])


if A.rt_w > 0:
    # 実測の保持時間を目的にする: 分析法ごとの表（lkw の格子 → 保持時間）を補間して、予測の保持時間を出す
    RT = np.load(os.path.join(A.data, "rt_table_%s.npz" % A.group))
    g0, dg, ng = float(RT["grid"][0]), float(RT["grid"][1] - RT["grid"][0]), len(RT["grid"])
    Ft = T(np.nan_to_num(RT["F"], nan=0.0))
    run_m = RT["run"]
    has_rt = np.zeros(len(ci), bool)
    has_rt[:n0] = np.isfinite(RT["rt"]) & np.isfinite(run_m[mi[:n0]]) & (run_m[mi[:n0]] > 0)
    rt_t = T(np.concatenate([np.nan_to_num(RT["rt"], nan=0.0), np.zeros(len(ci) - n0, np.float32)]))
    runrow_t = T(np.concatenate([np.nan_to_num(run_m[mi[:n0]], nan=1.0), np.ones(len(ci) - n0, np.float32)]))
    hasrt_t = T(has_rt, torch.bool)
    RT_SCALE = 0.03          # 分析時間の 3% の外れを、損失の 1 の大きさにする


def loss_fn(p, rows):
    yy, cc = y_t[rows], cens_t[rows]
    d = p - yy
    d = torch.where(cc == 1, torch.clamp(d, max=0.0), torch.where(cc == -1, torch.clamp(d, min=0.0), d))
    l = F.huber_loss(d, torch.zeros_like(d), delta=1.0, reduction="none")
    if A.rt_w > 0:
        x = ((p - g0) / dg).clamp(0.0, ng - 1.001)
        i0 = x.floor().long(); fr = x - i0
        mm = mi_t[rows].clamp(max=Ft.shape[0] - 1)          # 合成の条件（蒸留用の組）は表が無い。has_rt が偽なので使われない
        tr_p = Ft[mm, i0] * (1 - fr) + Ft[mm, i0 + 1] * fr
        dr = (tr_p - rt_t[rows]) / runrow_t[rows] / RT_SCALE
        lr = F.huber_loss(dr, torch.zeros_like(dr), delta=1.0, reduction="none")
        h = hasrt_t[rows]
        l = torch.where(h, A.rt_w * lr + A.lkw_w * l, l)
    return (w_t[rows] * l).sum() / w_t[rows].sum()


@torch.no_grad()
def predict(net, idx):
    net.eval(); out = np.zeros(len(idx), np.float32)
    for k in range(0, len(idx), 4096):
        out[k:k + 4096] = net(T(idx[k:k + 4096], torch.long)).float().cpu().numpy()
    return out


def mae_unc(p, idx):
    m = cens[idx] == 0
    return float(np.mean(np.abs(p[m] - y[idx][m])))


net = Net(A.model, A.hid, A.layers, A.drop).to(dev)
npar = sum(p.numel() for p in net.parameters())
if A.predict:
    # 教師の予測: 学習済みのモデルごとに、組の行を予測して平均とばらつきを書く
    ps = []
    for d in A.predict.split(","):
        net.load_state_dict(torch.load(os.path.join(d, "best.pt"), map_location=dev))
        ps.append(predict(net, np.arange(n0, n0 + na)))
        print("predicted", d, float(ps[-1].mean()), flush=True)
    ps = np.array(ps)
    np.save(os.path.join(A.out, "aug_pred.npy"), ps.mean(0)); np.save(os.path.join(A.out, "aug_pred_std.npy"), ps.std(0))
    print("aug rows", na, "mean", float(ps.mean()), "mean std between models", float(ps.std(0).mean()))
    sys.exit(0)
tri, vai = np.where(tr)[0], np.where(va)[0]
if A.aug_y:
    tri = np.concatenate([tri, n0 + np.where(aug_ok)[0]])
steps = A.epochs * (len(tri) // A.bs)
opt = torch.optim.AdamW(net.parameters(), lr=A.lr, weight_decay=1e-5)
sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=A.lr, total_steps=steps, pct_start=0.05)
print("device", dev, "| model", A.model, "params", npar, "| train rows", len(tri), "val rows", len(vai), "| steps", steps, flush=True)
log, best, t0, rng = [], (1e9, -1), time.time(), np.random.default_rng(A.seed)
for ep in range(A.epochs):
    net.train(); perm = rng.permutation(tri); tot = 0.0; nb = len(perm) // A.bs
    for b in range(nb):
        rows = T(perm[b * A.bs:(b + 1) * A.bs], torch.long)
        l = loss_fn(net(rows), rows)
        opt.zero_grad(set_to_none=True); l.backward(); nn.utils.clip_grad_norm_(net.parameters(), 5.0); opt.step(); sched.step()
        tot += float(l)
    v = mae_unc(predict(net, vai), vai)
    log.append(dict(epoch=ep, train_loss=tot / nb, val_mae=v, sec=time.time() - t0))
    if v < best[0] or A.select == "last":
        best = (v, ep); torch.save(net.state_dict(), os.path.join(A.out, "best.pt"))
    if ep % 2 == 0 or ep == A.epochs - 1 or A.aug_y:
        print("ep %d loss %.4f val_mae %.4f best %.4f@%d  %.0fs" % (ep, tot / nb, v, best[0], best[1], time.time() - t0), flush=True)
net.load_state_dict(torch.load(os.path.join(A.out, "best.pt")))
pred = predict(net, np.arange(n0))
np.save(os.path.join(A.out, "pred.npy"), pred)
te = in_tr & test_c & ok
if not te.any():
    te = va.copy()          # --all の時: 評価用の行は無い（学習に入っている行での値を、目安として書くだけ）
res = dict(args=vars(A), params=npar, n_train=int(len(tri)), n_val=int(len(vai)), best_val_mae=best[0], best_epoch=best[1], sec=time.time() - t0,
           test_mae_uncensored=mae_unc(pred[te], np.where(te)[0]), tree_test_mae_uncensored=mae_unc(R["tree_" + A.group][te], np.where(te)[0]),
           by_source={meta["sources"][s]: dict(n=int((te & (R["src"] == s)).sum()), model=mae_unc(pred[te & (R["src"] == s)], np.where(te & (R["src"] == s))[0]),
                                               tree=mae_unc(R["tree_" + A.group][te & (R["src"] == s)], np.where(te & (R["src"] == s))[0]))
                      for s in range(len(meta["sources"])) if (te & (R["src"] == s) & (cens[:n0] == 0)).sum() > 0},
           log=log, device=str(dev), gpu=torch.cuda.get_device_name(0) if torch.cuda.is_available() else None)
if A.save_z:
    net.eval(); Z = np.zeros((len(cks), net.zd), np.float32)
    with torch.no_grad():
        for k in range(0, len(cks), 2048):
            Z[k:k + 2048] = net.encode(torch.arange(k, min(k + 2048, len(cks)), device=dev)).cpu().numpy()
    np.save(os.path.join(A.out, "z.npy"), Z)
json.dump(res, open(os.path.join(A.out, "result.json"), "w"), indent=1)
print(json.dumps({k: res[k] for k in ["params", "best_val_mae", "best_epoch", "sec", "test_mae_uncensored", "tree_test_mae_uncensored", "by_source"]}, indent=1), flush=True)
if A.export:
    # ブラウザ用の書き出し: 重み（float32 を並べた 1 つのファイル）、並びと正規化の定数（manifest.json）、照合用の例（check.json）
    sd = {k: v.detach().cpu().numpy().astype(np.float32) for k, v in net.state_dict().items()}
    off, man, blobs = 0, {}, []
    for k, v in sd.items():
        man[k] = dict(shape=list(v.shape), offset=off); blobs.append(v.ravel()); off += int(v.size)
    np.concatenate(blobs).tofile(os.path.join(A.out, "weights.bin"))
    fl = lambda a: [float(x) for x in np.asarray(a).ravel()]
    manifest = dict(format=1, kind=A.model, group=A.group, hid=A.hid, layers=A.layers, head=A.head, n_params=int(off), atom_dims=A_DIMS, bond_dims=B_DIMS,
                    graph_version=int(G["version"]) if "version" in G.files else 0, feature_version=meta.get("fv"), n_desc=int(Xd.shape[1]),
                    method_features=meta["mnames"], ion_keep=meta["ion_keep"], ion_signs=meta["ion_signs"], jc=[int(x) for x in jc], ju=[int(x) for x in ju],
                    norm=dict(d_mu=fl(d_mu), d_sd=fl(d_sd), m_mu=fl(m_mu), m_sd=fl(m_sd), q_sd=fl(q_sd), u_mu=fl(u_mu), u_sd=fl(u_sd)),
                    y_floor=YF, y_ceil=YC, tensors=man, train_rows=int(len(tri)), epochs=A.epochs, seed=A.seed)
    json.dump(manifest, open(os.path.join(A.out, "manifest.json"), "w"))
    rs = np.random.default_rng(5)
    pick = np.concatenate([rs.choice(np.where(in_tr & ok & (R["src"] == k))[0], min(40, int((in_tr & ok & (R["src"] == k)).sum())), replace=False)
                           for k in range(len(meta["sources"])) if (in_tr & ok & (R["src"] == k)).any()])
    net.eval()
    with torch.no_grad():
        zz = net.encode(T(ci[pick], torch.long)).cpu().numpy()
    chk = dict(rows=[dict(smiles=smis[ci[i]], desc=fl(np.hstack([R["XC"][ci[i]], R["IC"][ci[i]]])), M=fl(np.nan_to_num(R["M"][mi[i]])), IF=fl(np.nan_to_num(R["IF"][i])),
                          lkw=float(pred[i]), z8=fl(zz[j][:8]), source=meta["sources"][int(R["src"][i])]) for j, i in enumerate(pick)])
    json.dump(chk, open(os.path.join(A.out, "check.json"), "w"))
    print("exported", off, "params ->", os.path.join(A.out, "weights.bin"), "| check rows", len(pick), flush=True)
