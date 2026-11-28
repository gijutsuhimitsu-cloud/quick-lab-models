// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE（MIT）。
// 溶解度のグラフのモデル（分子グラフ → log10 S）を、ブラウザの中で直接計算する。外部の推論ライブラリは使わない。
// 学習は 学習の道具（このリポジトリには未同梱）（L4 で教師を学習し、この大きさの生徒へ蒸留）。骨組みは 逆相 HPLC のモデルの gnn_model.js と同じ。
//
//   var M = SolGNN.open(json);                       // json: weights/net_<mode>.json（重みは float16 を base64 で）
//   var z = M.encode(graph, desc);                   // graph は SolGraph.fromSmiles、desc は features.js の desc・groups に Abraham の E S A B V を足した 102 個
//   var y = M.predict(z, desc, 'etoh', TmK, 298.15); // 溶媒の id（weights/base.json）。TmK は融点（K、無ければ null）
//
// 出口: log10 S(mol/L) = A + B·t + κ·ideal(T)       t = 1000/T − 1000/298.15、ideal(T) = −2.951·max(0, Tm/T − 1)
(function (root) {
  'use strict';

  function f16(b64) {
    var bin = typeof atob === 'function' ? atob(b64) : Buffer.from(b64, 'base64').toString('binary');
    var n = bin.length >> 1, out = new Float32Array(n);
    for (var i = 0; i < n; i++) {
      var h = bin.charCodeAt(2 * i) | (bin.charCodeAt(2 * i + 1) << 8);
      var s = (h & 0x8000) ? -1 : 1, e = (h >> 10) & 0x1f, f = h & 0x3ff;
      out[i] = e === 0 ? s * f * 5.960464477539063e-8 : e === 31 ? (f ? NaN : s * Infinity) : s * Math.pow(2, e - 15) * (1 + f / 1024);
    }
    return out;
  }

  function open(J) {
    var W = J.weights_f32 ? Float32Array.from(J.weights_f32) : f16(J.weights_f16);
    if (W.length !== J.n_params) throw new Error('グラフのモデルの重みの大きさが合わない（' + W.length + ' と ' + J.n_params + '）');
    var T = J.tensors, d = J.hid, L = J.layers, N = J.norm, AD = J.atom_dims, BD = J.bond_dims, nA = AD.length, nB = BD.length;
    function t(name) {
      var e = T[name]; if (!e) throw new Error('重みが無い: ' + name);
      var n = 1; for (var i = 0; i < e.shape.length; i++) n *= e.shape[i];
      return W.subarray(e.offset, e.offset + n);
    }
    function silu(x) { return x / (1 + Math.exp(-x)); }
    function softplus(x) { return x > 30 ? x : Math.log(1 + Math.exp(x)); }
    function linear(w, b, x, nin, nout, out, act) {
      for (var o = 0; o < nout; o++) {
        var s = b ? b[o] : 0, wo = o * nin;
        for (var i = 0; i < nin; i++) s += w[wo + i] * x[i];
        out[o] = act ? silu(s) : s;
      }
      return out;
    }
    var gnn = !!J.gnn, P = { ae: [], be: [], m0w: [], m0b: [], m2w: [], m2b: [], nw: [], nb: [] }, k, l;
    P.d0w = t('desc.0.weight'); P.d0b = t('desc.0.bias'); P.d3w = t('desc.3.weight'); P.d3b = t('desc.3.bias');
    if (gnn) {
      P.eps = t('enc.eps'); P.ow = t('enc.out.0.weight'); P.ob = t('enc.out.0.bias');
      for (k = 0; k < nA; k++) P.ae.push(t('enc.ae.' + k + '.weight'));
      for (l = 0; l < L; l++) {
        var bl = []; for (k = 0; k < nB; k++) bl.push(t('enc.be.' + l + '.' + k + '.weight')); P.be.push(bl);
        P.m0w.push(t('enc.mlp.' + l + '.0.weight')); P.m0b.push(t('enc.mlp.' + l + '.0.bias'));
        P.m2w.push(t('enc.mlp.' + l + '.2.weight')); P.m2b.push(t('enc.mlp.' + l + '.2.bias'));
        P.nw.push(t('enc.nrm.' + l + '.weight')); P.nb.push(t('enc.nrm.' + l + '.bias'));
      }
    }
    P.iw = t('inter.weight'); P.f0w = t('f.0.weight'); P.f0b = t('f.0.bias'); P.f3w = t('f.3.weight'); P.f3b = t('f.3.bias'); P.f5w = t('f.5.weight'); P.f5b = t('f.5.bias');
    P.kappa = softplus(t('kappa')[0]); P.sve = J.svemb ? t('svemb.weight') : null;
    var nd = J.n_desc, dh = T['desc.0.weight'].shape[0], dz = T['desc.3.weight'].shape[0], zd = (gnn ? d : 0) + dz;
    var h1 = T['f.0.weight'].shape[0], h2 = T['f.3.weight'].shape[0], nin = T['f.0.weight'].shape[1], ns = J.svemb || 0;
    if (nin !== 3 * zd + 21 + ns) throw new Error('グラフのモデルの出口の入力の数が合わない');

    function encode(graph, desc) {
      if (desc.length !== nd) throw new Error('記述子の数が合わない（' + desc.length + ' と ' + nd + '）');
      var i, j, c, kk, ll, x = new Float32Array(nd);
      for (i = 0; i < nd; i++) {
        var v = desc[i]; if (v === null || v === undefined || !isFinite(v)) v = 0;
        v = (v < 0 ? -1 : 1) * Math.log(1 + Math.abs(v));
        x[i] = (v - N.d_mu[i]) / N.d_sd[i];
      }
      var zdv = linear(P.d3w, P.d3b, linear(P.d0w, P.d0b, x, nd, dh, new Float32Array(dh), true), dh, dz, new Float32Array(dz), true);
      var z = new Float32Array(zd);
      if (gnn) {
        if (!graph || !graph.n) throw new Error('分子グラフが空');
        var n = graph.n, nb = graph.nb, A = graph.atom, B = graph.bond;
        var h = new Float32Array(n * d), m = new Float32Array(n * d), e = new Float32Array(d), tn = new Float32Array(d), u1 = new Float32Array(2 * d), u2 = new Float32Array(d);
        for (i = 0; i < n; i++) for (kk = 0; kk < nA; kk++) { var ea = P.ae[kk], o = A[i * nA + kk] * d; for (c = 0; c < d; c++) h[i * d + c] += ea[o + c]; }
        for (ll = 0; ll < L; ll++) {
          m.fill(0);
          for (j = 0; j < nb; j++) {
            e.fill(0);
            for (kk = 0; kk < nB; kk++) { var eb = P.be[ll][kk], ob = B[j * nB + kk] * d; for (c = 0; c < d; c++) e[c] += eb[ob + c]; }
            var p = graph.a[j] * d, q = graph.b[j] * d;
            for (c = 0; c < d; c++) { m[q + c] += silu(h[p + c] + e[c]); m[p + c] += silu(h[q + c] + e[c]); }
          }
          var ep1 = 1 + P.eps[ll], nw = P.nw[ll], nbb = P.nb[ll];
          for (i = 0; i < n; i++) {
            var mean = 0, vr = 0, off = i * d;
            for (c = 0; c < d; c++) { tn[c] = ep1 * h[off + c] + m[off + c]; mean += tn[c]; }
            mean /= d;
            for (c = 0; c < d; c++) { var dv = tn[c] - mean; vr += dv * dv; }
            var inv = 1 / Math.sqrt(vr / d + 1e-5);
            for (c = 0; c < d; c++) tn[c] = (tn[c] - mean) * inv * nw[c] + nbb[c];
            linear(P.m0w[ll], P.m0b[ll], tn, d, 2 * d, u1, true);
            linear(P.m2w[ll], P.m2b[ll], u1, 2 * d, d, u2, false);
            for (c = 0; c < d; c++) h[off + c] += u2[c];
          }
        }
        var pool = new Float32Array(3 * d);
        for (c = 0; c < d; c++) pool[2 * d + c] = -Infinity;
        for (i = 0; i < n; i++) for (c = 0; c < d; c++) { var hv = h[i * d + c]; pool[c] += hv; if (hv > pool[2 * d + c]) pool[2 * d + c] = hv; }
        for (c = 0; c < d; c++) { pool[d + c] = pool[c] / n; pool[c] = pool[c] / 8.0; }
        z.set(linear(P.ow, P.ob, pool, 3 * d, d, new Float32Array(d), true), 0);
      }
      z.set(zdv, gnn ? d : 0);
      return z;
    }

    var xin = new Float32Array(nin), a1 = new Float32Array(h1), a2 = new Float32Array(h2), a3 = new Float32Array(2), iz = new Float32Array(zd);
    function predict(zu, desc, svId, TmK, TK) {
      var S = J.solvents[svId]; if (!S) throw new Error('溶媒が無い: ' + svId);
      var zv = S.z, i, o = 0;
      linear(P.iw, null, zv, zd, zd, iz, false);
      xin.set(zu, 0); o += zd; xin.set(zv, o); o += zd;
      for (i = 0; i < zd; i++) xin[o + i] = zu[i] * iz[i]; o += zd;
      var has = S.ab ? 1 : 0, ab = [], sc = S.ab || [0, 0, 0, 0, 0, 0], lser = 0;
      for (i = 0; i < 5; i++) { var v = desc[nd - 5 + i]; ab.push(isFinite(v) && v !== null ? v : 0); }
      for (i = 0; i < 6; i++) xin[o + i] = has ? (sc[i] - N.s_mu[i]) / N.s_sd[i] : 0; o += 6;
      xin[o++] = has;
      lser = sc[0] * has;
      for (i = 0; i < 5; i++) { var pr = ab[i] * sc[i + 1] * has; lser += pr; xin[o + i] = pr / 3; } o += 5;
      xin[o++] = lser / 3;
      for (i = 0; i < 5; i++) xin[o + i] = ab[i] / 3; o += 5;
      var tk = (TmK !== null && TmK !== undefined && isFinite(TmK) && J.mode === 'tm') ? TmK : null, h = tk === null ? 0 : 1;
      xin[o++] = h; xin[o++] = h ? Math.max(0, tk - 298.15) / 100 : 0; xin[o++] = h ? -2.951 * Math.max(0, tk / 298.15 - 1) : 0;
      if (ns) { var se = P.sve, r = S.svid * ns; for (i = 0; i < ns; i++) xin[o + i] = se[r + i]; o += ns; }
      linear(P.f0w, P.f0b, xin, nin, h1, a1, true);
      linear(P.f3w, P.f3b, a1, h1, h2, a2, true);
      linear(P.f5w, P.f5b, a2, h2, 2, a3, false);
      var T = TK || 298.15, tt = 1000 / T - 1000 / 298.15, ideal = h ? -2.951 * Math.max(0, tk / T - 1) : 0;
      return a3[0] + a3[1] * tt + P.kappa * ideal;
    }
    return { json: J, zdim: zd, encode: encode, predict: predict };
  }

  var api = { VERSION: 1, open: open };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.SolGNN = api;
})(typeof self !== 'undefined' ? self : this);
