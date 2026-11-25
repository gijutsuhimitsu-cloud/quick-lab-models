// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE。
// 小さなグラフのモデル（分子グラフ → lkw）を、ブラウザの中で直接計算する。外部の推論ライブラリは使わない。
//
//   var G = ChromGNN.create(manifest, weights);          // manifest: gnn_<name>.json、weights: gnn_<name>.bin の ArrayBuffer（float32）
//   var z = G.encode(graph, desc127);                    // 分子ごとに 1 回だけ（graph は ChromGraph.fromSmiles、desc127 は features.js の desc・groups・ion をつないだ 127 個）
//   var lkw = G.lkw(z, method16, ion10);                 // カラム・条件を変えるたび（軽い）。method16・ion10 は retention_model.js が作る並びと同じ
//
// 出口の形: lkw = f0(z, c, u) + Σ_j s_j · softplus(g_j(z, c) − 2) · q_j
//   q_j は向きの決まったイオン化の量（s_j = ±1）。係数が必ず正なので、酸・塩基の pH の向きは必ず守られる。
(function (root) {
  'use strict';

  function create(manifest, weights) {
    var W = (weights instanceof Float32Array) ? weights : new Float32Array(weights);
    if (W.length !== manifest.n_params) throw new Error('グラフのモデルの重みの大きさが、manifest と合わない（' + W.length + ' と ' + manifest.n_params + '）');
    var T = manifest.tensors, d = manifest.hid, L = manifest.layers, N = manifest.norm;
    var AD = manifest.atom_dims, BD = manifest.bond_dims, nA = AD.length, nB = BD.length;
    function t(name) {
      var e = T[name]; if (!e) throw new Error('重みが無い: ' + name);
      var n = 1; for (var i = 0; i < e.shape.length; i++) n *= e.shape[i];
      return W.subarray(e.offset, e.offset + n);
    }
    function silu(x) { return x / (1 + Math.exp(-x)); }
    function softplus(x) { return x > 30 ? x : Math.log(1 + Math.exp(x)); }
    // y = Wx + b（W は [出力, 入力] を行ごとに並べた物）。x の off から nin 個を読む
    function linear(w, b, x, xoff, nin, nout, out, ooff, act) {
      for (var o = 0; o < nout; o++) {
        var s = b[o], wo = o * nin;
        for (var i = 0; i < nin; i++) s += w[wo + i] * x[xoff + i];
        out[ooff + o] = act ? silu(s) : s;
      }
    }
    var P = {
      d0w: t('desc.0.weight'), d0b: t('desc.0.bias'), d3w: t('desc.3.weight'), d3b: t('desc.3.bias'),
      eps: t('enc.eps'), ow: t('enc.out.0.weight'), ob: t('enc.out.0.bias'),
      f0w: t('f0.0.weight'), f0b: t('f0.0.bias'), f3w: t('f0.3.weight'), f3b: t('f0.3.bias'), f5w: t('f0.5.weight'), f5b: t('f0.5.bias'),
      g0w: t('g.0.weight'), g0b: t('g.0.bias'), g2w: t('g.2.weight'), g2b: t('g.2.bias'), ae: [], be: [], m0w: [], m0b: [], m2w: [], m2b: [], nw: [], nb: []
    };
    var k, l;
    for (k = 0; k < nA; k++) P.ae.push(t('enc.ae.' + k + '.weight'));
    for (l = 0; l < L; l++) {
      var bl = []; for (k = 0; k < nB; k++) bl.push(t('enc.be.' + l + '.' + k + '.weight')); P.be.push(bl);
      P.m0w.push(t('enc.mlp.' + l + '.0.weight')); P.m0b.push(t('enc.mlp.' + l + '.0.bias'));
      P.m2w.push(t('enc.mlp.' + l + '.2.weight')); P.m2b.push(t('enc.mlp.' + l + '.2.bias'));
      P.nw.push(t('enc.nrm.' + l + '.weight')); P.nb.push(t('enc.nrm.' + l + '.bias'));
    }
    var nd = manifest.n_desc, dh = T['desc.0.weight'].shape[0], dz = T['desc.3.weight'].shape[0], zd = d + dz;
    var nc = N.m_mu.length, jc = manifest.jc, ju = manifest.ju, nq = jc.length, nu = ju.length;
    var h1 = T['f0.0.weight'].shape[0], h2 = T['f0.3.weight'].shape[0], gh = T['g.0.weight'].shape[0];
    var sg = jc.map(function (j) { return manifest.ion_signs[j]; });

    function encode(graph, desc) {
      if (!graph || !graph.n) throw new Error('分子グラフが空');
      if (desc.length !== nd) throw new Error('記述子の数が合わない（' + desc.length + ' と ' + nd + '）');
      var n = graph.n, nb = graph.nb, A = graph.atom, B = graph.bond, i, j, c, kk, ll;
      // 記述子の枝
      var x = new Float32Array(nd);
      for (i = 0; i < nd; i++) {
        var v = desc[i]; if (v === null || v === undefined || !isFinite(v)) v = 0;
        v = (v < 0 ? -1 : 1) * Math.log(1 + Math.abs(v));
        x[i] = (v - N.d_mu[i]) / N.d_sd[i];
      }
      var t1 = new Float32Array(dh), zdv = new Float32Array(dz);
      linear(P.d0w, P.d0b, x, 0, nd, dh, t1, 0, true);
      linear(P.d3w, P.d3b, t1, 0, dh, dz, zdv, 0, true);
      // グラフの枝
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
          linear(P.m0w[ll], P.m0b[ll], tn, 0, d, 2 * d, u1, 0, true);
          linear(P.m2w[ll], P.m2b[ll], u1, 0, 2 * d, d, u2, 0, false);
          for (c = 0; c < d; c++) h[off + c] += u2[c];
        }
      }
      var pool = new Float32Array(3 * d);
      for (c = 0; c < d; c++) pool[2 * d + c] = -Infinity;
      for (i = 0; i < n; i++) for (c = 0; c < d; c++) { var hv = h[i * d + c]; pool[c] += hv; if (hv > pool[2 * d + c]) pool[2 * d + c] = hv; }
      for (c = 0; c < d; c++) { pool[d + c] = pool[c] / n; pool[c] = pool[c] / 8.0; }
      var z = new Float32Array(zd);
      linear(P.ow, P.ob, pool, 0, 3 * d, d, z, 0, true);
      z.set(zdv, d);
      return z;
    }

    var xin = new Float32Array(zd + nc + nu), a1 = new Float32Array(h1), a2 = new Float32Array(h2), a3 = new Float32Array(1), b1 = new Float32Array(gh), b2 = new Float32Array(nq);
    function lkw(z, method, ion) {
      var i;
      xin.set(z, 0);
      for (i = 0; i < nc; i++) { var mv = method[i]; if (!isFinite(mv)) mv = 0; xin[zd + i] = (mv - N.m_mu[i]) / N.m_sd[i]; }
      for (i = 0; i < nu; i++) { var uv = ion[ju[i]]; if (!isFinite(uv)) uv = 0; xin[zd + nc + i] = (uv - N.u_mu[i]) / N.u_sd[i]; }
      linear(P.f0w, P.f0b, xin, 0, zd + nc + nu, h1, a1, 0, true);
      linear(P.f3w, P.f3b, a1, 0, h1, h2, a2, 0, true);
      linear(P.f5w, P.f5b, a2, 0, h2, 1, a3, 0, false);
      linear(P.g0w, P.g0b, xin, 0, zd + nc, gh, b1, 0, true);
      linear(P.g2w, P.g2b, b1, 0, gh, nq, b2, 0, false);
      var out = a3[0];
      for (i = 0; i < nq; i++) { var qv = ion[jc[i]]; if (!isFinite(qv)) qv = 0; out += softplus(b2[i] - 2.0) * sg[i] * (qv / N.q_sd[i]); }
      return out;
    }

    return { manifest: manifest, zdim: zd, encode: encode, lkw: lkw };
  }

  var api = { VERSION: 1, create: create };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.ChromGNN = api;
})(typeof self !== 'undefined' ? self : this);
