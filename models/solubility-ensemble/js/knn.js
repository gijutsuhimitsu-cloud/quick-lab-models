// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE（MIT）。
// 溶けるか・混ざるか：似た分子の外れで直す補正（水のモデル用）。ページ（ブラウザ）と確かめ（Node）で同じファイルを使う。
// 材料は weights/knn_water.json（学習の道具（このリポジトリには未同梱） が書く）。計算は 学習の道具（このリポジトリには未同梱） の knn_corr と同じ：
//   指紋（Morgan、半径 2・1024 ビット）の Tanimoto で上位 K 個を選び、似ている度合いが smin 以上の物だけ、その pow 乗で重みづけした外れの和を、(重みの和 + lam) で割る。
(function (root) {
  'use strict';
  function b64(s) {
    if (typeof Buffer !== 'undefined') { var b = Buffer.from(s, 'base64'); return new Uint8Array(b.buffer, b.byteOffset, b.length); }
    var bin = atob(s), u = new Uint8Array(bin.length); for (var i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); return u;
  }
  // 材料を読める形にする
  function load(j) {
    var cnt = b64(j.cnt), bitsB = b64(j.bits), dv = new DataView(bitsB.buffer, bitsB.byteOffset, bitsB.byteLength);
    var bits = new Uint16Array(bitsB.byteLength / 2); for (var i = 0; i < bits.length; i++) bits[i] = dv.getUint16(2 * i, true);
    var off = new Uint32Array(j.n + 1); for (var k = 0; k < j.n; k++) off[k + 1] = off[k] + cnt[k];
    function res(s) { var u = b64(s), d = new DataView(u.buffer, u.byteOffset, u.byteLength), r = new Float64Array(j.n); for (var q = 0; q < j.n; q++) { var v = d.getInt16(2 * q, true); r[q] = v === -32768 ? NaN : v / 1000; } return r; }
    return { n: j.n, K: j.K, pow: j.pow, lam: j.lam, smin: j.smin, nBits: j.nBits, radius: j.radius, cnt: cnt, bits: bits, off: off, r: { tm: res(j.r_tm), notm: res(j.r_notm) } };
  }
  // RDKit.js の get_morgan_fp が返す '0101…' から、立っているビットの番号
  function onBits(str) { var a = []; for (var i = 0; i < str.length; i++) if (str.charCodeAt(i) === 49) a.push(i); return a; }
  // 補正の値と、使った近い溶質の数。key は 'tm' か 'notm'。exclude は除く溶質の番号（確かめ用）
  function corr(M, on, key, exclude) {
    var q = new Uint8Array(M.nBits), cq = on.length, R = M.r[key], top = [], i, k;
    for (i = 0; i < on.length; i++) q[on[i]] = 1;
    for (k = 0; k < M.n; k++) {
      if (k === exclude || isNaN(R[k])) continue;
      var inter = 0; for (i = M.off[k]; i < M.off[k + 1]; i++) inter += q[M.bits[i]];
      var s = inter / Math.max(1, cq + M.cnt[k] - inter);
      if (top.length < M.K || s > top[top.length - 1][0]) {          // 同じ値なら先の番号を残す（Python の安定な並べ替えと同じ）
        var p = top.length; while (p > 0 && top[p - 1][0] < s) p--;
        top.splice(p, 0, [s, k]); if (top.length > M.K) top.pop();
      }
    }
    var num = 0, den = 0, used = 0;
    for (i = 0; i < top.length; i++) { if (top[i][0] < M.smin) continue; var w = Math.pow(top[i][0], M.pow); num += w * R[top[i][1]]; den += w; used++; }
    return { c: num / (den + M.lam), used: used, best: top.length ? top[0][0] : 0 };
  }
  var api = { load: load, onBits: onBits, corr: corr };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.SolKNN = api;
})(typeof self !== 'undefined' ? self : this);
