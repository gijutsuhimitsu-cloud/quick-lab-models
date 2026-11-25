// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE。
// 混合物のクロマトグラム（ピークの並び・重なり・順番の確かさ・曲線）を、保持の予測から作る。
// retention_model.js の P（ChromRetention.create の戻り値）を使う。条件を動かすたびに呼んでよい（化合物 10 個で 2 ミリ秒ほど）。
//
//   var m = ChromMixture.run(P, compounds, columnId, cond, { areas: [1, 1, 0.5] });
//     compounds : 化合物の配列（stock_compounds.json の 1 件、または features.js で作った物。z や cal が付いていればそれも使われる）
//     戻り値 m = {
//       t0, tMax（横軸の終わりの目安。分）,
//       peaks: [ { index（compounds の中の番号）, name, tR, sigma, k（等組成の時）, lkw, model, calibrated, inRange,
//                  tR68:[早い側, 遅い側], tR90:[..]（予測の幅がある時）, area, height } ... ]   // 保持時間の早い順
//       pairs: [ { a, b（peaks の番号。隣どうし）, Rs（分離度）, orderProb（この順番で出る確からしさ。0.5〜1。幅が無い時は null）, bandsOverlap（7 割の幅が重なる） } ... ],
//       minRs, uncertainPairs（orderProb が 0.8 未満の隣どうしの数）
//     }
//   var y = ChromMixture.trace(m, { tMax: m.tMax, n: 600 });      // Float64Array。ピーク（ガウス型）を足し合わせた曲線。y[i] は t = i·tMax/(n − 1) の値
//   var yi = ChromMixture.peakTrace(m.peaks[0], { tMax: m.tMax, n: 600 });   // 1 本だけの曲線（色を分けて描く時）
//
// 注意（ページに書くこと）:
//   - 保持時間は予測で、未知の分子では分析時間の 2〜5% ほど外れる。近い 2 本の順番や分離は、予測では決められない。
//     orderProb は「2 つの予測の外れが独立で、幅が正規分布の 1σ にあたる」とした時の目安。
//   - ピークの幅は、カラムの長さと粒径から決めた段数の目安（換算段高 2.5）によるガウス型。テーリングや過負荷は入っていない。
//   - ピークの高さ・面積は、濃度や検出の感度を表していない（areas を渡さなければ、どの成分も同じ面積）。
(function (root) {
  'use strict';
  var SQ2PI = Math.sqrt(2 * Math.PI);

  // 正規分布の累積（誤差関数の近似。絶対誤差 1.5e-7）
  function phi(x) {
    var z = Math.abs(x) / Math.SQRT2, t = 1 / (1 + 0.3275911 * z);
    var e = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-z * z);
    return x >= 0 ? 0.5 * (1 + e) : 0.5 * (1 - e);
  }

  function run(P, compounds, columnId, cond, opts) {
    opts = opts || {};
    var peaks = [], i;
    for (i = 0; i < compounds.length; i++) {
      var c = compounds[i], r = P.predictRetention(c, columnId, cond);
      var area = opts.areas && opts.areas[i] !== undefined ? opts.areas[i] : 1;
      var p = { index: i, name: c.name || c.smiles || ('#' + (i + 1)), tR: r.tR, sigma: r.sigma, lkw: r.lkw, model: r.model, calibrated: !!r.calibration, inRange: r.inRange,
        area: area, height: area / (r.sigma * SQ2PI), t0: r.t0 };
      if (r.k !== undefined) p.k = r.k;
      if (r.uncertainty) { p.tR68 = r.uncertainty.tR68; p.tR90 = r.uncertainty.tR90; p.seen = r.uncertainty.seen; }
      peaks.push(p);
    }
    peaks.sort(function (a, b) { return a.tR - b.tR; });
    var pairs = [], minRs = Infinity, unsure = 0, tMax = 0;
    for (i = 0; i < peaks.length; i++) {
      tMax = Math.max(tMax, peaks[i].tR + 5 * peaks[i].sigma);
      if (i === 0) continue;
      var a = peaks[i - 1], b = peaks[i], dt = b.tR - a.tR;
      var q = { a: i - 1, b: i, Rs: dt / (2 * (a.sigma + b.sigma)), orderProb: null, bandsOverlap: null };
      if (a.tR68 && b.tR68) {
        var sa = (a.tR68[1] - a.tR68[0]) / 2, sb = (b.tR68[1] - b.tR68[0]) / 2, sd = Math.sqrt(sa * sa + sb * sb);
        q.orderProb = sd > 0 ? phi(dt / sd) : 1;
        q.bandsOverlap = a.tR68[1] > b.tR68[0];
        if (q.orderProb < 0.8) unsure++;
      }
      if (q.Rs < minRs) minRs = q.Rs;
      pairs.push(q);
    }
    var t0 = peaks.length ? peaks[0].t0 : null;
    if (cond && cond.gradient && cond.gradient.length) tMax = Math.max(tMax, cond.gradient[cond.gradient.length - 1][0] + (cond.dwell || 0));
    return { t0: t0, tMax: tMax * 1.05, peaks: peaks, pairs: pairs, minRs: pairs.length ? minRs : null, uncertainPairs: unsure };
  }

  function peakTrace(p, o, into) {
    var n = o.n || 600, tMax = o.tMax, y = into || new Float64Array(n), dt = tMax / (n - 1);
    var lo = Math.max(0, Math.floor((p.tR - 6 * p.sigma) / dt)), hi = Math.min(n - 1, Math.ceil((p.tR + 6 * p.sigma) / dt));
    for (var i = lo; i <= hi; i++) { var z = (i * dt - p.tR) / p.sigma; y[i] += p.height * Math.exp(-0.5 * z * z); }
    return y;
  }
  function trace(m, o) {
    o = o || {}; var oo = { n: o.n || 600, tMax: o.tMax || m.tMax }, y = new Float64Array(oo.n);
    for (var i = 0; i < m.peaks.length; i++) peakTrace(m.peaks[i], oo, y);
    return y;
  }

  var api = { run: run, trace: trace, peakTrace: peakTrace, normalCdf: phi };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.ChromMixture = api;
})(typeof self !== 'undefined' ? self : this);
