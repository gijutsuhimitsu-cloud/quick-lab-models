// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE（MIT）。
// 溶けるか・混ざるか：溶解度の見込みを出す式。学習・検証（Node）とページ（ブラウザ）で同じファイルを使う。
// 単位：溶解度は log10(mol L⁻¹)、温度は °C。どの値も「見込み」で、幅（band）と一緒に使う。
//
//   固体の項   crystal(Tm, T) = −ΔS_f (Tm − T) / (2.303 R T)     理想溶解度（モル分率）の対数。ΔS_f は一定値で近似（Walden の規則）
//   水         log S_w = 0.5 − log P_ow + crystal                 一般溶解度式（Yalkowsky）。25 °C では crystal ≈ −0.01 (Tm − 25)
//   有機溶媒   log S_s = k0 + kP(−log P_ow) + kC·crystal + kT·transfer   transfer = log P(水 → 溶媒)（Abraham の式、abraham.js）。係数 4 つは実測に合わせた値
//   混合溶媒   log S_mix = Σ f_i log S_i + Σ f_i f_j [A0 + A1 (f_i − f_j)]   対数線形＋溶媒の組ごとの補正（f は体積分率）
//   実測 1 点  log S_b = log S_a + (transfer_b − transfer_a)          固体の項が消える。同じ固体・中性の分子・薄い溶液が前提
(function (root) {
  'use strict';
  var R = 8.314462, DSF = 56.5;          // J mol⁻¹ K⁻¹。ΔS_f = 56.5 は剛直な有機分子の目安（Walden の規則）
  var K0 = 273.15;

  // 理想溶解度の対数（固体の項）。融点が測る温度より低い（液体）なら 0。
  function crystal(TmC, TC) {
    if (TmC === null || TmC === undefined || !isFinite(TmC)) return null;
    var Tm = TmC + K0, T = TC + K0;
    return Tm <= T ? 0 : -DSF * (Tm - T) / (2.302585 * R * T);
  }

  // 水への溶解度（中性の形）。logPow は 1-オクタノール／水の分配係数の対数。
  function logSwater(logPow, TmC, TC) {
    var c = crystal(TmC, TC);
    return c === null ? null : 0.5 - logPow + c;
  }

  // 純粋な溶質のモル濃度の対数（これより濃くはならない）。密度は 1.2 g cm⁻³ と置く。
  function logSmax(M) { return Math.log10(1200 / M); }

  // 1 つの溶媒での見込み（実測で合わせていない時）。solute = {M, logPow, TmC}、sv = {water:true} か {transfer: log P(水 → 溶媒)}。
  //   水      log S = 0.5 − log Pow + crystal                                   （一般溶解度式そのまま）
  //   有機溶媒 log S = k0 + kP·(−log Pow) + kC·crystal + kT·transfer             （K は weights/base.json。BigSolDB の実測に合わせた 4 つの係数）
  // 融点が分からない時は tmDefault（実測のある固体の融点の中央値）を使い、幅を広く取る。返す parts は内訳。
  function predict(solute, sv, TC, K, tmDefault) {
    var known = solute.TmC !== null && solute.TmC !== undefined && isFinite(solute.TmC);
    var c = crystal(known ? solute.TmC : tmDefault, TC), cap = logSmax(solute.M), parts, raw;
    if (sv.water) { parts = { base: 0.5, polarity: -solute.logPow, crystal: c, transfer: 0 }; }
    else { parts = { base: K.k0, polarity: -K.kP * solute.logPow, crystal: K.kC * c, transfer: K.kT * sv.transfer }; }
    raw = parts.base + parts.polarity + parts.crystal + parts.transfer;
    return { logS: Math.min(raw, cap), raw: raw, capped: raw > cap, tmKnown: known, parts: parts };
  }

  // 実測 1 点から、ほかの溶媒を見込む：log S_b = log S_a + (transfer_b − transfer_a)。水の transfer は 0。
  function fromReference(logSa, transferA, transferB, M) {
    var raw = logSa + (transferB - transferA), cap = logSmax(M);
    return { logS: Math.min(raw, cap), raw: raw, capped: raw > cap };
  }

  // 混合溶媒（対数線形）。items = [{logS, f}]、f は体積分率（和が 1 でなくても割って使う）。
  function mix(items) {
    var s = 0, w = 0;
    for (var i = 0; i < items.length; i++) { s += items[i].f * items[i].logS; w += items[i].f; }
    return w > 0 ? s / w : null;
  }

  // 混合溶媒の補正（対数線形からのずれ）。comp = {溶媒の id: 体積分率}、MIX = weights/base.json の {pairs, cls}。
  // 溶媒の組ごとに f_i f_j [A0 + A1 (f_i − f_j)] を足す（組の名前は id の辞書順で i|j）。3 成分は 2 成分の組の和で近似する。
  // 返す w は、ずれの大きさの目安（2 成分のまんなかで 1、純粋な溶媒で 0）。est は、実測が足りない組の見込み（MIX.est）を使ったかどうか。
  // 係数の選び方：実測の組（MIX.pairs）→ 溶媒の係数の似ていなさに比例させた見込み（MIX.est）→ 組の種類の平均（MIX.cls）。
  function mixCorr(comp, MIX) {
    var ids = Object.keys(comp).filter(function (k) { return comp[k] > 0; }).sort(), s = 0, w = 0, est = false;
    for (var i = 0; i < ids.length; i++) for (var j = i + 1; j < ids.length; j++) {
      var key = ids[i] + "|" + ids[j], fi = comp[ids[i]], fj = comp[ids[j]], A = MIX.pairs[key];
      if (!A && MIX.est && MIX.est[key]) { A = MIX.est[key]; est = true; }
      if (!A) A = MIX.cls[(ids[i] === "water" || ids[j] === "water") ? "wo" : "oo"];
      s += fi * fj * (A[0] + A[1] * (fi - fj)); w += 4 * fi * fj;
    }
    return { corr: s, w: Math.min(1, w), est: est };
  }

  // 原液を薄めた後の濃度と、薄めた後の溶媒の組成（体積は足し算で近似）。
  // stockUL の原液（濃度 stockConc）に、diluentUL の薄める液を足す。diluent = [{key, f}]（体積分率）。返す comp は溶媒ごとの体積分率。
  function dilute(stockKey, stockConc, stockUL, diluent, diluentUL) {
    var total = stockUL + diluentUL, comp = {}, w = 0, i;
    if (!(total > 0)) return null;
    comp[stockKey] = stockUL / total;
    for (i = 0; i < diluent.length; i++) w += diluent[i].f;
    for (i = 0; i < diluent.length; i++) comp[diluent[i].key] = (comp[diluent[i].key] || 0) + (diluentUL / total) * (diluent[i].f / w);
    return { conc: stockConc * stockUL / total, factor: total / stockUL, comp: comp };
  }

  // 目標の濃度と見込みの幅から、3 段の判定を返す。band は log10 の片側の幅（外れの約 7 割が入る幅）。
  //   'yes' 幅の下の端でも目標より上／'no' 幅の上の端でも目標より下／'edge' 目標が幅の中
  function verdict(logS, logTarget, band) {
    if (logS - band >= logTarget) return 'yes';
    if (logS + band < logTarget) return 'no';
    return 'edge';
  }

  var toMolL = function (mgPerML, M) { return mgPerML / M; };          // mg mL⁻¹ → mol L⁻¹
  var toMgML = function (molL, M) { return molL * M; };

  var api = { R: R, DSF: DSF, crystal: crystal, logSwater: logSwater, logSmax: logSmax, predict: predict, fromReference: fromReference, mix: mix, mixCorr: mixCorr, dilute: dilute, verdict: verdict, toMolL: toMolL, toMgML: toMgML };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.SolCore = api;
})(typeof self !== 'undefined' ? self : this);
