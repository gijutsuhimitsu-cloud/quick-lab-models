// SPDX-License-Identifier: MIT
// 溶けるか・混ざるか：学習済みモデル（勾配ブースティングの木）の入力と計算。学習の書き出し（Node）とページ（ブラウザ）で同じファイルを使う。
// 入力の並びを変えたら、モデルを学習し直すこと（train/build_tables.cjs → train/train_gbt.py）。
(function (root) {
  'use strict';
  var SC = (typeof module !== 'undefined' && module.exports) ? require('./solub.js') : root.SolCore;
  var AB = (typeof module !== 'undefined' && module.exports) ? require('./abraham.js') : root.SolAbraham;

  // 入力の名前。sx は features.js の desc 41 個 ＋ 官能基の数 56 個（names は features.js の names()）。
  function names(sxNames) {
    return sxNames.concat(['abr_E', 'abr_S', 'abr_A', 'abr_B', 'abr_V', 'logPow', 'crystal', 'sv_c', 'sv_e', 'sv_s', 'sv_a', 'sv_b', 'sv_v', 'sv_water', 'transfer', 'formula', 'formula0']);
  }

  // 1 行の入力。solute = {x:[97], abr:{E,S,A,B,V}, logPow, M, TmC（なければ null）}、sv = 溶媒の係数 {c,e,s,a,b,v, water}、K = 式の係数、tmDefault = 融点が分からない時に置く値。
  // crystal と formula は融点を使う値（融点がなければ 0 と formula0）、formula0 は融点を tmDefault と置いた式の値。
  function build(solute, sv, K, tmDefault, TC) {
    var tr = sv.water ? 0 : AB.logP(sv, solute.abr), has = solute.TmC !== null && solute.TmC !== undefined && isFinite(solute.TmC);
    var arg = sv.water ? { water: true } : { transfer: tr };
    var f0 = SC.predict({ M: solute.M, logPow: solute.logPow, TmC: null }, arg, TC, K, tmDefault).raw;
    var f1 = has ? SC.predict(solute, arg, TC, K, tmDefault).raw : f0;
    return solute.x.concat([solute.abr.E, solute.abr.S, solute.abr.A, solute.abr.B, solute.abr.V, solute.logPow, has ? SC.crystal(solute.TmC, TC) : 0,
      sv.c, sv.e, sv.s, sv.a, sv.b, sv.v, sv.water ? 1 : 0, tr, f1, f0]);
  }

  // 木のモデルの計算。model = {init, lr, cols:[入力の列の番号], trees:[{f,t,l,r,v}]}（f は cols の中の番号、葉は l = -1）。
  function evalModel(model, row) {
    var y = model.init, T = model.trees, cols = model.cols;
    for (var i = 0; i < T.length; i++) {
      var t = T[i], n = 0;
      if (t.l) { while (t.l[n] >= 0) n = row[cols[t.f[n]]] <= t.t[n] ? t.l[n] : t.r[n]; }
      else { while (t.r[n] >= 0) n = row[cols[t.f[n]]] <= t.t[n] ? n + 1 : t.r[n]; }   // l を省いた形：左の子は次の番号
      y += model.lr * t.v[n];
    }
    return y;
  }

  var api = { names: names, build: build, evalModel: evalModel };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.SolFeat = api;
})(typeof self !== 'undefined' ? self : this);
