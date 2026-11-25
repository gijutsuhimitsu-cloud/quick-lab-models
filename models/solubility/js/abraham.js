// SPDX-License-Identifier: MIT
// 溶けるか・混ざるか：Abraham の溶媒和の式（線形自由エネルギー関係）の部品。学習（Node）とページ（ブラウザ）で同じファイルを使う。
//   log10 P(水 → 溶媒) = c + e·E + s·S + a·A + b·B + v·V
//   E 過剰モル屈折、S 双極性・分極率、A 水素結合の酸性度、B 水素結合の塩基性度、V McGowan の特性体積（cm³ mol⁻¹ ÷ 100）
// 溶媒の係数 c・e・s・a・b・v は weights/solvents.json（CC0 の表から）。
// 溶質の E・S・A・B は weights/abraham_solute_model.json の一次式で推算する（CC BY 4.0 の表に当てはめた物）。V は構造から計算する。
(function (root) {
  'use strict';

  // McGowan の原子の体積（cm³ mol⁻¹）。結合 1 本につき 6.56 を引く（結合次数に関わらず 1 本と数える。H への結合も数える）。
  // Abraham, McGowan, Chromatographia 23, 243 (1987)。値は tools/fit_abraham.cjs が表の V の列と照合する。
  var MCG = { C: 16.35, H: 8.71, N: 14.39, O: 12.43, F: 10.48, Si: 26.83, P: 24.87, S: 22.91, Cl: 20.95, B: 18.32, Br: 26.21, I: 34.53, Se: 27.81, As: 29.42, Ge: 31.02, Sn: 39.35 };
  var BOND = 6.56;

  // H を付けた構造から、元素ごとの数と結合の数を数える。読めなければ null。
  function composition(RDKit, smiles) {
    var m = null;
    try { m = RDKit.get_mol(smiles); } catch (e) { m = null; }
    if (!m || (m.is_valid && !m.is_valid())) { if (m) m.delete(); return null; }
    var L;
    try { m.add_hs_in_place(); L = m.get_molblock().split('\n'); } finally { m.delete(); }
    var na = +L[3].slice(0, 3), nb = +L[3].slice(3, 6), el = {};
    if (!(na > 0)) return null;
    for (var i = 0; i < na; i++) { var s = L[4 + i].slice(31, 34).trim(); el[s] = (el[s] || 0) + 1; }
    return { atoms: el, nAtoms: na, nBonds: nb };
  }

  // McGowan の特性体積 V（cm³ mol⁻¹ ÷ 100）。表にない元素を含む時は null。
  function mcgowanV(RDKit, smiles) {
    var c = composition(RDKit, smiles);
    if (!c) return null;
    var v = 0;
    for (var k in c.atoms) { if (!(k in MCG)) return null; v += MCG[k] * c.atoms[k]; }
    return (v - BOND * c.nBonds) / 100;
  }

  // 一次式のモデル（標準化した特徴量 × 重み ＋ 切片）で E・S・A・B を推算する。x は ChromFeatures.compute の desc と groups を並べた配列。
  // A と B は負にならない量なので、0 で切る。
  function predictSolute(model, x) {
    var out = {}, T = model.targets;
    for (var t in T) {
      var w = T[t].w, y = T[t].b;
      for (var j = 0; j < w.length; j++) { var sd = model.sd[j]; if (sd > 0) y += w[j] * (x[j] - model.mean[j]) / sd; }
      out[t] = y;
    }
    if (out.A < 0) out.A = 0;
    if (out.B < 0) out.B = 0;
    return out;
  }

  // log10 P(水 → 溶媒)。solv は {c,e,s,a,b,v}、d は {E,S,A,B,V}。
  function logP(solv, d) {
    return solv.c + solv.e * d.E + solv.s * d.S + solv.a * d.A + solv.b * d.B + solv.v * d.V;
  }

  var api = { MCG: MCG, BOND: BOND, composition: composition, mcgowanV: mcgowanV, predictSolute: predictSolute, logP: logP };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.SolAbraham = api;
})(typeof self !== 'undefined' ? self : this);
