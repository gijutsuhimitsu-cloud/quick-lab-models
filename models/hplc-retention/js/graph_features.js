// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE。
// 分子グラフ（原子と結合の表）を、RDKit.js の分子から作る。グラフのモデル（gnn_model.js）の入力。
// 学習データのグラフも、この同じ関数を Node で動かして作っている（ブラウザと学習で、特徴量の作り方を完全に同じにするため）。
//
//   var g = ChromGraph.fromSmiles(RDKit, 'CC(C)Cc1ccc(cc1)C(C)C(=O)O');
//   g = { n, atom: Uint8Array(n * 8), nb, bond: Uint8Array(nb * 3), a: Int32Array(nb), b: Int32Array(nb) }   // 結合は向きなしで 1 本ずつ
//   失敗（読めない構造）なら null。
//
// 原子の特徴量（8 個。どれも小さな整数の番号）:
//   0 元素（C N O S F Cl Br I P B Si Se H Na K As の順に 0〜15、その他 16）  1 結合の数（0〜6）  2 電荷（−2〜+2 → 0〜4）
//   3 付いている水素の数（0〜4）  4 芳香族（0/1）  5 入っている環のうち最小の大きさ（無し 0、3〜7 → 1〜5、8 以上 6）
//   6 入っている環の数（0〜3）  7 立体中心の表示（無し 0、R 1、S 2）
// 結合の特徴量（3 個）:
//   0 種類（単 0、二重 1、三重 2、芳香族 3、その他 4）  1 環の中（0/1）  2 二重結合の立体の表示（無し 0、cis 1、trans 2）
// 立体の表示が入力の書き方で変わらないよう、正規化した SMILES に直してから読む。
(function (root) {
  'use strict';
  var VERSION = 1;
  var ELEM = [6, 7, 8, 16, 9, 17, 35, 53, 15, 5, 14, 34, 1, 11, 19, 33];
  var ATOM_DIMS = [17, 7, 5, 5, 2, 7, 4, 3], BOND_DIMS = [5, 2, 3];

  function fromMol(RDKit, mol) {
    var can = null;
    try {
      can = RDKit.get_mol(mol.get_smiles());
      if (!can || (can.is_valid && !can.is_valid())) { if (can) can.delete(); return null; }
      can.remove_hs_in_place();
      var J = JSON.parse(can.get_json());
      var M = J.molecules[0], dA = J.defaults.atom, dB = J.defaults.bond;
      var ext = null, i, k;
      for (i = 0; i < (M.extensions || []).length; i++) if (M.extensions[i].name === 'rdkitRepresentation') ext = M.extensions[i];
      ext = ext || {};
      var atoms = M.atoms || [], bonds = M.bonds || [], n = atoms.length, nb = bonds.length;
      var deg = new Int32Array(n), aro = new Uint8Array(n), nring = new Int32Array(n), minring = new Int32Array(n), cip = new Uint8Array(n);
      for (i = 0; i < (ext.aromaticAtoms || []).length; i++) aro[ext.aromaticAtoms[i]] = 1;
      var rings = ext.atomRings || [], ringPair = {};
      for (i = 0; i < rings.length; i++) {
        var r = rings[i];
        for (k = 0; k < r.length; k++) {
          var x = r[k], y = r[(k + 1) % r.length];
          nring[x]++; if (minring[x] === 0 || r.length < minring[x]) minring[x] = r.length;
          ringPair[x < y ? x + '_' + y : y + '_' + x] = 1;
        }
      }
      for (i = 0; i < (ext.cipCodes || []).length; i++) {
        var c = ext.cipCodes[i];
        if (c[0] < n) cip[c[0]] = (String(c[1]).toUpperCase() === 'R') ? 1 : (String(c[1]).toUpperCase() === 'S') ? 2 : 0;
      }
      var aroB = new Uint8Array(nb);
      for (i = 0; i < (ext.aromaticBonds || []).length; i++) aroB[ext.aromaticBonds[i]] = 1;
      var a = new Int32Array(nb), b = new Int32Array(nb), bond = new Uint8Array(nb * 3);
      for (i = 0; i < nb; i++) {
        var B = bonds[i], p = B.atoms[0], q = B.atoms[1];
        a[i] = p; b[i] = q; deg[p]++; deg[q]++;
        var bo = (B.bo === undefined) ? dB.bo : B.bo, st = (B.stereo === undefined) ? dB.stereo : B.stereo;
        bond[i * 3] = aroB[i] ? 3 : (bo === 1 ? 0 : bo === 2 ? 1 : bo === 3 ? 2 : 4);
        bond[i * 3 + 1] = ringPair[p < q ? p + '_' + q : q + '_' + p] ? 1 : 0;
        bond[i * 3 + 2] = st === 'cis' ? 1 : st === 'trans' ? 2 : 0;
      }
      var atom = new Uint8Array(n * 8);
      for (i = 0; i < n; i++) {
        var A = atoms[i];
        var z = (A.z === undefined) ? dA.z : A.z, chg = (A.chg === undefined) ? dA.chg : A.chg, h = (A.impHs === undefined) ? dA.impHs : A.impHs;
        var e = ELEM.indexOf(z);
        atom[i * 8] = e < 0 ? 16 : e;
        atom[i * 8 + 1] = Math.min(deg[i], 6);
        atom[i * 8 + 2] = Math.max(-2, Math.min(2, chg)) + 2;
        atom[i * 8 + 3] = Math.max(0, Math.min(h, 4));
        atom[i * 8 + 4] = aro[i];
        atom[i * 8 + 5] = minring[i] === 0 ? 0 : (minring[i] >= 8 ? 6 : Math.max(1, minring[i] - 2));
        atom[i * 8 + 6] = Math.min(nring[i], 3);
        atom[i * 8 + 7] = cip[i];
      }
      return { n: n, atom: atom, nb: nb, bond: bond, a: a, b: b };
    } catch (err) {
      return null;
    } finally {
      if (can) { try { can.delete(); } catch (e2) { /* 何もしない */ } }
    }
  }

  function fromSmiles(RDKit, smiles) {
    var m = null;
    try {
      m = RDKit.get_mol(smiles);
      if (!m || (m.is_valid && !m.is_valid())) return null;
      return fromMol(RDKit, m);
    } catch (err) { return null; } finally { if (m) { try { m.delete(); } catch (e) { /* 何もしない */ } } }
  }

  var api = { VERSION: VERSION, ATOM_DIMS: ATOM_DIMS, BOND_DIMS: BOND_DIMS, ELEM: ELEM, fromMol: fromMol, fromSmiles: fromSmiles };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.ChromGraph = api;
})(typeof self !== 'undefined' ? self : this);
