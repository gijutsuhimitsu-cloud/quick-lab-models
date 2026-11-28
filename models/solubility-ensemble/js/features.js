// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE（MIT）。
// 化合物の特徴量。学習（Node）とページ（ブラウザ）で同じファイルを使う。
// 入力は RDKit.js（@rdkit/rdkit 2025.3.4-1.0.0、BSD-3-Clause）の分子オブジェクト。
// 並びを変えたら、モデルを学習し直すこと（FEATURE_VERSION を上げる）。
(function (root) {
  'use strict';
  var FEATURE_VERSION = 2;

  // RDKit.js の get_descriptors() から使う 41 個（43 個のうち exactmw と NumUnspecifiedAtomStereoCenters を外した物。NumAtoms は最後に入っている）
  var DESC = ['amw', 'lipinskiHBA', 'lipinskiHBD', 'NumRotatableBonds', 'NumHBD', 'NumHBA', 'NumHeavyAtoms',
    'NumHeteroatoms', 'NumAmideBonds', 'FractionCSP3', 'NumRings', 'NumAromaticRings', 'NumAliphaticRings',
    'NumSaturatedRings', 'NumHeterocycles', 'NumAromaticHeterocycles', 'NumSaturatedHeterocycles',
    'NumAliphaticHeterocycles', 'NumSpiroAtoms', 'NumBridgeheadAtoms', 'NumAtomStereoCenters',
    'labuteASA', 'tpsa', 'CrippenClogP', 'CrippenMR', 'chi0v', 'chi1v', 'chi2v', 'chi3v', 'chi4v',
    'chi0n', 'chi1n', 'chi2n', 'chi3n', 'chi4n', 'hallKierAlpha', 'kappa1', 'kappa2', 'kappa3', 'Phi', 'NumAtoms'];

  // 官能基の数（SMARTS）。名前, SMARTS
  var GROUPS = [
    ['carboxylic_acid', '[CX3](=O)[OX2H1,OX1-]'],
    ['phenol', '[c][OX2H1,OX1-]'],
    ['alcohol', '[CX4][OX2H1]'],
    ['amine_primary_aliph', '[NX3;H2;!$(N-[a]);!$(N-[C,S,P]=[O,S,N])]'],
    ['amine_secondary_aliph', '[NX3;H1;!$(N-[a]);!$(N-[C,S,P]=[O,S,N]);!$(N-[N,O])]([#6])[#6]'],
    ['amine_tertiary_aliph', '[NX3;H0;!$(N-[a]);!$(N-[C,S,P]=[O,S,N]);!$(N-[N,O])]([#6])([#6])[#6]'],
    ['aniline_n', '[NX3;!$(N-[C,S,P]=[O,S,N])]-[c]'],
    ['quaternary_n', '[NX4+]'],
    ['pyridine_n', '[nX2;r6]'],
    ['azole_n', '[nX2;r5]'],
    ['aromatic_nh', '[nX3;H1]'],
    ['amidine_guanidine', '[NX3][CX3]=[NX2]'],
    ['amide', '[CX3](=O)[NX3]'],
    ['sulfonamide', '[SX4](=O)(=O)[NX3]'],
    ['sulfonic_acid', '[SX4](=O)(=O)[OX2H1,OX1-]'],
    ['sulfate_phosphate', '[S,P](=O)([OX2H1,OX1-])[OX2]'],
    ['phosphorus', '[P]'],
    ['nitro', '[$([NX3](=O)=O),$([NX3+](=O)[O-])]'],
    ['nitrile', '[CX2]#[NX1]'],
    ['ketone', '[#6][CX3](=O)[#6]'],
    ['aldehyde', '[CX3H1](=O)[#6]'],
    ['ester', '[#6][CX3](=O)[OX2][#6]'],
    ['ether', '[OD2]([#6;!$(C=O)])[#6;!$(C=O)]'],
    ['carbamate_urea', '[NX3][CX3](=O)[OX2,NX3]'],
    ['thiol', '[SX2H1]'],
    ['thioether', '[SX2]([#6])[#6]'],
    ['sulfone_sulfoxide', '[$([SX4](=O)(=O)([#6])[#6]),$([SX3](=O)([#6])[#6])]'],
    ['thiocarbonyl', '[CX3]=[SX1]'],
    ['fluorine', '[F]'],
    ['chlorine', '[Cl]'],
    ['bromine', '[Br]'],
    ['iodine', '[I]'],
    ['cf3', '[CX4](F)(F)F'],
    ['aromatic_c', '[c]'],
    ['aromatic_o_s', '[o,s]'],
    ['sp3_c', '[CX4]'],
    ['methyl', '[CX4H3]'],
    ['methylene', '[CX4H2]'],
    ['alkene', '[CX3]=[CX3]'],
    ['alkyne', '[CX2]#[CX2]'],
    ['carbonyl_any', '[CX3]=[OX1]'],
    ['hydroxyl_any', '[OX2H1]'],
    ['nh_any', '[#7;!H0]'],
    ['n_any', '[#7]'],
    ['o_any', '[#8]'],
    ['s_any', '[#16]'],
    ['positive', '[+;!$([N+](=O)[O-]);!$([n+][O-]);!$([N+][O-])]'],
    ['negative', '[-;!$([O-][N+]=O);!$([O-][n+]);!$([O-][N+])]'],
    ['tetrazole', 'c1nnn[nX3H1,nX2]1'],
    ['enol_acidic_ch', '[OX2H1][CX3]=[CX3]'],
    ['n_oxide', '[$([#7+][OX1-]),$([#7]=[OX1]);!$([NX3](=O)=O);!$([NX3+](=O)[O-])]'],
    ['fused_aromatic_atom', '[c;R2]'],
    ['long_chain', '[CX4H2][CX4H2][CX4H2][CX4H2]'],
    ['tert_butyl_like', '[CX4]([CH3])([CH3])[CH3]'],
    ['ortho_oh_carbonyl', '[OX2H1]c:c[CX3]=O'],
    ['halo_aromatic', 'c[F,Cl,Br,I]']
  ];

  // ---- イオン化する官能基（特徴量の第 2 版）----
  // 規則ごとの部位の数を数える。その pH での電荷への換算は retention_model.js（規則の種類と pKa は係数ファイルの meta.ion_rules）。
  // [名前, SMARTS, 'acid'|'base'|'perm+', 代表的な pKa（塩基は共役酸の pKa）]
  var NOAMIDE = '!$(N-[C,S,P]=[O,S,N]);!$(N-[N,O])';
  var ION = [
    ['sulfonic_acid', '[SX4](=O)(=O)[OX2H1,OX1-]', 'acid', -1.0],
    ['phosphoric_acid', '[PX4](=O)([OX2H1,OX1-])', 'acid', 1.8],
    ['phosphoric_acid_2', '[PX4](=O)([OX2H1,OX1-])[OX2H1,OX1-]', 'acid', 6.5],          // 2 段目の解離（P 原子ごとに 1 回）
    ['cooh_alpha_amino', '[CX3](=O)([OX2H1,OX1-])[CX4][NX3,NX4+;!$(N-C=O)]', 'acid', 2.3],
    ['cooh_conjugated', '[CX3](=O)([OX2H1,OX1-])[c,$([CX3]=[CX3])]', 'acid', 4.1],
    ['cooh_any', '[CX3](=O)[OX2H1,OX1-]', 'acid', 4.6],                  // 上の 2 つを引いた数を使う
    ['tetrazole', '[nX3H1,nX2-]1nnnc1', 'acid', 4.9],
    ['acyl_sulfonamide', '[SX4](=O)(=O)[NX3H1]C(=O)', 'acid', 4.5],
    ['enol_vinylogous_acid', '[OX2H1,OX1-][CX3]=[CX3][CX3,SX4]=O', 'acid', 4.6],
    ['aryl_sulfonamide_nh', '[SX4](=O)(=O)[NX3H1]c', 'acid', 7.2],
    ['nitrophenol', '[OX2H1,OX1-;$(Oc:c[$([NX3](=O)=O),$([NX3+](=O)[O-])]),$(Oc1ccc([$([NX3](=O)=O),$([NX3+](=O)[O-])])cc1)]', 'acid', 7.2],   // オルトかパラにニトロ基（O 原子ごとに 1 回）
    ['imide', '[CX3](=O)[NX3H1][CX3](=O)', 'acid', 8.8],
    ['thiol', '[SX2H1,SX1-]', 'acid', 9.5],
    ['phenol', '[OX2H1,OX1-]c', 'acid', 10.0],                             // ニトロフェノールを引いた数を使う
    ['sulfonamide_nh2', '[SX4](=O)(=O)[NX3H2]', 'acid', 10.1],
    ['guanidine', '[CX3](=[NX2,NX3+])([NX3])[NX3]', 'base', 13.0],                 // 中心の C ごとに 1 回。グアニジン 13.6・アルギニン 12.5 の間
    ['amidine', '[CX3;!$(C([#7])([#7])[#7]);!r5;!r6](=[NX2;!$(N-[O,N])])[NX3;!$(N-C=O)]', 'base', 11.0],
    ['amine_primary', '[NX3;H2;!$(N-[a]);' + NOAMIDE + '][CX4]', 'base', 9.8],
    ['amine_secondary', '[NX3;H1;!$(N-[a]);' + NOAMIDE + ']([CX4])[CX4]', 'base', 10.3],
    ['amine_tertiary', '[NX3;H0;!$(N-[a]);' + NOAMIDE + ']([CX4])([CX4])[CX4]', 'base', 9.3],
    ['amine_protonated', '[NX4+;!H0;!$([N+]-[O-])]', 'base', 10.0],                        // プロトン化した形で描かれたアミン（永久電荷ではない）
    ['amine_beta_hetero', '[NX3;!$(N-[a]);' + NOAMIDE + ']([CX4])[CX4][CX4][OX2,F,$([NX3]),$(C#N)]', 'mod', -1.5],  // 近くの O・N・F で pKa が下がる分
    ['aniline', '[NX3;' + NOAMIDE + ';!$(N(c)c)]c', 'base', 4.6],
    ['aminopyridine_n', '[nX2;r6]:c-[NX3;' + NOAMIDE + ']', 'base', 7.0],
    ['diazine_n', '[nX2;r6;$(n:c:n),$(n:n),$(n1ccncc1)]', 'base', 1.5],
    ['pyridine_n', '[nX2,$([nX3+;!H0]);r6]', 'base', 5.2],                            // 上の 2 つを引いた数を使う
    ['imidazole_n', '[nX2;r5]:c:[nX3;r5]', 'base', 6.8],
    ['azole_n_other', '[nX2;r5]', 'base', 2.0],                         // イミダゾールを引いた数を使う
    ['quaternary_n', '[NX4+;H0;!$([N+]-[O-])]', 'perm+', 0],                              // 水素を持たない物だけ（本当の第四級）
    ['pyridinium', '[nX3+;H0;!$([n+]-[O-])]', 'perm+', 0]
  ];
  var _qi = null;
  function ionQueries(RDKit) {
    if (_qi && _qi.R === RDKit) return _qi.list;
    var list = ION.map(function (g) {
      var q = RDKit.get_qmol(g[1]);
      if (!q || !q.is_valid()) throw new Error('SMARTS が読めない: ' + g[0]);
      return q;
    });
    _qi = { R: RDKit, list: list };
    return list;
  }
  function countSites(mol, q) {
    var s = mol.get_substruct_matches(q);
    if (!s || s === '{}' || s === '[]') return 0;
    var a = JSON.parse(s);
    if (!Array.isArray(a)) return a && a.atoms ? 1 : 0;
    var seen = {};                       // 同じ先頭の原子は 1 回だけ数える
    a.forEach(function (m) { seen[m.atoms[0]] = 1; });
    return Object.keys(seen).length;
  }
  function ionCounts(RDKit, mol) {
    var qs = ionQueries(RDKit), c = {};
    ION.forEach(function (g, i) { c[g[0]] = countSites(mol, qs[i]); });
    c.cooh_any = Math.max(0, c.cooh_any - c.cooh_alpha_amino - c.cooh_conjugated);
    c.phenol = Math.max(0, c.phenol - c.nitrophenol);
    c.aryl_sulfonamide_nh = Math.max(0, c.aryl_sulfonamide_nh - c.acyl_sulfonamide);
    c.pyridine_n = Math.max(0, c.pyridine_n - c.aminopyridine_n - c.diazine_n);
    c.azole_n_other = Math.max(0, c.azole_n_other - c.imidazole_n);
    c.amidine = Math.max(0, c.amidine - c.guanidine);            // グアニジンはアミジンの検索にも当たるので引く
    return ION.map(function (g) { return c[g[0]]; });
  }


  var FP_BITS = 512, FP_RADIUS = 2;
  var _q = null;

  function queries(RDKit) {
    if (_q && _q.R === RDKit) return _q.list;
    var list = GROUPS.map(function (g) {
      var q = RDKit.get_qmol(g[1]);
      if (!q || !q.is_valid()) throw new Error('SMARTS が読めない: ' + g[0]);
      return q;
    });
    _q = { R: RDKit, list: list };
    return list;
  }

  function countMatches(mol, q) {
    var s = mol.get_substruct_matches(q);
    if (!s || s === '{}' || s === '[]') return 0;
    var a = JSON.parse(s);
    return Array.isArray(a) ? a.length : (a && a.atoms ? 1 : 0);
  }

  // mol -> {desc:[41], groups:[56], ion:[規則ごとの数], fp:"0101..."(512)}
  // 呼ぶ前に mol.remove_hs_in_place() で明示の水素を外しておくこと（学習側と同じ）。
  function compute(RDKit, mol) {
    var d = JSON.parse(mol.get_descriptors());
    // 小数 4 桁に丸める（学習データの特徴量と同じ桁。丸めないと、木のしきい値の近くで、ストックの化合物と同じ分子でも予測がずれる）
    var desc = DESC.map(function (k) { var v = d[k]; return (typeof v === 'number' && isFinite(v)) ? +v.toFixed(4) : 0; });
    var qs = queries(RDKit);
    var groups = qs.map(function (q) { return countMatches(mol, q); });
    var fp = mol.get_morgan_fp(JSON.stringify({ radius: FP_RADIUS, nBits: FP_BITS }));
    return { desc: desc, groups: groups, ion: ionCounts(RDKit, mol), fp: fp };
  }

  function names() {
    return DESC.concat(GROUPS.map(function (g) { return 'g_' + g[0]; }));
  }

  var api = { FEATURE_VERSION: FEATURE_VERSION, DESC: DESC, GROUPS: GROUPS, ION: ION, ionCounts: ionCounts, FP_BITS: FP_BITS, FP_RADIUS: FP_RADIUS, compute: compute, names: names };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.ChromFeatures = api;
})(typeof self !== 'undefined' ? self : this);
