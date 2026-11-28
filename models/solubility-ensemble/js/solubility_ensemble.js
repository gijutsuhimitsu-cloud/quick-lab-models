// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE（MIT）。
// 溶けるか・混ざるか（組み合わせ版）：溶解度の見込みを出す入口。ブラウザ（<script>）でも Node（require）でも、ビルドなしで読める。
// 計算はサイトのページ https://gijutsu-himitsu.com/lab/solubility/detail/ （2026-11-27 の版）と同じ。
// 数はどれも 25 °C の「見込み」で、幅（誤差の約 7 割が入る、log10 の片側の幅）と一緒に返す。
//
//   有機溶媒：物理式の残差を学ぶ木 と 溶質・溶媒の積カーネル回帰 の平均に、蒸留した全結合のネットワークを混ぜる
//   水      ：水だけの木のモデル（BigSolDB の水と AqSolDB）＋ 表にない分子は、似た分子の誤差で直す
//
//   ブラウザ: js/ の features → abraham → solub → solfeat → knn → orgmodel → solgnn → solubility_ensemble の順に <script> で読み、
//             const S = await SolubilityEnsemble.load("weights/");
//   Node:     const SolubilityEnsemble = require("./js/solubility_ensemble.js");
//             const S = SolubilityEnsemble.loadNode();                    // 既定は ../weights
//
//   const c = S.fromExample("カフェイン")                    例の固体（RDKit.js なしで使える）
//   const c = S.prepare(RDKit, "CC(=O)Nc1ccc(O)cc1")           自分の分子（RDKit.js が要る）
//   S.estimate(c, "meoh", { target: 10 })                        1 つの溶媒。target は目標の濃度 mg/mL（任意）
//   S.estimate(c, "acn", { references: [{ solvent: "meoh", mgPerML: 25 }] })   手元の実測（3 溶媒まで）から換算
//   S.estimateAll(c, { target: 10 })                             14 溶媒
//   S.mix(c, "dmso", "water", 0.3)                               混合溶媒（2 つ目の溶媒の体積分率 0.3）
//   S.dilution(c, { stockSolvent: "dmso", stockMgPerML: 10, stockUL: 10, diluentUL: 990, diluent: { water: 0.5, acn: 0.5 } })
//   S.estimateTm(c)                                              融点の推算（目安。溶解度の計算には使わない）
(function (root) {
  'use strict';
  var isNode = typeof module !== 'undefined' && module.exports;
  var SC = isNode ? require('./solub.js') : root.SolCore;
  var AB = isNode ? require('./abraham.js') : root.SolAbraham;
  var SF = isNode ? require('./solfeat.js') : root.SolFeat;
  var FE = isNode ? require('./features.js') : root.ChromFeatures;
  var KNN = isNode ? require('./knn.js') : root.SolKNN;
  var ORGM = isNode ? require('./orgmodel.js') : root.SolOrg;
  var GNN = isNode ? require('./solgnn.js') : root.SolGNN;

  var VERSION = '0.2';
  var FILES = { base: 'base.json', gbt: 'gbt.json', calib: 'calib.json', orgTm: 'org_tm.json', orgNotm: 'org_notm.json', netTm: 'net_tm.json', netNotm: 'net_notm.json',
    knn: 'knn_water.json', tm: 'tm_model.json', lookup: 'lookup_25C.json', tmLookup: 'tm_lookup.json' };
  var REQUIRED = ['base', 'gbt', 'calib', 'orgTm', 'orgNotm', 'netTm', 'netNotm', 'knn'];
  var MEAS_BAND = 0.1;
  // 25 °C で混ざらない組（"x"）と、両端だけ混ざる組（"p"、[辞書順で前の溶媒に後ろが混ざる上限, 逆]、体積分率でおよそ）。値は溶媒の混和の一般的な表から
  var MISC = { 'toluene|water': 'x', 'dcm|water': 'x', 'chcl3|water': 'x', 'etoac|water': ['p', 0.03, 0.08] };
  var ERR = { empty: 'SMILES がない', multi: '塩と混合物（「.」を含む構造）は扱わない', parse: 'SMILES として読めない', charged: '電荷を持つ構造（イオン）は扱わない。中性の形で入れる',
    small: '炭素を含まない分子と、重い原子が 1 つだけの分子は扱わない', big: '重い原子が 80 個を超える分子は扱わない', elem: '扱えない元素を含む' };
  var OOTXT = { mw: '分子量が 500 以上なので、有機溶媒の見込みは出さない', amph: '酸と塩基の両方を持つ両性の分子なので、有機溶媒の見込みは出さない' };

  function fail(code) { var e = new Error(ERR[code] || code); e.code = code; return e; }
  var isNum = function (v) { return typeof v === 'number' && isFinite(v); };
  var copy = function (o) { var r = {}; for (var k in o) r[k] = o[k]; return r; };
  var toMg = function (logS, M) { return logS + Math.log10(M); };                                   // log10(mol/L) → log10(mg/mL)
  // 学習の表（BigSolDB）の有機溶媒の値は溶媒 1 L あたり。返す値は溶液 1 L あたり（溶質の密度を 1.2 g cm⁻³ と置く）
  var toSoln = function (logS, M) { return logS - Math.log10(1 + Math.pow(10, logS) * M / 1200); };
  var toSolv = function (logS, M) { var c = Math.pow(10, logS), v = Math.min(0.98, c * M / 1200); return Math.log10(c / (1 - v)); };
  // 文献の実測と目標を比べる時の幅（log）。AqSolDB の水の値は BigSolDB とのずれが大きいので広く取る
  var litTol = function (m) { return String(m[2]).indexOf('10.7910') === 0 ? 0.3 : 0.1; };
  function inverse(A) {                                          // 小さな行列の逆（掃き出し）
    var n = A.length, M = A.map(function (r, i) { return r.concat(A.map(function (_, j) { return +(i === j); })); });
    for (var c = 0; c < n; c++) {
      var p = c; for (var r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      var t = M[p]; M[p] = M[c]; M[c] = t; var d = M[c][c];
      for (var j = 0; j < 2 * n; j++) M[c][j] /= d;
      for (r = 0; r < n; r++) if (r !== c) { var f = M[r][c]; for (j = 0; j < 2 * n; j++) M[r][j] -= f * M[c][j]; }
    }
    return M.map(function (row) { return row.slice(n); });
  }

  function create(w) {
    REQUIRED.forEach(function (k) { if (!w || !w[k]) throw new Error('重みが足りない: ' + FILES[k]); });
    var D = w.base, G = w.gbt, CAL = w.calib, K = D.fit.K, TMDEF = D.fit.tmDefaultC, T = D.T_C, MIX = D.mix;
    var SV = D.solvents, byKey = {}, byId = {};
    SV.forEach(function (s) { byId[s.id] = s; [s.id, s.name, s.en, s.short].forEach(function (k) { if (k) byKey[String(k).toLowerCase()] = s; }); });
    var CI = {}; CAL.svs.forEach(function (id, i) { CI[id] = i; });
    var MWMAX = CAL.mw_max || 500;
    var ORG = { tm: ORGM.openModel(w.orgTm), notm: ORGM.openModel(w.orgNotm) };
    var NET = { tm: GNN.open(w.netTm), notm: GNN.open(w.netNotm) };
    var KM = KNN.load(w.knn);
    var lookup = w.lookup || null, tmLookup = w.tmLookup || null;
    var GN = FE.GROUPS.map(function (g) { return g[0]; });
    var simCache = new Map(), zCache = new Map();

    function solvent(x) {
      if (x && typeof x === 'object' && byId[x.id]) return byId[x.id];
      var s = byId[x] || byKey[String(x).toLowerCase()];
      if (!s) throw new Error('溶媒が表にない: ' + x + '（' + SV.map(function (v) { return v.id; }).join(', ') + '）');
      return s;
    }

    /* ---------- 対象外（有機溶媒）：分子量 500 以上と、酸と塩基の両方を持つ分子 ---------- */
    function gv(x, n) { var i = GN.indexOf(n); return x && i >= 0 ? (x[x.length - GN.length + i] || 0) : 0; }
    function amphoteric(x) {
      if (!x || gv(x, 'carboxylic_acid') + gv(x, 'sulfonic_acid') === 0) return false;
      var amd = gv(x, 'nitro') === 0 && gv(x, 'amide') === 0 ? gv(x, 'amidine_guanidine') : 0;
      return gv(x, 'amine_primary_aliph') + gv(x, 'amine_secondary_aliph') + gv(x, 'amine_tertiary_aliph') + amd > 0;
    }
    function outOfScope(sol) { return sol.M >= MWMAX ? 'mw' : amphoteric(sol.x) ? 'amph' : null; }

    /* ---------- モデル ---------- */
    function orgRow(sol, sv) { var v = SF.build(sol, sv, K, TMDEF, T), row = {}; G.names.forEach(function (n, i) { row[n] = +v[i].toFixed(5); }); return row; }
    function transfer(sv, sol) { return sv.water ? 0 : AB.logP(sv, sol.abr); }
    function pickTree(has, water) { return water && G.tm_w ? (has ? G.tm_w : G.notm_w) : (has ? G.tm : G.notm); }
    function netDesc(sol) { return (sol.x && sol.abr) ? sol.x.concat([sol.abr.E, sol.abr.S, sol.abr.A, sol.abr.B, sol.abr.V]) : null; }
    // 学習済みモデルの値（log10 mol/L。有機溶媒は溶媒 1 L あたり）と、式の内訳。融点があれば融点ありのモデル
    function model(sol, sv) {
      var has = isNum(sol.TmC), mode = has ? 'tm' : 'notm', O = ORG[mode], y, org = false, net = false;
      if (!sv.water && sol.words && (has || sol.gr)) {
        var row = orgRow(sol, sv), key = sol.smiles + '|' + mode + '|' + sol.TmC, sim = simCache.get(key);
        if (!sim) { if (simCache.size > 2000) simCache.clear(); sim = ORGM.soluteSim(O, row, sol.words); simCache.set(key, sim); }
        y = ORGM.predict(O, row, sol.words, sv.en, sol.gr, sim); org = true;
        var N = NET[mode], dsc = netDesc(sol);
        if (N && dsc) {
          var zk = sol.smiles + '|' + mode, z = zCache.get(zk);
          if (!z) { if (zCache.size > 2000) zCache.clear(); z = N.encode(null, dsc); zCache.set(zk, z); }
          var wm = N.json.mix_w; y = (1 - wm) * y + wm * N.predict(z, dsc, sv.id, has ? sol.TmC + 273.15 : null, 298.15); net = true;
        }
      } else {
        y = SF.evalModel(pickTree(has, sv.water), SF.build(sol, sv, K, TMDEF, T));
        if (sv.water && sol.knn) y += sol.knn[has ? 'tm' : 'notm'];     // 表にない分子：似た分子の誤差で直す（水だけ）
      }
      var p = SC.predict(sol, sv.water ? { water: true } : { transfer: transfer(sv, sol) }, T, K, TMDEF), parts = copy(p.parts);
      parts.corr = y - p.raw;
      return { y: y, has: has, org: org, net: net, parts: parts };
    }
    // 手元の実測（3 溶媒まで）。同じ溶媒を 2 回入れた時は対数で平均する
    function reference(sol, refs) {
      if (!refs || !refs.length) return null;
      var by = {};
      refs.forEach(function (r) {
        var v = +r.mgPerML; if (!(isNum(v) && v > 0)) throw new Error('references の mgPerML は正の数: ' + r.mgPerML);
        var id = solvent(r.solvent).id; (by[id] = by[id] || []).push(Math.log10(v / sol.M));
      });
      var ids = Object.keys(by); if (ids.length > 3) throw new Error('references の溶媒は 3 つまで');
      var P = CAL[isNum(sol.TmC) ? 'tm' : 'notm'];
      var list = ids.map(function (id) {
        var sv = byId[id], logS = by[id].reduce(function (a, b) { return a + b; }, 0) / by[id].length;
        return { sv: sv, logS: logS, r: (sv.water ? logS : toSolv(logS, sol.M)) - model(sol, sv).y };
      });
      var R = list.map(function (x) { return CI[x.sv.id]; }), Ai = inverse(R.map(function (i) { return R.map(function (j) { return P.C[i][j]; }); }));
      // C_RR⁻¹ (r_R − mu_R)
      var wv = Ai.map(function (row) { return row.reduce(function (s, a, j) { return s + a * (list[j].r - P.mu[R[j]]); }, 0); });
      var own = {}; list.forEach(function (x) { own[x.sv.id] = x; });
      return { R: R, Ai: Ai, w: wv, P: P, by: own };
    }
    function core(sol, sv, ref) {                                  // logS は log10(mol/L)（溶液 1 L あたり）
      var own = ref && ref.by[sv.id]; if (own) return { logS: own.logS, band: 0, ref: true, same: true };
      var oo = outOfScope(sol);
      if (!sv.water && oo) return { logS: NaN, band: 0, oo: oo };
      var m = model(sol, sv), cv = G.cv[m.has ? 'tm' : 'notm'], cap = SC.logSmax(sol.M), BS = CAL[m.has ? 'tm' : 'notm'].band;
      // 換算なしの幅：有機溶媒は溶媒ごと（学習に使っていない溶質で測った値）
      var base = sv.water ? (sol.knn && cv.water_knn ? cv.water_knn.p68 : cv.water.p68) : ((m.net || m.org) && BS && BS[sv.id]) ? BS[sv.id] : (m.org ? ORG[m.has ? 'tm' : 'notm'].cv.organic.p68 : cv.organic.p68);
      // どの表にも融点がない分子は、データの表に融点がない溶質で測った（広い）幅を使う
      var XU = D.extra && D.extra.notm_unknown_org;
      if (!sv.water && !m.has && sol.noTmData && XU) base = Math.max(base, XU.p68);
      if (ref) {
        // 入れた実測との差を、溶媒どうしの誤差の共分散で重みづけして移す：r_t = C_tR C_RR⁻¹ (r_R − mu_R) + mu_t
        var t = CI[sv.id], P = ref.P, c = ref.R.map(function (j) { return P.C[t][j]; });
        var shift = c.reduce(function (s, a, j) { return s + a * ref.w[j]; }, 0) + P.mu[t];
        var v = P.C[t][t] - c.reduce(function (s, a, i) { return s + a * ref.Ai[i].reduce(function (q, b, j) { return q + b * c[j]; }, 0); }, 0);
        var raw = m.y + shift, yR = sv.water ? raw : toSoln(raw, sol.M);
        return { logS: Math.min(yR, cap), band: base * Math.sqrt(Math.max(1e-4, v) / P.C[t][t]), ref: true, capped: yR > cap, shift: shift, has: m.has };
      }
      var yS = sv.water ? m.y : toSoln(m.y, sol.M);
      return { logS: Math.min(yS, cap), band: base, parts: m.parts, capped: yS > cap, has: m.has };
    }
    function measOf(sol, id) { var m = sol.meas && sol.meas[id]; return m ? { logS: m[0], mgPerML: Math.pow(10, toMg(m[0], sol.M)), n: m[1], doi: m[2] } : null; }
    function range(lmg, band) { return [Math.pow(10, lmg - band), Math.pow(10, lmg + band)]; }
    function capMg(sol) { return toMg(SC.logSmax(sol.M), sol.M); }

    // 1 つの溶媒での見込み。判定は、手元の実測 → 文献の実測 → 見込み の順
    function estimate(sol, svx, opt) {
      opt = opt || {};
      var sv = solvent(svx), ref = reference(sol, opt.references), e = core(sol, sv, ref), ms = measOf(sol, sv.id), lmg = toMg(e.logS, sol.M);
      var oo = !!e.oo;
      var r = { solvent: sv.id, logS: oo ? null : e.logS, mgPerML: oo ? null : Math.pow(10, lmg), band: oo ? null : e.band, factor: oo ? null : Math.pow(10, e.band),
        range: oo ? null : range(lmg, e.band), tmUsed: isNum(sol.TmC), capped: !!e.capped, parts: e.parts || null, fromReference: !!e.ref, shift: isNum(e.shift) ? e.shift : null,
        outOfScope: e.oo || null, reason: e.oo ? OOTXT[e.oo] : null, measured: ms };
      var tg = opt.target;
      if (isNum(tg) && tg > 0) {
        var lt = Math.log10(tg);
        if (e.same) { r.verdict = SC.verdict(lmg, lt, MEAS_BAND); r.basis = 'reference'; }
        else if (ms) { r.verdict = SC.verdict(toMg(ms.logS, sol.M), lt, litTol(sol.meas[sv.id])); r.basis = 'measured'; }
        else if (oo) { r.verdict = null; r.basis = null; }
        else { r.verdict = SC.verdict(lmg, lt, e.band); r.basis = e.ref ? 'reference' : 'model'; }
      }
      return r;
    }
    function estimateAll(sol, opt) {
      var list = (opt && opt.solvents) ? opt.solvents.map(solvent) : SV;
      return list.map(function (sv) { return estimate(sol, sv, opt); });
    }

    /* ---------- 混合溶媒と希釈 ---------- */
    // 端の値：手元の実測 → 文献の実測 → モデルの見込み（換算を含む）の順に使う
    function best(sol, sv, ref, useMeas) {
      var own = ref && ref.by[sv.id]; if (own) return { logS: own.logS, band: MEAS_BAND, own: true };
      var m = useMeas !== false && sol.meas && sol.meas[sv.id];
      return m ? { logS: m[0], band: litTol(m), meas: true } : core(sol, sv, ref);
    }
    function misc(a, b) { return MISC[[a, b].sort().join('|')] || null; }
    function oneLiquid(a, b, f) {
      var m = misc(a, b); if (!m || a === b) return true;
      if (f <= 0 || f >= 1) return true;
      if (m === 'x') return false;
      var ka = [a, b].sort()[0], fKb = a === ka ? f : 1 - f;      // 辞書順で後ろの溶媒の割合
      return fKb <= m[1] || 1 - fKb <= m[2];
    }
    // 混合溶媒（2 成分）。fB は溶媒 B の体積分率。対数線形の内挿に、溶媒の組ごとの補正を足す
    function mix(sol, ax, bx, fB, opt) {
      opt = opt || {};
      if (!(fB >= 0 && fB <= 1)) throw new Error('fB は 0〜1');
      var a = solvent(ax), b = solvent(bx), ref = reference(sol, opt.references), ea = best(sol, a, ref, opt.useMeasured), eb = best(sol, b, ref, opt.useMeasured);
      var r = { solvents: [a.id, b.id], fB: fB, mgPerML: null, band: null, factor: null, range: null, onePhase: oneLiquid(a.id, b.id, fB), outOfScope: null, reason: null, ends: { a: ea, b: eb } };
      if (!isFinite(ea.logS) || !isFinite(eb.logS)) { r.outOfScope = outOfScope(sol); r.reason = OOTXT[r.outOfScope]; return r; }
      if (!r.onePhase) { r.reason = 'この割合では 2 層に分かれるので、見込みを出さない'; return r; }
      var ya = toMg(ea.logS, sol.M), yb = toMg(eb.logS, sol.M), y = ya, band = ea.band;
      if (a.id !== b.id) {
        var comp = {}; comp[a.id] = 1 - fB; comp[b.id] = fB;
        var c = SC.mixCorr(comp, MIX), wo = a.water || b.water, mb = c.est ? (wo ? MIX.band.est_wo : MIX.band.est_oo) : (wo ? MIX.band.corrected_wo : MIX.band.corrected_oo), eb2 = (1 - fB) * ea.band + fB * eb.band;
        y = Math.min(ya + (yb - ya) * fB + c.corr, capMg(sol)); band = Math.sqrt(eb2 * eb2 + mb * mb * c.w * c.w);
        r.correction = c.corr;
      }
      r.mgPerML = Math.pow(10, y); r.band = band; r.factor = Math.pow(10, band); r.range = range(y, band);
      if (isNum(opt.target) && opt.target > 0) r.verdict = SC.verdict(y, Math.log10(opt.target), band);
      return r;
    }
    function mixCurve(sol, ax, bx, opt) {
      var n = (opt && opt.points) || 51, out = [];
      for (var i = 0; i < n; i++) out.push(mix(sol, ax, bx, i / (n - 1), opt));
      return out;
    }
    // 原液を薄めた時に析出するか。3 成分は 2 成分の組の和で近似し、補正の幅を 2 倍に取る（3 成分の実測との照合はしていない）
    function dilution(sol, o) {
      var sv = solvent(o.stockSolvent), ref = reference(sol, o.references), dil = [], D0 = o.diluent || { water: 1 };
      if (Array.isArray(D0)) D0.forEach(function (d) { dil.push({ key: solvent(d.solvent).id, f: +d.f }); });
      else for (var k in D0) dil.push({ key: solvent(k).id, f: +D0[k] });
      var conc = +o.stockMgPerML, take = +o.stockUL, vol = o.diluentUL === undefined || o.diluentUL === null ? 0 : +o.diluentUL;
      if (!(isNum(conc) && conc > 0)) throw new Error('stockMgPerML（原液の濃度）は正の数: ' + o.stockMgPerML);
      if (!(isNum(take) && take > 0)) throw new Error('stockUL（原液の量）は正の数: ' + o.stockUL);
      if (!(isNum(vol) && vol >= 0)) throw new Error('diluentUL（薄める液の量）は 0 以上の数: ' + o.diluentUL);
      var fsum = 0;
      dil.forEach(function (x) { if (!(isNum(x.f) && x.f >= 0)) throw new Error('diluent の成分の比は 0 以上の数: ' + x.key + ' = ' + x.f); fsum += x.f; });
      if (vol > 0 && !(fsum > 0)) throw new Error('diluent の成分の比の合計が 0（薄める液の組成が決まらない）');
      var d = SC.dilute(sv.id, conc, take, dil, vol);
      if (!d) throw new Error('量が正しくない');
      var items = Object.keys(d.comp).filter(function (k) { return d.comp[k] > 0; }).map(function (k) { var e = best(sol, byId[k], ref, o.useMeasured); return { id: k, f: d.comp[k], logS: toMg(e.logS, sol.M), band: e.band }; });
      var r = { factor: d.factor, concMgPerML: d.conc, comp: d.comp, solubilityMgPerML: null, band: null, range: null, verdict: null, threeComponent: items.length > 2, outOfScope: null, reason: null, warnings: [] };
      if (items.some(function (i) { return !isFinite(i.logS); })) { r.outOfScope = outOfScope(sol); r.reason = OOTXT[r.outOfScope]; return r; }
      var wo = d.comp.water > 0, c = SC.mixCorr(d.comp, MIX);
      var mb = (c.est ? (wo ? MIX.band.est_wo : MIX.band.est_oo) : (wo ? MIX.band.corrected_wo : MIX.band.corrected_oo)) * (items.length > 2 ? 2 : 1), eb2 = items.reduce(function (t, i) { return t + i.f * i.band; }, 0);
      var y = Math.min(SC.mix(items) + c.corr, capMg(sol)), band = Math.sqrt(eb2 * eb2 + mb * mb * c.w * c.w), v = SC.verdict(y, Math.log10(d.conc), band);
      r.solubilityMgPerML = Math.pow(10, y); r.band = band; r.range = range(y, band);
      r.verdict = v === 'yes' ? 'stays_dissolved' : v === 'no' ? 'precipitates' : 'edge';
      // 原液の段階で溶け残る見込み／原液の溶媒が水と混ざりにくい
      var es = best(sol, sv, ref, o.useMeasured), ys = toMg(es.logS, sol.M);
      if (isFinite(ys) && SC.verdict(ys, Math.log10(conc), es.band) === 'no') r.warnings.push('stock_not_dissolved');
      var mw = misc(sv.id, 'water'), fS = d.comp[sv.id] || 0;
      if (wo && mw && (mw === 'x' ? fS >= 0.01 : fS > mw[2])) r.warnings.push('may_separate_with_water');
      return r;
    }

    /* ---------- 構造から溶質を作る ---------- */
    function flatSmiles(R, can) {                      // 立体の指定を除いた形
      var t = can.split('@').join('').split('/').join('').split(String.fromCharCode(92)).join('');
      if (t === can) return can;
      var m = null; try { m = R.get_mol(t); } catch (e) { m = null; }
      if (!m || !m.is_valid()) { if (m) m.delete(); return null; }
      var c = m.get_smiles(); m.delete(); return c;
    }
    function lookUp(tab, field, can, flat) {
      if (!tab) return null;
      var T0 = tab[field] || {};
      if (T0[can] !== undefined) return { v: T0[can], flat: false };
      var k = flat && tab.flat ? tab.flat[flat] : null;
      if (k && T0[k] !== undefined) return { v: T0[k], flat: true };
      if (flat && flat !== can && T0[flat] !== undefined) return { v: T0[flat], flat: true };
      return null;
    }
    // RDKit は RDKit.js の initRDKitModule() の結果。融点は、入れた値（opt.TmC）＞ 溶解度の表 ＞ 融点の表。
    // どれもなければ null で、tmHint に推算値（5 °C 刻み。溶解度の計算には使わない）
    function prepare(R, smiles, opt) {
      opt = opt || {};
      var text = String(smiles || '').trim();
      if (!text) throw fail('empty');
      if (text.indexOf('.') >= 0) throw fail('multi');
      var m = null; try { m = R.get_mol(text); } catch (e) { m = null; }
      if (!m || !m.is_valid()) { if (m) m.delete(); throw fail('parse'); }
      var can, f, d, js, fp = null, gr = null;
      try {
        m.remove_hs_in_place(); can = m.get_smiles(); js = JSON.parse(m.get_json()).molecules[0]; f = FE.compute(R, m); d = JSON.parse(m.get_descriptors());
        try { fp = m.get_morgan_fp(JSON.stringify({ radius: 2, nBits: 1024 })); } catch (e) { fp = null; }
        try { gr = ORGM.packing(js); } catch (e) { gr = null; }
      } finally { m.delete(); }
      var atoms = js.atoms || [];
      if (atoms.reduce(function (t, a) { return t + (a.chg || 0); }, 0) !== 0) throw fail('charged');
      if (d.NumHeavyAtoms < 2 || !atoms.some(function (a) { return (a.z === undefined ? 6 : a.z) === 6; })) throw fail('small');
      if (d.NumHeavyAtoms > 80) throw fail('big');
      var V = AB.mcgowanV(R, can); if (V === null) throw fail('elem');
      var abr = AB.predictSolute(D.model, f.desc.concat(f.groups)); abr.V = V;
      var x = f.desc.concat(f.groups), flat = flatSmiles(R, can);
      var hl = lookUp(lookup, 'solutes', can, flat), hit = hl ? hl.v : null;
      var TmC = hit && isNum(hit.t) ? hit.t : null, tmSrc = TmC !== null ? (hit.s || 'sol') : null;
      if (TmC === null) { var tl = lookUp(tmLookup, 'mp', can, flat); if (tl && isNum(tl.v)) { TmC = tl.v; tmSrc = 'mp'; } }
      var noTmData = TmC === null;
      if (isNum(opt.TmC)) { TmC = opt.TmC; tmSrc = 'user'; }
      var tmHint = TmC === null && w.tm ? Math.round(SF.evalModel(w.tm.model, x) / 5) * 5 : null;
      var knn = null;
      if (!hit && fp) { var on = KNN.onBits(fp), ka = KNN.corr(KM, on, 'tm'), kb = KNN.corr(KM, on, 'notm'); knn = { tm: ka.c, notm: kb.c, used: kb.used }; }
      var sol = { name: opt.name || can, smiles: can, M: +d.amw.toFixed(2), logPow: d.CrippenClogP, TmC: TmC, tmSrc: tmSrc, tmHint: tmHint, noTmData: noTmData,
        abr: { E: abr.E, S: abr.S, A: abr.A, B: abr.B, V: abr.V }, x: x, words: fp ? ORGM.wordsFromRDKit(fp) : null, gr: gr, knn: knn, meas: hit ? hit.m || {} : {}, inLookup: !!hit };
      sol.outOfScope = outOfScope(sol);
      return sol;
    }
    var examples = D.presets;
    function fromExample(key) {
      var k = String(key).toLowerCase();
      for (var i = 0; i < examples.length; i++) {
        var e = examples[i];
        if (e.id === key || e.name.toLowerCase() === k || e.smiles === key) {
          var c = copy(e); c.words = ORGM.wordsFromOn(e.on); delete c.on; c.inLookup = true; c.knn = null; c.outOfScope = outOfScope(c); return c;
        }
      }
      throw new Error('例の固体にない: ' + key + '（weights/base.json の presets）');
    }
    function withTm(sol, TmC) { var c = copy(sol); c.TmC = isNum(TmC) ? TmC : null; c.tmSrc = c.TmC === null ? null : 'user'; return c; }
    function estimateTm(sol) {
      if (!w.tm) throw new Error('融点のモデル（tm_model.json）を読んでいない');
      return { TmC: SF.evalModel(w.tm.model, sol.x), band: w.tm.cv ? w.tm.cv.p68 : null, note: '目安。溶解度の計算には使わない（推算した融点では幅が狭くならなかった）' };
    }

    return { version: VERSION, T_C: T, solvents: SV, examples: examples, solvent: solvent, prepare: prepare, fromExample: fromExample, withTm: withTm,
      estimate: estimate, estimateAll: estimateAll, mix: mix, mixCurve: mixCurve, dilution: dilution, estimateTm: estimateTm };
  }

  // 重みを読む（fetch）。opt.lookup / opt.tmLookup / opt.tm を false にすると、その表を読まない（その時は融点や文献の実測を引かない）
  function load(base, opt) {
    opt = opt || {}; base = base || 'weights/'; if (base.slice(-1) !== '/') base += '/';
    var keys = REQUIRED.concat(['tm', 'lookup', 'tmLookup'].filter(function (k) { return opt[k] !== false; }));
    // 通信が切れた時は、少し待って 2 回までやり直す
    function get(file, n) {
      return fetch(base + file).then(function (r) { if (!r.ok) throw new Error(file + ': ' + r.status); return r.json(); })
        .catch(function (e) { if (n <= 0) throw e; return new Promise(function (ok) { setTimeout(ok, 600); }).then(function () { return get(file, n - 1); }); });
    }
    return Promise.all(keys.map(function (k) { return get(FILES[k], 2); })).then(function (vals) { var w = {}; keys.forEach(function (k, i) { w[k] = vals[i]; }); return create(w); });
  }
  // Node だけ：フォルダから同期で読む（既定はこのファイルの隣の ../weights）
  function loadNode(dir, opt) {
    opt = opt || {};
    var fs = require('fs'), path = require('path'), d = dir || path.join(__dirname, '..', 'weights'), w = {};
    Object.keys(FILES).forEach(function (k) {
      if (REQUIRED.indexOf(k) < 0 && opt[k] === false) return;
      var f = path.join(d, FILES[k]);
      if (fs.existsSync(f)) w[k] = JSON.parse(fs.readFileSync(f, 'utf8'));
    });
    return create(w);
  }

  var api = { VERSION: VERSION, FILES: FILES, create: create, load: load, loadNode: loadNode };
  if (isNode) module.exports = api; else root.SolubilityEnsemble = api;
})(typeof self !== 'undefined' ? self : this);
