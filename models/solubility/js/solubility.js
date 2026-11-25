// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE（MIT）。
// 溶けるか・混ざるか：溶解度の見込みを出す入口。ブラウザ（<script>）でも Node（require）でも、ビルドなしで読める。
// 計算はサイトのページ（第 5 版）と同じ。数はどれも 25 °C の「見込み」で、幅（外れの約 7 割が入る、log10 の片側の幅）と一緒に返す。
//
//   ブラウザ: <script src="js/features.js"></script> <script src="js/abraham.js"></script> <script src="js/solub.js"></script>
//             <script src="js/solfeat.js"></script> <script src="js/solubility.js"></script>
//             const S = await SolubilityModel.load("weights/");          // 重みの JSON を fetch で読む
//   Node:     const SolubilityModel = require("./js/solubility.js");
//             const S = SolubilityModel.loadNode();                      // 既定は ../weights
//
//   const c = S.fromExample("カフェイン")                 例の固体（RDKit.js なしで使える）
//   const c = S.prepare(RDKit, "CC(=O)Nc1ccc(O)cc1")        自分の分子（RDKit.js が要る。照合の表に実測があれば重ねる）
//   S.estimate(c, "meoh", { target: 10 })                     1 つの溶媒。target は目標の濃度 mg/mL（任意）
//   S.estimate(c, "acn", { reference: { solvent: "meoh", mgPerML: 25 } })   手元の実測 1 点から換算
//   S.estimateAll(c, { target: 10 })                          主な 16 溶媒
//   S.mix(c, "dmso", "water", 0.3)                            混合溶媒（2 つ目の溶媒の体積分率 0.3）
//   S.dilution(c, { stockSolvent: "dmso", stockMgPerML: 10, stockUL: 10, diluentUL: 990, diluent: { water: 0.5, acn: 0.5 } })
//   S.estimateTm(c)                                           融点の推算（目安。計算には使わない）
(function (root) {
  'use strict';
  var isNode = typeof module !== 'undefined' && module.exports;
  var SC = isNode ? require('./solub.js') : root.SolCore;
  var AB = isNode ? require('./abraham.js') : root.SolAbraham;
  var SF = isNode ? require('./solfeat.js') : root.SolFeat;
  var FE = isNode ? require('./features.js') : root.ChromFeatures;

  var VERSION = '0.1';
  var FILES = { gbt: 'solubility_gbt.json', formula: 'formula.json', mix: 'mix_fit.json', solvents: 'solvents.json', abraham: 'abraham_solute_model.json',
    tm: 'tm_model.json', lookup: 'lookup_25C.json', tmLookup: 'tm_lookup.json', examples: 'example_compounds.json' };
  var REQUIRED = ['gbt', 'formula', 'mix', 'solvents', 'abraham'];
  var MEAS_BAND = 0.1;                 // 実測の値を使う時の幅（log10）。ページと同じ
  var ERR = { empty: 'SMILES がない', multi: '塩と混合物（「.」を含む構造）は扱わない', parse: 'SMILES として読めない', big: '重い原子が 80 個を超える分子は扱わない', elem: '扱えない元素を含む（McGowan の体積が出ない）' };

  function fail(code) { var e = new Error(ERR[code] || code); e.code = code; return e; }
  var isNum = function (v) { return typeof v === 'number' && isFinite(v); };
  var log10mg = function (logS, M) { return logS + Math.log10(M); };         // log10(mol/L) → log10(mg/mL)
  function copy(o) { var r = {}; for (var k in o) r[k] = o[k]; return r; }

  function create(w) {
    REQUIRED.forEach(function (k) { if (!w || !w[k]) throw new Error('重みが足りない: ' + FILES[k]); });
    var G = w.gbt, F = w.formula, K = F.K, TMDEF = F.tm_default_C, T = F.T_C || 25, MIX = w.mix;
    if (G.meta && G.meta.version && String(G.meta.version).split('+')[0] !== VERSION) console.warn('solubility: 重みの版 ' + G.meta.version + ' とプログラムの版 ' + VERSION + ' が違う');
    var solvents = w.solvents.solvents.map(function (s) { var o = copy(s); o.water = !!s.water; return o; });
    var byId = {}, byKey = {};
    solvents.forEach(function (s) {
      byId[s.id] = s;
      [s.id, s.name, s.ja, s.short].forEach(function (k) { if (k) { var key = String(k).toLowerCase(); if (!(key in byKey)) byKey[key] = s; } });
    });
    var main = solvents.filter(function (s) { return s.main; });
    var lookup = w.lookup ? w.lookup.solutes : null, tmLookup = w.tmLookup ? w.tmLookup.mp : null;
    var examples = w.examples ? w.examples.compounds : [];

    function solvent(x) {
      if (x && typeof x === 'object' && x.id && byId[x.id]) return byId[x.id];
      var s = byId[x] || byKey[String(x).toLowerCase()];
      if (!s) throw new Error('溶媒が表にない: ' + x + '（weights/solvents.json の id・name・ja）');
      return s;
    }
    function tmOf(sol) { return isNum(sol.TmC) ? sol.TmC : null; }
    function norm(sol) { return { M: sol.M, logPow: sol.logPow, TmC: tmOf(sol), abr: sol.abr, x: sol.x }; }
    // 使うモデルの名前：融点あり／なし × 水／有機溶媒。学習し直した重みに無いモデルは、ある方を使う（自社のデータだけで学習した時など）
    function modelKey(has, water) {
      var sfx = water && (G.tm_w || G.notm_w) ? '_w' : '', k = (has ? 'tm' : 'notm') + sfx;
      return G[k] ? k : (G[(has ? 'notm' : 'tm') + sfx] ? (has ? 'notm' : 'tm') + sfx : (has ? 'tm' : 'notm'));
    }
    function pickModel(has, water) { return G[modelKey(has, water)]; }
    function cvOf(has) { return G.cv[has ? 'tm' : 'notm'] && G.cv[has ? 'tm' : 'notm'].organic ? G.cv[has ? 'tm' : 'notm'] : G.cv[has ? 'notm' : 'tm']; }
    function transfer(sv, sol) { return sv.water ? 0 : AB.logP(sv, sol.abr); }
    // 幅（外れの 68%）。学習し直した重みで、その種類の検証の行が 0 件の時（自社のデータが有機溶媒だけ、など）は、
    // 検証できた幅のうちいちばん広い物（どれも無ければ 1.0 = 10 倍）を使い、bandFallback の印を付ける
    function p68(o) { return o && isNum(o.p68) ? o.p68 : null; }
    var FALLBACK_BAND = (function () {
      var v = [];
      ['tm', 'notm'].forEach(function (k) { var c = G.cv && G.cv[k]; if (c) ['organic', 'water', 'one_point_organic'].forEach(function (s) { var x = p68(c[s]); if (x !== null) v.push(x); }); });
      return v.length ? Math.max.apply(null, v) : 1.0;
    })();
    function bandOf(cv, kind) { var x = p68(cv && cv[kind]); return x !== null ? { b: x, fb: false } : { b: FALLBACK_BAND, fb: true }; }

    // 学習済みモデルの値（log10 mol/L）と、式の内訳
    function model(sol, sv) {
      var has = sol.TmC !== null, y = SF.evalModel(pickModel(has, sv.water), SF.build(sol, sv, K, TMDEF, T));
      var p = SC.predict(sol, sv.water ? { water: true } : { transfer: transfer(sv, sol) }, T, K, TMDEF), parts = copy(p.parts);
      parts.corr = y - p.raw;
      return { y: y, has: has, parts: parts };
    }
    // ページの estimate と同じ。返す logS は log10(mol/L)、band は log10 の片側の幅
    function core(sol, sv, ref) {
      var m = model(sol, sv), cv = cvOf(m.has), cap = SC.logSmax(sol.M);
      if (ref) {
        if (ref.sv.id === sv.id) return { logS: ref.logS, band: 0, ref: true, same: true, has: m.has };
        var delta = m.y - model(sol, ref.sv).y, raw = ref.logS + delta;      // 実測 1 点に、モデルが見込む溶媒どうしの差を足す
        var bo = bandOf(cv, (ref.sv.water || sv.water) ? 'water' : 'one_point_organic');
        return { logS: Math.min(raw, cap), band: bo.b, bandFallback: bo.fb, ref: true, capped: raw > cap, delta: delta, has: m.has };
      }
      var bm = bandOf(cv, sv.water ? 'water' : 'organic');
      return { logS: Math.min(m.y, cap), band: bm.b, bandFallback: bm.fb, parts: m.parts, capped: m.y > cap, has: m.has };
    }
    function refOf(sol, r) {
      if (!r) return null;
      var v = +r.mgPerML; if (!(isNum(v) && v > 0)) throw new Error('reference.mgPerML は正の数');
      return { sv: solvent(r.solvent), logS: Math.log10(v / sol.M) };
    }
    function measOf(sol, id) { var m = sol.meas && sol.meas[id]; return m ? { logS: m[0], mgPerML: Math.pow(10, log10mg(m[0], sol.M)), n: m[1], doi: m[2] } : null; }
    function best(sol, sv, ref, useMeas) {                     // 混合溶媒・希釈の端の値：実測があれば実測（幅 0.1）
      var m = useMeas !== false ? measOf(sol, sv.id) : null;
      return m ? { logS: m.logS, band: MEAS_BAND, meas: true } : core(norm(sol), sv, ref);
    }
    function verdict(lmg, band, target) { return isNum(target) && target > 0 ? SC.verdict(lmg, Math.log10(target), band) : null; }
    function range(lmg, band) { return [Math.pow(10, lmg - band), Math.pow(10, lmg + band)]; }

    // 1 つの溶媒での見込み
    function estimate(sol, svx, opt) {
      opt = opt || {};
      var sv = solvent(svx), s = norm(sol), ref = refOf(s, opt.reference), e = core(s, sv, ref), lmg = log10mg(e.logS, s.M), ms = measOf(sol, sv.id);
      var r = { solvent: sv.id, logS: e.logS, mgPerML: Math.pow(10, lmg), band: e.band, factor: Math.pow(10, e.band), range: range(lmg, e.band),
        model: modelKey(e.has, sv.water), tmUsed: e.has, capped: !!e.capped, parts: e.parts || null, fromReference: !!e.ref, delta: isNum(e.delta) ? e.delta : null,
        measured: ms, extrapolated: !(sv.n_rows > 0), bandFallback: !!e.bandFallback };
      if (isNum(opt.target) && opt.target > 0) {   // 目標の濃度は正の数の時だけ判定する（0 や負の数では判定を付けない）
        if (ms) { r.verdict = SC.verdict(log10mg(ms.logS, s.M), Math.log10(opt.target), MEAS_BAND); r.basis = 'measured'; }
        else { r.verdict = verdict(lmg, e.same ? MEAS_BAND : e.band, opt.target); r.basis = e.ref ? 'reference' : 'model'; }
      }
      return r;
    }
    function estimateAll(sol, opt) {
      var list = (opt && opt.solvents) ? opt.solvents.map(solvent) : main;
      return list.map(function (sv) { return estimate(sol, sv, opt); });
    }
    // 混合溶媒（2 成分）。fB は溶媒 B の体積分率。対数線形の内挿に、溶媒の組ごとの補正を足す（ページの mixAt と同じ）
    function mixAt(sol, a, b, ea, eb, f) {
      var ya = log10mg(ea.logS, sol.M), yb = log10mg(eb.logS, sol.M);
      if (a.id === b.id) return { y: ya, band: ea.band };
      var comp = {}; comp[a.id] = 1 - f; comp[b.id] = f;
      var c = SC.mixCorr(comp, MIX), mb = (a.water || b.water) ? MIX.band.corrected_wo : MIX.band.corrected_oo, eb2 = (1 - f) * ea.band + f * eb.band;
      return { y: Math.min(ya + (yb - ya) * f + c.corr, log10mg(SC.logSmax(sol.M), sol.M)), band: Math.sqrt(eb2 * eb2 + mb * mb * c.w * c.w), corr: c.corr };
    }
    function mix(sol, ax, bx, fB, opt) {
      opt = opt || {};
      if (!(fB >= 0 && fB <= 1)) throw new Error('fB は 0〜1');
      var a = solvent(ax), b = solvent(bx), s = norm(sol), ref = refOf(s, opt.reference);
      var ea = best(sol, a, ref, opt.useMeasured), eb = best(sol, b, ref, opt.useMeasured), q = mixAt(s, a, b, ea, eb, fB);
      var r = { solvents: [a.id, b.id], fB: fB, mgPerML: Math.pow(10, q.y), band: q.band, factor: Math.pow(10, q.band), range: range(q.y, q.band), correction: q.corr || 0, ends: { a: ea, b: eb } };
      if (isNum(opt.target)) r.verdict = verdict(q.y, q.band, opt.target);
      return r;
    }
    function mixCurve(sol, ax, bx, opt) {
      var n = (opt && opt.points) || 51, out = [];
      for (var i = 0; i < n; i++) out.push(mix(sol, ax, bx, i / (n - 1), opt));
      return out;
    }
    // 原液を薄めた時に析出するか（ページの図 3 と同じ。3 成分は 2 成分の組の和で近似し、補正の幅を 2 倍に取る。3 成分の実測との照合はしていない）
    function dilution(sol, o) {
      var sv = solvent(o.stockSolvent), s = norm(sol), ref = refOf(s, o.reference), dil = [];
      var D = o.diluent || { water: 1 };
      if (Array.isArray(D)) D.forEach(function (d) { dil.push({ key: solvent(d.solvent).id, f: +d.f }); });
      else for (var k in D) dil.push({ key: solvent(k).id, f: +D[k] });
      // 入力の確かめ：原液の濃度と量は正の数、薄める液の量と成分の比は 0 以上、薄める液を足すなら成分の比の合計は正
      var conc = +o.stockMgPerML, take = +o.stockUL, vol = o.diluentUL === undefined || o.diluentUL === null ? 0 : +o.diluentUL;
      if (!(isNum(conc) && conc > 0)) throw new Error('stockMgPerML（原液の濃度）は正の数: ' + o.stockMgPerML);
      if (!(isNum(take) && take > 0)) throw new Error('stockUL（原液の量）は正の数: ' + o.stockUL);
      if (!(isNum(vol) && vol >= 0)) throw new Error('diluentUL（薄める液の量）は 0 以上の数: ' + o.diluentUL);
      var fsum = 0;
      dil.forEach(function (x) { if (!(isNum(x.f) && x.f >= 0)) throw new Error('diluent の成分の比は 0 以上の数: ' + x.key + ' = ' + x.f); fsum += x.f; });
      if (vol > 0 && !(fsum > 0)) throw new Error('diluent の成分の比の合計が 0（薄める液の組成が決まらない）');
      var d = SC.dilute(sv.id, conc, take, dil, vol);
      if (!d) throw new Error('量が正しくない');
      var items = Object.keys(d.comp).filter(function (k) { return d.comp[k] > 0; }).map(function (k) { var e = best(sol, byId[k], ref, o.useMeasured); return { id: k, f: d.comp[k], logS: log10mg(e.logS, s.M), band: e.band }; });
      var c = SC.mixCorr(d.comp, MIX), mb = (d.comp.water > 0 ? MIX.band.corrected_wo : MIX.band.corrected_oo) * (items.length > 2 ? 2 : 1), eb2 = items.reduce(function (t, i) { return t + i.f * i.band; }, 0);
      var y = Math.min(SC.mix(items) + c.corr, log10mg(SC.logSmax(s.M), s.M)), band = Math.sqrt(eb2 * eb2 + mb * mb * c.w * c.w), v = SC.verdict(y, Math.log10(d.conc), band);
      return { factor: d.factor, concMgPerML: d.conc, comp: d.comp, solubilityMgPerML: Math.pow(10, y), band: band, range: range(y, band),
        verdict: v === 'yes' ? 'stays_dissolved' : v === 'no' ? 'precipitates' : 'edge', threeComponent: items.length > 2 };
    }

    // 構造（SMILES）から溶質を作る。RDKit は RDKit.js の initRDKitModule() の結果。融点は ページと同じ順：入れた値（opt.TmC）＞ 溶解度の照合の表（lookup_25C）
    // ＞ 融点の照合の表（tm_lookup）。どれもなければ null で、tmHint に推算値（5 °C 刻み。計算には使わない）。名前は opt.name か、なければ SMILES
    function prepare(RDKit, smiles, opt) {
      opt = opt || {};
      var text = String(smiles || '').trim();
      if (!text) throw fail('empty');
      if (text.indexOf('.') >= 0) throw fail('multi');
      var m = null; try { m = RDKit.get_mol(text); } catch (e) { m = null; }
      if (!m || !m.is_valid()) { if (m) m.delete(); throw fail('parse'); }
      var can, f, d;
      try { m.remove_hs_in_place(); can = m.get_smiles(); f = FE.compute(RDKit, m); d = JSON.parse(m.get_descriptors()); } finally { m.delete(); }
      if (d.NumHeavyAtoms > 80) throw fail('big');
      var V = AB.mcgowanV(RDKit, can); if (V === null) throw fail('elem');
      var x = f.desc.concat(f.groups), abr = AB.predictSolute(w.abraham, x); abr.V = V;
      var hit = lookup ? lookup[can] || null : null;
      var TmC = hit && isNum(hit.t) ? hit.t : null, tmSrc = TmC !== null ? (hit.s || 'sol') : null;
      if (TmC === null && tmLookup && isNum(tmLookup[can])) { TmC = tmLookup[can]; tmSrc = 'mp'; }
      if (isNum(opt.TmC)) { TmC = opt.TmC; tmSrc = 'user'; }
      var tmHint = null;
      if (TmC === null && w.tm) tmHint = Math.round(SF.evalModel(w.tm.model, x) / 5) * 5;
      return { name: opt.name || can, smiles: can, M: +d.amw.toFixed(2), logPow: d.CrippenClogP, TmC: TmC, tmSrc: tmSrc, tmHint: tmHint, abr: { E: abr.E, S: abr.S, A: abr.A, B: abr.B, V: abr.V }, x: x, meas: hit ? hit.m : {}, inLookup: !!hit };
    }
    function fromExample(key) {
      var k = String(key).toLowerCase();
      for (var i = 0; i < examples.length; i++) { var e = examples[i]; if (e.id === key || e.name.toLowerCase() === k || e.smiles === key) { var c = copy(e); c.inLookup = true; return c; } }
      throw new Error('例の固体にない: ' + key + '（weights/example_compounds.json）');
    }
    function withTm(sol, TmC) { var c = copy(sol); c.TmC = isNum(TmC) ? TmC : null; c.tmSrc = c.TmC === null ? null : 'user'; return c; }
    function estimateTm(sol) {
      if (!w.tm) throw new Error('融点のモデル（tm_model.json）を読んでいない');
      return { TmC: SF.evalModel(w.tm.model, sol.x), band: w.tm.cv ? w.tm.cv.p68 : null, note: '目安。溶解度の計算には使わない（推算した融点では幅が狭くならなかった）' };
    }

    return { version: VERSION, weightsVersion: G.meta && G.meta.version, T_C: T, solvents: solvents, mainSolvents: main, examples: examples,
      solvent: solvent, prepare: prepare, fromExample: fromExample, withTm: withTm, estimate: estimate, estimateAll: estimateAll, mix: mix, mixCurve: mixCurve,
      dilution: dilution, estimateTm: estimateTm, lookup: function (can) { return lookup ? lookup[can] || null : null; },
      lookupTm: function (can) { return tmLookup && isNum(tmLookup[can]) ? tmLookup[can] : null; },
      verdict: function (mgPerML, band, target) { return verdict(Math.log10(mgPerML), band, target); },
      _raw: { model: model, core: core, pickModel: pickModel } };
  }

  // 重みを読む（ブラウザ・Deno・新しい Node の fetch）。opt.lookup / opt.tmLookup / opt.tm / opt.examples を false にすると、その表を読まない
  function load(base, opt) {
    opt = opt || {}; base = base || 'weights/'; if (base.slice(-1) !== '/') base += '/';
    var keys = REQUIRED.concat(['tm', 'lookup', 'tmLookup', 'examples'].filter(function (k) { return opt[k] !== false; }));
    return Promise.all(keys.map(function (k) {
      return fetch(base + FILES[k]).then(function (r) { if (!r.ok) throw new Error(FILES[k] + ': ' + r.status); return r.json(); })
        .catch(function (e) { if (REQUIRED.indexOf(k) >= 0) throw e; return null; });
    })).then(function (vals) { var w = {}; keys.forEach(function (k, i) { w[k] = vals[i]; }); return create(w); });
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
  if (isNode) module.exports = api; else root.SolubilityModel = api;
})(typeof self !== 'undefined' ? self : this);
