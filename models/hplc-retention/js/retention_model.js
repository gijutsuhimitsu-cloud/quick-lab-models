// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE。
// 保持の予測（ページから呼ぶ入口）。第 2 版の式（2026-11-24）。係数ファイルは meta.version が 2 の物を使う。
// 特徴量の第 2 版（meta.feature_version が 2）では、化合物にイオン化の規則ごとの数（compound.ion。ChromFeatures.compute が返す）が要る。
// 学習済みの木のモデル（model_*.json）とカラムの表（columns.json）を読んで使う。
//
//   var P = ChromRetention.create(modelJson, columnsJson);
//   var r = P.predictRetention(compound, columnId, cond);
//     compound : { desc:[41], groups:[56] }  … ChromFeatures.compute(RDKit, mol) の結果、または stock_compounds.json の 1 件
//     columnId : columns.json の id（文字列）
//     cond     : { organic:'acn'|'meoh', phi:0〜1（等組成）, gradient:[[分, phi], ...]（勾配。あれば phi より優先）,
//                  t0:分（省略時はカラム寸法と流量から）, dwell:分（省略時 0）, tempC, pH,
//                  additive:'formic'|'acetic'|'tfa'|'ammonium'|'basic'|'phosphate'|'none',
//                  flow:mL/min, length_mm, id_mm, particle_um }
//   戻り値 : { lkw, S, k（等組成のとき）, tR（分）, t0, N（段数の目安）, sigma（ピークの標準偏差・分）, kElution, inRange }
//   予測の幅（任意）: create の 3 つ目に uncertainty_*.json を渡すと、戻り値に uncertainty が付く。
//     uncertainty : { seen（その化合物が学習データに入っているか。compound.seen による）, s68, s90（lkw の外れの見積もり。外れの 68%・90% がこの値より小さい）,
//                     tR68:[早い側, 遅い側], tR90:[..], k68:[小, 大], k90:[..]（k は等組成のとき） }
//
// グラフのモデル（任意）: create の 4 つ目に ChromGNN.create(manifest, weights) を渡すと、compound.z（分子ごとに 1 回 P.embed で作る）を持つ
//   化合物では、lkw をグラフのモデルで計算する（戻り値の model が 'gnn'。z の無い化合物は木のモデルのまま＝'tree'）。
//     var P = ChromRetention.create(modelJson, columnsJson, uncJson, ChromGNN.create(gnnManifest, gnnWeights));
//     P.embed(compound, ChromGraph.fromSmiles(RDKit, compound.smiles));      // compound.z が入る。構造を変えた時だけ
//
// 測定値 1 点での合わせ直し（任意）: 利用者が、その分子を自分の装置で 1 回測った値を入れると、同じ分子のほかの条件の予測を合わせ直す。
//     var c = P.calibrate(compound, columnId, cond, { tR: 分 });      // 等組成なら { k: 保持係数 } でもよい。cond は測った時の条件
//       戻り値 { ok, dlkw（実測 − 予測。lkw の単位）, lkwMeasured, lkwPredicted, ratio90（|dlkw| ÷ 90% の幅。1 を大きく超えたら入力の間違いを疑う）, reason（ok が false の時） }
//       compound.cal に 1 件足される（何件でも足せる。予測の時は、当てる条件にいちばん近い 1 件を使う）。消すのは P.clearCalibration(compound)。
//     そのあとの predictRetention は、lkw に w·dlkw を足す。w は 0〜0.85 で、測った条件と当てる条件が近いほど大きい
//       （同じカラム・同じ有機溶媒・同じイオン化の状態で 0.85。カラムが別なら ×0.89、有機溶媒が別なら ×0.73、イオン化の状態が違えば ×0.70 まで下がる）。
//     戻り値に calibration: { w, dlkw, shift（実際に足した量）, lkwRaw（合わせ直す前の lkw）, point（使った測定の番号） } が付き、予測の幅は sqrt(1 − w²) 倍に狭くなる。
//     cond.calibrate を false にすると、その 1 回だけ合わせ直しを使わない（比べたい時用）。
//     係数は、公開の重みと同じ学習データ（継承のないデータだけ）での交差検証の予測から決めた（train/scripts/fit_calibration.py）。
//     uncertainty の meta.calibration があれば、そちらを使う（自社データで学習し直して当て直した時用）。
//
// 保持の式（model.meta.par に定数 a, b, B。有機溶媒ごと）:
//   u = phi / (1 + B·phi),  S = a + b·lkw,  log10 k = lkw − S·u
// モデルが出すのは lkw（見かけの log kw）だけ。勾配は ∫ dt / (t0·k(phi(t))) = 1 を数値積分で解く。
(function (root) {
  'use strict';
  var HSM = ['H', 'A', 'B', 'C', 'V', 'D'];
  var PHASES = ['C18', 'C8', 'Phenyl', 'CN', 'EP', 'F', 'short'];
  var ADDS = ['formic', 'acetic', 'tfa', 'ammonium', 'basic', 'phosphate'];
  // 測定値 1 点での合わせ直しの係数（継承のないデータだけで当て直した値。train/scripts/fit_calibration.py の rounded）
  var CAL_DEFAULT = { w0: 0.85, col: 0.89, solv: 0.73, floor: 0.70, tau: 0.02 };

  function evalTrees(model, x) {
    var s = model.baseline, trees = model.trees;
    for (var t = 0; t < trees.length; t++) {
      var tr = trees[t], n = 0;
      for (;;) {
        var f = tr.f[n];
        if (f < 0) { s += tr.v[n]; break; }
        var xv = x[f];
        var left = (xv !== xv) ? (tr.m[n] === 1) : (Math.fround(xv) <= tr.t[n]);
        n = left ? tr.l[n] : tr.r[n];
      }
    }
    return s;
  }

  function uOf(phi, B) { return phi / (1 + B * phi); }
  function logK(lkw, phi, par) { var u = uOf(phi, par.B); return lkw * (1 - par.b * u) - par.a * u; }

  function phiAt(prog, t, dwell) {
    var d = dwell || 0;
    if (t <= prog[0][0] + d) return prog[0][1];
    for (var i = 0; i + 1 < prog.length; i++) {
      var ta = prog[i][0] + d, tb = prog[i + 1][0] + d;
      if (t <= tb) return tb > ta ? prog[i][1] + (prog[i + 1][1] - prog[i][1]) * (t - ta) / (tb - ta) : prog[i + 1][1];
    }
    return prog[prog.length - 1][1];
  }

  // 保持時間（分）。prog は時刻の順に並んだ [[分, phi], ...]。等組成は 1 点だけ。
  function retentionTime(lkw, t0, prog, par, dwell, nSteps) {
    if (prog.length === 1) return t0 * (1 + Math.pow(10, logK(lkw, prog[0][1], par)));
    var tEnd = prog[prog.length - 1][0] + (dwell || 0);
    var n = nSteps || 1500, dt = tEnd / n, cum = 0, fPrev = Math.pow(10, -logK(lkw, phiAt(prog, 0, dwell), par)) / t0;
    for (var i = 1; i <= n; i++) {
      var f = Math.pow(10, -logK(lkw, phiAt(prog, i * dt, dwell), par)) / t0;
      var step = (f + fPrev) * 0.5 * dt;
      if (cum + step >= 1) return t0 + (i - 1) * dt + (1 - cum) / step * dt;
      cum += step; fPrev = f;
    }
    return t0 + tEnd + (1 - cum) / fPrev;      // 表の終わりのあとは最後の割合のまま
  }

  // 保持時間（分）→ lkw（retentionTime の逆。lkw について単調なので二分法）。範囲（−3〜14）の端に当たったら edge に 'early'・'late'
  function lkwFromRetentionTime(tR, t0, prog, par, dwell) {
    var lo = -3, hi = 14;
    if (!(tR > retentionTime(lo, t0, prog, par, dwell))) return { lkw: lo, edge: 'early' };
    if (!(tR < retentionTime(hi, t0, prog, par, dwell))) return { lkw: hi, edge: 'late' };
    for (var i = 0; i < 48; i++) { var mid = (lo + hi) / 2; if (retentionTime(mid, t0, prog, par, dwell) < tR) lo = mid; else hi = mid; }
    return { lkw: (lo + hi) / 2, edge: null };
  }

  // 死時間の目安（分）: 空隙率 0.65 の充填カラム
  function deadTime(lengthMm, idMm, flow) {
    return Math.PI * Math.pow(idMm / 20, 2) * (lengthMm / 10) * 0.65 / flow;
  }

  // 段数の目安: 換算段高 2.5（粗い目安。学習した値ではない）
  function plates(lengthMm, particleUm) { return lengthMm * 1000 / (2.5 * particleUm); }

  // イオン化の特徴量（学習側の scripts/ionfeat.py と同じ計算）。counts は規則ごとの部位の数。
  function ionFeatures(meta, counts, pH, raw) {
    var R = meta.ion_rules, K = meta.ion_const, i, rep;
    var ph = (pH === undefined || pH === null || pH !== pH) ? K.ph_default : pH;
    var names = R.map(function (r) { return r[0]; });
    var beta = counts[names.indexOf('amine_beta_hetero')] > 0;
    var pk = [], ty = [], perm = 0;
    for (i = 0; i < R.length; i++) {
      var t = R[i][1], c = counts[i] | 0;
      if (t === 'perm+') { perm += c; continue; }
      if (t !== 'acid' && t !== 'base') continue;
      for (rep = 1; rep <= K.max_rep; rep++) {
        if (c >= rep && pk.length < K.max_groups) {
          pk.push(R[i][2] + ((K.amines.indexOf(R[i][0]) >= 0 && beta) ? K.beta_shift : 0)); ty.push(t === 'acid' ? -1 : 1);
        }
      }
    }
    var qneg = 0, qpos = perm, logneu = 0, na = 0, nb = 0, pa = Infinity, pb = -Infinity;
    for (i = 0; i < pk.length; i++) {
      var f;
      if (ty[i] < 0) { f = 1 / (1 + Math.pow(10, pk[i] - ph)); qneg += f; na++; pa = Math.min(pa, pk[i]); }
      else { f = 1 / (1 + Math.pow(10, ph - pk[i])); qpos += f; nb++; pb = Math.max(pb, pk[i]); }
      logneu += Math.log10(Math.max(1 - f, 1e-6));
    }
    if (perm > 0) logneu = -6;
    if (pa === Infinity) pa = 14;
    if (pb === -Infinity) pb = 0;
    var da = pa - ph, db = ph - pb, cap = K.cap_no_group;
    if (cap !== undefined && cap !== null) { if (na === 0) da = cap; if (nb === 0) db = cap; }
    var all = [qneg, qpos, qpos - qneg, qpos + qneg, logneu, na, nb, pa, pb, da, db, qpos * qneg];
    return (K.keep && !raw) ? K.keep.map(function (j) { return all[j]; }) : all;
  }

  function create(model, columns, unc, gnn) {
    if (!model.meta || model.meta.version !== 2) throw new Error('係数ファイルの版が合わない（meta.version が 2 の model_*.json が要る）');
    var colById = {};
    columns.columns.forEach(function (c) { colById[c.id] = c; });
    var nComp = model.n_comp, PAR = model.meta.par;
    if (unc && (unc.meta.feature_version !== model.meta.feature_version || unc.meta.n_features !== model.n_features || unc.meta.name !== model.meta.name))
      throw new Error('予測の幅の係数ファイルが、モデルの係数ファイルと合わない（同じ版・同じ名前の uncertainty_*.json が要る）');
    var nStruct = model.n_comp + (model.n_ion || 0), nMethod = model.method_features.length;
    var CAL = (unc && unc.meta.calibration) || CAL_DEFAULT;
    if (gnn) {
      var gm = gnn.manifest;
      if (gm.n_desc !== nStruct || gm.feature_version !== model.meta.feature_version || gm.method_features.join() !== model.method_features.join() || gm.group !== model.meta.name)
        throw new Error('グラフのモデルが、木のモデルの係数ファイルと合わない（同じ版・同じ名前・同じ特徴量の物が要る）');
    }
    function embed(compound, graph) {
      if (!gnn) throw new Error('グラフのモデルが渡されていない');
      if (!graph) { compound.z = null; return null; }
      var d = [].concat(Array.prototype.slice.call(compound.desc), Array.prototype.slice.call(compound.groups), Array.prototype.slice.call(compound.ion));
      compound.z = gnn.encode(graph, d);
      return compound.z;
    }

    function vector(compound, col, cond) {
      var x = new Array(model.n_features), i, k = 0;
      for (i = 0; i < compound.desc.length; i++) x[k++] = compound.desc[i];
      for (i = 0; i < compound.groups.length; i++) x[k++] = compound.groups[i];
      if (k !== nComp) throw new Error('化合物の特徴量の数が合わない: ' + k + ' / ' + nComp);
      var fv2 = model.meta.feature_version >= 2;
      if (fv2) {
        if (!compound.ion || compound.ion.length !== model.n_ion) throw new Error('イオン化の規則ごとの数（compound.ion）が無い・数が合わない');
        for (i = 0; i < compound.ion.length; i++) x[k++] = compound.ion[i];
      }
      var add = cond.additive || 'formic', MF = model.method_features;
      for (i = 0; i < MF.length; i++) {
        var n = MF[i], v;
        if (n === 'logkEB') v = col.logkEB;
        else if (HSM.indexOf(n) >= 0) v = col[n];
        else if (n.indexOf('phase_') === 0) v = col.phase === n.substr(6) ? 1 : 0;
        else if (n === 'frac_acn') v = cond.organic === 'meoh' ? 0 : 1;
        else if (n === 'pH_weak') v = (cond.pH === undefined || cond.pH === null) ? NaN : cond.pH;
        else if (n === 'temp_C') v = (cond.tempC === undefined || cond.tempC === null) ? NaN : cond.tempC;
        else if (n.indexOf('add_') === 0) v = add === n.substr(4) ? 1 : 0;
        else throw new Error('知らないメソッドの特徴量: ' + n);
        x[k++] = v;
      }
      if (fv2) { var f = ionFeatures(model.meta, compound.ion, cond.pH); for (i = 0; i < f.length; i++) x[k++] = f[i]; }
      return x;
    }

    // モデルそのままの lkw（合わせ直しなし）
    function rawLkw(compound, col, cond) {
      var xvec = vector(compound, col, cond);
      var lkw = evalTrees(model, xvec), lkwTree = lkw, used = 'tree';
      if (gnn && compound.z) {
        lkw = Math.max(-3, Math.min(14, gnn.lkw(compound.z, xvec.slice(nStruct, nStruct + nMethod), xvec.slice(nStruct + nMethod))));
        used = 'gnn';
      }
      return { lkw: lkw, lkwTree: lkwTree, used: used, xvec: xvec };
    }
    function widths(xvec, compound) {
      var s68 = Math.max(evalTrees(unc.s68, xvec), unc.meta.floor) * unc.s68.scale;
      var s90 = Math.max(evalTrees(unc.s90, xvec), unc.meta.floor) * unc.s90.scale;
      var seen = compound.seen === true || (compound.seen && compound.seen[model.meta.name] === true) ||
        !!(unc.meta.seen_stock && compound.inchikey && unc.meta.seen_stock.indexOf(String(compound.inchikey).slice(0, 14)) >= 0);
      return { s68: s68, s90: Math.max(s90, s68), seen: !!seen };
    }
    // 割合は 0〜1（% ではない）。範囲の外や数でない値は、黙って計算せずに止める
    function fraction(v, what) {
      var x = +v;
      if (!(x >= 0 && x <= 1)) throw new Error(what + ' は 0〜1 の割合で書く（50% なら 0.5）: ' + v);
      return x;
    }
    function program(cond) {
      if (!(cond.gradient && cond.gradient.length)) return [[0, cond.phi === undefined ? 0.5 : fraction(cond.phi, 'phi')]];
      var prev = 0;
      return cond.gradient.map(function (p, i) {
        var t = +p[0];
        if (!(t >= prev)) throw new Error('gradient の時間は 0 以上で、増える順に書く（' + (i + 1) + ' 点目）: ' + p[0]);
        prev = t;
        return [t, fraction(p[1], 'gradient の ' + (i + 1) + ' 点目の割合')];
      });
    }
    // 寸法・流量・時間：書いてなければ既定の値。書いてあれば正の数（遅れ時間 dwell は 0 以上）
    function positive(v, def, what, zeroOk) {
      if (v === undefined || v === null || v === '') return def;
      var x = +v;
      if (!(isFinite(x) && (zeroOk ? x >= 0 : x > 0))) throw new Error(what + ' は' + (zeroOk ? ' 0 以上' : '正') + 'の数で書く: ' + v);
      return x;
    }
    function setup(cond) {
      if (cond.organic !== undefined && cond.organic !== 'acn' && cond.organic !== 'meoh') throw new Error("organic は 'acn' か 'meoh': " + cond.organic);
      var L = positive(cond.length_mm, 100, 'length_mm'), id = positive(cond.id_mm, 2.1, 'id_mm'), dp = positive(cond.particle_um, 2.7, 'particle_um'),
        flow = positive(cond.flow, 0.3, 'flow');
      positive(cond.dwell, 0, 'dwell', true);
      var t0 = positive(cond.t0, null, 't0');
      return { L: L, dp: dp, t0: t0 !== null ? t0 : deadTime(L, id, flow), par: PAR[cond.organic === 'meoh' ? 'meoh' : 'acn'], prog: program(cond) };
    }

    // ---- 測定値 1 点での合わせ直し ----
    function charges(compound, cond) {
      if (!(model.meta.feature_version >= 2) || !compound.ion) return [0, 0];
      var a = ionFeatures(model.meta, compound.ion, cond.pH, true);
      return [a[0], a[1]];
    }
    function sameColumn(a, b) {
      if (a.id === b.id) return true;
      var ks = HSM.concat(['logkEB']);
      for (var i = 0; i < ks.length; i++) if (Math.abs(a[ks[i]] - b[ks[i]]) > 5e-5) return false;
      return true;
    }
    // 合わせ直しの割合 w（測った条件 pt と、当てる条件がどのくらい近いか）
    function calWeight(compound, col, cond, pt) {
      var pc = colById[pt.columnId];
      if (!pc) return -1;
      var w = CAL.w0;
      if (!sameColumn(col, pc)) w *= CAL.col;
      if ((cond.organic === 'meoh') !== (pt.cond.organic === 'meoh')) w *= CAL.solv;
      var q = charges(compound, cond), dq = Math.abs(q[0] - pt.q[0]) + Math.abs(q[1] - pt.q[1]);
      return w * (CAL.floor + (1 - CAL.floor) * Math.exp(-dq / CAL.tau));
    }
    function calibrate(compound, columnId, cond, measured) {
      var col = colById[columnId];
      if (!col) throw new Error('カラムが見つからない: ' + columnId);
      cond = cond || {}; measured = measured || {};
      var st = setup(cond), lkwM, hasK = measured.k !== undefined && measured.k !== null;
      if (hasK) {
        if (st.prog.length !== 1) throw new Error('保持係数 k で合わせ直せるのは等組成の時だけ（勾配は保持時間 tR で）');
        if (!(measured.k > 0)) return { ok: false, reason: 'not_retained' };
        var u = uOf(st.prog[0][1], st.par.B);
        lkwM = (Math.log(measured.k) / Math.LN10 + st.par.a * u) / (1 - st.par.b * u);
        if (!(lkwM > -3)) return { ok: false, reason: 'too_early' };
        if (!(lkwM < 14)) return { ok: false, reason: 'too_late' };
      } else if (measured.tR !== undefined && measured.tR !== null) {
        if (!(measured.tR > st.t0)) return { ok: false, reason: 'not_retained' };          // 死時間より早い（保持されていない）
        var inv = lkwFromRetentionTime(measured.tR, st.t0, st.prog, st.par, cond.dwell || 0);
        if (inv.edge) return { ok: false, reason: inv.edge === 'early' ? 'too_early' : 'too_late' };
        lkwM = inv.lkw;
      } else throw new Error('測定値が無い（{ tR: 分 } か { k: 保持係数 } を渡す）');
      var r0 = rawLkw(compound, col, cond), c2 = {}, k;
      for (k in cond) if (k !== 'calibrate') c2[k] = cond[k];
      var pt = { columnId: columnId, cond: c2, measured: hasK ? { k: measured.k } : { tR: measured.tR },
        lkwMeasured: lkwM, lkwPredicted: r0.lkw, dlkw: lkwM - r0.lkw, model: r0.used, q: charges(compound, cond) };
      compound.cal = (compound.cal || []).concat([pt]);
      var res = { ok: true, dlkw: pt.dlkw, lkwMeasured: lkwM, lkwPredicted: r0.lkw, model: r0.used, point: compound.cal.length - 1 };
      if (unc) { var wd = widths(r0.xvec, compound); res.ratio90 = Math.abs(pt.dlkw) / wd.s90; res.within68 = Math.abs(pt.dlkw) <= wd.s68; res.within90 = Math.abs(pt.dlkw) <= wd.s90; }
      return res;
    }
    function clearCalibration(compound) { compound.cal = null; }

    function predictRetention(compound, columnId, cond) {
      var col = colById[columnId];
      if (!col) throw new Error('カラムが見つからない: ' + columnId);
      cond = cond || {};
      var st = setup(cond), par = st.par;
      var r0 = rawLkw(compound, col, cond), xvec = r0.xvec, lkw = r0.lkw, lkwTree = r0.lkwTree, used = r0.used, calInfo = null;
      if (compound.cal && compound.cal.length && cond.calibrate !== false) {
        var bi = -1, bw = -1;
        for (var ci = 0; ci < compound.cal.length; ci++) { var wi = calWeight(compound, col, cond, compound.cal[ci]); if (wi > bw) { bw = wi; bi = ci; } }
        if (bi >= 0 && bw > 0) {
          var pt = compound.cal[bi];
          if (pt.model !== used) {      // 合わせ直した時と使うモデルが違う（z を足した・外した）: その条件の予測を取り直す
            var rr = rawLkw(compound, colById[pt.columnId], pt.cond);
            pt.lkwPredicted = rr.lkw; pt.model = rr.used; pt.dlkw = pt.lkwMeasured - rr.lkw;
          }
          lkw = Math.max(-3, Math.min(14, r0.lkw + bw * pt.dlkw));
          calInfo = { w: bw, dlkw: pt.dlkw, shift: lkw - r0.lkw, lkwRaw: r0.lkw, point: bi };
        }
      }
      var L = st.L, dp = st.dp, t0 = st.t0, prog = st.prog;
      var out = { lkw: lkw, logkw: lkw, S: par.a + par.b * lkw, t0: t0, column: col.name, model: used, lkwTree: lkwTree };
      if (calInfo) out.calibration = calInfo;
      out.tR = retentionTime(lkw, t0, prog, par, cond.dwell || 0);
      var phiE = prog.length === 1 ? prog[0][1] : phiAt(prog, Math.max(out.tR - t0, 0), cond.dwell || 0);
      out.kElution = Math.pow(10, logK(lkw, phiE, par));
      if (prog.length === 1) out.k = out.kElution;
      out.N = plates(L, dp);
      out.sigma = t0 * (1 + out.kElution) / Math.sqrt(out.N);     // 等組成は tR/√N。勾配は溶出時の k での近似
      if (unc) {
        var isIso = prog.length === 1, dw = cond.dwell || 0;
        var at = function (l) { return retentionTime(l, t0, prog, par, dw); };
        var kAt = function (l) { return Math.pow(10, logK(l, prog[0][1], par)); };
        var wd = widths(xvec, compound), s68 = wd.s68, s90 = wd.s90, seen = wd.seen;
        var a68 = s68, a90 = s90;
        if (seen && unc.meta.seen_scale) { a68 = s68 * unc.meta.seen_scale.s68; a90 = s90 * unc.meta.seen_scale.s90; }      // 学習データに入っている化合物は、外れが小さい
        if (calInfo) {      // 合わせ直したあとの幅は sqrt(1 − w²) 倍（交差検証で確かめた）。学習済みの化合物では、2 つの狭め方を掛け合わせず、狭いほうを使う
          var fc = Math.sqrt(1 - calInfo.w * calInfo.w);
          a68 = Math.min(a68, Math.max(s68 * fc, unc.meta.floor)); a90 = Math.min(a90, Math.max(s90 * fc, unc.meta.floor));
        }
        s68 = a68; s90 = Math.max(a90, a68);
        out.uncertainty = { seen: !!seen, s68: s68, s90: s90, tR68: [at(lkw - s68), at(lkw + s68)], tR90: [at(lkw - s90), at(lkw + s90)] };
        if (isIso) { out.uncertainty.k68 = [kAt(lkw - s68), kAt(lkw + s68)]; out.uncertainty.k90 = [kAt(lkw - s90), kAt(lkw + s90)]; }
      }
      var mw = compound.desc[0], r = model.range || {};
      out.inRange = !(r.mw && (mw < r.mw[0] || mw > r.mw[1]));
      return out;
    }

    return { predictRetention: predictRetention, embed: embed, calibrate: calibrate, clearCalibration: clearCalibration, calibrationRule: CAL,
      hasGraphModel: !!gnn, columns: columns.columns, meta: model.meta };
  }

  var api = { create: create, ionFeatures: ionFeatures, retentionTime: retentionTime, lkwFromRetentionTime: lkwFromRetentionTime, logK: logK, phiAt: phiAt, deadTime: deadTime, plates: plates, evalTrees: evalTrees };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.ChromRetention = api;
})(typeof self !== 'undefined' ? self : this);
