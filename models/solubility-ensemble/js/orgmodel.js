// SPDX-License-Identifier: MIT
// Copyright (c) 2026 技術の秘密 (gijutsu-himitsu.com)。ライセンス文はリポジトリの直下の LICENSE（MIT）。
/* 有機溶媒の見込み（log10 mol/L）を、学習の道具（このリポジトリには未同梱） が書いた weights/org_tm.json・org_notm.json で計算する。
   融点あり：物理式の残差を学ぶ木 と 溶質・溶媒の積カーネル回帰 の平均。融点なし：その平均 と 融解熱を補助学習した木 を半分ずつ。
   元の形は 開発時の試作（Python との一致を確かめた物）。ページと Node の両方で動く。MIT。 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.SolOrg = factory();
}(typeof self !== "undefined" ? self : this, function () {
  "use strict";
  const b64 = typeof atob === "function" ? atob : s => Buffer.from(s, "base64").toString("binary");
  function unpack(encoded, Type){ const raw = Uint8Array.from(b64(encoded), c => c.charCodeAt(0)); return new Type(raw.buffer); }
  function openTrees(t){ t.right = unpack(t.right, Uint16Array); t.features = unpack(t.features, Int16Array); t.cuts = unpack(t.cuts, Float64Array); t.leaves = unpack(t.leaves, Float32Array); }
  function openModel(model){
    for(const key of ["refs", "solventRefs", "weights"]) model[key] = unpack(model[key], Float32Array);
    model.bits = unpack(model.bits, Uint32Array); model.counts = unpack(model.counts, Uint16Array);
    for(const t of model.fusion ? [model.trees, model.fusion.trees, model.fusion.student] : [model.trees]) openTrees(t);
    model.solventIndex = Object.fromEntries(model.solventNames.map((s, j) => [s, j]));
    return model;
  }
  function countBits(x){ x -= (x >>> 1) & 0x55555555; x = (x & 0x33333333) + ((x >>> 2) & 0x33333333); return (((x + (x >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24; }
  // 構造の指紋（Morgan、半径 2・1024 ビット）を 32 ビットの語 32 個に
  function wordsFromRDKit(bits){ const w = new Uint32Array(32); for(let j = 0; j < bits.length; j++) if(bits[j] === "1") w[j >>> 5] |= 1 << (j & 31); return w; }
  function wordsFromOn(on){ const w = new Uint32Array(32); for(const j of on) w[j >>> 5] |= 1 << (j & 31); return w; }
  const abc = ["abr_E", "abr_S", "abr_A", "abr_B", "abr_V"], svc = ["sv_c", "sv_e", "sv_s", "sv_a", "sv_b", "sv_v"];
  function distance(values, table, row, scales){ let sq = 0; for(let j = 0; j < values.length; j++) sq += ((values[j] - table[row * values.length + j]) / scales[j]) ** 2; return sq / values.length; }
  function matern(d){ const l = Math.sqrt(3 * d); return (1 + l) * Math.exp(-l); }
  function trees(t, vector){
    const x = vector.map(Math.fround); let sum = t.init;
    for(const root of t.roots){ let n = root; while(t.features[n] >= 0) n = x[t.features[n]] <= t.cuts[n] ? n + 1 : n + t.right[n]; sum += t.lr * t.leaves[n]; }
    return sum;
  }
  // 溶質ごとに 1 回だけ計算する部分（訓練の溶質との似かた）。溶媒が変わっても同じ
  function soluteSim(model, row, words){
    const tm = model.mode === "tm", q = model.kernelCols.map(c => row[c]), nb = words.reduce((v, w) => v + countBits(w), 0), n = model.counts.length, sim = new Float64Array(n);
    for(let i = 0; i < n; i++){
      let hits = 0; for(let j = 0; j < 32; j++) hits += countBits(words[j] & model.bits[i * 32 + j]);
      const tani = hits / Math.max(1, nb + model.counts[i] - hits);
      let dist = distance(q, model.refs, i, model.scales);
      if(tm) dist += 0.28 * (row.crystal - model.refs[(i + 1) * q.length - 1]) ** 2;
      sim[i] = 0.53 * matern(dist) + 0.47 * tani;
    }
    return sim;
  }
  // row：入力の列の名前 → 値（solfeat.js の build と同じ値を小数 5 桁に丸めた物）。solvent：学習の溶媒の名前（weights/base.json の溶媒の en）。
  // graph：packing() の値（融点なしだけ使う）。sim：soluteSim() の値（省くと中で計算する）
  function predict(model, row, words, solvent, graph, sim){
    const tm = model.mode === "tm";
    const solv = svc.map(c => row[c]), w = row.sv_water, lp = row.logPow;
    const p = [w, lp, w * lp, Math.log10(row.amw), ...solv, ...abc.map((c, j) => row[c] * solv[j + 1])];
    if(tm) p.push(row.crystal);
    const offset = p.reduce((o, x, j) => o + x * model.ridge[j], model.intercept);
    const resid = trees(model.trees, [...model.treeCols.map(c => row[c]), ...p, offset]);
    const kv = model.solventNames.map((s, j) => 0.45 + 0.46 * matern(2 * distance(solv, model.solventRefs, j, model.solventScales)) + 0.09 * (s === solvent));
    const S = sim || soluteSim(model, row, words), m = kv.length;
    let kernel = 0;
    for(let i = 0; i < S.length; i++){ let wsum = 0; for(let j = 0; j < m; j++) wsum += model.weights[i * m + j] * kv[j]; kernel += S[i] * wsum; }
    const two = offset + (kernel + resid) / 2;
    if(!model.fusion) return two;
    const f = model.fusion, gv = f.graphCols.map(k => graph[k]);
    const aux = trees(f.student, [...f.studentCols.map(k => row[k]), ...gv]);
    const c = f.formulaK, formula = w ? 0.5 - lp + f.defaultCrystal : c[0] - c[1] * lp + c[2] * f.defaultCrystal + c[3] * row.transfer;
    const vector = f.cols.map(k => k === "formula0" ? +formula.toFixed(5) : row[k]);
    return (two + trees(f.trees, [...vector, ...gv, aux])) / 2;
  }
  // 重原子のグラフから、原子環境の同等さ・径・分岐・近い極性原子の対など 10 個の値（RDKit.js の get_json の分子）
  function packing(graph){
    const n = graph.atoms.length, ext = (graph.extensions || []).find(x => x.name === "rdkitRepresentation") || {}, aromatic = new Set(ext.aromaticBonds || []);
    const atoms = graph.atoms.map(a => ({ z: a.z ?? 6, h: a.impHs ?? 0, q: a.chg ?? 0 })), adj = atoms.map(() => []);
    (graph.bonds || []).forEach((b, i) => { const [a, c] = b.atoms, order = aromatic.has(i) ? 1.5 : b.bo ?? 1; adj[a].push([c, order]); adj[c].push([a, order]); });
    const classes = keys => { const ranks = new Map(); return keys.map(k => { if(!ranks.has(k)) ranks.set(k, ranks.size); return ranks.get(k); }); };
    let ranks = classes(atoms.map(a => [a.z, a.h, a.q].join(",")));
    for(let step = 0; step < n; step++){
      const next = classes(adj.map((nb, i) => ranks[i] + "/" + nb.map(([j, b]) => b + ":" + ranks[j]).sort().join(",")));
      if(next.every((r, i) => r === ranks[i])) break; ranks = next;
    }
    const counts = new Map(); ranks.forEach(r => counts.set(r, (counts.get(r) ?? 0) + 1));
    let similar = 0, crowding = 0, diameter = 0, sum = 0, radius = n;
    for(const count of counts.values()) similar += count * Math.log2(count) / n;
    const isPolar = i => [7, 8].includes(atoms[i].z);
    for(let i = 0; i < n; i++){
      const dist = Array(n).fill(-1), queue = [i]; dist[i] = 0;
      for(let k = 0; k < queue.length; k++){ const a = queue[k]; for(const [j] of adj[a]) if(dist[j] === -1){ dist[j] = dist[a] + 1; queue.push(j); } }
      const far = Math.max(...dist); diameter = Math.max(diameter, far); radius = Math.min(radius, far);
      for(let j = i + 1; j < n; j++){ sum += dist[j]; if(dist[j] >= 3 && dist[j] <= 5 && isPolar(i) && isPolar(j) && (atoms[i].h > 0 || atoms[j].h > 0)) crowding++; }
    }
    return { distinctFraction: counts.size / n, equivalentBits: similar, largestClass: Math.max(...counts.values()), diameter, radius, diameterFraction: diameter / Math.max(1, n - 1),
      averageDistance: 2 * sum / Math.max(1, n * (n - 1)), nearbyPolarPairs: crowding, aromaticFraction: (ext.aromaticAtoms || []).length / n, branchFraction: adj.filter(x => x.length >= 3).length / n };
  }
  return { openModel, wordsFromRDKit, wordsFromOn, soluteSim, predict, packing };
}));
