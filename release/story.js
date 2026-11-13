/* How it works: three animations of real samples, drawn on a canvas. The whole resource, from the records of the
   searches to the final counts by class; one cohort group close up, in eight steps from the candidate studies to the
   training and test sides of one benchmark fold (steps 2 to 5 and 8 follow the group of GSE26440); and twins and splits,
   the records of one benchmark task under a random split, leave one study out and the patient-safe split. Each step is
   a pure function of the time since it started, so that any step can be replayed or opened directly (#story/group/3;
   ?t=1200 holds that moment of the step). The data are in data/story.js (story() in build_data.py), and the positions
   before and after the correction are those of the start page (D.hero, in data.js). */
'use strict';
(() => {
  const GROUP_SCENES = [
    { label: 'Screening', dur: 4700 }, { label: 'One group', dur: 2900 }, { label: 'Repeats', dur: 4500 },
    { label: 'Labels', dur: 3600 }, { label: 'Patients', dur: 3800 }, { label: 'Sex', dur: 4000 },
    { label: 'Correction', dur: 5400 }, { label: 'Benchmark', dur: 4400 }];
  const SIM = { anim: 'whole', count: null, scene: 0, t: 0, running: false, auto: false, started: false, raf: 0, last: 0, hold: 7000, L: null, ctx: null, f: 1,
                hi: null, hit: null, near: -1, fixedT: null, epoch: -1 };
  let S = null, M = null;
  const N = D.numbers, H = D.hero;

  // ---------- timing and drawing helpers (virtual units: 880 across, or 400 on narrow screens) ----------
  const prog = (t, t0, d) => d <= 0 ? (t >= t0 ? 1 : 0) : clamp((t - t0) / d, 0, 1);
  const eio = ease, eout = p => 1 - Math.pow(1 - p, 3), lerp = (a, b, e) => a + (b - a) * e;
  const mixC = (a, b, e) => e <= 0 ? a : e >= 1 ? b : mixHex(a, b, e);
  const font = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim() || 'system-ui, sans-serif';
  const SC = { study: [], sex: {} };
  function readColors() {
    const cs = getComputedStyle($('story')), v = n => cs.getPropertyValue(n).trim();
    SC.study = [1, 2, 3, 4, 5, 6, 7, 8].map(k => v(`--sy-study-${k}`));
    ['female', 'male', 'not_recorded', 'ambiguous'].forEach(k => { SC.sex[k] = v(`--sy-sex-${k}`); });
    SC.ghost = v('--sy-ghost'); SC.sans = font('--sans'); SC.serif = font('--serif');
    SC.twinArray = v('--twin-array'); SC.twinPatient = v('--twin-patient'); SC.twinBoth = v('--twin-both'); SC.exposed = v('--twin-exposed');
  }
  const studyColor = k => SC.study[k] || C.muted;
  function tx(s, x, y, o = {}) {
    const c = SIM.ctx; if (o.a === 0) return;
    c.save(); c.globalAlpha = o.a == null ? 1 : o.a;
    c.font = `${o.w || 400} ${o.size || 12}px ${o.serif ? SC.serif : SC.sans}`;
    c.fillStyle = o.c || C.ink2; c.textAlign = o.align || 'left'; c.textBaseline = o.base || 'alphabetic';
    let str = String(s);
    if (o.max) { while (str.length > 3 && c.measureText(str).width > o.max) str = str.slice(0, -2); if (str !== String(s)) str = str.replace(/\s*$/, '') + '…'; }
    c.fillText(str, x, y); c.restore();
  }
  // a text broken into lines of at most width w, by words, three lines at most
  function wrap(s, w, size) {
    const out = [''];
    String(s).split(' ').forEach(word => { const t = out[out.length - 1] ? `${out[out.length - 1]} ${word}` : word; if (textW(t, size) <= w || !out[out.length - 1]) out[out.length - 1] = t; else out.push(word); });
    return out.length > 3 ? [...out.slice(0, 2), out.slice(2).join(' ')] : out;
  }
  function textW(s, size, w = 400) { const c = SIM.ctx; c.save(); c.font = `${w} ${size}px ${SC.sans}`; const r = c.measureText(String(s)).width; c.restore(); return r; }
  function dot(x, y, r, fill, a = 1, ring = null, rw = 1) {
    if (a <= 0 || r <= 0) return;
    const c = SIM.ctx; c.globalAlpha = Math.min(1, a); c.beginPath(); c.arc(x, y, r, 0, 2 * Math.PI);
    if (fill) { c.fillStyle = fill; c.fill(); }
    if (ring) { c.strokeStyle = ring; c.lineWidth = rw; c.stroke(); }
    c.globalAlpha = 1;
  }
  function ln(x0, y0, x1, y1, col, a = 1, w = 1, dashed = null) {
    if (a <= 0) return;
    const c = SIM.ctx; c.save(); c.globalAlpha = a; c.strokeStyle = col; c.lineWidth = w; if (dashed) c.setLineDash(dashed);
    c.beginPath(); c.moveTo(x0, y0); c.lineTo(x1, y1); c.stroke(); c.restore();
  }
  function box(x, y, w, h, col, a = 1, r = 3, stroke = null) {
    if (a <= 0 || w <= 0 || h <= 0) return;
    const c = SIM.ctx; c.save(); c.globalAlpha = a; c.beginPath();
    if (c.roundRect) c.roundRect(x, y, w, h, Math.min(r, w / 2, h / 2)); else c.rect(x, y, w, h);
    if (col) { c.fillStyle = col; c.fill(); }
    if (stroke) { c.strokeStyle = stroke; c.lineWidth = 1; c.stroke(); }
    c.restore();
  }
  // dots drawn once and stamped (for the thousands of the correction step), crisp at the canvas resolution
  let SPR = {};
  function sprite(col, r, alpha, rim) {
    const key = `${col}|${r}|${alpha}|${rim}`;
    if (SPR[key]) return SPR[key];
    const f = SIM.f, pad = rim ? 1.6 : 0.8, m = r + pad, size = Math.ceil(2 * m * f), cv = document.createElement('canvas');
    cv.width = cv.height = size;
    const c = cv.getContext('2d'); c.scale(f, f);
    c.beginPath(); c.arc(m, m, r, 0, 2 * Math.PI); c.fillStyle = rgba(col, alpha); c.fill();
    c.lineWidth = rim ? 1.1 : Math.max(0.35, r / 4);
    c.strokeStyle = rim || rgba(mixHex(col, C.dark ? '#ffffff' : '#000000', C.dark ? 0.25 : 0.3), Math.min(1, alpha + 0.15)); c.stroke();
    cv.m = m; cv.css = size / f;
    return (SPR[key] = cv);
  }
  const stampV = (sp, x, y) => SIM.ctx.drawImage(sp, x - sp.m, y - sp.m, sp.css, sp.css);
  // a hovered key picks out its category; the rest fade
  const hiA = key => SIM.hi == null || SIM.hi === String(key) ? 1 : 0.13;
  // things under the pointer: points (x, y, r) and rectangles, each with the function that writes its tooltip
  function hitPoint(x, y, r, html) { const h = SIM.hit; h.x.push(x); h.y.push(y); h.r.push(r); h.f.push(html); }
  function hitRect(x, y, w, h, html) { SIM.hit.rects.push([x, y, w, h, html]); }
  const q = s => `“${esc(s)}”`, qq = s => `“${s}”`;
  const day = s => s.unit ? `${s.unit} ${s.time}` : `time point ${s.time}`;

  // ---------- the model: every item of every step, built once from the data ----------
  function buildModel() {
    const sc = S.screen, g = S.group;
    const rows = sc.reasons.map((r, i) => ({ ...r, i, items: [] }));
    const studies = sc.studies.map((s, i) => { const d = { ...s, i, j: rows[s.r].items.length }; rows[s.r].items.push(i); return d; });
    const r = rng(574), perm = [...studies.keys()];
    for (let i = perm.length - 1; i > 0; i--) { const j = Math.floor(r() * (i + 1)); [perm[i], perm[j]] = [perm[j], perm[i]]; }
    studies.forEach((d, i) => { d.g = perm[i]; });
    // the cohort group: columns (study platforms), records and distinct samples
    const cols = g.columns.map((c, k) => ({ ...c, k, year: String(c.date).slice(0, 4), recs: [], samples: [] }));
    const samples = g.samples.map((s, i) => ({ ...s, i, recs: [] }));
    const recs = g.rec_col.map((c, i) => ({ i, col: c, acc: g.rec_acc[i], s: g.rec_sample[i] }));
    recs.forEach(x => { x.slot = cols[x.col].recs.length; cols[x.col].recs.push(x.i); samples[x.s].recs.push(x.i); });
    recs.forEach(x => { x.primary = x.slot < cols[x.col].primary; });
    samples.forEach(s => { s.slot = cols[s.col].samples.length; cols[s.col].samples.push(s.i); s.prim = s.recs.find(ri => recs[ri].primary); });
    cols.forEach(c => {
      let k = 0; c.recs.forEach(ri => { if (!recs[ri].primary) recs[ri].q = k++; });
      c.samples.slice().sort((a, b) => samples[a].cls - samples[b].cls || a - b).forEach((s, j) => { samples[s].cslot = j; });
    });
    // patients: several samples, one sample, no identifier
    const pats = g.patients.map((id, p) => ({ p, id, samples: [] }));
    samples.forEach(s => { if (s.pat >= 0) pats[s.pat].samples.push(s.i); });
    const ord = s => samples[s].order == null ? 1e9 : samples[s].order;
    pats.forEach(p => p.samples.sort((a, b) => ord(a) - ord(b) || a - b));
    const exP = g.example_patient, several = pats.filter(p => p.samples.length > 1)
      .sort((a, b) => b.samples.length - a.samples.length || (a.p === exP) - (b.p === exP) || a.samples[0] - b.samples[0]);
    const sections = [
      { key: 'several', tokens: several.map(p => ({ pat: p.p, samples: p.samples })) },
      { key: 'single', tokens: pats.filter(p => p.samples.length === 1).map(p => ({ pat: p.p, samples: p.samples })) },
      { key: 'none', tokens: samples.filter(s => s.pat < 0).map(s => ({ pat: -1, samples: [s.i] })) }].filter(x => x.tokens.length);
    let ti = 0; sections.forEach(sec => sec.tokens.forEach(tk => { tk.i = ti++; tk.samples.forEach(s => { samples[s].tok = tk; }); }));
    // sex in the study of step 6
    const sx = S.sex, sexKeys = ['female', 'male', 'not_recorded'], callKeys = ['female', 'male', 'ambiguous'];
    const jitter = rng(802);
    const sexPts = sx.x.map((x, i) => {
      const rec = sexKeys[+sx.rec[i]], call = callKeys[+sx.call[i]];
      return { i, x, y: sx.y[i], rec, call, acc: sx.acc[i], bad: rec !== 'not_recorded' && call !== 'ambiguous' && rec !== call,
               jx: 0.5 + 0.09 * gauss(jitter), jy: 0.5 + 0.09 * gauss(jitter), delay: jitter() };
    });
    // the correction: the samples of the group among those of the start page
    const want = new Set(S.correction.samples), mine = new Uint8Array(H.n);
    for (let i = 0; i < H.n; i++) if (want.has(`${H.studies[H.study[i]]}:${H.sample[i]}`)) mine[i] = 1;
    const STAGGER = 0.3, delay = (() => { const rr = rng(H.n); return H.studies.map(() => rr() * STAGGER); })();
    // the benchmark fold
    const blocks = S.bench.groups.map((b, i) => ({ ...b, i, n: b.pos + b.neg, test: b.role !== 'train' }));
    const held = S.bench.held, ghosts = held.records - held.pos - held.neg;
    return { rows, studies, cols, samples, recs, pats, sections, sexPts, mine, delay, STAGGER, blocks, ghosts,
             test: blocks.find(b => b.test), noOwn: cols.filter(c => c.primary === 0) };
  }

  // ---------- layout of every step in virtual units ----------
  function layout(cssW = $('storyStage').clientWidth || 820) {
    const n = cssW < 620;
    const L = { n, W: n ? 400 : 880, H: n ? 660 : 560, x0: n ? 14 : 24, top: n ? 68 : 84, fs: n ? 12 : 13, fsS: n ? 11 : 11.5 };
    // 1. screening: a grid of all candidates, then one row per reason
    const gc = n ? 22 : 41, gsp = n ? 14 : 12.5, grows = Math.ceil(M.studies.length / gc);
    L.grid = { x0: (L.W - (gc - 1) * gsp) / 2, y0: L.top + Math.max(10, (L.H - L.top - 30 - (grows - 1) * gsp) / 2), sp: gsp, cols: gc, r: gsp * 0.34 };
    L.labX = n ? 184 : 300; L.cntX = n ? 212 : 342; L.dotX = n ? 220 : 354;
    const dotsW = L.W - L.x0 - L.dotX, hh = n ? 18 : 22, lineH = n ? 13 : 16, groups = [...new Set(M.rows.map(r => r.group))];
    for (L.sp = n ? 7 : 10; L.sp > 3; L.sp -= 0.5) {
      L.per = Math.floor(dotsW / L.sp);
      const h = M.rows.reduce((t, r) => t + Math.max(lineH, Math.ceil(r.n / L.per) * L.sp) + (n ? 3 : 4), 0) + (groups.length - 1) * hh;
      if (L.top + h <= L.H - 12) break;
    }
    let y = L.top; L.heads = []; let lastGroup = null;
    M.rows.forEach(r => {
      if (r.group !== lastGroup) { if (r.group !== 'included') { L.heads.push([r.group, y + hh - (n ? 7 : 9)]); y += hh; } lastGroup = r.group; }
      const lines = Math.ceil(r.n / L.per); r.y = y; r.ly = y + Math.min(lines * L.sp, lineH) / 2;
      y += Math.max(lineH, lines * L.sp) + (n ? 3 : 4);
    });
    // 2-4. the cohort group as columns, in one band or two
    const bands = n ? 2 : 1, perBand = Math.ceil(M.cols.length / bands), colW = (L.W - 2 * L.x0) / perBand, chh = n ? 32 : 40;
    const bandCols = [...Array(bands)].map((_, b) => M.cols.slice(b * perBand, (b + 1) * perBand));
    for (L.csp = n ? 8 : 12; L.csp > 3; L.csp -= 0.5) {
      L.cper = Math.floor((colW - 8) / L.csp);
      const h = bandCols.reduce((t, cs) => t + chh + Math.max(...cs.map(c => Math.ceil(c.records / L.cper))) * L.csp + 14, 0);
      if (L.top + h <= L.H - (n ? 96 : 70)) break;
    }
    L.cr = L.csp * 0.38; L.chh = chh; L.colW = colW;
    y = L.top;
    bandCols.forEach(cs => {
      cs.forEach((c, k) => { c.x = L.x0 + k * colW; c.hy = y; c.dx = c.x + (colW - L.cper * L.csp) / 2; c.dy = y + chh; });
      y += chh + Math.max(...cs.map(c => Math.ceil(c.records / L.cper))) * L.csp + 14;
    });
    L.colsBottom = y;
    L.primBottom = Math.max(...M.cols.map(c => c.dy + Math.ceil(c.primary / L.cper) * L.csp));
    // 5. patients: tokens of one or more samples in rows, by section
    const secH = n ? 18 : 22;
    for (L.psp = n ? 11 : 15; L.psp > 4; L.psp -= 0.5) {
      const rowH = L.psp * 1.3, gap = L.psp * 0.75, W = L.W - 2 * L.x0;
      let yy = L.top;
      M.sections.forEach(sec => {
        sec.y = yy; yy += secH; let x = 0, rowsUsed = 1;
        sec.tokens.forEach(tk => { const w = tk.samples.length * L.psp; if (x > 0 && x + w > W) { x = 0; rowsUsed++; } tk.px = L.x0 + x; tk.row = rowsUsed - 1; x += w + gap; });
        sec.tokens.forEach(tk => { tk.py = yy + tk.row * rowH + L.psp / 2; });
        yy += rowsUsed * rowH + (n ? 6 : 10);
        if (sec.tokens.some(tk => tk.pat === S.group.example_patient)) { L.noteY = yy; yy += n ? 50 : 54; }
      });
      if (yy <= L.H - 16) break;
    }
    L.pr = L.psp * 0.36;
    // 6. the plane of sex: XIST across, the Y genes up, one scale for both
    const xs = M.sexPts.map(p => p.x), ys = M.sexPts.map(p => p.y);
    const lo = Math.min(-0.15, ...xs, ...ys) - 0.06, hi = Math.max(1.15, ...xs, ...ys) + 0.06;
    const side = n ? L.W - 2 * L.x0 - 34 : Math.min(L.H - L.top - 64, 440);
    L.sex = { x0: L.x0 + (n ? 34 : 50), y0: L.top + 8, side, lo, hi };
    L.sex.sx = lin(lo, hi, L.sex.x0, L.sex.x0 + side); L.sex.sy = lin(lo, hi, L.sex.y0 + side, L.sex.y0);
    L.sex.cardX = n ? L.x0 : L.sex.x0 + side + 46; L.sex.cardY = n ? L.sex.y0 + side + 52 : L.top + 40;
    // 7. the principal components before and after, with equal units on both axes
    const all = k => [...H[k]].map(v => v / H.scale);
    const X = [...all('x0'), ...all('x1')], Y = [...all('y0'), ...all('y1')];
    const ext = a => [Math.min(...a), Math.max(...a)], [xa, xb] = ext(X), [ya, yb] = ext(Y);
    const pw = n ? L.W - L.x0 - 34 : 560, ph = n ? 330 : L.H - L.top - 56;
    const k = Math.min(pw / (xb - xa), ph / (yb - ya));
    const px0 = (n ? 34 : 56), py0 = L.top + 10, cx = px0 + pw / 2, cy = py0 + ph / 2;
    L.pca = { k, x0: px0, y0: py0, w: pw, h: ph, sx: v => cx + (v - (xa + xb) / 2) * k, sy: v => cy - (v - (ya + yb) / 2) * k, xa, xb, ya, yb };
    L.pca.panelX = n ? L.x0 : px0 + pw + 30; L.pca.panelY = n ? py0 + ph + 40 : L.top + 34;
    // 8. the blocks of the benchmark: all together, then training and test apart
    L.bench = benchLayout(L);
    L.w = layoutWhole(L);
    L.tw = layoutTwins(L);
    return L;
  }
  function benchLayout(L) {
    const n = L.n, test = M.test, allCells = test.n + M.ghosts;
    const areas = n
      ? { all: [L.x0, L.top + 20, L.W - 2 * L.x0], train: [L.x0, L.top + 20, L.W - 2 * L.x0], test: [L.x0, 0, L.W - 2 * L.x0] }
      : { all: [L.x0, L.top + 24, L.W - 2 * L.x0], train: [L.x0, L.top + 24, 560], test: [L.x0 + 600, L.top + 24, L.W - 2 * L.x0 - 600] };
    const shelf = (blocks, [x0, y0, w], pitch) => {
      let x = 0, y = 0, rowH = 0; const gap = Math.max(3, pitch * 1.4), pos = {};
      blocks.forEach(b => {
        const bw = b.cols * pitch, bh = Math.ceil(b.cells / b.cols) * pitch;
        if (x > 0 && x + bw > w) { x = 0; y += rowH + gap; rowH = 0; }
        pos[b.i] = [x0 + x, y0 + y]; x += bw + gap; rowH = Math.max(rowH, bh);
      });
      return { pos, h: y + rowH };
    };
    let out = null;
    for (let pitch = n ? 6 : 8; pitch >= 2; pitch -= 0.5) {
      M.blocks.forEach(b => { b.cells = b.n; b.cols = b === test ? Math.ceil(Math.sqrt(allCells * 1.4)) : Math.max(1, Math.ceil(Math.sqrt(b.n))); });
      const together = shelf(M.blocks, areas.all, pitch);
      const train = shelf(M.blocks.filter(b => !b.test), areas.train, pitch);
      const testH = Math.ceil(allCells / test.cols) * pitch;
      if (n) areas.test[1] = areas.train[1] + train.h + 46;
      const bottom = Math.max(together.h + areas.all[1], n ? areas.test[1] + testH : Math.max(areas.train[1] + train.h, areas.test[1] + testH));
      const fits = test.cols * pitch <= areas.test[2] && bottom <= L.H - (n ? 56 : 64);
      out = { pitch, together: together.pos, split: { ...train.pos, [test.i]: [areas.test[0], areas.test[1]] }, areas, trainH: train.h, testH, togetherH: together.h };
      if (fits) break;
    }
    return out;
  }
  function setupCanvas() {
    readColors(); SPR = {};
    const L = SIM.L = layout(), cv = $('storyCanvas'), dpr = window.devicePixelRatio || 1, cssW = $('storyStage').clientWidth || 820;
    cv.width = Math.round(cssW * dpr); cv.height = Math.round(cssW * L.H / L.W * dpr);
    cv.style.height = `${cssW * L.H / L.W}px`;
    SIM.f = cv.width / L.W;
    SIM.ctx = cv.getContext('2d'); SIM.ctx.setTransform(SIM.f, 0, 0, SIM.f, 0, 0);
    SIM.epoch = themeEpoch;
  }

  // ---------- one frame ----------
  function simDraw() {
    if (SIM.dry) return;                                   // the self-test of playback runs without drawing
    if (!SIM.ctx || SIM.epoch !== themeEpoch) setupCanvas();
    const c = SIM.ctx, L = SIM.L, s = SIM.scene, t = SIM.started ? SIM.t : 1e9;
    c.save(); c.setTransform(1, 0, 0, 1, 0, 0); c.fillStyle = C.raised; c.fillRect(0, 0, c.canvas.width, c.canvas.height); c.restore();
    SIM.hit = { x: [], y: [], r: [], f: [], rects: [] }; SIM.count = null;
    A().draw[s](t);
    const [title, sub, short] = A().header(s, t), cnt = SIM.count, big = cnt && (typeof cnt[0] === 'string' ? cnt[0] : fmt(cnt[0]));
    let room = L.W - 2 * L.x0;
    if (cnt && !L.n) {
      tx(big, L.W - L.x0, 36, { size: 26, w: 600, serif: true, align: 'right', c: C.ink });
      tx(cnt[1], L.W - L.x0, 54, { size: L.fsS, align: 'right', c: C.muted });
      room -= Math.max(textW(cnt[1], L.fsS), textW(big, 26, 600)) + 24;
    }
    tx(L.n && short ? short : title, L.x0, L.n ? 26 : 34, { size: L.n ? 16 : 22, w: 600, serif: true, c: C.ink, max: room });
    if (cnt && L.n) tx(`${big} ${cnt[1]}`, L.x0, 45, { size: L.fs, w: 600, c: C.ink2, max: room });
    else tx(sub, L.x0, L.n ? 45 : 58, { size: L.fs, c: C.muted, max: room });
    const h = SIM.hit;
    if (SIM.near >= 0 && SIM.near < h.x.length) dot(h.x[SIM.near], h.y[SIM.near], h.r[SIM.near] + 3, null, 1, C.ink, 1.5);
    if (!SIM.started) { c.save(); c.setTransform(1, 0, 0, 1, 0, 0); c.fillStyle = C.raised; c.globalAlpha = 0.55; c.fillRect(0, 0, c.canvas.width, c.canvas.height); c.restore(); }
  }
  function headerText(s, t) {
    const g = S.group, st = g.stated, sx = S.sex, B = S.bench, merged = t >= GROUP_SCENES[2].dur;
    return [
      [`${fmt(S.screen.candidates)} candidate studies`, `one dot per study · ${fmt(N.studies)} included after screening`],
      [`Cohort group ${g.id}: ${fmt(st.records)} sample records`, `in ${spell(M.cols.length)} studies; one dot per record, colored by the study it was deposited in`],
      [`${fmt(st.records)} records, ${fmt(st.distinct)} distinct samples`, merged ? 'each array once, in the column of its first deposit' : 'records of one array join its first deposit'],
      [`${fmt(st.distinct)} samples, each with a class`, 'one dot per distinct sample, colored by the class of its record'],
      [`${fmt(st.patients)} patients and ${fmt(st.without_patient)} samples without one`, 'samples of one patient joined; a hollow dot is not the marked sample'],
      [`Sex read from the genes in ${sx.study}`, `${fmt(M.sexPts.length)} arrays, colored by the sex in their records`],
      [`${fmt(H.n)} bulk-blood samples, ${t < 1200 + 0.5 * 3600 ? 'before' : 'after'} the correction`, `first two principal components; ringed: the ${fmt(S.correction.samples.length)} samples of ${g.id}`,
        `${fmt(H.n)} samples, ${t < 1200 + 0.5 * 3600 ? 'before' : 'after'} correction`],
      [`${taskName()}: ${B.fold} held out`, 'one square per sample of the task, one block per patient-safe group', `${B.fold} held out as the test set`]][s];
  }
  const taskName = () => `${lab('class', S.bench.positive)} against ${labLower('class', S.bench.negative)}s`;

  // step 1: candidates in a grid, then in rows by the reason that excluded them
  function drawScreen(t) {
    const L = SIM.L, G0 = L.grid, rows = M.rows;
    const rowCol = r => r.group === 'included' ? C.accent : r.group === 'access' ? C.muted : mixHex(C.muted, C.raised, 0.4);
    L.heads.forEach(([g, y], k) => {
      if (g === 'included') return;
      tx(`Excluded: ${labLowerText(S.screen.groups.find(x => x[0] === g)[1])}`, L.x0, y, { size: L.fsS, w: 600, c: C.muted, a: eout(prog(t, 1700 + k * 120, 500)) });
    });
    rows.forEach(r => {
      const a = eout(prog(t, 1800 + r.i * 90, 500)) * hiA(r.group === 'included' ? 'included' : 'excluded');
      tx(r.label, L.labX, r.ly + 4, { size: L.fs, align: 'right', c: r.group === 'included' ? C.accentInk : C.ink2, w: r.group === 'included' ? 600 : 400, a, max: L.labX - L.x0 - 4 });
      tx(fmt(r.n), L.cntX, r.ly + 4, { size: L.fs, align: 'right', c: C.ink, w: 600, a });
    });
    const gridA = 1 - eout(prog(t, 4200, 400));
    tx(`${fmt(S.screen.census_records)} records found, ${fmt(S.screen.candidates)} candidate studies after removing copies`, L.W / 2,
      G0.y0 + (Math.ceil(M.studies.length / G0.cols) - 1) * G0.sp + 30, { size: L.fs, align: 'center', c: C.muted, a: eout(prog(t, 300, 400)) * (1 - prog(t, 1300, 300)) * gridA });
    const neutral = C.muted;
    M.studies.forEach(d => {
      const r = rows[d.r], gx = G0.x0 + (d.g % G0.cols) * G0.sp, gy = G0.y0 + Math.floor(d.g / G0.cols) * G0.sp;
      const rx = L.dotX + (d.j % L.per + 0.5) * L.sp, ry = r.y + (Math.floor(d.j / L.per) + 0.5) * L.sp;
      const pop = prog(t, 80 + d.g * 1.5, 380), t0 = 1500 + r.i * 85 + d.j * 3, e = eio(prog(t, t0, 900));
      const x = lerp(gx, rx, e), y = lerp(gy, ry, e) - Math.sin(Math.PI * e) * 18, rad = lerp(G0.r, L.sp * 0.36, e) * (pop < 1 ? eout(pop) : 1);
      const col = mixC(neutral, rowCol(r), prog(t, t0 + 300, 500));
      const a = (pop < 1 ? pop : 1) * (e > 0 ? 1 : 0.7) * hiA(r.group === 'included' ? 'included' : 'excluded');
      dot(x, y, rad, col, a);
      if (e >= 1) hitPoint(x, y, rad, () => screenTip(d));
    });
  }
  function screenTip(d) {
    const r = M.rows[d.r], s = S.screen;
    const assay = { array: 'array', 'rna-seq': 'RNA-seq', 'array+rna-seq': 'array and RNA-seq' }[d.k] || d.k;
    return `<b>${esc(d.a)}</b> <span class="m">· ${esc(s.sources[d.s])}${d.y ? ` · ${d.y}` : ''}${assay ? ` · ${esc(assay)}` : ''}</span><br>${esc(d.t)}<br>` +
      (r.group === 'included' ? '<b>Included</b>' : `<b>Excluded:</b> ${esc(labLowerText(r.label))}<br><span class="m">${esc(r.rule)}</span>`);
  }
  const labLowerText = s => /^[A-Z][a-z]/.test(s) ? s[0].toLowerCase() + s.slice(1) : s;

  // steps 2 to 4: the columns of the cohort group
  const slotXY = (c, k) => { const L = SIM.L; return [c.dx + (k % L.cper + 0.5) * L.csp, c.dy + (Math.floor(k / L.cper) + 0.5) * L.csp]; };
  function drawColumnHeads(line2, a = 1) {
    const L = SIM.L;
    M.cols.forEach(c => {
      const ca = a * hiA(c.k);
      box(c.x + 4, c.hy + 1, 3, L.n ? 11 : 13, studyColor(c.k), ca, 1.5);
      tx(c.study, c.x + 11, c.hy + (L.n ? 11 : 13), { size: L.n ? 10.5 : 13, w: 600, c: C.ink, a: ca, max: L.colW - 14 });
      tx(line2(c), c.x + 11, c.hy + (L.n ? 24 : 29), { size: L.fsS, c: C.muted, a: ca, max: L.colW - 14 });
    });
  }
  const colYear = c => SIM.L.n ? '' : `${c.year} · `;
  function recordTip(x) {
    const r = M.recs[x], s = M.samples[r.s], pc = M.cols[s.col], c = M.cols[r.col];
    return `<b>${esc(c.study)} · ${esc(r.acc)}</b><br>` + (r.primary
      ? `First deposit of its array (primary record)${s.nrec > 1 ? `; ${spell(s.nrec)} records in all` : ''}`
      : `Repeat deposit: the array was first deposited as ${esc(pc.study)}:${esc(s.acc)}${pc.date ? `, released ${fdate(pc.date)}` : ''}`);
  }
  function sampleTip(i, withPatient = false) {
    const s = M.samples[i], c = M.cols[s.col], cls = D.order.class[s.cls];
    const others = s.recs.filter(ri => !M.recs[ri].primary).map(ri => M.cols[M.recs[ri].col].study);
    let html = `<b>${esc(c.study)}:${esc(s.acc)}</b><br>${sw(C.cls[cls])}${lab('class', cls)} <span class="m">· ${s.field && s.field !== '(no key)' ? `${esc(s.field)}: ` : ''}${q(s.said)}</span>`;
    if (others.length) html += `<br><span class="m">Also deposited in ${andList([...new Set(others)].map(esc))}${s.differs ? '; its records differ in class' : ''}</span>`;
    if (withPatient) {
      html += s.pat >= 0 ? `<br>Patient ${esc(M.pats[s.pat].id)}${s.time ? `, ${esc(day(s))}` : ''} · ${s.marked ? '<b>marked</b> as the sample of this patient' : 'not the marked sample'}`
        : '<br>No patient identifier: marked as its own patient';
      if (s.flag) html += `<br><span class="m">Quality flag: ${esc(s.flag.replace(/_/g, ' '))}</span>`;
    }
    return html;
  }
  function drawGroup(t) {
    const L = SIM.L;
    drawColumnHeads(c => `${colYear(c)}${fmt(c.records)} records`, 1);
    M.cols.forEach(c => c.recs.forEach((ri, k) => {
      const e = prog(t, 150 + c.k * 220 + k * 3, 380); if (e <= 0) return;
      const [x, y] = slotXY(c, k), r = L.cr * eout(e) * (1 + 0.3 * Math.sin(Math.PI * e));
      dot(x, y, r, studyColor(c.k), Math.min(1, e * 2) * hiA(c.k));
      if (e >= 1) hitPoint(x, y, L.cr, () => recordTip(ri));
    }));
  }
  // step 3: each repeat deposit flies to the primary record of its array and merges with it
  function drawMerge(t) {
    const L = SIM.L;
    const start = r => 300 + M.recs[r].col * 300 + M.recs[r].q * 3.6, D_ = 950;
    const left = c => c.primary + c.recs.filter(ri => !M.recs[ri].primary && prog(t, start(ri), D_) < 1).length;
    const done = t >= GROUP_SCENES[2].dur;
    drawColumnHeads(c => `${colYear(c)}${fmt(left(c))} ${done ? (left(c) === 1 ? 'sample' : 'samples') : 'records'}`, 1);
    M.cols.forEach(c => c.recs.forEach((ri, k) => {
      const r = M.recs[ri]; if (!r.primary) return;
      const [x, y] = slotXY(c, k);
      dot(x, y, L.cr, studyColor(c.k), hiA(c.k));
      hitPoint(x, y, L.cr, () => sampleTip(r.s));
    }));
    M.cols.forEach(c => c.recs.forEach((ri, k) => {
      const r = M.recs[ri]; if (r.primary) return;
      const e = prog(t, start(ri), D_);
      if (e > 0) { const [gx, gy] = slotXY(c, k); dot(gx, gy, L.cr * 0.8, null, 0.35 * prog(e, 0, 0.3) * hiA(c.k), studyColor(c.k), 0.8); }
      if (e >= 1) return;
      const s = M.samples[r.s], [x0, y0] = slotXY(c, k), [x1, y1] = slotXY(M.cols[s.col], s.slot), ee = eio(e);
      const x = lerp(x0, x1, ee), y = lerp(y0, y1, ee) - Math.sin(Math.PI * ee) * (L.n ? 10 : 18);
      dot(x, y, L.cr * (1 - 0.35 * prog(e, 0.75, 0.25)), studyColor(c.k), (1 - prog(e, 0.8, 0.2)) * hiA(c.k), e > 0 ? C.raised : null, 0.8);
    }));
    // the count, falling from records to distinct samples
    const now = S.group.stated.distinct + M.recs.filter(r => !r.primary && prog(t, start(r.i), D_) < 1).length;
    const big = L.n ? 26 : 34, bx = L.W - L.x0, by = L.H - (L.n ? 30 : 30);
    tx(fmt(now), bx, by, { size: big, w: 600, serif: true, align: 'right', c: C.ink });
    tx(now === S.group.stated.distinct ? 'distinct samples' : 'records', bx, by + (L.n ? 16 : 18), { size: L.fsS, align: 'right', c: C.muted });
  }
  // step 4: every sample takes the color of its class, then the classes gather within each column
  function drawLabels(t) {
    const L = SIM.L, g = S.group;
    drawColumnHeads(c => `${colYear(c)}${fmt(c.primary)} ${c.primary === 1 ? 'sample' : 'samples'}`, 1);
    M.samples.forEach(s => {
      const c = M.cols[s.col], e = prog(t, 100 + c.k * 150 + s.slot * 2.4, 450), m = eio(prog(t, 1750 + c.k * 70, 900));
      const [x0, y0] = slotXY(c, s.slot), [x1, y1] = slotXY(c, s.cslot), x = lerp(x0, x1, m), y = lerp(y0, y1, m);
      const cls = D.order.class[s.cls];
      dot(x, y, L.cr, mixC(studyColor(c.k), C.cls[cls], e), e < 1 ? 1 : hiA(cls));
      hitPoint(x, y, L.cr, () => sampleTip(s.i));
      s.lx = x; s.ly = y;
    });
    // two worked examples, with the words of their records
    const a = eout(prog(t, 2700, 500)); if (a <= 0) return;
    const ex = [M.samples[g.example_sepsis], M.samples[g.example_twin]];
    const twinHere = S.group.twin_records.find(r => M.cols[r.col].study === g.study && !r.primary);
    const notes = [
      [`${M.cols[ex[0].col].study}:${ex[0].acc}`, `${ex[0].field}: ${qq(ex[0].said)} → ${labLower('class', D.order.class[ex[0].cls])}`],
      [`${M.cols[ex[1].col].study}:${ex[1].acc}`, `${qq(ex[1].said)} → ${labLower('class', D.order.class[ex[1].cls])}` + (twinHere ? `; in ${g.study}: ${qq(twinHere.said)} → ${labLower('class', D.order.class[twinHere.cls])}` : '')]];
    callouts(ex.map((s, k) => ({ x: s.lx, y: s.ly, head: notes[k][0], body: notes[k][1] })), a);
  }
  // boxes with a leader line to a dot, in the free space under the columns
  function callouts(items, a) {
    const L = SIM.L, size = L.fsS, pad = 7;
    items.forEach((it, k) => {
      const head = it.head, body = it.body;
      const w = Math.min(L.W - 2 * L.x0, Math.max(textW(head, size, 600), textW(body, size)) + 2 * pad), h = 2 * size + 3 * pad - 2;
      const bx = L.n ? L.x0 : clamp(it.x - w / 2 + (k ? 60 : -60), L.x0, L.W - L.x0 - w), by = L.n ? L.H - 12 - (items.length - k) * (h + 8) : L.primBottom + 26 + k * (h + 12);
      ln(it.x, it.y, clamp(it.x, bx + 8, bx + w - 8), by, C.ink, a, 1);
      box(bx, by, w, h, C.bg, a * 0.95, 6, C.lineStrong);
      tx(head, bx + pad, by + pad + size - 1, { size, w: 600, c: C.ink, a, max: w - 2 * pad });
      tx(body, bx + pad, by + 2 * pad + 2 * size - 3, { size, c: C.ink2, a, max: w - 2 * pad });
      dot(it.x, it.y, L.cr + 2.5, null, a, C.ink, 1.5);
    });
  }
  // step 5: samples of one patient joined, one of them marked
  function drawPatients(t) {
    const L = SIM.L, hd = 1 - prog(t, 0, 400), g = S.group;
    if (hd > 0) drawColumnHeads(c => `${colYear(c)}${fmt(c.primary)} samples`, hd);
    const names = { several: n => `Patients with two or more samples · ${fmt(n)}`, single: n => `Patients with one sample · ${fmt(n)}`, none: n => `No patient identifier, each its own patient · ${fmt(n)}` };
    const ha = eout(prog(t, 1600, 500));
    M.sections.forEach(sec => tx(names[sec.key](sec.tokens.length), L.x0, sec.y + (L.n ? 12 : 15), { size: L.fsS, w: 600, c: C.muted, a: ha }));
    const pos = s => { const tk = s.tok, k = tk.samples.indexOf(s.i); return [tk.px + (k + 0.5) * L.psp, tk.py]; };
    const linkA = eout(prog(t, 2100, 500)), hollow = eio(prog(t, 2600, 600));
    // links first, under the dots
    if (linkA > 0) M.sections.forEach(sec => sec.tokens.forEach(tk => {
      if (tk.samples.length < 2) return;
      const [xa, ya] = pos(M.samples[tk.samples[0]]), [xb] = pos(M.samples[tk.samples[tk.samples.length - 1]]);
      ln(xa, ya, xb, ya, C.ink2, linkA * hiA(D.order.class[M.samples[tk.samples[0]].cls]) * 0.9, L.n ? 1.5 : 2);
    }));
    M.samples.forEach(s => {
      const c = M.cols[s.col], e = eio(prog(t, 100 + s.tok.i * 4, 1100)), [x0, y0] = slotXY(c, s.cslot), [x1, y1] = pos(s);
      const x = lerp(x0, x1, e), y = lerp(y0, y1, e), r = lerp(L.cr, L.pr, e), cls = D.order.class[s.cls], col = C.cls[cls], a = hiA(cls);
      if (!s.marked && hollow > 0) { dot(x, y, r, mixC(col, C.raised, hollow), a, col, 1.4); }
      else dot(x, y, r, col, a);
      hitPoint(x, y, r, () => sampleTip(s.i, true));
      s.px = x; s.py = y;
    });
    // the patient of the worked example
    const a = eout(prog(t, 3000, 500)), p = M.pats[g.example_patient];
    if (a > 0 && p) {
      const first = M.samples[p.samples[0]], last = M.samples[p.samples[p.samples.length - 1]];
      const head = `Patient ${p.id}`, body = p.samples.map(i => { const s = M.samples[i]; return `${day(s)}: ${labLower('class', D.order.class[s.cls])}${s.marked ? ' (marked)' : ''}`; }).join(' · ');
      const size = L.fsS, pad = 7, w = Math.max(textW(head, size, 600), textW(body, size)) + 2 * pad, h = 2 * size + 3 * pad - 2;
      const bx = clamp((first.px + last.px) / 2 - w / 2, L.x0, L.W - L.x0 - w), by = L.noteY + 2;
      ln((first.px + last.px) / 2, first.py + L.pr + 2, clamp((first.px + last.px) / 2, bx + 8, bx + w - 8), by, C.ink, a, 1);
      box(bx, by, w, h, C.bg, a * 0.95, 6, C.lineStrong);
      tx(head, bx + pad, by + pad + size - 1, { size, w: 600, c: C.ink, a });
      tx(body, bx + pad, by + 2 * pad + 2 * size - 3, { size, c: C.ink2, a, max: w - 2 * pad });
      box(first.px - L.pr - 3, first.py - L.pr - 3, last.px - first.px + 2 * L.pr + 6, 2 * L.pr + 6, null, a, L.pr + 3, C.ink);
    }
  }
  // step 6: the arrays of one study placed by XIST and the Y-chromosome genes
  function drawSex(t) {
    const L = SIM.L, P = L.sex, sx = S.sex, aAx = eout(prog(t, 0, 500)), n = L.n;
    const x0 = P.sx(P.lo), x1 = P.sx(P.hi), y0 = P.sy(P.hi), y1 = P.sy(P.lo);
    box(x0, y0, x1 - x0, y1 - y0, null, aAx, 0, C.line);
    ticks(P.lo, P.hi, n ? 4 : 6).forEach(v => {
      ln(P.sx(v), y0, P.sx(v), y1, v === 0 || v === 1 ? C.lineStrong : C.grid, aAx, 1, v === 0 || v === 1 ? [3, 3] : null);
      ln(x0, P.sy(v), x1, P.sy(v), v === 0 || v === 1 ? C.lineStrong : C.grid, aAx, 1, v === 0 || v === 1 ? [3, 3] : null);
      tx(fx(v, 1), P.sx(v), y1 + 14, { size: L.fsS - 0.5, align: 'center', c: C.muted, a: aAx });
      tx(fx(v, 1), x0 - 5, P.sy(v) + 4, { size: L.fsS - 0.5, align: 'right', c: C.muted, a: aAx });
    });
    tx(`${sx.female_gene}: 0 = average of the arrays read as male, 1 = female`, (x0 + x1) / 2, y1 + (n ? 30 : 34), { size: L.fsS, align: 'center', c: C.ink2, a: aAx, max: x1 - x0 + 30 });
    const c = SIM.ctx; c.save(); c.translate(x0 - (n ? 26 : 34), (y0 + y1) / 2); c.rotate(-Math.PI / 2);
    tx('Y-chromosome genes: 0 = female average, 1 = male', 0, 0, { size: L.fsS, align: 'center', c: C.ink2, a: aAx, max: y1 - y0 }); c.restore();
    // the groups: where the arrays read as female and as male gather
    const la = eout(prog(t, 2600, 500));
    tx('read as female', P.sx(1), P.sy(0.3), { size: L.fsS, align: 'center', c: C.ink2, w: 600, a: la });
    tx('read as male', P.sx(0.16), P.sy(1.25), { size: L.fsS, c: C.ink2, w: 600, a: la });
    const r = n ? 2.3 : 2.9, ring = eout(prog(t, 3100, 450));
    const order = M.sexPts.slice().sort((a, b) => (a.bad || a.call === 'ambiguous') - (b.bad || b.call === 'ambiguous'));
    order.forEach(p => {
      const e = eio(prog(t, 300 + p.delay * 1500, 1100)), x = lerp(P.sx(p.jx), P.sx(p.x), e), y = lerp(P.sy(p.jy), P.sy(p.y), e);
      const a = (0.25 + 0.75 * Math.min(1, prog(t, 0, 300) * 2)) * hiA(p.rec);
      const special = p.bad ? C.ink : p.call === 'ambiguous' ? SC.sex.ambiguous : null;
      dot(x, y, special ? r + 0.6 * ring : r, rgba(SC.sex[p.rec], special ? 0.95 : 0.6), a, special && ring > 0 ? special : null, special ? 1.6 * ring : 1);
      hitPoint(x, y, r + 1, () => sexTip(p));
      p.vx = x; p.vy = y;
    });
    // the two worked examples of the guide
    const ca = eout(prog(t, 3300, 500)); if (ca <= 0) return;
    const pad = 8, size = L.fsS, cw = n ? (L.W - 2 * L.x0 - 10) / 2 : L.W - L.x0 - P.cardX, ch = n ? 3 * size + 4 * pad : 4 * size + 5 * pad;
    sx.examples.forEach((ex, k) => {
      const p = M.sexPts.find(d => d.acc === ex.acc); if (!p) return;
      const bx = n ? L.x0 + k * (cw + 10) : P.cardX, by = n ? P.cardY : P.cardY + k * (ch + 16);
      ln(p.vx, p.vy, n ? bx + cw / 2 : bx, n ? by : by + ch / 2, C.ink, ca * 0.8, 1);
      box(bx, by, cw, ch, C.bg, ca * 0.95, 6, C.lineStrong);
      dot(p.vx, p.vy, r + 3, null, ca, C.ink, 1.5);
      const lines = [[`${ex.acc} · recorded ${ex.sex}`, 600, C.ink], [`Y genes ${fx(ex.y, 2)} − ${sx.female_gene} ${fx(ex.xist, 2)} = ${fmt(ex.d, 1)}`, 400, C.ink2],
        [`read as ${ex.sex}`, 600, SC.sex[ex.sex]]];
      if (!n) lines.splice(1, 0, [`mean of the Y genes minus ${sx.female_gene}, on the log₂ scale`, 400, C.muted]);
      lines.forEach(([s, w, col], j) => tx(s, bx + pad, by + pad + size + j * (size + pad) - 1, { size, w, c: col, a: ca, max: cw - 2 * pad }));
    });
  }
  function sexTip(p) {
    const sx = S.sex;
    return `<b>${esc(sx.study)} · ${esc(p.acc)}</b><br>${sw(SC.sex[p.rec])}Recorded: ${p.rec === 'not_recorded' ? 'not recorded' : p.rec}<br>Read from the genes: ${p.call}` +
      (p.bad ? '<br><b>Disagrees with its record</b>' : p.call === 'ambiguous' ? '<br><span class="m">Too close to the middle to call</span>' : '');
  }
  // step 7: every sample of the bulk-blood table, from before the correction to after it
  function drawCorrection(t) {
    const L = SIM.L, P = L.pca, n = L.n, fade = eout(prog(t, 0, 600)), p = prog(t, 1200, 3600);
    const x0 = P.sx(P.xa), x1 = P.sx(P.xb), y0 = P.sy(P.yb), y1 = P.sy(P.ya), state = p < 0.5 ? 'before' : 'after';
    ticks(P.xa, P.xb, n ? 4 : 6).forEach(v => { ln(P.sx(v), y0, P.sx(v), y1, v === 0 ? C.lineStrong : C.grid, fade); tx(fmt(v), P.sx(v), y1 + 14, { size: L.fsS - 0.5, align: 'center', c: C.muted, a: fade }); });
    ticks(P.ya, P.yb, n ? 4 : 5).forEach(v => { ln(x0, P.sy(v), x1, P.sy(v), v === 0 ? C.lineStrong : C.grid, fade); tx(fmt(v), x0 - 5, P.sy(v) + 4, { size: L.fsS - 0.5, align: 'right', c: C.muted, a: fade }); });
    const share = H.variance[state];
    tx(`Component 1 (${pct(share[0], 0)} of the variation)`, (x0 + x1) / 2, y1 + (n ? 30 : 34), { size: L.fsS, align: 'center', c: C.ink2, a: fade });
    const c = SIM.ctx; c.save(); c.translate(x0 - (n ? 24 : 34), (y0 + y1) / 2); c.rotate(-Math.PI / 2);
    tx(`Component 2 (${pct(share[1], 0)})`, 0, 0, { size: L.fsS, align: 'center', c: C.ink2, a: fade }); c.restore();
    tx(state === 'before' ? 'Before the correction' : 'After the correction', x0 + 8, y0 + (n ? 16 : 20), { size: n ? 13 : 16, w: 600, serif: true, c: C.ink, a: fade });
    const fams = D.order.family, r = n ? 1.4 : 1.8, R = n ? 2.4 : 3.1;
    const sprites = fams.map(f => sprite(C.fam[f], r, 0.55, null)), mineSp = fams.map(f => sprite(C.fam[f], R, 0.95, C.ink));
    const at = i => { const e = ease(clamp((p - M.delay[H.study[i]]) / (1 - M.STAGGER), 0, 1)); return [lerp(H.x0[i], H.x1[i], e) / H.scale, lerp(H.y0[i], H.y1[i], e) / H.scale]; };
    c.globalAlpha = fade;
    for (let i = 0; i < H.n; i++) {
      if (M.mine[i]) continue;
      const f = +H.fam[i], a = hiA(fams[f]); if (a < 1) c.globalAlpha = fade * a;
      const [x, y] = at(i), px = P.sx(x), py = P.sy(y); stampV(sprites[f], px, py); if (a < 1) c.globalAlpha = fade;
      hitPoint(px, py, r, () => heroTip(i));
    }
    for (let i = 0; i < H.n; i++) {
      if (!M.mine[i]) continue;
      const [x, y] = at(i), px = P.sx(x), py = P.sy(y), f = +H.fam[i]; c.globalAlpha = fade * (SIM.hi == null || SIM.hi === 'mine' || SIM.hi === fams[f] ? 1 : 0.13);
      stampV(mineSp[f], px, py); hitPoint(px, py, R, () => heroTip(i));
    }
    c.globalAlpha = 1;
    // what happened to the samples of the group, on the right (below on narrow screens)
    const g = S.group, rows = M.corr, pa = eout(prog(t, 4300, 600)), size = L.fsS;
    let y = P.panelY;
    tx(`The ${fmt(S.correction.samples.length)} samples of ${g.id}`, P.panelX, y, { size: L.fs, w: 600, c: C.ink, a: pa }); y += size + 10;
    rows.forEach(([method, studies, k]) => {
      tx(`${fmt(k)} ${method === 'healthy_anchored' ? 'on healthy controls' : 'by matched groups'}`, P.panelX, y, { size, c: C.ink2, w: 600, a: pa, max: L.W - L.x0 - P.panelX });
      y += size + 4; tx(andList(studies), P.panelX, y, { size, c: C.muted, a: pa, max: L.W - L.x0 - P.panelX }); y += size + 10;
    });
  }
  function heroTip(i) {
    const c = D.order.class[+H.cls[i]], f = D.order.family[+H.fam[i]];
    return `<b>${esc(H.studies[H.study[i]])}:${esc(H.sample[i])}</b>${M.mine[i] ? ` <span class="m">· ${esc(S.group.id)}</span>` : ''}<br>${sw(C.fam[f])}${lab('family', f)} &nbsp;${sw(C.cls[c])}${lab('class', c)}`;
  }
  // step 8: patient-safe groups as blocks, one of them held out as the test set
  function drawBench(t) {
    const L = SIM.L, B = L.bench, n = L.n, P = B.pitch, cell = P >= 3 ? P - 1 : P * 0.8, test = M.test, bench = S.bench;
    const neg = C.cls[bench.negative], pos = C.cls[bench.positive];
    const appear = b => eout(prog(t, b.i * 12, 500)), move = b => eio(prog(t, b.test ? 900 : 1150 + b.i * 10, 1100));
    const la = eout(prog(t, 2300, 500)), ga = prog(t, 2500, 900);
    // the two sides
    const [trX, trY, trW] = B.areas.train, [teX, teY, teW] = B.areas.test;
    tx(`Training · ${fmt(M.blocks.length - 1)} groups`, trX, trY - 8, { size: L.fs, w: 600, c: C.ink, a: la });
    tx(`Test · ${test.g}`, teX, teY - 8, { size: L.fs, w: 600, c: C.ink, a: la });
    const sumTrain = k => sum(M.blocks.filter(b => !b.test).map(b => b[k]));
    tx(`${fmt(sumTrain('pos'))} ${labLower('class', bench.positive)} · ${fmt(sumTrain('neg'))} ${labLower('class', bench.negative)}s`, trX, trY + B.trainH + 18, { size: L.fsS, c: C.ink2, a: la });
    tx(`${fmt(test.pos)} ${labLower('class', bench.positive)} · ${fmt(test.neg)} ${labLower('class', bench.negative)}s`, teX, teY + B.testH + 18, { size: L.fsS, c: C.ink2, a: la });
    if (ga > 0) {
      tx(`and the other ${fmt(M.ghosts)} records of ${test.g},`, teX, teY + B.testH + 18 + L.fsS + 6, { size: L.fsS, c: C.muted, a: la, max: teW });
      tx('which are not samples of the task', teX, teY + B.testH + 18 + 2 * (L.fsS + 6), { size: L.fsS, c: C.muted, a: la, max: teW });
    }
    const c = SIM.ctx;
    M.blocks.forEach(b => {
      const a = appear(b); if (a <= 0) return;
      const e = move(b), [ax, ay] = B.together[b.i], [bx, by] = B.split[b.i], x = lerp(ax, bx, e), y = lerp(ay, by, e) - (b.test ? Math.sin(Math.PI * e) * 26 : 0);
      const total = b.test ? b.n + Math.round(M.ghosts * ga) : b.n, rows = Math.ceil(total / b.cols);
      c.globalAlpha = a;
      for (let k = 0; k < total; k++) {
        const cx = x + (k % b.cols) * P, cy = y + Math.floor(k / b.cols) * P, col = k < b.pos ? pos : k < b.n ? neg : SC.ghost;
        const hk = k < b.pos ? bench.positive : k < b.n ? bench.negative : 'ghost';
        c.globalAlpha = a * hiA(hk); c.fillStyle = col; c.fillRect(cx, cy, cell, cell);
      }
      c.globalAlpha = 1;
      if (b.test) box(x - 3, y - 3, b.cols * P + 5, rows * P + 5, null, a, 3, C.ink);
      hitRect(x - 1, y - 1, b.cols * P + 1, rows * P + 1, () => blockTip(b));
    });
  }
  function blockTip(b) {
    const bench = S.bench, st = b.studies, list = st.length > 6 ? `${st.slice(0, 6).join(', ')} and ${fmt(st.length - 6)} more` : andList(st);
    return `<b>${esc(b.g)}</b> <span class="m">· ${b.tech === 'rna_seq' ? 'RNA-seq' : 'arrays'} · ${b.test ? 'test' : 'training'}</span><br>${esc(list)}<br>` +
      `${sw(C.cls[bench.positive])}${fmt(b.pos)} ${labLower('class', bench.positive)} &nbsp;${sw(C.cls[bench.negative])}${fmt(b.neg)} ${labLower('class', bench.negative)}s` +
      (b.test ? `<br><span class="m">All ${fmt(bench.held.records)} records of the group are on this side</span>` : '');
  }

  // =====================================================================================================================
  // The whole resource: every record of the searches, then every sample record of the included studies, down to the
  // classes of the distinct samples with expression values (the flow of Figure 1). One dot per record or sample.
  const WHOLE_SCENES = [
    { label: 'Searches', dur: 4900 }, { label: 'Records', dur: 2700 }, { label: 'No values', dur: 3300 },
    { label: 'Repeats', dur: 4000 }, { label: 'Classes', dur: 4600 }, { label: 'Patients', dur: 3000 }];
  function buildWhole() {
    const w = S.whole, n = w.code.length, units = w.units.map((u, i) => ({ ...u, i }));
    const unit = new Uint16Array(n); let k = 0;
    units.forEach(u => { u.start = k; for (let j = 0; j < u.n; j++) unit[k++] = u.i; });
    const cls = new Uint8Array(n), prim = new Uint8Array(n), tab = new Uint8Array(n), mark = new Uint8Array(n), why = new Uint8Array(n);
    w.code.forEach((c, i) => { cls[i] = c & 7; prim[i] = c >> 3 & 1; tab[i] = c >> 4 & 1; mark[i] = c >> 5 & 1; why[i] = c >> 6; });
    const target = new Int32Array(n).fill(-1), repRank = new Int32Array(n).fill(-1);
    for (let p = 0, r = 0; p < w.rep.length; p += 2, r++) { target[w.rep[p]] = w.rep[p + 1]; repRank[w.rep[p]] = r; }
    const fam = new Uint8Array(n); for (let i = 0; i < n; i++) fam[i] = units[unit[i]].fam;
    // the records without values, in the order of their reason; the distinct samples with values, by class, marked ones first
    const gone = [], byClass = D.order.class.map(() => []);
    for (let i = 0; i < n; i++) if (!tab[i]) gone.push(i);
    gone.sort((a, b) => why[a] - why[b] || a - b);
    const goneRank = new Int32Array(n).fill(-1), reasonRank = new Int32Array(n).fill(-1), perReason = w.no_values.map(() => 0);
    gone.forEach((i, r) => { goneRank[i] = r; reasonRank[i] = perReason[why[i] - 1]++; });
    for (let i = 0; i < n; i++) if (prim[i] && tab[i]) byClass[cls[i]].push(i);
    const classRank = new Int32Array(n).fill(-1);
    byClass.forEach(list => { list.sort((a, b) => mark[b] - mark[a] || a - b); list.forEach((i, r) => { classRank[i] = r; }); });
    // the records of the searches: one dot each, in the order of the flow, on a shuffled grid
    const cats = [...w.not_candidate.map(([key, label, cnt]) => ({ col: 0, key, label, n: cnt })),
      ...w.excluded.map(([key, label, cnt]) => ({ col: 1, key, label, n: cnt })),
      { col: 2, key: 'included', label: 'Included', n: w.flow.included_studies }].filter(c => c.n > 0);
    const nS = sum(cats.map(c => c.n)), scat = new Uint8Array(nS), sj = new Int32Array(nS), scol = new Int32Array(nS);
    const colN = [0, 0, 0];
    let d = 0;
    cats.forEach((c, ci) => { c.start = colN[c.col]; for (let j = 0; j < c.n; j++, d++) { scat[d] = ci; sj[d] = j; scol[d] = colN[c.col] + j; } colN[c.col] += c.n; });
    const rr = rng(5454), perm = Int32Array.from({ length: nS }, (_, i) => i);
    for (let i = nS - 1; i > 0; i--) { const j = Math.floor(rr() * (i + 1)); const x = perm[i]; perm[i] = perm[j]; perm[j] = x; }
    // per technology: all records (n), the records with values (rec) and the distinct samples with values (tab)
    const families = D.order.family.map((f, fi) => ({ key: f, fi, n: 0, rec: 0, tab: 0 }));
    for (let i = 0; i < n; i++) { const f = families[fam[i]]; f.n++; if (tab[i]) { f.rec++; if (prim[i]) f.tab++; } }
    const classes = D.order.class.map((c, ci) => ({ key: c, ci, n: byClass[ci].length, marked: byClass[ci].filter(i => mark[i]).length }));
    return { n, units, unit, cls, prim, tab, mark, why, target, repRank, fam, gone, goneRank, reasonRank, byClass, classRank,
             cats, nS, scat, sj, scol, colN, perm, families, classes, reps: w.rep.length / 2 };
  }
  // stacked bands, one per group, filled row by row; labels on the left, or above each band on narrow screens
  function bands(L, groups, y0, bottom, opt) {
    const X0 = L.x0, Wd = L.W - 2 * X0, n = L.n, labelW = n ? 0 : opt.labelW, labelH = n ? 16 : 0, gap = n ? 8 : opt.gap;
    for (let sp = opt.sp; sp >= 2; sp -= 0.1) {
      const cols = Math.floor((Wd - labelW) / sp); let y = y0; const out = [];
      groups.forEach(g => {
        y += labelH; const rows = Math.ceil(g.n / cols), h = Math.max(rows * sp, n ? 0 : opt.minH);
        out.push({ ...g, x: X0 + labelW, y, cols, rows, h, ly: n ? y - 4 : y }); y += h + gap;
      });
      if (y - gap <= bottom) return { sp, cols, bands: out, bottom: y - gap };
    }
    throw new Error('the bands of the whole resource do not fit the drawing area');
  }
  function layoutWhole(L) {
    const W = M.W, n = L.n, X0 = L.x0, Wd = L.W - 2 * X0, P = {};
    // 1. the records of the searches: a grid, then blocks by their fate
    const gsp = n ? 3.9 : 4.6, gcols = Math.floor(Wd / gsp), grows = Math.ceil(W.nS / gcols);
    P.gsp = gsp; P.gcols = gcols; P.gx = X0 + (Wd - (gcols - 1) * gsp) / 2; P.gy = L.top + 26 + Math.max(0, (L.H - L.top - 90 - grows * gsp) / 2);
    const colsX = n ? [[X0, Wd], [X0, Wd], [X0, Wd]] : [[X0, 390], [X0 + 414, 268], [X0 + 706, Wd - 706]];
    const isp = n ? 7 : 8, hd = n ? 16 : 18, lab = n ? 14 : 16, gap = n ? 8 : 12;
    let y = L.top + 4;
    P.heads = [];
    W.cats.forEach((c, ci) => {
      const [x, w] = colsX[c.col], sp = c.key === 'included' ? isp : gsp, per = Math.floor(w / sp);
      if (ci === 0 || W.cats[ci - 1].col !== c.col) { if (!n) y = L.top + 4; else if (ci) y += 6; P.heads.push([c.col, x, y + hd - 5]); y += hd; }
      c.x = x; c.w = w; c.ly = y + lab - 4; c.bx = x; c.by = y + lab; c.sp = sp; c.per = per; c.rows = Math.ceil(c.n / per);
      y = c.by + c.rows * sp + gap;
    });
    P.bottom1 = y;
    // 2-4. the sample records in bands by technology, with room below for the records without values
    const fams = W.families.filter(f => f.n > 0);
    const nv = S.whole.no_values.filter(r => r[2] > 0), asideH = 28 + nv.length * (n ? 15 : 16);
    P.fam = bands(L, fams, L.top + 4, L.H - asideH - (n ? 14 : 18), { labelW: 168, sp: 5.4, minH: 30, gap: 8 });
    P.r2 = P.fam.sp * 0.4;
    P.aside = { y: P.fam.bottom + (n ? 22 : 26), labX: n ? 214 : 280, cntX: n ? 240 : 312, dotX: n ? 248 : 322, sp: n ? 5 : 6, rowH: n ? 15 : 16 };
    P.aside.rows = nv.map(([key, label, cnt], k) => ({ key, label, n: cnt, y: P.aside.y + 18 + k * P.aside.rowH }));
    // 5-6. the samples with values in bands by class
    const cl = W.classes.filter(c => c.n > 0);
    P.cls = bands(L, cl, L.top + 4, L.H - (n ? 12 : 14), { labelW: 176, sp: 5.6, minH: 42, gap: 10 });
    P.r5 = P.cls.sp * 0.4;
    // where every record sits in each arrangement
    const N_ = W.n, x2 = new Float32Array(N_), y2 = new Float32Array(N_), x5 = new Float32Array(N_), y5 = new Float32Array(N_);
    const famBand = {}; P.fam.bands.forEach(b => { famBand[b.fi] = b; });
    const seen = new Int32Array(D.order.family.length);
    for (let i = 0; i < N_; i++) {
      const b = famBand[W.fam[i]], k = seen[W.fam[i]]++, sp = P.fam.sp;
      x2[i] = b.x + (k % b.cols + 0.5) * sp; y2[i] = b.y + (Math.floor(k / b.cols) + 0.5) * sp;
    }
    const clsBand = {}; P.cls.bands.forEach(b => { clsBand[b.ci] = b; });
    for (let i = 0; i < N_; i++) {
      if (W.classRank[i] < 0) continue;
      const b = clsBand[W.cls[i]], k = W.classRank[i], sp = P.cls.sp;
      x5[i] = b.x + (k % b.cols + 0.5) * sp; y5[i] = b.y + (Math.floor(k / b.cols) + 0.5) * sp;
    }
    Object.assign(P, { x2, y2, x5, y5 });
    return P;
  }
  const WHOLE_COLS = ['Not candidate studies', 'Excluded at screening', 'Included'];
  const unitName = u => `${u.study}${u.platform ? ` · ${u.platform}` : ''}`;
  function wholeTip(i, kind) {
    const W = M.W, u = W.units[W.unit[i]], f = D.order.family[u.fam], c = D.order.class[W.cls[i]];
    let html = `<b>${esc(unitName(u))}</b> <span class="m">· ${lab('family', f)} · ${fmt(u.n)} ${u.n === 1 ? 'record' : 'records'}</span>`;
    if (kind === 'record') html += W.tab[i] ? `<br>${W.prim[i] ? 'First deposit of its array or library' : 'Repeat deposit of a sample first deposited in another study'}`
      : `<br>No expression values: ${esc(labLowerText(S.whole.no_values[W.why[i] - 1][1]))}`;
    if (kind === 'gone') html += `<br>No expression values: ${esc(labLowerText(S.whole.no_values[W.why[i] - 1][1]))}`;
    if (kind === 'class') html += `<br>${sw(C.cls[c])}${lab('class', c)}${W.mark[i] ? '' : ' <span class="m">· not the marked sample of its patient</span>'}`;
    return html;
  }
  function counter(value, label) { SIM.count = [value, label]; }
  // step 1: 5,454 records found, sorted by what became of them
  function drawW1(t) {
    const L = SIM.L, P = L.w, W = M.W, n = L.n, starts = [700, 2500, 3600], durs = [1500, 800, 700];
    const pale = mixHex(C.muted, C.raised, 0.55), colors = [pale, C.muted, C.accent];
    P.heads.forEach(([col, x, y]) => tx(`${WHOLE_COLS[col]} · ${fmt(W.colN[col])}`, x, y, { size: L.fs, w: 600, c: col === 2 ? C.accentInk : C.ink, a: eout(prog(t, starts[col], 400)) * hiA(['nc', 'ex', 'in'][col]) }));
    W.cats.forEach(c => {
      const a = eout(prog(t, starts[c.col] + 200, 400)) * hiA(['nc', 'ex', 'in'][c.col]);
      if (c.key !== 'included') tx(`${c.label} · ${fmt(c.n)}`, c.x, c.ly, { size: L.fsS, c: C.ink2, a, max: c.w });
      const rows = Math.ceil(c.n / c.per);
      hitRect(c.bx, c.ly - 10, c.w, rows * c.sp + 12, () => `<b>${esc(c.label)}</b> · ${fmt(c.n)}<br><span class="m">${esc(WHOLE_COLS[c.col])}</span>`);
    });
    const spr = [sprite(C.muted, P.gsp * 0.36, 0.55, null), sprite(pale, P.gsp * 0.36, 0.9, null), sprite(C.muted, P.gsp * 0.36, 0.9, null)];
    const done = [0, 0, 0], c2 = SIM.ctx;
    for (let d = 0; d < W.nS; d++) {
      const c = W.cats[W.scat[d]], col = c.col, g = W.perm[d];
      const pop = prog(t, (g / W.nS) * 700, 300);
      const t0 = starts[col] + (W.scol[d] / W.colN[col]) * durs[col], e = eio(prog(t, t0, 650));
      if (e >= 1) done[col]++;
      const gx = P.gx + (g % P.gcols) * P.gsp, gy = P.gy + Math.floor(g / P.gcols) * P.gsp, j = W.sj[d];
      const bx = c.bx + (j % c.per + 0.5) * c.sp, by = c.by + (Math.floor(j / c.per) + 0.5) * c.sp;
      const x = lerp(gx, bx, e), y = lerp(gy, by, e) - Math.sin(Math.PI * e) * 14;
      const a = (pop < 1 ? pop : 1) * hiA(['nc', 'ex', 'in'][col]);
      if (a <= 0) continue;
      c2.globalAlpha = a;
      if (col === 2 && e > 0) { c2.globalAlpha = 1; dot(x, y, lerp(P.gsp * 0.36, c.sp * 0.38, e), colors[2], a); }
      else stampV(e > 0 ? spr[col + 1] || spr[0] : spr[0], x, y);
    }
    c2.globalAlpha = 1;
    const left = W.nS - done[0] - done[1];
    counter(left, done[0] < W.colN[0] ? 'search records' : done[1] < W.colN[1] ? 'candidate studies' : 'included studies');
  }
  // the sample records in their technology bands, with the label and count of each band
  function famLabels(t, a, count, unitWord) {
    const L = SIM.L, P = L.w;
    P.fam.bands.forEach(b => {
      const k = count(b), ba = a * hiA(b.key);
      if (L.n) { tx(`${lab('family', b.key)} · ${fmt(k)} ${unitWord(k)}`, b.x, b.ly, { size: L.fsS, w: 600, c: C.ink2, a: ba }); return; }
      tx(lab('family', b.key), L.x0, b.y + 11, { size: L.fs, w: 600, c: C.ink, a: ba, max: 160 });
      tx(`${fmt(k)} ${unitWord(k)}`, L.x0, b.y + 27, { size: L.fsS, c: C.muted, a: ba, max: 160 });
    });
  }
  const famSprites = r => D.order.family.map(f => [sprite(C.fam[f], r, 0.95, null), sprite(mixHex(C.fam[f], C.raised, 0.35), r, 0.95, null)]);
  function drawW2(t) {
    const L = SIM.L, P = L.w, W = M.W, c = SIM.ctx, spr = famSprites(P.r2);
    famLabels(t, 1, b => b.n, k => k === 1 ? 'record' : 'records');
    let shown = 0;
    for (let i = 0; i < W.n; i++) {
      const e = prog(t, (i / W.n) * 1900, 350); if (e <= 0) continue;
      shown++;
      const f = W.fam[i]; c.globalAlpha = Math.min(1, e * 1.5) * hiA(D.order.family[f]);
      stampV(spr[f][W.unit[i] % 2], P.x2[i], P.y2[i]);
      hitPoint(P.x2[i], P.y2[i], P.r2, () => wholeTip(i, 'record'));
    }
    c.globalAlpha = 1;
    counter(shown, 'sample records');
  }
  // the records without values in their rows below the bands, as step 3 leaves them
  function asideStill(a) {
    const L = SIM.L, P = L.w, A_ = P.aside, W = M.W;
    W.gone.forEach(i => {
      const row = A_.rows[W.why[i] - 1], x = A_.dotX + (W.reasonRank[i] + 0.5) * A_.sp, y = row.y - 4;
      dot(x, y, A_.sp * 0.38, C.fam[D.order.family[W.fam[i]]], a);
      if (a >= 1) hitPoint(x, y, Math.max(P.r2, 2.5), () => wholeTip(i, 'gone'));
    });
    tx(`No expression values · ${fmt(W.gone.length)}`, L.x0, A_.y, { size: L.fs, w: 600, c: C.ink, a });
    A_.rows.forEach(row => {
      tx(row.label, A_.labX, row.y, { size: L.fsS, align: 'right', c: C.ink2, a, max: A_.labX - L.x0 });
      tx(fmt(row.n), A_.cntX, row.y, { size: L.fsS, align: 'right', c: C.ink, w: 600, a });
    });
  }
  // step 4: each repeat deposit with values joins the first deposit of its array or library
  function drawW3(t) {
    const L = SIM.L, P = L.w, W = M.W, c = SIM.ctx, spr = famSprites(P.r2), D_ = 800, span = 2800;
    const start = i => 300 + (W.repRank[i] / Math.max(1, W.reps)) * span;
    const left = new Int32Array(D.order.family.length); let arrived = 0;
    for (let i = 0; i < W.n; i++) {
      if (!W.tab[i]) continue;
      const f = W.fam[i];
      if (W.prim[i]) { left[f]++; c.globalAlpha = hiA(D.order.family[f]); stampV(spr[f][W.unit[i] % 2], P.x2[i], P.y2[i]); hitPoint(P.x2[i], P.y2[i], P.r2, () => wholeTip(i, 'record')); continue; }
      const e = prog(t, start(i), D_);
      if (e <= 0) { left[f]++; c.globalAlpha = hiA(D.order.family[f]); stampV(spr[f][W.unit[i] % 2], P.x2[i], P.y2[i]); continue; }
      if (e >= 1) { arrived++; continue; }
      const j = W.target[i], ee = eio(e), x = lerp(P.x2[i], P.x2[j], ee), y = lerp(P.y2[i], P.y2[j], ee) - Math.sin(Math.PI * ee) * (L.n ? 10 : 16);
      c.globalAlpha = (1 - prog(e, 0.8, 0.2)) * hiA(D.order.family[f]); stampV(spr[f][0], x, y);
    }
    c.globalAlpha = 1;
    asideStill(1);
    const done = arrived === W.reps;
    famLabels(t, 1, b => left[b.fi], k => done ? (k === 1 ? 'sample' : 'samples') : (k === 1 ? 'record' : 'records'));
    counter(W.n - W.gone.length - arrived, done ? 'distinct samples' : 'records with values');
  }
  // step 3: the records without expression values step aside, one row per reason
  function drawW4(t) {
    const L = SIM.L, P = L.w, A_ = P.aside, W = M.W, c = SIM.ctx, spr = famSprites(P.r2);
    famLabels(t, 1, b => b.rec, k => k === 1 ? 'record with values' : 'records with values');
    for (let i = 0; i < W.n; i++) {
      if (!W.tab[i]) continue;
      const f = W.fam[i]; c.globalAlpha = hiA(D.order.family[f]) * (SIM.hi === 'gone' ? 0.13 : 1);
      stampV(spr[f][W.unit[i] % 2], P.x2[i], P.y2[i]); hitPoint(P.x2[i], P.y2[i], P.r2, () => wholeTip(i, 'record'));
    }
    c.globalAlpha = 1;
    const ha = eout(prog(t, 700, 400)), cnt = A_.rows.map(() => 0);
    tx(`No expression values · ${fmt(W.gone.length)}`, L.x0, A_.y, { size: L.fs, w: 600, c: C.ink, a: ha });
    let arrived = 0;
    W.gone.forEach(i => {
      const r = W.goneRank[i], e = eio(prog(t, 900 + r * 14, 800)), row = A_.rows[W.why[i] - 1], k = W.reasonRank[i];
      const ax = A_.dotX + (k + 0.5) * A_.sp, ay = row.y - 4, x = lerp(P.x2[i], ax, e), y = lerp(P.y2[i], ay, e);
      const ring = 1 - prog(t, 2600, 500), pulse = 1 + 0.5 * Math.sin(Math.min(1, t / 800) * Math.PI * 3) ** 2;
      if (e >= 1) { arrived++; cnt[W.why[i] - 1]++; }
      const col = C.fam[D.order.family[W.fam[i]]];
      dot(x, y, lerp(P.r2, A_.sp * 0.38, e), col, SIM.hi == null || SIM.hi === 'gone' ? 1 : 0.13, ring > 0 ? C.ink : null, 1.2 * ring * pulse);
      hitPoint(x, y, Math.max(P.r2, 2.5), () => wholeTip(i, 'gone'));
    });
    A_.rows.forEach((row, k) => {
      tx(row.label, A_.labX, row.y, { size: L.fsS, align: 'right', c: C.ink2, a: ha, max: A_.labX - L.x0 });
      tx(fmt(cnt[k]), A_.cntX, row.y, { size: L.fsS, align: 'right', c: C.ink, w: 600, a: ha });
    });
    counter(W.n - arrived, arrived === W.gone.length ? 'sample records with values' : 'sample records');
  }
  // the class bands: name, count and, in the last step, the marked samples
  function classLabels(a, counts, marked) {
    const L = SIM.L, P = L.w;
    P.cls.bands.forEach(b => {
      const k = counts[b.ci], ba = a * hiA(b.key);
      if (L.n) { tx(`${lab('class', b.key)} · ${fmt(k)}${marked ? ` · ${fmt(marked[b.ci])} marked` : ''}`, b.x, b.ly, { size: L.fsS, w: 600, c: C.ink, a: ba }); return; }
      tx(lab('class', b.key), L.x0, b.y + 12, { size: L.fs, w: 600, c: C.ink2, a: ba, max: 168 });
      tx(fmt(k), L.x0, b.y + 34, { size: 22, w: 600, serif: true, c: C.ink, a: ba });
      if (marked) tx(`${fmt(marked[b.ci])} marked`, L.x0 + 84, b.y + 33, { size: L.fsS, c: C.muted, a: ba, max: 90 });
    });
  }
  // step 5: the samples with values move into their classes
  function drawW5(t) {
    const L = SIM.L, P = L.w, W = M.W, c = SIM.ctx, spr = famSprites(P.r2), out = 1 - prog(t, 0, 500);
    const cspr = D.order.class.map(k => sprite(C.cls[k], P.r5, 0.95, null));
    if (out > 0) {
      famLabels(t, out, b => b.tab, k => k === 1 ? 'sample with values' : 'samples with values');
      asideStill(out);
    }
    const arrived = D.order.class.map(() => 0), ncls = W.classes.length;
    for (let i = 0; i < W.n; i++) {
      if (W.classRank[i] < 0) continue;
      const k = W.cls[i], list = W.byClass[k], t0 = 400 + k * (1400 / ncls) + (W.classRank[i] / list.length) * 1800, e = eio(prog(t, t0, 900));
      if (e >= 1) arrived[k]++;
      const x = lerp(P.x2[i], P.x5[i], e), y = lerp(P.y2[i], P.y5[i], e), key = D.order.class[k];
      c.globalAlpha = hiA(key);
      stampV(e < 0.5 ? spr[W.fam[i]][W.unit[i] % 2] : cspr[k], x, y);
      hitPoint(x, y, P.r5, () => wholeTip(i, e < 0.5 ? 'record' : 'class'));
    }
    c.globalAlpha = 1;
    classLabels(eout(prog(t, 500, 500)), arrived, null);
    counter(sum(arrived), 'distinct samples with values');
  }
  // step 6: one sample for each patient; the further samples of a patient fade
  function drawW6(t) {
    const L = SIM.L, P = L.w, W = M.W, c = SIM.ctx, cspr = D.order.class.map(k => sprite(C.cls[k], P.r5, 0.95, null));
    const fade = D.order.class.map((_, k) => eio(prog(t, 300 + k * 260, 1000)));
    let faded = 0;
    for (let i = 0; i < W.n; i++) {
      if (W.classRank[i] < 0) continue;
      const k = W.cls[i], f = W.mark[i] ? 0 : fade[k];
      if (!W.mark[i]) faded += f;
      c.globalAlpha = (1 - 0.84 * f) * hiA(D.order.class[k]);
      stampV(cspr[k], P.x5[i], P.y5[i]); hitPoint(P.x5[i], P.y5[i], P.r5, () => wholeTip(i, 'class'));
    }
    c.globalAlpha = 1;
    const ma = eout(prog(t, 400, 500));
    classLabels(1, W.classes.map(x => x.n), ma > 0 ? W.classes.map(x => x.marked) : null);
    counter(Math.round(W.classes.reduce((s_, x) => s_ + x.n, 0) - faded), faded > 0.5 ? 'marked samples, one per patient' : 'distinct samples with values');
  }
  function wholeHeader(s, t) {
    const F = S.whole.flow;
    return [
      [`From ${fmt(F.records_found)} search records to ${fmt(F.included_studies)} studies`, 'one dot per record found, then per study', `${fmt(F.records_found)} records, ${fmt(F.included_studies)} studies`],
      [`${fmt(F.included_sample_records)} sample records of the ${fmt(F.included_studies)} studies`, `one dot per record, in bands by technology; shades alternate between the ${fmt(F.included_study_platforms)} study platforms`, `${fmt(F.included_sample_records)} sample records`],
      ['Records without expression values', 'listed apart with their reason; the resource is the records with values', 'Records without values'],
      ['Repeat deposits join their first deposit', 'records of one array or library become one distinct sample', 'Repeat deposits merge'],
      [`The ${fmt(F.distinct_samples)} distinct samples, by class`, 'one dot per distinct sample with expression values', 'The final counts by class'],
      ['One sample marked for each patient', 'faded: further samples of a patient', 'One sample per patient']][s];
  }
  function wholeCaption(s) {
    const w = S.whole, F = w.flow, W = M.W, b = v => `<b>${v}</b>`;
    const nc = w.not_candidate, ncTotal = sum(nc.map(r => r[2])), ex = w.excluded, exTotal = sum(ex.map(r => r[2]));
    const src = N.studies_by_source, fams = W.families.filter(f => f.n > 0).sort((a, c) => c.n - a.n);
    const cls = W.classes.filter(c => c.n > 0);
    const nv = w.no_values.filter(r => r[2] > 0).sort((a, c) => c[2] - a[2]);
    const marked = sum(W.classes.map(c => c.marked));
    const plural = { sepsis: 'sepsis samples', healthy_control: 'healthy controls', non_septic_control: 'non-septic controls', other_condition: 'samples of other conditions', not_recorded: 'samples whose class is not recorded' };
    const list = (get, k) => `${fmt(get(k))} ${plural[k] || labLower('class', k)}`;
    return [
      ['From the searches to the included studies',
        `Searches of public repositories, mainly GEO, ArrayExpress and the European Nucleotide Archive, returned ${fmt(F.records_found)} records, one dot each. ${fmt(ncTotal)} of them were not <i>candidate studies</i>: ${andList(nc.map(r => `${fmt(r[2])} ${labLowerText(r[1])}`))}. Screening checked the ${fmt(F.candidate_studies)} candidates against the eligibility rules and excluded ${fmt(exTotal)}, each for the first rule it failed. ${fmt(F.included_studies)} studies were included.`,
        `${fmt(F.records_found)} − ${fmt(ncTotal)} = ${b(fmt(F.candidate_studies))} candidate studies<br>${fmt(F.candidate_studies)} − ${fmt(exTotal)} = ${b(fmt(F.included_studies))} included studies (${andList(Object.entries(src).sort((a, c) => c[1] - a[1]).map(([k, v]) => `${fmt(v)} from ${{ geo: 'GEO', arrayexpress: 'ArrayExpress', ena: 'ENA' }[k] || k}`))})<br>Excluded: ${ex.map(r => `${labLowerText(r[1])} ${fmt(r[2])}`).join('; ')}`],
      ['The sample records of the included studies',
        `The ${fmt(F.included_studies)} included studies hold ${fmt(F.included_sample_records)} sample records, one dot each, on ${fmt(F.included_study_platforms)} <i>study platforms</i> (the part of a study measured on one platform). Each band is one technology, and within a band the shade alternates from one study platform to the next.`,
        fams.map(f => `${lab('family', f.key)} ${b(fmt(f.n))}`).join(' · ') + `<br>${fmt(F.sample_records_with_values)} of the ${fmt(F.included_sample_records)} records have expression values`],
      ['Records without expression values',
        `${fmt(W.gone.length)} records have no expression values in any table: ${andList(nv.map(r => `${labLowerText(r[1])} (${fmt(r[2])})`))}. They are not in the sample table; they are listed apart, each with its reason and its recorded fields. The resource is the ${fmt(F.sample_records_with_values)} records with values, on ${fmt(F.study_platforms_with_values)} study platforms.`,
        `${fmt(F.included_sample_records)} − ${fmt(W.gone.length)} = ${b(fmt(F.sample_records_with_values))} sample records with expression values`],
      ['Repeat deposits join their first deposit',
        `Some arrays and libraries were deposited in more than one study. Their records were linked by what the raw files record about the physical array or library, or by identical counts, and each repeat deposit joins the first deposit of its array or library, its <i>primary record</i>. The sample is then counted once and keeps all its records. The second animation, one group close up, shows this step for seven studies.`,
        `${fmt(F.sample_records_with_values)} records − ${fmt(W.reps)} repeat deposits = ${b(fmt(F.distinct_samples))} distinct samples<br>${fmt(N.repeated_samples)} samples have one or more repeat deposits`],
      ['The final counts by class',
        `Each sample takes the class of its primary record, read from the words of the record: sepsis only where the record or a sentence of its study says sepsis, septic shock or septic, and a healthy control only where it is called healthy or normal. The ${fmt(F.distinct_samples)} distinct samples are ${andList(cls.map(c => list(k => F[`class_${k}`], c.key)))}.`,
        cls.map(c => `${lab('class', c.key)} ${b(fmt(F[`class_${c.key}`]))}`).join(' · ') + `<br>${cls.map(c => fmt(F[`class_${c.key}`])).join(' + ')} = ${b(fmt(F.distinct_samples))}`],
      ['One sample for each patient',
        `Some patients gave several samples, on several days or on two kinds of array. One sample is marked for each patient, and a sample without a patient identifier counts as its own patient. Of the ${fmt(F.distinct_samples)} distinct samples, ${fmt(F.one_sample_per_patient)} are marked: one for each of ${fmt(F.identified_patients)} identified patients, and ${fmt(F.samples_without_identifier)} samples without an identifier. Since some of those may come from one person, ${fmt(F.one_sample_per_patient)} is an upper bound on the number of patients.`,
        cls.map(c => `${lab('class', c.key)} ${b(fmt(c.marked))}`).join(' · ') + `<br>${fmt(F.identified_patients)} + ${fmt(F.samples_without_identifier)} = ${b(fmt(marked))} marked samples`]][s];
  }
  const WHOLE_HOWTO = [
    'Each small dot is one record that the searches returned. The blocks show what became of them: not a candidate study, excluded at screening, or included. Hover over a block to see its reason.',
    'Each dot is one sample record. The bands are technologies, and the shades within a band alternate between study platforms. Hover over a dot to see its study.',
    'The ringed dots are the records without expression values; they move to the rows below, one row per reason, and leave the count.',
    'Dots that leave their place are repeat deposits; they join the first deposit of their array, in whatever band it lies. The counts of the bands and the number at the top fall as they merge.',
    'Each band is one class, and each dot one distinct sample with expression values. The numbers count the dots as they arrive.',
    'Solid dots are the samples marked as one per patient; faded dots are further samples of a patient.'];
  function wholeKeys() {
    const s = SIM.scene, W = M.W, el = $('storyKeys'), hov = k => { if (SIM.hi !== k) { SIM.hi = k; if (!SIM.running) simDraw(); } };
    let items = [], extra = '';
    if (s === 0) items = [{ key: 'nc', label: WHOLE_COLS[0], color: mixHex(C.muted, C.raised, 0.55), count: W.colN[0] }, { key: 'ex', label: WHOLE_COLS[1], color: C.muted, count: W.colN[1] },
      { key: 'in', label: WHOLE_COLS[2], color: C.accent, count: W.colN[2] }];
    if (s >= 1 && s <= 3) items = W.families.filter(f => f.n > 0).map(f => ({ key: f.key, label: lab('family', f.key), color: C.fam[f.key], count: s === 1 ? f.n : s === 2 ? f.rec : f.tab }));
    if (s === 2) extra = keys([[`no expression values (${fmt(W.gone.length)})`, C.ink, 'ring']]);
    if (s >= 4) items = W.classes.filter(c => c.n > 0).map(c => ({ key: c.key, label: lab('class', c.key), color: C.cls[c.key], count: s === 4 ? c.n : c.marked }));
    if (s === 5) extra = keys([['marked sample', C.ink2, 'dot'], ['further sample of a patient (faded)', mixHex(C.ink2, C.raised, 0.7), 'dot']]);
    legend(el, items, hov);
    if (extra) el.insertAdjacentHTML('beforeend', extra);
    $('storyHow').lastElementChild.textContent = WHOLE_HOWTO[s];
  }

  // =====================================================================================================================
  // Twins and splits: every record of one benchmark task as a naive analysis pools it, in blocks by patient-safe group
  // and, inside a block, by study; then one real random split, leave one study out and the folds of leave one group out
  // (the leakage audit of the benchmark stage). Every count on screen is recomputed here from the twins of the records.
  const TWIN_SCENES = [{ label: 'Twins', dur: 3800 }, { label: 'Random split', dur: 5600 }, { label: 'One study out', dur: 6400 },
    { label: 'Patient-safe', dur: 7600 }];
  function buildTwins() {
    const K = S.leak, n = K.rec_band.length, st = K.studies;
    const band = Int16Array.from(K.rec_band), study = Int16Array.from(K.rec_band, b => K.bands[b].s), group = Int16Array.from(K.rec_band, b => K.bands[b].g);
    const arr = Int32Array.from(K.rec_array), pat = Int32Array.from(K.rec_patient), elig = new Uint8Array(n), test = new Uint8Array(n);
    for (let i = 0; i < n; i++) { const f = +K.rec_flag[i]; elig[i] = f & 1; test[i] = f >> 1 & 1; }
    const members = codes => { const m = new Map(); codes.forEach((c, i) => { if (c >= 0) { if (!m.has(c)) m.set(c, []); m.get(c).push(i); } }); return m; };
    const arrays = members(arr), patients = members(pat);
    const akey = i => arr[i] >= 0 ? arr[i] : -1 - i;           // the physical array of a record
    // the twins of every record: 1 another record of its array, 2 another array of its patient, 3 both
    const kind = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      if (arr[i] >= 0) kind[i] |= 1;
      if (pat[i] >= 0 && patients.get(pat[i]).some(j => akey(j) !== akey(i))) kind[i] |= 2;
    }
    // the test records of a split that have a twin in training, by the rule of the audit, and one such twin for the thread
    const mark = new Uint8Array(n), twin = new Int32Array(n).fill(-1);
    function expose(isTest, isTrain, into = null) {
      const c = { test: 0, exposed: 0, array: 0, patient: 0 };
      for (let i = 0; i < n; i++) {
        if (!isTest(i)) continue;
        c.test++;
        let m = 0, tw = -1;
        if (arr[i] >= 0) for (const j of arrays.get(arr[i])) if (j !== i && isTrain(j)) { m |= 1; tw = j; break; }
        if (pat[i] >= 0) for (const j of patients.get(pat[i])) if (akey(j) !== akey(i) && isTrain(j)) { m |= 2; if (tw < 0) tw = j; break; }
        if (into) { into.mark[i] = m; into.twin[i] = tw; }
        if (m) c.exposed++; if (m & 1) c.array++; if (m & 2) c.patient++;
      }
      return c;
    }
    const random = { mark: new Uint8Array(n), twin: new Int32Array(n).fill(-1) };
    random.c = expose(i => test[i] === 1, i => test[i] === 0, random);
    // leave one study out: every study once, in the order of the blocks; each record is tested once
    const loso = { mark, twin, folds: st.map((_, s) => expose(i => study[i] === s, i => study[i] !== s, { mark, twin })) };
    loso.c = ['test', 'exposed', 'array', 'patient'].reduce((o, k) => ({ ...o, [k]: sum(loso.folds.map(f => f[k])) }), {});
    // the folds of leave one group out: the samples of the task in the held-out groups against every record of the training groups
    const folds = K.folds.map(f => { const out = new Set(f.out), tr = new Set(f.train); return { ...f, c: expose(i => elig[i] === 1 && out.has(group[i]), i => tr.has(group[i])) }; });
    // the records of every study in every group (a band of its block) and the bands of every group, in the order of the blocks
    const bands = K.bands.map(b => ({ s: b.s, g: b.g, start: -1, n: 0 }));
    for (let i = 0; i < n; i++) { const b = bands[band[i]]; if (b.start < 0) b.start = i; b.n++; }
    const groups = K.groups.map((name, g) => ({ g, name, bands: bands.filter(b => b.g === g), n: 0 }));
    groups.forEach(G => { G.n = sum(G.bands.map(b => b.n)); });
    const testList = []; for (let i = 0; i < n; i++) if (test[i]) testList.push(i);
    const lit = testList.filter(i => random.mark[i]);
    const rank = new Int32Array(n).fill(-1); testList.forEach((i, k) => { rank[i] = k; });
    const litRank = new Int32Array(n).fill(-1); lit.forEach((i, k) => { litRank[i] = k; });
    return { n, band, study, group, arr, pat, elig, test, arrays, patients, akey, kind, random, loso, folds, bands, groups, testList, rank, litRank, lit: lit.length,
             samples: sum([...elig]), exArray: arrays.get(K.example_array) || [], exPatient: patients.get(K.example_patient) || [] };
  }
  // the blocks of the groups in the training area, the test area beside them (below on narrow screens), and every place
  function layoutTwins(L) {
    const T = M.T, n = L.n, area = { x: L.x0, y: L.top + (n ? 22 : 28), w: n ? L.W - 2 * L.x0 : 596 };
    const test = { x: n ? L.x0 : L.x0 + 628, y: area.y, w: n ? L.W - 2 * L.x0 : L.W - 2 * L.x0 - 628 };
    let out = null;
    for (let p = n ? 6 : 8; p >= 2; p -= 0.25) {
      const maxCols = Math.floor(Math.min(area.w, test.w) / p), gap = Math.max(4, p * 1.25);
      const blocks = T.groups.map(G => {
        const cols = clamp(Math.ceil(Math.sqrt(G.n * 1.25)), 1, maxCols);
        return { g: G.g, cols, rows: sum(G.bands.map(b => Math.ceil(b.n / cols))), G };
      });
      let x = 0, y = 0, rowH = 0;
      blocks.forEach(b => {
        const bw = b.cols * p, bh = b.rows * p;
        if (x > 0 && x + bw > area.w) { x = 0; y += rowH + gap; rowH = 0; }
        b.x = area.x + x; b.y = area.y + y; b.w = bw; b.h = bh; x += bw + gap; rowH = Math.max(rowH, bh);
      });
      const homeBottom = area.y + y + rowH;
      if (n) test.y = homeBottom + 30;
      const tcols = Math.floor(test.w / p), randomH = Math.ceil(T.testList.length / tcols) * p;
      const foldH = Math.max(...T.folds.map(f => sum(f.out.map(g => blocks[g].h)) + (f.out.length - 1) * gap));
      const room = L.H - 14 - test.y;
      out = { p, r: p * 0.36, gap, area, test: { ...test, cols: tcols, h: Math.max(randomH, foldH) }, blocks, homeBottom, randomH, foldH };
      if (homeBottom <= L.H - 14 && randomH <= room && foldH <= room) break;
    }
    // home places: each study of a block starts a new row
    const hx = new Float32Array(T.n), hy = new Float32Array(T.n), rx = new Float32Array(T.n), ry = new Float32Array(T.n), p = out.p;
    out.blocks.forEach(b => {
      let row = 0;
      b.G.bands.forEach(band => {
        band.x = b.x; band.y = b.y + row * p; band.rows = Math.ceil(band.n / b.cols); band.w = Math.min(band.n, b.cols) * p;
        for (let k = 0; k < band.n; k++) { const i = band.start + k; hx[i] = b.x + (k % b.cols + 0.5) * p; hy[i] = b.y + (row + Math.floor(k / b.cols) + 0.5) * p; }
        row += band.rows;
      });
    });
    // the test places of the random split, and where the held-out groups of every fold go
    T.testList.forEach((i, k) => { rx[i] = out.test.x + (k % out.test.cols + 0.5) * p; ry[i] = out.test.y + (Math.floor(k / out.test.cols) + 0.5) * p; });
    out.shift = T.folds.map(f => { let y = out.test.y; const m = new Map(); f.out.forEach(g => { const b = out.blocks[g]; m.set(g, [out.test.x - b.x, y - b.y]); y += b.h + out.gap; }); return m; });
    out.divX = n ? null : L.x0 + 612; out.divY = n ? test.y - 16 : null;
    Object.assign(out, { hx, hy, rx, ry, cx: new Float32Array(T.n), cy: new Float32Array(T.n) });
    return out;
  }
  const twinPale = () => mixHex(C.muted, C.raised, 0.5);
  const TWIN_KEYS = ['none', 'array', 'patient', 'both'];
  const twinShare = (a, b) => pct(b ? a / b : 0);
  // the dots of one frame: every record at (x, y) in a color, under the pointer as one hit per record, in record order
  function twinDots(at, colorOf, keyOf, ring = () => false) {
    const P = SIM.L.tw, T = M.T, c = SIM.ctx, spr = {};
    const get = col => spr[col] || (spr[col] = sprite(col, P.r, 0.95, null));
    for (let i = 0; i < T.n; i++) {
      const [x, y, a] = at(i); P.cx[i] = x; P.cy[i] = y;
      if (a <= 0) { hitPoint(-1e4, -1e4, 0, null); continue; }
      const k = keyOf(i); c.globalAlpha = a * hiA(k);
      if (ring(i)) { c.globalAlpha = 1; dot(x, y, P.r * 0.9, null, a * hiA(k), colorOf(i), Math.max(0.8, P.r * 0.45)); }
      else stampV(get(colorOf(i)), x, y);
      hitPoint(x, y, P.r, () => twinTip(i));
    }
    c.globalAlpha = 1;
  }
  // threads from test records to a twin of theirs in training, drawn as one path
  function threads(list, twinOf, a) {
    if (a <= 0 || !list.length) return;
    const P = SIM.L.tw, c = SIM.ctx; c.save(); c.globalAlpha = a * (SIM.hi == null || SIM.hi === 'exp' ? 1 : 0.15); c.strokeStyle = SC.exposed; c.lineWidth = SIM.L.n ? 0.6 : 0.8;
    c.beginPath(); list.forEach(i => { const j = twinOf(i); if (j < 0) return; c.moveTo(P.cx[i], P.cy[i]); c.lineTo(P.cx[j], P.cy[j]); }); c.stroke(); c.restore();
  }
  // the twins of the record under the pointer, joined to it
  function twinHover() {
    const T = M.T, i = SIM.near, P = SIM.L.tw;
    if (SIM.anim !== 'twins' || i < 0 || i >= T.n) return;
    const others = [...(T.arr[i] >= 0 ? T.arrays.get(T.arr[i]) : []), ...(T.pat[i] >= 0 ? T.patients.get(T.pat[i]) : [])].filter(j => j !== i);
    others.forEach(j => { ln(P.cx[i], P.cy[i], P.cx[j], P.cy[j], C.ink, 0.8, 1); dot(P.cx[j], P.cy[j], P.r + 2, null, 1, C.ink, 1.2); });
  }
  function twinTip(i) {
    const T = M.T, K = S.leak, name = j => `${K.studies[T.study[j]]}:${K.rec_acc[j]}`, s = SIM.scene;
    const same = T.arr[i] >= 0 ? T.arrays.get(T.arr[i]).filter(j => j !== i) : [];
    const mine = T.pat[i] >= 0 ? T.patients.get(T.pat[i]).filter(j => T.akey(j) !== T.akey(i)) : [];
    const list = a => a.length > 4 ? `${a.slice(0, 4).map(esc).join(', ')} and ${fmt(a.length - 4)} more` : andList(a.map(esc));
    let h = `<b>${esc(name(i))}</b> <span class="m">· ${esc(K.groups[T.group[i]])}</span>`;
    if (same.length) h += `<br>${sw(SC.twinArray)}The same array as ${list(same.map(name))}`;
    if (mine.length) h += `<br>${sw(SC.twinPatient)}${plural(new Set(mine.map(T.akey)).size, 'other array')} of its patient: ${list(mine.map(name))}`;
    if (!same.length && !mine.length) h += '<br><span class="m">No twin among these records</span>';
    if (s === 1) h += `<br>${T.test[i] ? (T.random.mark[i] ? '<b>Test record with a twin in training</b>' : 'Test record, no twin in training') : 'Training record'}`;
    if (s === 2) h += `<br>${T.loso.mark[i] ? '<b>A twin was in training when its study was tested</b>' : 'No twin in training when its study was tested'}`;
    if (s === 3) h += `<br>${T.elig[i] ? 'A sample of the task' : '<span class="m">Not a sample of the task: it moves with its group, but neither trains nor is tested</span>'}`;
    return h;
  }
  // step 1: the records appear block by block, and those with a twin take the color of its kind
  function drawT1(t) {
    const L = SIM.L, P = L.tw, T = M.T, pale = twinPale(), cols = [pale, SC.twinArray, SC.twinPatient, SC.twinBoth];
    let shown = 0;
    twinDots(i => {
      const g = T.group[i], pop = eout(prog(t, 100 + g * 22, 420)), e = prog(t, 1300 + g * 16, 500);
      if (T.kind[i] && e >= 1) shown++;
      return [P.hx[i], P.hy[i], pop];
    }, i => { const g = T.group[i], e = prog(t, 1300 + g * 16, 500); return T.kind[i] ? mixC(pale, cols[T.kind[i]], e) : pale; }, i => TWIN_KEYS[T.kind[i]]);
    counter(shown, shown === 1 ? 'record with a twin' : 'records with a twin');
    // the two examples: one array deposited in several studies, and the patient with the most arrays
    const a = eout(prog(t, 2500, 500)); if (a <= 0) { twinHover(); return; }
    const K = S.leak, name = j => K.studies[T.study[j]];
    const exA = T.exArray, exP = T.exPatient, studiesA = [...new Set(exA.map(name))], studiesP = [...new Set(exP.map(name))];
    const items = [
      { list: exA, col: SC.twinArray, head: `One array, ${plural(exA.length, 'record')} in ${plural(studiesA.length, 'study', 'studies')}`, body: andList(studiesA) },
      { list: exP, col: SC.twinPatient, head: `${K.example_patient_id ? `Patient ${K.example_patient_id}` : 'One patient'}: ${plural(new Set(exP.map(T.akey)).size, 'array')}`,
        body: `${plural(exP.length, 'record')} in ${andList(studiesP)}` }].filter(x => x.list.length);
    const size = L.fsS, pad = 7, bw = P.test.w;
    let by = L.n ? P.test.y : P.test.y + 30;
    items.forEach(it => {
      const lines = wrap(it.body, bw - 2 * pad, size), bh = (1 + lines.length) * size + (2 + lines.length) * pad - 2, bx = P.test.x;
      const pts = it.list.map(i => [P.hx[i], P.hy[i]]);
      SIM.ctx.save(); SIM.ctx.globalAlpha = a; SIM.ctx.strokeStyle = C.ink; SIM.ctx.lineWidth = 1; SIM.ctx.beginPath();
      pts.forEach(([x, y], j) => { if (j) SIM.ctx.lineTo(x, y); else SIM.ctx.moveTo(x, y); }); SIM.ctx.stroke(); SIM.ctx.restore();
      pts.forEach(([x, y]) => dot(x, y, P.r + 1.8, null, a, C.ink, 1.2));
      const [lx, ly] = pts.reduce((best, q) => Math.hypot(q[0] - bx, q[1] - by) < Math.hypot(best[0] - bx, best[1] - by) ? q : best, pts[0]);
      ln(lx, ly, L.n ? clamp(lx, bx + 8, bx + bw - 8) : bx, L.n ? by : by + bh / 2, C.ink, a * 0.8, 1);
      box(bx, by, bw, bh, C.bg, a * 0.95, 6, C.lineStrong);
      box(bx + pad, by + pad + 1, 3, size + 1, it.col, a, 1.5);
      tx(it.head, bx + pad + 8, by + pad + size - 1, { size, w: 600, c: C.ink, a, max: bw - 2 * pad - 8 });
      lines.forEach((line, j) => tx(line, bx + pad, by + (2 + j) * pad + (2 + j) * size - 3, { size, c: C.ink2, a, max: bw - 2 * pad }));
      by += bh + 14;
    });
    twinHover();
  }
  // the line between training and test, with the name of each side
  function twinSides(a, testName = 'Test') {
    const L = SIM.L, P = L.tw;
    if (L.n) ln(L.x0, P.divY, L.W - L.x0, P.divY, C.lineStrong, a, 1, [4, 4]);
    else ln(P.divX, P.area.y - 22, P.divX, L.H - 10, C.lineStrong, a, 1, [4, 4]);
    tx('Training', P.area.x, P.area.y - 9, { size: L.fs, w: 600, c: C.ink, a });
    tx(testName, P.test.x, P.test.y - 9, { size: L.fs, w: 600, c: C.ink, a, max: P.test.w });
  }
  // step 2: one real random split; the test records fly to their side, and those with a twin in training light up
  function drawT2(t) {
    const L = SIM.L, P = L.tw, T = M.T, R = T.random, pale = twinPale(), nt = T.testList.length;
    twinSides(eout(prog(t, 0, 500)));
    const move = i => eio(prog(t, 500 + T.rank[i] / nt * 1300, 900)), light = i => prog(t, 2900 + T.litRank[i] / Math.max(1, T.lit) * 1700, 250);
    let on = 0;
    T.testList.forEach(i => { if (R.mark[i] && light(i) >= 1) on++; });
    twinDots(i => {
      if (!T.test[i]) return [P.hx[i], P.hy[i], 1];
      const e = move(i); return [lerp(P.hx[i], P.rx[i], e), lerp(P.hy[i], P.ry[i], e) - Math.sin(Math.PI * e) * (L.n ? 10 : 16), 1];
    }, i => !T.test[i] ? pale : R.mark[i] ? mixC(C.ink2, SC.exposed, light(i)) : C.ink2,
    i => !T.test[i] ? 'train' : R.mark[i] && light(i) >= 1 ? 'exp' : 'test');
    const shown = T.testList.filter(i => R.mark[i] && light(i) > 0);
    threads(shown, i => R.twin[i], (C.dark ? 0.2 : 0.26) * Math.min(1, prog(t, 2900, 400)));
    shown.forEach(i => { const j = R.twin[i]; if (j >= 0) dot(P.cx[j], P.cy[j], P.r * 1.1, null, light(i) * hiA('exp'), SC.exposed, 1); });
    counter(twinShare(on, nt), `of ${fmt(nt)} test records have a twin in training`);
    twinHover();
  }
  // step 3: leave one study out; each study in turn is outlined and tested, and its records with a twin in training stay lit
  function drawT3(t) {
    const L = SIM.L, P = L.tw, T = M.T, R = T.random, O = T.loso, pale = twinPale(), nst = O.folds.length;
    twinSides(1 - eio(prog(t, 0, 700)));
    const back = 1 - eio(prog(t, 0, 700)), t0 = 900, dt = (TWIN_SCENES[2].dur - t0 - 500) / nst;
    const k = Math.floor((t - t0) / dt), cur = k >= 0 && k < nst ? k : -1;
    const tested = s => t >= t0 + (s + 0.5) * dt;
    twinDots(i => {
      if (!T.test[i] || back <= 0) return [P.hx[i], P.hy[i], 1];
      return [lerp(P.hx[i], P.rx[i], back), lerp(P.hy[i], P.ry[i], back), 1];
    }, i => {
      const s = T.study[i];
      if (back > 0 && T.test[i]) return R.mark[i] ? mixC(pale, SC.exposed, back) : mixC(pale, C.ink2, back);
      if (s === cur) return O.mark[i] && tested(s) ? SC.exposed : C.ink2;
      return O.mark[i] && tested(s) ? SC.exposed : pale;
    }, i => O.mark[i] && tested(T.study[i]) ? 'exp' : T.study[i] === cur ? 'test' : 'train');
    threads(T.testList.filter(i => R.mark[i]), i => R.twin[i], (C.dark ? 0.2 : 0.26) * back);
    if (cur >= 0) {
      const pad = 2, mine = [];
      T.bands.filter(b => b.s === cur).forEach(b => {
        box(b.x - pad, b.y - pad, b.w + 2 * pad, b.rows * P.p + 2 * pad, null, 1, 3, C.ink);
        for (let i = b.start; i < b.start + b.n; i++) if (O.mark[i]) mine.push(i);
      });
      if (tested(cur)) threads(mine, i => O.twin[i], 0.55);
      const f = O.folds[cur], size = L.fsS, a = eout(prog(t, t0, 300));
      tx(`Testing ${S.leak.studies[cur]}`, P.test.x, P.test.y + size, { size: L.fs, w: 600, c: C.ink, a, max: P.test.w });
      tx(`study ${fmt(cur + 1)} of ${fmt(nst)}, ${plural(f.test, 'record')}`, P.test.x, P.test.y + 2 * size + 6, { size, c: C.ink2, a, max: P.test.w });
      tx(`${fmt(f.exposed)} with a twin in training`, P.test.x, P.test.y + 3 * size + 12, { size, w: 600, c: f.exposed ? SC.exposed : C.ink2, a, max: P.test.w });
    }
    const done = O.folds.filter((_, s) => tested(s)), e = sum(done.map(f => f.exposed)), n_ = sum(done.map(f => f.test));
    counter(twinShare(e, n_), `of ${fmt(n_)} records had a twin in training when their study was tested`);
    twinHover();
  }
  // step 4: the patient-safe split; each held-out group goes to the test side with all its records, twins included
  function drawT4(t) {
    const L = SIM.L, P = L.tw, T = M.T, O = T.loso, pale = twinPale(), F = T.folds, nf = F.length;
    const fade = eio(prog(t, 0, 700)), t0 = 900, dt = (TWIN_SCENES[3].dur - t0 - 400) / nf;
    const k = clamp(Math.floor((t - t0) / dt), -1, nf - 1), ph = k >= 0 ? (t - t0 - k * dt) / dt : 0;
    twinSides(eout(prog(t, 300, 500)), k >= 0 ? `Test · ${F[k].fold}, fold ${k + 1} of ${nf}` : 'Test');
    const e = k < 0 ? 0 : ph < 0.35 ? eio(ph / 0.35) : ph < 0.75 || k === nf - 1 ? 1 : 1 - eio((ph - 0.75) / 0.25);
    const shift = k >= 0 ? P.shift[k] : null;
    twinDots(i => {
      const d = shift && shift.get(T.group[i]);
      if (!d) return [P.hx[i], P.hy[i], 1];
      return [P.hx[i] + d[0] * e, P.hy[i] + d[1] * e - Math.sin(Math.PI * e) * (L.n ? 8 : 14), 1];
    }, i => {
      const out = shift && shift.has(T.group[i]) && e > 0.5;
      if (!T.elig[i]) return fade > 0.5 ? mixHex(C.muted, C.raised, 0.2) : pale;
      if (out) return C.ink2;
      return O.mark[i] ? mixC(SC.exposed, pale, fade) : pale;
    }, i => !T.elig[i] ? 'ghost' : shift && shift.has(T.group[i]) && e > 0.5 ? 'test' : 'train', i => !T.elig[i] && fade > 0.5);
    const started = F.filter((_, j) => j < k || (j === k && ph >= 0.35)), n_ = sum(started.map(f => f.c.test)), x = sum(started.map(f => f.c.exposed));
    if (t < t0) counter(pct(lerp(O.c.exposed / O.c.test, 0, fade)), 'of the records have a twin in training');
    else counter(twinShare(x, n_), `of ${fmt(n_)} test samples have a twin in training`);
    twinHover();
  }
  function twinHeader(s) {
    const K = S.leak, task = `${lab('class', K.positive)} against ${labLower('class', K.negative)}s`, n = M.T.n;
    return [
      [`${fmt(n)} records of the task ${task.toLowerCase()}`, 'one dot per record, in blocks by patient-safe group; colored: records with a twin', `${fmt(n)} records of one task`],
      ['A random split', `${pct(K.test_fraction, 0)} of the records drawn by chance for the test side`, 'A random split'],
      ['Leave one study out', 'each study in turn is the test side, all the others train', 'One study out'],
      ['The patient-safe split', 'each patient-safe group in turn is the test side, with all its records', 'Patient-safe split']][s];
  }
  function twinCaption(s) {
    const K = S.leak, A = K.audit, AR = A.all_records, T = M.T, R = T.random.c, RA = A.random_all, O = T.loso.c, G = A.gts, b = v => `<b>${v}</b>`;
    const task = `${lab('class', K.positive)} against ${labLower('class', K.negative)}s`, sc = K.scores, models = Object.keys(K.models);
    const exA = T.exArray, exP = T.exPatient, name = j => K.studies[T.study[j]];
    return [
      ['Records that have a twin',
        `A classifier is judged on records it has not seen. In public data many records are not independent: the same array was often deposited in several studies, and many patients gave several samples. This animation takes the task ${task.toLowerCase()} as a naive analysis would pool it: the ${fmt(T.n)} records of the two classes and of the blood sample types of the benchmark, one dot each, in blocks by patient-safe group and, inside a block, by study. A record has a <i>twin</i> when another record is the same array, deposited again, or another sample of the same patient. ${fmt(A.linked)} of the ${fmt(T.n)} records have at least one.` +
          (exA.length ? ` One array here has ${spell(exA.length)} records, in ${andList([...new Set(exA.map(name))])}.` : ''),
        `${b(fmt(A.linked))} of ${fmt(T.n)} records have a twin (${pct(A.linked / T.n)})<br>The same array ${fmt(A.linked_array)} · another sample of the patient ${fmt(A.linked_patient)}<br>In the whole resource ${fmt(AR.with_a_linked_record)} of ${fmt(AR.total)} records (${pct(AR.share_with_a_linked_record)}); ${fmt(AR.arrays_in_more_than_one_study)} arrays were deposited in more than one study, in ${fmt(AR.studies_sharing_an_array)} studies`],
      ['A random split leaves twins on both sides',
        `A random split draws the test records by chance, here ${pct(K.test_fraction, 0)} of them: ${fmt(R.test)} records go to the test side and the other ${fmt(T.n - R.test)} train. A test record whose twin is in training is no longer unseen, because the classifier has learned from the same array or the same patient. In this split ${fmt(R.exposed)} of the ${fmt(R.test)} test records have a twin in training, about the share that chance predicts from the number of twins of each record (${pct(A.random.expected)}). Each thread joins one of them to its twin.`,
        `${b(`${fmt(R.exposed)} of ${fmt(R.test)}`)} test records have a twin in training (${pct(R.exposed / R.test)})<br>By the same array ${fmt(R.array)}, by the same patient ${fmt(R.patient)}<br>All ${spell(RA.folds)} random splits of the benchmark: ${fmt(RA.exposed)} of ${fmt(RA.test)} (${pct(RA.exposed / RA.test)})`],
      ['Holding out whole studies is not enough',
        `Leaving one study out tests each of the ${fmt(T.loso.folds.length)} studies in turn, with all the other studies as training data. Twins inside one study now stay together, but an array deposited in two studies, or a patient sampled in two, still has a copy in training. Over all the folds, ${fmt(O.exposed)} of the ${fmt(O.test)} records had a twin in training when their study was tested; they stay lit.`,
        `${b(`${fmt(O.exposed)} of ${fmt(O.test)}`)} records (${pct(O.exposed / O.test)}) over ${fmt(T.loso.folds.length)} folds<br>By the same array ${fmt(O.array)}, by the same patient ${fmt(O.patient)}`],
      ['Patient-safe groups keep every twin on one side',
        `The benchmark of the resource keeps one record of every array and of every patient, the samples of the task, and splits by patient-safe group: studies joined wherever they share an array or a patient, so that all the twins of a record lie in its group. Leave one group out tests each group with enough samples of both classes in turn, and every record of the group moves with it. None of the ${fmt(G.test)} test samples of its ${fmt(G.folds)} folds has a twin in training, and neither has any test sample of the ${fmt(A.patient_safe_folds)} patient-safe folds of the three tasks.`,
        `${fmt(A.samples)} samples of the task among the ${fmt(T.n)} records; the hollow dots are the others<br>${fmt(G.folds)} folds, ${fmt(G.test)} test samples, ${b('0')} with a twin in training<br>Mean AUROC, random split against leave one group out: ${models.map(m => `${esc(K.models[m])} ${fx(sc[m].naive_record_split, 3)} and ${fx(sc[m][K.scheme], 3)}`).join('; ')}`]][s];
  }
  const TWIN_HOWTO = [
    'Each dot is one record, and each block one patient-safe group, with its studies stacked inside it. Colored dots have a twin among these records. Hover over a dot to join it to its twins.',
    'Dots on the test side are the test records. A lit dot has a twin in training, and its thread crosses the split to that twin.',
    'The outlined band is the study being tested. Its records that still have a twin in training, in another study of the same group, light up and stay lit.',
    'Each patient-safe group in turn moves to the test side with all its records, its twins included. Hollow dots are records that are not samples of the task.'];
  function twinKeys() {
    const s = SIM.scene, T = M.T, el = $('storyKeys'), hov = k => { if (SIM.hi !== k) { SIM.hi = k; if (!SIM.running) simDraw(); } };
    const pale = twinPale(), cnt = k => T.kind.filter(v => v === k).length;
    let items = [];
    if (s === 0) items = [{ key: 'array', label: 'Twin: the same array', color: SC.twinArray, count: cnt(1) }, { key: 'patient', label: 'Twin: another sample of the patient', color: SC.twinPatient, count: cnt(2) },
      { key: 'both', label: 'Both', color: SC.twinBoth, count: cnt(3) }, { key: 'none', label: 'No twin', color: pale, count: cnt(0) }];
    if (s === 1) items = [{ key: 'train', label: 'Training', color: pale, count: T.n - T.testList.length }, { key: 'test', label: 'Test, no twin in training', color: C.ink2, count: T.random.c.test - T.random.c.exposed },
      { key: 'exp', label: 'Test, with a twin in training', color: SC.exposed, count: T.random.c.exposed }];
    if (s === 2) items = [{ key: 'train', label: 'Training', color: pale, count: null }, { key: 'test', label: 'The study being tested', color: C.ink2, count: null },
      { key: 'exp', label: 'Had a twin in training when tested', color: SC.exposed, count: T.loso.c.exposed }];
    if (s === 3) items = [{ key: 'train', label: 'Training or not yet tested', color: pale, count: null }, { key: 'test', label: 'Test samples of the held-out group', color: C.ink2, count: null },
      { key: 'ghost', label: 'Not a sample of the task (hollow)', color: pale, count: T.n - T.samples }];
    legend(el, items, hov);
    $('storyHow').lastElementChild.textContent = TWIN_HOWTO[s];
  }

  // ---------- captions, keys and the note under the stage (one group close up) ----------
  function caption(s) {
    const g = S.group, st = g.stated, sc = S.screen, sx = S.sex, B = S.bench, n = v => `<b>${v}</b>`;
    const cols = M.cols, focusCol = cols.find(c => c.study === g.study);
    const excluded = sc.candidates - N.studies, rows = M.rows.filter(r => r.group !== 'included');
    const top = rows.slice().sort((a, b) => b.n - a.n).slice(0, 3);
    const access = sum(rows.filter(r => r.group === 'access').map(r => r.n));
    const ex = M.samples[g.example_sepsis], twin = M.samples[g.example_twin], twinHere = g.twin_records.find(r => cols[r.col].study === g.study && !r.primary);
    const pat = M.pats[g.example_patient], several = M.pats.filter(p => p.samples.length > 1).length;
    const differs = M.samples.filter(s => s.differs).length;
    const recSex = sx.matrix, groupSex = sx.group_sex.female + sx.group_sex.male;
    const corr = M.corr, E = H.explained;
    const train = M.blocks.filter(b => !b.test), trP = sum(train.map(b => b.pos)), trN = sum(train.map(b => b.neg)), te = M.test;
    const firstIn = Object.entries(g.focus_first_in).map(([c, k]) => `${cols[+c].study} (${fmt(k)})`);
    const linked = g.linked_by === 'chip_type_scan_time_and_median_intensity'
      ? 'The raw file of an Affymetrix array records its chip type, the moment it was scanned and its median intensity; records that share all three are one array.'
      : 'Records with the same evidence in their raw files are one array or library.';
    const platforms = [...new Set(cols.map(c => c.platform))];
    const capsFor = [
      ['Searching and screening',
        `Searches of public repositories, mainly GEO, ArrayExpress and the European Nucleotide Archive, returned ${fmt(sc.census_records)} records. Copies of one study in several repositories were counted once, which left ${fmt(sc.candidates)} <i>candidate studies</i>, one dot each. Every candidate was checked against the eligibility rules in a fixed order, and the first rule it failed is its reason for exclusion. The ${fmt(N.studies)} included studies measured RNA in whole blood, leukocytes or PBMC of patients whom the study itself diagnosed with sepsis.`,
        `${fmt(sc.candidates)} candidates − ${fmt(excluded)} excluded = ${n(fmt(N.studies))} included<br>Most common reasons: ${andList(top.map(r => `${labLowerText(r.label)} (${fmt(r.n)})`))}<br>${fmt(access)} met the criteria, but no expression values could be obtained`],
      ['One cohort group, close up',
        `Investigators sometimes deposit the same arrays in several studies. ${g.study}, a study of children with septic shock, shares arrays with ${spell(cols.length - 1)} other studies of the same investigators, released from ${cols[0].year} to ${cols[cols.length - 1].year}. Study platforms that share a sample form a <i>cohort group</i>; this one is ${g.id}, and the next steps follow it. Each dot is one sample record, colored by the study it was deposited in.`,
        `${g.id}: ${n(fmt(st.records))} sample records in ${spell(cols.length)} studies${platforms.length === 1 && platforms[0] ? `, all on ${platforms[0]}` : ''}<br>${g.study}: ${fmt(focusCol.records)} records`],
      ['Repeat deposits join their first deposit',
        `${linked} Each repeat deposit joins the first deposit of its array, its <i>primary record</i>, so that the array counts once and keeps all its records. ${M.noOwn.length ? `${andList(M.noOwn.map(c => c.study))} hold${M.noOwn.length === 1 ? 's' : ''} no array of ${M.noOwn.length === 1 ? 'its' : 'their'} own: every one of ${M.noOwn.length === 1 ? 'its' : 'their'} records is a repeat deposit.` : ''}`,
        `${fmt(st.records)} records → ${n(fmt(st.distinct))} distinct samples, ${fx(st.records / st.distinct, 1)} records for every array<br>${g.study}: ${fmt(focusCol.records)} records, ${fmt(focusCol.primary)} of them primary; the other ${fmt(focusCol.records - focusCol.primary)} were first deposited in ${andList(firstIn)}<br>${Spell(g.most_records)} records of one array at most, for ${fmt(g.samples_with_most)} arrays`],
      ['A class from the words of the record',
        `Every sample takes its class from what its record says. In ${cols[ex.col].study}, the field ${q(ex.field)} of ${ex.acc} reads ${q(ex.said)}, which makes it a ${labLower('class', D.order.class[ex.cls])} sample. A control is a healthy control only where its record calls it healthy or normal. The array ${cols[twin.col].study}:${twin.acc} is recorded as ${q(twin.said)} in ${cols[twin.col].study}${twin.stmt >= 0 ? `, whose description says ${q(g.statements[twin.stmt])}` : ''}; nothing there calls the controls healthy, so that record is a ${labLower('class', D.order.class[twin.cls])}.` +
          (twinHere ? ` Its record in ${g.study}, ${twinHere.acc}, says ${q(twinHere.said)} and is a ${labLower('class', D.order.class[twinHere.cls])}. Each record keeps its own label, and the disagreement is marked.` : '') + ' A dot takes the class of its primary record.',
        D.order.class.filter(k => st.classes[k] > 0).map(k => `${lab('class', k)} ${n(fmt(st.classes[k]))}`).join(' · ') + `<br>${fmt(differs)} samples of ${g.id} have records that differ in class`],
      ['Samples of one patient, one marked for each',
        `Some patients were sampled on more than one day. Their records name the patient, and samples with one identifier are joined.` +
          (pat ? ` Patient ${pat.id}, for example, gave ${andList(pat.samples.map(i => { const s = M.samples[i]; return `${s.acc} on ${day(s)}, recorded as ${q(s.said)}`; }))}.` : '') +
          ` One sample is marked for each patient: preferably one without a quality flag, and the earliest. A sample without a patient identifier counts as its own patient, so the number of marked samples is an upper bound on the number of patients.`,
        `${fmt(st.distinct)} samples · ${n(fmt(st.patients))} identified patients, ${fmt(several)} of them with two or more samples<br>${Spell(st.without_patient)} samples without a patient identifier<br>${n(fmt(st.one_per_patient))} marked samples, one for every patient`],
      ['Sex read from the genes',
        `Genes of the Y chromosome are active in the blood cells of men, and ${sx.female_gene} in those of women. In every table, the samples were split into two groups by these genes alone, and only then compared with the recorded sex. The records of ${g.id} give the sex of only ${fmt(groupSex)} of its ${fmt(st.distinct)} samples, so this step uses ${sx.study}, whose ${fmt(recSex.samples)} arrays all have one. Each dot is placed by its ${sx.female_gene} and its Y genes, scaled so that the two groups of the table sit at 0 and 1, and colored by its recorded sex; ringed dots disagree with their record.`,
        `Read as female ${fmt(recSex.called_female)}, as male ${fmt(recSex.called_male)}, too close to call ${fmt(recSex.ambiguous)}<br>Agree with the record: ${n(`${fmt(recSex.concordant)} of ${fmt(recSex.concordant + recSex.discordant)}`)} (${pct(recSex.concordant / (recSex.concordant + recSex.discordant))})<br>All tables: ${pct(N.sex_concordance)} of ${fmt(N.sex_compared)} samples with a recorded sex`],
      ['Removing the differences between studies',
        `Values measured in different studies are not on one scale. Each dot is one of the ${fmt(H.n)} samples of the bulk-blood table, placed by the first two principal components of its ${fmt(H.genes)} most variable genes. Before the correction the samples fall into islands by technology and study, and the ${fmt(S.correction.samples.length)} samples of ${g.id} (ringed) form one of them. Studies with at least ${spell(S.correction.minimum_healthy)} healthy controls were corrected on those controls; the others by matched groups, which give each class the average and spread of the same class in the first. After the correction, the samples of ${g.id} spread through the common cloud. Because matched groups use the classes, part of the separation of classes afterwards comes from the correction itself.`,
        `Technology explains ${fx(E.before.family[0], 2)} of the variation along component 1 before and ${n(fx(E.after.family[0], 2))} after; study ${fx(E.before.study_accession[0], 2)} and ${n(fx(E.after.study_accession[0], 2))}<br>` +
          corr.map(([m, studies, k]) => `${g.id}: ${fmt(k)} ${m === 'healthy_anchored' ? 'on healthy controls' : 'by matched groups'} (${andList(studies)})`).join('<br>')],
      ['Training and test data that share no sample',
        `A classifier is judged on samples it has not seen. In the task ${taskName().toLowerCase()}, each block is one <i>patient-safe group</i> and each square one sample of the task: the sample marked for its patient, in a table, without a flag, of bulk blood, and not taken before sepsis began or during recovery. Leave one group out holds out, in turn, every group with at least ${spell(B.minimum_per_class)} samples of each class; in this fold ${te.g} is the test set and the ${fmt(train.length)} other groups train. All ${fmt(B.held.records)} records of ${te.g} go to the test side: its ${fmt(B.held.reasons.not_primary_record || 0)} repeat deposits and the ${fmt(B.held.reasons.not_one_per_patient || 0)} other samples of its patients are not samples of the task, and none of them can reach the training data.`,
        `Training: ${fmt(train.length)} groups, ${fmt(trP)} ${labLower('class', B.positive)} and ${fmt(trN)} ${labLower('class', B.negative)}s<br>Test: ${te.g}, ${fmt(te.pos)} ${labLower('class', B.positive)} and ${fmt(te.neg)} ${labLower('class', B.negative)}s<br>${fmt(trP)} + ${fmt(te.pos)} = ${n(fmt(B.totals.pos))} and ${fmt(trN)} + ${fmt(te.neg)} = ${n(fmt(B.totals.neg))}, the samples of the task<br>${fmt(B.folds)} folds of this task hold out one group each`]];
    return capsFor[s];
  }
  const GROUP_HOWTO = [
    'Each dot is one candidate study. The rows are the reasons for exclusion, in groups, with the number of studies; hover over a dot to see its study and the rule that excluded it.',
    'Each column is one study of the group, in the order of release; each dot is one sample record. Hover over a dot to see whether it is the first deposit of its array.',
    'Dots that leave a column are repeat deposits; they fly to the first deposit of their array and merge with it. The number under each study counts what is left in its column.',
    'Each dot is one distinct sample, now colored by its class; within each study the classes gather. Hover over a dot to read the words its class was taken from.',
    'Dots joined by a line are samples of one patient, in the order of sampling. A filled dot is the sample marked for its patient; a hollow one is not.',
    'Dots to the lower right have high XIST and low Y genes, as in women; dots to the upper left the reverse. A dot whose color does not match its corner disagrees with its record.',
    'Dots close together have similar expression profiles. Watch the islands of technologies merge, and the ringed samples of the group spread out. Hover over a key to pick out a technology.',
    'Each block is one patient-safe group and each square one sample of the task. The outlined block moves to the test side with every record of its group, including those that are not samples of the task.'];
  function groupKeys() {
    const s = SIM.scene, g = S.group, el = $('storyKeys');
    const hov = k => { if (SIM.hi !== k) { SIM.hi = k; if (!SIM.running) simDraw(); } };
    let items = [], extra = '';
    if (s === 0) items = [{ key: 'included', label: 'Included', color: C.accent, count: N.studies }, { key: 'excluded', label: 'Excluded', color: C.lineStrong, count: S.screen.candidates - N.studies }];
    if (s === 1 || s === 2) items = M.cols.map(c => ({ key: String(c.k), label: c.study, color: studyColor(c.k), count: s === 1 ? c.records : c.primary }));
    if (s === 3 || s === 4) items = D.order.class.filter(k => g.stated.classes[k] > 0).map(k => ({ key: k, label: lab('class', k), color: C.cls[k], count: g.stated.classes[k] }));
    if (s === 4) extra = keys([['marked sample', C.ink2, 'dot'], ['another sample of the patient', C.ink2, 'ring']]);
    if (s === 5) {
      const cnt = k => M.sexPts.filter(p => p.rec === k).length;
      items = ['female', 'male', 'not_recorded'].map(k => ({ key: k, label: k === 'not_recorded' ? 'Sex not recorded' : `Recorded ${k}`, color: SC.sex[k], count: cnt(k) })).filter(x => x.count > 0);
      extra = keys([[`disagrees with its record (${fmt(M.sexPts.filter(p => p.bad).length)})`, C.ink, 'ring'], [`too close to call (${fmt(M.sexPts.filter(p => p.call === 'ambiguous').length)})`, SC.sex.ambiguous, 'ring']]);
    }
    if (s === 6) {
      const cnt = D.order.family.map(() => 0); for (let i = 0; i < H.n; i++) cnt[+H.fam[i]]++;
      items = D.order.family.map((f, i) => ({ key: f, label: lab('family', f), color: C.fam[f], count: cnt[i] })).filter(x => x.count > 0);
      items.push({ key: 'mine', label: `Samples of ${g.id} (ringed)`, color: C.ink, count: S.correction.samples.length });
    }
    if (s === 7) {
      const B = S.bench;
      items = [{ key: B.positive, label: lab('class', B.positive), color: C.cls[B.positive], count: B.totals.pos }, { key: B.negative, label: lab('class', B.negative), color: C.cls[B.negative], count: B.totals.neg },
        { key: 'ghost', label: `Other records of ${M.test.g}`, color: SC.ghost, count: M.ghosts }];
    }
    legend(el, items, hov);
    if (extra) el.insertAdjacentHTML('beforeend', extra);
    $('storyHow').lastElementChild.textContent = GROUP_HOWTO[s];
  }
  // the three animations: their steps, drawings, captions and keys
  const ANIM = {
    whole: { scenes: WHOLE_SCENES, draw: [drawW1, drawW2, drawW4, drawW3, drawW5, drawW6], header: wholeHeader, caption: wholeCaption, keys: wholeKeys,
             end: '<button class="btn primary" type="button" data-go="story/group/1">Watch one group close up</button><button class="btn" type="button" data-go="resource">Explore the resource</button>' },
    group: { scenes: GROUP_SCENES, draw: [drawScreen, drawGroup, drawMerge, drawLabels, drawPatients, drawSex, drawCorrection, drawBench], header: headerText, caption, keys: groupKeys,
             end: '<button class="btn primary" type="button" data-go="story/twins/1">Watch twins and splits</button><button class="btn" type="button" data-go="resource">Explore the resource</button>' },
    twins: { scenes: TWIN_SCENES, draw: [drawT1, drawT2, drawT3, drawT4], header: twinHeader, caption: twinCaption, keys: twinKeys,
             end: '<button class="btn primary" type="button" data-go="resource/r-evaluation">See what the splits change</button><button class="btn" type="button" data-go="methods/m-bench">Read the method in detail</button>' } };
  const ANIM_NAME = { whole: 'The whole resource', group: 'One group close up', twins: 'Twins and splits' };
  const A = () => ANIM[SIM.anim];
  function renderKeys() { A().keys(); }
  function renderCaption() {
    const s = SIM.scene, an = A(), [title, text, nums] = an.caption(s), n = an.scenes.length, last = s === n - 1;
    $('storyCap').innerHTML = `<p class="story-cap-k">${ANIM_NAME[SIM.anim]} · step ${s + 1} of ${n}</p><h3>${title}</h3><p>${text}</p><div class="story-cap-e">${nums}</div>` +
      (last ? `<div class="story-cap-go">${an.end}</div>` : '');
    const plain = $('storyCap').querySelector('p:not(.story-cap-k)').textContent;
    $('storyCanvas').setAttribute('aria-label', `Animation, step ${s + 1} of ${n}: ${title}. ${plain}`);
    SIM.hold = clamp((title + ' ' + plain).split(/\s+/).length * 140, 5000, 11000);
    renderKeys();
  }

  // ---------- playback ----------
  function renderScrub() {
    const s = SIM.scene, d = A().scenes[s].dur, p = !SIM.started ? 0 : SIM.auto ? Math.min(1, SIM.t / (d + SIM.hold)) : Math.min(1, SIM.t / d);
    $('storyScrub').querySelectorAll('button').forEach((b, i) => {
      b.classList.toggle('on', i === s); b.setAttribute('aria-selected', i === s);
      b.querySelector('b').style.width = `${i === s ? 100 * p : SIM.started && i < s ? 100 : 0}%`;
    });
    const pl = $('storyPlay');
    pl.innerHTML = SIM.auto ? ICON_PAUSE : ICON_PLAY;
    pl.setAttribute('aria-label', SIM.auto ? 'Pause' : 'Play'); pl.title = SIM.auto ? 'Pause (space)' : 'Play (space)';
    $('storyBack').disabled = SIM.started && s === 0; $('storyNext').disabled = SIM.started && s === A().scenes.length - 1;
  }
  function loop(now) {
    SIM.raf = 0;
    if ($('story').hidden) { SIM.auto = false; SIM.running = false; return; }
    const dt = Math.min(100, now - (SIM.last || now)); SIM.last = now;
    if (SIM.running) {
      const d = A().scenes[SIM.scene].dur;
      SIM.t = (reducedMotion() ? Math.max(SIM.t, d) : SIM.t) + dt;
      if (SIM.auto && SIM.t >= d + SIM.hold) {
        if (SIM.scene < A().scenes.length - 1) { go(SIM.scene + 1); return; }
        SIM.auto = false;
      }
      if (!SIM.auto && SIM.t >= d) SIM.running = false;
    }
    simDraw(); renderScrub();
    if (SIM.running) SIM.raf = requestAnimationFrame(loop);
  }
  function kick() { if (!SIM.raf) { SIM.last = performance.now(); SIM.raf = requestAnimationFrame(loop); } }
  function go(s, animateIt = true) {
    const changed = s !== SIM.scene;
    SIM.scene = clamp(s, 0, A().scenes.length - 1); SIM.started = true; $('storyBigPlay').hidden = true;
    SIM.t = animateIt ? 0 : 1e9; SIM.running = animateIt; SIM.near = -1; untip();
    if (changed) SIM.hi = null;
    renderCaption();
    if (ROUTE.view === 'story' && !$('story').hidden) { const h = `story/${SIM.anim}/${SIM.scene + 1}`; try { if (location.hash.slice(1) !== h) history.replaceState(null, '', `#${h}`); } catch (e) { } }
    renderScrub();
    if (animateIt) kick(); else simDraw();
  }
  function toggle() {
    if (!SIM.started) { SIM.auto = true; go(0); return; }
    SIM.auto = !SIM.auto;
    if (SIM.auto) {
      if (SIM.scene === A().scenes.length - 1 && SIM.t >= A().scenes[SIM.scene].dur) { go(0); return; }
      SIM.running = true; kick();
    }
    renderScrub();
  }
  function step(dir) { if (!SIM.started) { go(0); return; } go(SIM.scene + dir); }

  // ---------- the pointer ----------
  function nearest(e) {
    const cv = $('storyCanvas'), r = cv.getBoundingClientRect(), L = SIM.L, h = SIM.hit;
    if (!h || !r.width) return { i: -1 };
    const x = (e.clientX - r.left) * L.W / r.width, y = (e.clientY - r.top) * L.W / r.width;
    let best = -1, bd = Infinity;
    for (let i = 0; i < h.x.length; i++) {
      const d = (h.x[i] - x) ** 2 + (h.y[i] - y) ** 2, lim = Math.max(5, h.r[i] + 2.5) ** 2;
      if (d < lim && d < bd) { bd = d; best = i; }
    }
    if (best >= 0) return { i: best, html: h.f[best] };
    const rect = h.rects.find(([rx, ry, rw, rh]) => x >= rx && x <= rx + rw && y >= ry && y <= ry + rh);
    return rect ? { i: -1, html: rect[4] } : { i: -1 };
  }
  function initPointer() {
    const cv = $('storyCanvas');
    const onMove = e => {
      if (!SIM.started) return;
      const hit = nearest(e);
      if (hit.i !== SIM.near) { SIM.near = hit.i; if (!SIM.running) simDraw(); }
      if (hit.html) { tip(e, hit.html()); cv.style.cursor = 'pointer'; } else { untip(); cv.style.cursor = 'default'; }
    };
    cv.addEventListener('pointermove', e => { if (e.pointerType === 'mouse') onMove(e); });
    cv.addEventListener('pointerdown', e => { if (e.pointerType !== 'mouse') onMove(e); });
    cv.addEventListener('mouseleave', () => { untip(); if (SIM.near >= 0) { SIM.near = -1; if (!SIM.running) simDraw(); } });
  }

  // ---------- the view ----------
  function init() {
    S = part('story'); M = buildModel(); M.W = buildWhole(); M.T = buildTwins();
    M.corr = ['healthy_anchored', 'matched_groups'].map(m => {
      const rows = S.correction.by_study.filter(r => r.method === m);
      return [m, rows.map(r => r.study).sort((a, b) => M.cols.findIndex(c => c.study === a) - M.cols.findIndex(c => c.study === b)), sum(rows.map(r => r.n))];
    }).filter(r => r[2] > 0);
    buildScrub();
    seg('storyAnim', v => { setAnim(v); try { history.replaceState(null, '', `#story/${v}`); } catch (e) { } });
    $('storyScrub').addEventListener('click', e => { const b = e.target.closest('button'); if (b) go([...$('storyScrub').children].indexOf(b)); });
    $('storyBack').onclick = () => step(-1);
    $('storyNext').onclick = () => step(1);
    $('storyPlay').onclick = toggle;
    $('storyBigPlay').onclick = () => { SIM.auto = true; go(0); };
    fullscreenButton('storyFull', 'storyGrid', () => { SIM.ctx = null; if (!$('story').hidden) simDraw(); });
    initPointer();
    const q = new URLSearchParams(location.search).get('t');
    if (q != null && q !== '' && isFinite(+q)) SIM.fixedT = Math.max(0, +q);
    let w = 0;
    new ResizeObserver(() => requestAnimationFrame(() => { const nw = $('storyStage').clientWidth; if (nw && Math.abs(nw - w) > 2) { w = nw; SIM.ctx = null; if (!$('story').hidden) { simDraw(); } } })).observe($('storyStage'));
    renderCaption(); renderScrub();
  }
  function buildScrub() {
    const sc = A().scenes, box = $('storyScrub');
    box.style.gridTemplateColumns = `repeat(${sc.length}, minmax(0, 1fr))`;
    box.innerHTML = sc.map((x, i) => `<button type="button" role="tab" aria-label="Step ${i + 1}: ${x.label}"><i><b></b></i><span>${x.label}</span></button>`).join('');
  }
  // switch to an animation, at its first step and not yet playing
  function setAnim(v) {
    if (!ANIM[v]) return;
    SIM.anim = v; setSeg('storyAnim', v); buildScrub();
    SIM.scene = 0; SIM.auto = false; SIM.running = false; SIM.started = false; SIM.t = 1e9; SIM.hi = null; SIM.near = -1; untip();
    $('storyBigPlay').hidden = false;
    renderCaption(); renderScrub(); if (SIM.ctx) simDraw();
  }
  // #story opens the whole resource; #story/whole/3, #story/group/3 and #story/twins/3 a step of one; #story/3 a step of the group
  function show(route) {
    const a = route.args, which = ANIM[a[0]] ? a[0] : /^\d+$/.test(a[0] || '') ? 'group' : null;
    const n = parseInt(which && a[0] === which ? a[1] : a[0], 10);
    if (which && which !== SIM.anim) setAnim(which);
    if (n >= 1 && n <= A().scenes.length) {
      go(n - 1, false);
      if (SIM.fixedT != null) { SIM.t = SIM.fixedT; SIM.fixedT = null; }
    }
    SIM.ctx = null; simDraw(); renderScrub(); renderKeys();
  }
  registerView('story', {
    part: 'story', init, show,
    redraw() { SIM.ctx = null; renderKeys(); simDraw(); },
    routes: () => Object.keys(ANIM).flatMap(k => ANIM[k].scenes.map((_, i) => `story/${k}/${i + 1}`)),
    onKey(e) {
      if (e.key === 'ArrowRight') { step(1); return true; }
      if (e.key === 'ArrowLeft') { step(-1); return true; }
      if (e.key === ' ' && !e.target.closest('button')) { toggle(); return true; }
      if ((e.key === 'f' || e.key === 'F') && canFullscreen()) { toggleFullscreen($('storyGrid')); return true; }
      return false;
    },
  });

  // ---------- checks (?selftest=1) ----------
  const ok = (cond, detail) => ({ ok: !!cond, detail });
  selftest('story', 'step 1: every candidate is one dot, and the dots of each reason equal its count', () => {
    const bad = M.rows.filter(r => r.items.length !== r.n).map(r => `${r.code} ${r.items.length} against ${r.n}`);
    const total = sum(M.rows.map(r => r.n)), incl = M.rows.find(r => r.code === 'eligible');
    return ok(!bad.length && M.studies.length === S.screen.candidates && total === N.candidates && incl.n === N.studies,
      bad.join('; ') || `${fmt(M.studies.length)} dots, ${fmt(total)} in the rows, ${fmt(incl.n)} included`);
  });
  selftest('story', 'step 2: the columns hold every record of the cohort group', () => {
    const bad = M.cols.filter(c => c.recs.length !== c.records).map(c => `${c.study} ${c.recs.length} against ${c.records}`);
    return ok(!bad.length && M.recs.length === S.group.stated.records, bad.join('; ') || `${fmt(M.recs.length)} records`);
  });
  selftest('story', 'step 3: the records collapse to the stated number of distinct samples', () => {
    const distinct = new Set(M.recs.map(r => r.s)).size, primaries = M.recs.filter(r => r.primary).length;
    const wrong = M.samples.filter(s => s.recs.length !== s.nrec || s.prim == null || M.recs[s.prim].col !== s.col || M.recs[s.prim].slot !== s.slot);
    return ok(distinct === S.group.stated.distinct && primaries === distinct && M.samples.length === distinct && !wrong.length,
      `${fmt(M.recs.length)} records, ${fmt(distinct)} samples, ${fmt(primaries)} primary records, stated ${fmt(S.group.stated.distinct)}; ${wrong.length} samples whose records do not add up`);
  });
  selftest('story', 'step 4: the classes of the samples add up to the stated counts', () => {
    const cnt = D.order.class.map((_, k) => M.samples.filter(s => s.cls === k).length);
    const bad = D.order.class.filter((k, i) => cnt[i] !== S.group.stated.classes[k]);
    return ok(!bad.length && sum(cnt) === S.group.stated.distinct, bad.map(k => `${k} ${cnt[D.order.class.indexOf(k)]} against ${S.group.stated.classes[k]}`).join('; '));
  });
  selftest('story', 'step 5: each identified patient has one marked sample, and the marks add up', () => {
    const st = S.group.stated, marks = M.samples.filter(s => s.marked).length;
    const bad = M.pats.filter(p => p.samples.filter(s => M.samples[s].marked).length !== 1).map(p => p.id);
    const without = M.samples.filter(s => s.pat < 0);
    return ok(!bad.length && M.pats.length === st.patients && without.length === st.without_patient && without.every(s => s.marked) &&
      marks === st.one_per_patient && marks === st.patients + st.without_patient, bad.length ? `patients ${bad.slice(0, 5).join(', ')}` : `${fmt(marks)} marked samples`);
  });
  selftest('story', 'step 6: the calls and agreements recomputed from the points equal those of the table', () => {
    const m = S.sex.matrix, P = M.sexPts, call = k => P.filter(p => p.call === k).length;
    const agree = P.filter(p => p.rec !== 'not_recorded' && p.call === p.rec).length, bad = P.filter(p => p.bad).length;
    return ok(P.length === m.samples && call('female') === m.called_female && call('male') === m.called_male && call('ambiguous') === m.ambiguous &&
      agree === m.concordant && bad === m.discordant, `${fmt(agree)} agree and ${fmt(bad)} disagree, against ${fmt(m.concordant)} and ${fmt(m.discordant)}`);
  });
  selftest('story', 'step 6: the worked examples recompute (Y genes minus XIST)', () => {
    const bad = S.sex.examples.filter(e => Math.abs(e.y - e.xist - e.d) > 0.06 || (e.d > 0) !== (e.sex === 'male'));
    return ok(!bad.length, bad.map(e => e.acc).join(', '));
  });
  selftest('story', 'step 7: every sample of the group in the bulk-blood table is a point of the start page', () => {
    const found = sum([...M.mine]), byMethod = sum(S.correction.by_study.map(r => r.n));
    return ok(found === S.correction.samples.length && byMethod === found, `${fmt(found)} found of ${fmt(S.correction.samples.length)}`);
  });
  selftest('story', 'step 8: the training and test sides add up to the samples of the task, and one group is held out', () => {
    const B = S.bench, test = M.blocks.filter(b => b.test), P = sum(M.blocks.map(b => b.pos)), Q = sum(M.blocks.map(b => b.neg));
    const r = B.held.reasons, held = B.held.pos + B.held.neg + sum(Object.values(r));
    return ok(test.length === 1 && test[0].g === B.fold && P === B.totals.pos && Q === B.totals.neg && test[0].pos === B.held.pos && test[0].neg === B.held.neg &&
      held === B.held.records && new Set(M.blocks.map(b => b.g)).size === M.blocks.length, `${fmt(P)} and ${fmt(Q)} against ${fmt(B.totals.pos)} and ${fmt(B.totals.neg)}; ${fmt(held)} records of the held-out group`);
  });
  selftest('story', 'every step fits its drawing area at both widths', () => {
    const out = [];
    [false, true].forEach(narrow => {
      const L = layout(narrow ? 380 : 820);
      const bottoms = { screen: M.rows[M.rows.length - 1].y + Math.ceil(M.rows[M.rows.length - 1].n / L.per) * L.sp, columns: L.colsBottom,
        bench: Math.max(L.bench.areas.test[1] + L.bench.testH, L.bench.areas.train[1] + L.bench.trainH),
        patients: Math.max(...M.samples.map(s => s.tok.py + L.psp)), sex: narrow ? L.sex.cardY + 70 : L.sex.y0 + L.sex.side + 40,
        'whole searches': L.w.bottom1, 'whole records': L.w.aside.rows[L.w.aside.rows.length - 1].y + 4, 'whole classes': L.w.cls.bottom,
        'twins training side': L.tw.homeBottom, 'twins test side': L.tw.test.y + L.tw.test.h };
      if (L.tw.blocks.some(b => b.w > L.tw.test.w)) out.push(`${narrow ? 'narrow' : 'wide'}: a patient-safe group is wider than the test side`);
      Object.entries(bottoms).forEach(([k, v]) => { if (!(v <= L.H)) out.push(`${narrow ? 'narrow' : 'wide'} ${k} ends at ${fmt(v)} of ${fmt(L.H)}`); });
    });
    SIM.ctx = null; if (!$('story').hidden) simDraw();     // the test moved every item; lay out again for the screen
    return ok(!out.length, out.join('; '));
  });
  // the whole resource: every step against the flow of Figure 1
  const FL = () => S.whole.flow;
  selftest('story', 'whole resource, step 1: the records of the searches add up to the candidate and included studies', () => {
    const W = M.W, F = FL();
    return ok(W.nS === F.records_found && W.colN[0] === F.records_found - F.candidate_studies && W.colN[1] === F.candidate_studies - F.included_studies &&
      W.colN[2] === F.included_studies && F.included_studies === N.studies, `${fmt(W.nS)} = ${W.colN.map(v => fmt(v)).join(' + ')}`);
  });
  selftest('story', 'whole resource, step 2: one dot per sample record of the included studies, on every included study platform', () => {
    const W = M.W, F = FL(), withValues = W.tab.reduce((t, v) => t + v, 0);
    return ok(W.n === F.included_sample_records && sum(W.units.map(u => u.n)) === W.n && W.units.length === F.included_study_platforms && withValues === F.sample_records_with_values,
      `${fmt(W.n)} records on ${fmt(W.units.length)} study platforms, ${fmt(withValues)} with values`);
  });
  selftest('story', 'whole resource, step 3: the records without values, by reason, leave the records with values of Figure 1', () => {
    const W = M.W, F = FL(), by = S.whole.no_values.map((_, k) => W.gone.filter(i => W.why[i] === k + 1).length);
    return ok(W.gone.length === F.records_without_values && W.n - W.gone.length === F.sample_records_with_values && by.every((v, k) => v === S.whole.no_values[k][2]),
      `${fmt(W.n)} − ${fmt(W.gone.length)} = ${fmt(W.n - W.gone.length)}`);
  });
  selftest('story', 'whole resource, step 4: the records with values collapse to the distinct samples of Figure 1', () => {
    const W = M.W, F = FL();
    let prims = 0, bad = 0;
    for (let i = 0; i < W.n; i++) {
      if (!W.tab[i]) continue;
      if (W.prim[i]) prims++; else if (W.target[i] < 0 || !W.prim[W.target[i]] || !W.tab[W.target[i]]) bad++;
    }
    return ok(prims === F.distinct_samples && W.reps === F.repeat_deposits && F.sample_records_with_values - W.reps === prims && !bad,
      `${fmt(F.sample_records_with_values)} − ${fmt(W.reps)} = ${fmt(F.sample_records_with_values - W.reps)}; ${bad} repeat deposits without a primary record`);
  });
  selftest('story', 'whole resource, step 5: the five classes equal Figure 1 and add up to the distinct samples', () => {
    const W = M.W, F = FL(), bad = W.classes.filter(c => c.n !== F[`class_${c.key}`]).map(c => `${c.key} ${c.n}`);
    return ok(!bad.length && sum(W.classes.map(c => c.n)) === F.distinct_samples, bad.join('; ') || `${W.classes.map(c => fmt(c.n)).join(' + ')} = ${fmt(F.distinct_samples)}`);
  });
  selftest('story', 'whole resource, step 6: the marked samples equal the identified patients plus the samples without an identifier', () => {
    const W = M.W, F = FL(), m = sum(W.classes.map(c => c.marked));
    return ok(m === F.one_sample_per_patient && m === F.identified_patients + F.samples_without_identifier, `${fmt(m)} marked`);
  });
  // twins and splits: every count against the leakage audit of the benchmark stage
  const AU = () => S.leak.audit, same4 = (c, a) => c.test === a.test && c.exposed === a.exposed && c.array === a.array && c.patient === a.patient;
  selftest('story', 'twins, step 1: the twins recomputed from the arrays and patients of the records equal those of the leakage audit', () => {
    const T = M.T, A = AU(), cnt = k => T.kind.filter(v => (v & k) === k).length, any = T.kind.filter(v => v > 0).length;
    const keyed = [0, 1, 2, 3].map(k => T.kind.filter(v => v === k).length);
    return ok(T.n === A.records && any === A.linked && cnt(1) === A.linked_array && cnt(2) === A.linked_patient && sum(keyed) === T.n,
      `${fmt(any)} of ${fmt(T.n)} records with a twin (audit ${fmt(A.linked)}); ${keyed.map(v => fmt(v)).join(' + ')} = ${fmt(sum(keyed))}`);
  });
  selftest('story', 'twins, step 2: the test records with a twin in training, recomputed for the random split, equal the audit of that split', () => {
    const T = M.T, c = T.random.c, A = AU().random;
    return ok(same4(c, A) && c.test === T.testList.length && Math.abs(c.exposed / c.test - A.share) < 1e-6,
      `${fmt(c.exposed)} of ${fmt(c.test)} (${fmt(c.array)} same array, ${fmt(c.patient)} same patient); audit ${fmt(A.exposed)} of ${fmt(A.test)}`);
  });
  selftest('story', 'twins, step 3: leave one study out, recomputed study by study, equals the audit and tests every record once', () => {
    const T = M.T, c = T.loso.c, A = AU().loso;
    return ok(same4(c, A) && c.test === T.n && T.loso.folds.length === A.folds && T.loso.mark.filter(v => v > 0).length === c.exposed,
      `${fmt(c.exposed)} of ${fmt(c.test)} over ${fmt(T.loso.folds.length)} studies; audit ${fmt(A.exposed)} of ${fmt(A.test)} over ${fmt(A.folds)}`);
  });
  selftest('story', 'twins, step 4: no test sample of a patient-safe fold has a twin in training, and the test samples add up to the audit', () => {
    const T = M.T, A = AU(), F = T.folds, tested = sum(F.map(f => f.c.test)), exposed = sum(F.map(f => f.c.exposed));
    const twinsApart = [...T.arrays.values(), ...T.patients.values()].filter(m => new Set(m.map(i => T.group[i])).size > 1).length;
    return ok(exposed === 0 && tested === A.gts.test && F.length === A.gts.folds && T.samples === A.samples && !twinsApart,
      `${fmt(tested)} test samples in ${fmt(F.length)} folds, ${fmt(exposed)} with a twin in training; ${twinsApart} sets of twins in two groups`);
  });
  selftest('story', 'twins: the number on screen at the end of each step is that of the audit', () => {
    const A = AU(), want = [fmt(A.linked), pct(A.random.exposed / A.random.test), pct(A.loso.exposed / A.loso.test), pct(0)];
    if (!SIM.ctx) setupCanvas();
    const got = ANIM.twins.draw.map(f => { SIM.hit = { x: [], y: [], r: [], f: [], rects: [] }; SIM.count = null; f(1e9); return SIM.count ? (typeof SIM.count[0] === 'string' ? SIM.count[0] : fmt(SIM.count[0])) : null; });
    if (!$('story').hidden) simDraw();
    return ok(got.every((v, k) => v === want[k]), got.map((v, k) => `${v} (${want[k]})`).join(', '));
  });
  selftest('story', 'twins: the colors of the kinds of twin and of a twin in training are set in both themes and stand out from the stage', () => {
    const before = themePref, out = [];
    try {
      ['light', 'dark'].forEach(theme => {
        applyTheme(theme, { persist: false, redraw: false }); readColors();
        const cols = [SC.twinArray, SC.twinPatient, SC.twinBoth, SC.exposed];
        if (cols.some(c => !/^#[0-9a-f]{6}$/i.test(c)) || new Set(cols).size !== cols.length) out.push(`${theme}: ${cols.join(' ')}`);
        cols.forEach(c => { if (/^#/.test(c) && Math.hypot(...hexRGB(c).map((v, k) => v - hexRGB(C.raised)[k])) < 90) out.push(`${theme}: ${c} is close to the stage`); });
      });
    } finally { applyTheme(before, { persist: false, redraw: false }); readColors(); SIM.ctx = null; }
    return ok(!out.length, out.join('; ') || 'four colors in each theme');
  });
  // Playback on a simulated clock of 100 ms frames, so that it does not depend on the browser running animation frames:
  // each animation must pass through every step in order, holding each long enough to read, and stop at its last step.
  selftest('story', 'autoplay plays every step of both animations in order, then stops at the last', async () => {
    const back = location.hash.slice(1) || 'start', raf = window.requestAnimationFrame, out = [];
    let good = true;
    try {
      await goTo('story/whole', false);
      window.requestAnimationFrame = () => 0; SIM.dry = true;
      for (const k of Object.keys(ANIM)) {
        setAnim(k); SIM.auto = true; go(0); SIM.raf = 0;
        const seen = [0], held = []; let now = performance.now(), frames = 0, since = 0;
        SIM.last = now;
        while (SIM.auto && frames < 6000) {
          now += 100; loop(now); SIM.raf = 0; frames++; since++;
          if (seen[seen.length - 1] !== SIM.scene) { seen.push(SIM.scene); held.push(since / 10); since = 0; }
        }
        const n = ANIM[k].scenes.length, inOrder = seen.length === n && seen.every((v, i) => v === i) && SIM.scene === n - 1 && !SIM.auto;
        const shortest = Math.min(...held);
        good = good && inOrder && shortest >= 7;
        out.push(`${k}: steps ${seen.map(v => v + 1).join(' ')} in ${fmt(frames / 10)} s, each shown ${fmt(shortest)} to ${fmt(Math.max(...held))} s`);
      }
    } finally {
      window.requestAnimationFrame = raf; SIM.dry = false; SIM.auto = false; SIM.running = false; setAnim('whole');
      await goTo(back, false);
    }
    return ok(good, out.join('; '));
  });
  selftest('story', 'whole resource: the number on screen at the end of each step is that of Figure 1', () => {
    const F = FL(), want = [F.included_studies, F.included_sample_records, F.sample_records_with_values, F.distinct_samples, F.distinct_samples, F.one_sample_per_patient];
    if (!SIM.ctx) setupCanvas();
    const got = ANIM.whole.draw.map(f => { SIM.hit = { x: [], y: [], r: [], f: [], rects: [] }; SIM.count = null; f(1e9); return SIM.count ? SIM.count[0] : null; });
    if (!$('story').hidden) simDraw();
    return ok(got.every((v, k) => v === want[k]), got.map((v, k) => `${fmt(v)} (${fmt(want[k])})`).join(', '));
  });
})();
