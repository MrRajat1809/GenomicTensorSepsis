/* The start page and the start of the whole page. The hero shows every sample of the corrected bulk-blood table on the
   first two principal components, moving from before the correction to after it (the points of Figure 6 A-C); then the
   main numbers, the four ways through the page and the glossary. Loaded last: it calls boot() in common.js. */
'use strict';
(() => {
  const N = D.numbers, H = D.hero, BULK = N.corrected.bulk_blood, PBMC = N.corrected.pbmc;

  // ---------- the hero: positions before and after, on one scale with equal units on both axes ----------
  const n = H.n, X0 = H.x0.map(v => v / H.scale), Y0 = H.y0.map(v => v / H.scale), X1 = H.x1.map(v => v / H.scale), Y1 = H.y1.map(v => v / H.scale);
  const SCHEMES = {
    cls: { codes: H.cls, levels: D.order.class, color: k => C.cls[k], label: k => lab('class', k) },
    fam: { codes: H.fam, levels: D.order.family, color: k => C.fam[k], label: k => lab('family', k) },
  };
  // Studies set off one after another, each as a block, so that the eye can follow them: study s starts at DELAY[s]
  // (as a share of the whole move) and each point takes the rest of the time, eased.
  const STAGGER = 0.3, DELAY = (() => { const r = rng(n); return H.studies.map(() => r() * STAGGER); })();
  const progress = (i, t) => ease(clamp((t - DELAY[H.study[i]]) / (1 - STAGGER), 0, 1));
  const posAt = (i, t) => { const p = progress(i, t); return [X0[i] + (X1[i] - X0[i]) * p, Y0[i] + (Y1[i] - Y0[i]) * p]; };
  const ext = a => [Math.min(...a), Math.max(...a)];
  const [xa, xb] = ext([...X0, ...X1]), [ya, yb] = ext([...Y0, ...Y1]);
  const PAD = 0.035, XD = [xa - PAD * (xb - xa), xb + PAD * (xb - xa)], YD = [ya - PAD * (yb - ya), yb + PAD * (yb - ya)];
  const hero = { t: 0, by: 'cls', hi: null, near: -1, timer: null, autoplay: null, geo: null, px: new Float32Array(n), py: new Float32Array(n) };

  // In full screen the plot is as tall as the screen less the controls and keys under it (FULL_ROOM), and the figure as
  // wide as that height needs for equal units on both axes; the caption below it scrolls.
  const FULL_ROOM = 170;
  function layout() {
    const cv = $('heroCanvas'), full = fullscreenEl() === $('heroFig'), maxH = full ? Math.max(280, innerHeight - FULL_ROOM) : 700;
    if (full) $('heroFig').style.setProperty('--hero-full-w', `${Math.round((maxH - 64 - 44) * (XD[1] - XD[0]) / (YD[1] - YD[0]) + 48 + 12)}px`);
    const w = cv.clientWidth || cv.parentNode.clientWidth;
    if (!w) return null;
    const small = w < 520, m = { l: small ? 42 : 48, r: 12, t: small ? 54 : 64, b: small ? 38 : 44 };
    // as tall as the plot needs for equal units on both axes
    const h = Math.round(clamp((w - m.l - m.r) * (YD[1] - YD[0]) / (XD[1] - XD[0]) + m.t + m.b, 280, maxH));
    const fit = canvasFit(cv, { height: h });
    const pw = fit.w - m.l - m.r, ph = fit.h - m.t - m.b, k = Math.min(pw / (XD[1] - XD[0]), ph / (YD[1] - YD[0]));
    const cx = m.l + pw / 2, cy = m.t + ph / 2, xm = (XD[0] + XD[1]) / 2, ym = (YD[0] + YD[1]) / 2;
    const rad = clamp(fit.w / 290, 1.5, 2.5);
    const sprites = {};
    Object.entries(SCHEMES).forEach(([s, sc]) => { sprites[s] = sc.levels.map(lv => dotSprite(sc.color(lv), rad)); });
    hero.geo = { ...fit, m, k, sx: v => cx + (v - xm) * k, sy: v => cy - (v - ym) * k, small, rad, sprites, theme: themeEpoch };
    return hero.geo;
  }
  function draw() {
    let g = hero.geo;
    if (!g || g.theme !== themeEpoch || Math.abs(($('heroCanvas').clientWidth || g.w) - g.w) > 1) g = layout();
    if (!g) return;
    const { ctx, w, h, m, sx, sy, small } = g, t = hero.t, sc = SCHEMES[hero.by];
    ctx.clearRect(0, 0, w, h);
    // grid with ticks (one scale for both states), and the axis titles of the state nearer to the slider
    ctx.font = `11px ${getComputedStyle(document.body).fontFamily}`; ctx.lineWidth = 1;
    const x0 = sx(XD[0]), x1 = sx(XD[1]), y0 = sy(YD[1]), y1 = sy(YD[0]);
    ticks(XD[0], XD[1], small ? 4 : 6).forEach(v => {
      ctx.strokeStyle = v === 0 ? C.lineStrong : C.grid; ctx.beginPath(); ctx.moveTo(Math.round(sx(v)) + 0.5, y0); ctx.lineTo(Math.round(sx(v)) + 0.5, y1); ctx.stroke();
      ctx.fillStyle = C.muted; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(fmt(v), sx(v), y1 + 5);
    });
    ticks(YD[0], YD[1], small ? 4 : 5).forEach(v => {
      ctx.strokeStyle = v === 0 ? C.lineStrong : C.grid; ctx.beginPath(); ctx.moveTo(x0, Math.round(sy(v)) + 0.5); ctx.lineTo(x1, Math.round(sy(v)) + 0.5); ctx.stroke();
      ctx.fillStyle = C.muted; ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(fmt(v), x0 - 5, sy(v));
    });
    const state = t < 0.5 ? 'before' : 'after', share = H.variance[state];
    ctx.font = `12px ${getComputedStyle(document.body).fontFamily}`; ctx.fillStyle = C.ink2; ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
    ctx.fillText(`Component 1 (${pct(share[0], 0)} of the variation)`, (x0 + x1) / 2, h - 8);
    ctx.save(); ctx.translate(small ? 3 : 8, (y0 + y1) / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top';
    ctx.fillText(`Component 2 (${pct(share[1], 0)} of the variation)`, 0, 0); ctx.restore();
    // the samples, in a fixed shuffled order; a highlighted category is drawn last, the others faded
    const hiIdx = hero.hi == null ? -1 : sc.levels.indexOf(hero.hi), sp = g.sprites[hero.by];
    for (let i = 0; i < n; i++) { const [x, y] = posAt(i, t); hero.px[i] = sx(x); hero.py[i] = sy(y); }
    if (hiIdx >= 0) ctx.globalAlpha = 0.1;
    for (let i = 0; i < n; i++) { const c = +sc.codes[i]; if (c !== hiIdx) stamp(ctx, sp[c], hero.px[i], hero.py[i]); }
    ctx.globalAlpha = 1;
    if (hiIdx >= 0) for (let i = 0; i < n; i++) if (+sc.codes[i] === hiIdx) stamp(ctx, sp[hiIdx], hero.px[i], hero.py[i]);
    if (hero.near >= 0) {
      ctx.beginPath(); ctx.arc(hero.px[hero.near], hero.py[hero.near], g.rad + 3.5, 0, 2 * Math.PI);
      ctx.strokeStyle = C.ink; ctx.lineWidth = 1.5; ctx.stroke();
    }
    renderState();
  }
  function renderState() {
    const after = hero.t >= 0.5;
    $('heroState').innerHTML = after ? '<b>After correction</b><span>Differences between studies removed</span>'
      : '<b>Before correction</b><span>Every sample on one distribution of values, not yet corrected</span>';
    $('heroEnd0').classList.toggle('on', !after); $('heroEnd1').classList.toggle('on', after);
    const slider = $('heroT'); if (+slider.value !== Math.round(hero.t * 1000)) slider.value = Math.round(hero.t * 1000);
    slider.setAttribute('aria-valuetext', hero.t <= 0 ? 'before the correction' : hero.t >= 1 ? 'after the correction' : `${fmt(100 * hero.t)}% of the way to after the correction`);
    $('heroCanvas').setAttribute('aria-label', `Scatter plot of the ${fmt(n)} samples of the corrected bulk-blood table on their first two principal components, ${after ? 'after' : 'before'} the correction, colored by ${hero.by === 'cls' ? 'class' : 'technology'}.`);
    const b = $('heroPlay'), playing = !!hero.timer;
    b.innerHTML = playing ? ICON_PAUSE : hero.t >= 1 ? ICON_REPLAY : ICON_PLAY;
    b.setAttribute('aria-label', playing ? 'Pause' : hero.t >= 1 ? 'Play the correction again' : 'Play the correction');
    b.title = b.getAttribute('aria-label');
  }

  // ---------- playing, the slider, the colors and the tooltip ----------
  const MOVE_MS = 3800;
  function stop() { if (hero.timer) cancelAnimationFrame(hero.timer); hero.timer = null; clearTimeout(hero.autoplay); hero.autoplay = null; renderState(); }
  function play() {
    stop();
    if (hero.t >= 1) hero.t = 0;
    if (reducedMotion()) { hero.t = 1; draw(); return; }
    let last = performance.now();
    const step = now => {
      if ($('start').hidden) { stop(); return; }
      hero.t = Math.min(1, hero.t + (now - last) / MOVE_MS); last = now; draw();
      if (hero.t < 1) hero.timer = requestAnimationFrame(step); else stop();
    };
    hero.timer = requestAnimationFrame(step); renderState();
  }
  function nearest(e) {
    const r = $('heroCanvas').getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
    let best = -1, bd = 64;
    for (let i = 0; i < n; i++) { const d = (hero.px[i] - x) ** 2 + (hero.py[i] - y) ** 2; if (d < bd) { bd = d; best = i; } }
    return best;
  }
  function pointTip(i) {
    const c = D.order.class[+H.cls[i]], f = D.order.family[+H.fam[i]], tis = D.order.tissue[+H.tis[i]], tier = D.order.tier[+H.tier[i]];
    const platform = H.platforms[H.platform[i]];
    return `<b>${esc(H.sample[i])}</b> <span class="m">· ${esc(H.studies[H.study[i]])}${platform ? ` · ${esc(platform)}` : ''}</span><br>` +
      `${sw(C.cls[c])}${lab('class', c)} &nbsp;${sw(C.fam[f])}${lab('family', f)}<br>` +
      `<span class="m">${lab('tissue', tis)} · ${esc(labLower('tier', tier))}</span>`;
  }
  function renderKeys() {
    const sc = SCHEMES[hero.by], counts = sc.levels.map((_, k) => 0);
    for (let i = 0; i < n; i++) counts[+sc.codes[i]]++;
    legend('heroKeys', sc.levels.map((k, i) => ({ key: k, label: sc.label(k), color: sc.color(k), count: counts[i] })).filter(x => x.count > 0),
      k => { if (k !== hero.hi) { hero.hi = k; draw(); } });
  }
  function initHero() {
    $('heroPlay').onclick = () => { if (hero.timer) stop(); else play(); };
    $('heroT').addEventListener('input', e => { stop(); hero.t = +e.target.value / 1000; draw(); });
    seg('heroBy', v => { hero.by = v; hero.hi = null; renderKeys(); draw(); });
    fullscreenButton('heroFull', 'heroFig', () => { hero.geo = null; draw(); });
    const cv = $('heroCanvas');
    const onMove = e => {
      const i = nearest(e);
      if (i !== hero.near) { hero.near = i; if (!hero.timer) draw(); }
      if (i >= 0) { tip(e, pointTip(i)); cv.style.cursor = 'pointer'; } else { untip(); cv.style.cursor = 'crosshair'; }
    };
    cv.addEventListener('pointermove', e => { if (e.pointerType === 'mouse') onMove(e); });
    cv.addEventListener('pointerdown', e => { if (e.pointerType !== 'mouse') onMove(e); });
    cv.addEventListener('mouseleave', () => { untip(); if (hero.near >= 0) { hero.near = -1; if (!hero.timer) draw(); } });
    cv.addEventListener('click', e => { if (e.pointerType && e.pointerType !== 'mouse') return; const i = nearest(e); if (i >= 0) { untip(); goTo(`studies/${H.studies[H.study[i]]}`); } });
    // for screenshots: ?morph=0 (before), ?morph=1 (after) or anything between holds the hero still; ?by=technology
    const q = new URLSearchParams(location.search), at = q.get('morph');
    if (at != null && at !== '' && isFinite(+at)) { hero.t = clamp(+at, 0, 1); hero.fixed = true; }
    if (q.get('by') === 'technology') { hero.by = 'fam'; setSeg('heroBy', 'fam'); }
  }

  // ---------- the text of the start page ----------
  const CLASS_PLURAL = { sepsis: 'sepsis', healthy_control: 'healthy controls', non_septic_control: 'non-septic controls', other_condition: 'other conditions', not_recorded: 'not recorded' };
  function renderText() {
    $('heroEyebrow').textContent = `${fmt(N.studies)} public studies · ${N.first_year}–${N.last_year}`;
    const classes = andList(D.order.class.filter(k => N.classes[k] > 0).map(k => `${fmt(N.classes[k])} ${CLASS_PLURAL[k] || labLower('class', k)}`));
    const tech = N.studies_by_technology;
    $('facts').innerHTML = `
      <div class="fact"><div class="big">${fmt(N.studies)}</div><div class="lab">studies included, of ${fmt(N.candidates)} candidate studies found by searches of GEO, ArrayExpress and the European Nucleotide Archive</div><div class="src">${fmt(tech.array)} microarray and ${fmt(tech.rna_seq)} RNA-seq studies, released ${N.first_year}–${N.last_year}</div></div>
      <div class="fact"><div class="big">${fmt(N.distinct)}</div><div class="lab">distinct samples in ${fmt(N.records)} sample records: ${classes}</div><div class="src">An array or library deposited in several studies counts once</div></div>
      <div class="fact"><div class="big">${fmt(BULK.samples)}</div><div class="lab">bulk-blood samples of ${fmt(BULK.studies)} studies on one scale after the correction, and ${fmt(PBMC.samples)} PBMC samples of ${fmt(PBMC.studies)} studies in a second table</div><div class="src">The per-study tables are kept as processed</div></div>
      <div class="fact"><div class="big">${pct(N.sex_concordance)}</div><div class="lab">of ${fmt(N.sex_compared)} samples with a recorded sex: the sex read from their genes agreed with the record</div><div class="src">${pct(N.sex_concordance_reliable)} outside the ${spell(N.sex_unreliable_tables)} tables whose recorded sex is marked as unreliable</div></div>`;
    $('pathStudies').textContent = `Any of the ${fmt(N.studies)} included studies on its own: its samples, platform, genes, checks, correction and benchmark folds, and the excluded candidates with the reason for each.`;
    // which side of the first component the healthy controls take after the correction (its sign is arbitrary)
    const side = k => { const xs = []; for (let i = 0; i < n; i++) if (D.order.class[+H.cls[i]] === k) xs.push(X1[i]); return mean(xs); };
    const left = side('healthy_control') < side('sepsis');
    const E = H.explained, f = (s, k) => fx(E[s][k][0], 2);
    const fig = H.figure.replace(/^Fig/, 'Figure ');
    $('heroCap').innerHTML = `Each dot is one of the ${fmt(n)} samples of the corrected bulk-blood table (whole blood, blood not further specified and leukocytes, from ${fmt(BULK.studies)} studies), placed by the first two principal components of its ${fmt(H.genes)} most variable genes, before the correction and after it (${fig}, panels A to C, of the paper).`;
    $('heroHow').lastElementChild.textContent = `Dots close together have similar expression profiles. Before the correction the samples fall into islands by technology and study, and the classes mix within each island. After it the technologies share one cloud, healthy controls gather on the ${left ? 'left' : 'right'} and sepsis spreads to the ${left ? 'right' : 'left'}. Hover over a dot to see its sample; click it to open its study.`;
    $('heroNote').textContent = `Technology explains ${f('before', 'family')} of the variation along the first component before the correction and ${f('after', 'family')} after it; class explains ${f('before', 'sample_class')} and ${f('after', 'sample_class')}. Studies without healthy controls were corrected by lining up the averages of their classes with those of the other studies, so part of the separation of classes comes from the correction itself; in the studies corrected on their healthy controls, class explains ${f('before', 'sample_class_anchored')} before and ${f('after', 'sample_class_anchored')} after.`;
  }

  let started = false;
  registerView('start', {
    init() { initHero(); renderText(); renderKeys(); },
    show() {
      renderKeys(); hero.geo = null; draw();
      if (!started && !hero.fixed && !reducedMotion()) { started = true; hero.autoplay = setTimeout(play, 900); }
      started = true;
    },
    redraw() { renderKeys(); hero.geo = null; draw(); },
    onKey(e) {
      if ((e.key === 'f' || e.key === 'F') && canFullscreen()) { toggleFullscreen($('heroFig')); return true; }
      return false;
    },
  });

  // ---------- checks of the core (?selftest=1) ----------
  const counts = (codes, levels) => levels.map((_, k) => [...codes].filter(c => +c === k).length);
  selftest('core', 'the class counts add up to the distinct-sample total', () => {
    const total = sum(D.order.class.map(k => N.classes[k]));
    return { ok: total === N.distinct, detail: `${fmt(total)} against ${fmt(N.distinct)}` };
  });
  selftest('core', 'the hero holds every sample and study of the corrected bulk-blood table', () => {
    const lens = [H.x0, H.y0, H.x1, H.y1, H.cls, H.fam, H.tis, H.tier, H.study, H.platform, H.sample].map(a => a.length);
    return { ok: n === BULK.samples && lens.every(l => l === n) && H.studies.length === BULK.studies && new Set(H.sample.map((s, i) => `${H.study[i]}:${s}`)).size === n,
      detail: `${fmt(n)} points, table ${fmt(BULK.samples)}; ${fmt(H.studies.length)} studies, table ${fmt(BULK.studies)}` };
  });
  selftest('core', 'the hero\'s classes, sample types and correction methods match those of the corrected table', () => {
    const cls = counts(H.cls, D.order.class), tis = counts(H.tis, D.order.tissue), tier = counts(H.tier, D.order.tier);
    const METHOD = { anchored: 'healthy_anchored', unanchored: 'matched_groups' };
    const bad = [...D.order.class.map((k, i) => [k, cls[i], BULK.by_class[k] || 0]), ...D.order.tissue.map((k, i) => [k, tis[i], BULK.by_tissue[k] || 0]),
      ...D.order.tier.map((k, i) => [k, tier[i], BULK.by_method[METHOD[k]] || 0])].filter(([, a, b]) => a !== b);
    return { ok: !bad.length, detail: bad.map(([k, a, b]) => `${k} ${a} against ${b}`).join('; ') };
  });
  selftest('core', 'the hero starts at the positions before the correction and ends at those after it', () => {
    let worst = 0;
    for (let i = 0; i < n; i++) { const a = posAt(i, 0), b = posAt(i, 1); worst = Math.max(worst, Math.abs(a[0] - X0[i]), Math.abs(a[1] - Y0[i]), Math.abs(b[0] - X1[i]), Math.abs(b[1] - Y1[i])); }
    return { ok: worst < 1e-9, detail: `largest difference ${worst}` };
  });
  selftest('core', 'the page\'s data colors are those of the figures (light theme)', () => {
    const before = themePref; applyTheme('light', { persist: false, redraw: false });
    // "not recorded" is darker on the page than in the figures, so that it shows on the paper background (index.html)
    const bad = [];
    [['cls', 'class'], ['fam', 'family'], ['tier', 'tier']].forEach(([c, k]) => Object.entries(D.colors[k]).forEach(([key, hex]) => {
      if (!(c === 'cls' && key === 'not_recorded') && C[c][key].toLowerCase() !== hex.toLowerCase()) bad.push(`${key} ${C[c][key]} against ${hex}`);
    }));
    applyTheme(before, { persist: false, redraw: false });
    return { ok: !bad.length, detail: bad.join('; ') };
  });
  selftest('core', 'the facts row shows the main numbers', () => {
    const txt = $('facts').textContent, want = [fmt(N.studies), fmt(N.candidates), fmt(N.distinct), fmt(N.records), fmt(BULK.samples), pct(N.sex_concordance)];
    const missing = want.filter(w => !txt.includes(w));
    return { ok: !missing.length, detail: missing.length ? `missing ${missing.join(', ')}` : '' };
  });
  selftest('core', 'each animation has a full-screen button, and what it fills holds the animation and its controls', () => {
    const sets = [['heroFull', 'heroFig', 'heroCanvas', 'heroPlay', 'heroT'], ['storyFull', 'storyGrid', 'storyCanvas', 'storyPlay', 'storyScrub', 'storyCap']];
    const bad = sets.filter(([b, t, ...inside]) => !$(b) || !$(t) || ![b, ...inside].every(id => $(t).contains($(id))) || (canFullscreen() && !$(b).getAttribute('aria-label')));
    return { ok: !bad.length, detail: bad.length ? `wrong: ${bad.map(s => s[0]).join(', ')}` : `full screen ${canFullscreen() ? 'available' : 'not available in this browser'}` };
  });
  selftest('core', 'the glossary holds every term of the loaded parts once, each with a definition', () => {
    const e = glossaryEntries(), terms = new Set(e.map(x => x[0].toLowerCase()));
    const all = Object.keys(PARTS).filter(loaded).flatMap(k => (PARTS[k].glossary || []).map(x => x[0].toLowerCase()));
    const shown = $('glossStart').querySelectorAll('dt').length;
    return { ok: terms.size === e.length && e.every(x => x[1] && x[1].length > 20) && all.every(t => terms.has(t)) && shown === e.length,
      detail: `${fmt(e.length)} terms, ${fmt(shown)} shown` };
  });

  boot();
})();
