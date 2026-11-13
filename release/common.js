/* Shared helpers of the results page: data parts, formatting, colors, charts, keys, controls, the tooltip, views and
   chapters with their addresses, the glossary and the self-test. Loaded after data.js and before the scripts of the views.

   For the writer of a view (story.js, resource.js, studies.js, methods.js):
   - registerView(name, {part, init, show, redraw, resolve, hashFor, groupLabel, navName, routes, onKey}). The part (or a
     list of parts) is loaded before init() runs, once; show(route) runs every time the view or one of its chapters is
     shown, and again after a change of theme or window width unless redraw() is given. route = {view, chapter, args}.
   - A book view has <section class="sec" id="..." data-nav data-sub [data-group]> chapters; #view/<id> opens one.
     resolve(args) may turn other addresses into a chapter id (#studies/GSE65682/checks), and hashFor(section) builds
     them back; groupLabel(group) names a group of chapters in the rail and navName(section) a chapter in the pager.
     Call refreshBook(view) after adding or renaming chapters.
   - routes() may list the addresses that the self-test opens (by default every chapter); onKey(event) returns true when
     the view handled a key (by default the arrow keys turn the chapters of a book).
   - need('methods') loads data/methods.js and resolves when window.DATA_PARTS.methods is set; part('methods') returns it.
     need('studies/GSE65682') loads data/studies/GSE65682.js.
   - A part may hold a glossary: [[term, definition], ...]; its entries join the glossary once the part is loaded.
   - selftest(part, name, fn): fn returns true, false or {ok, detail}, or throws; it runs with ?selftest=1 after every
     part has loaded and every view and chapter has been drawn in both themes.
   - fullscreenButton(buttonId, targetId, onChange) lets a button put an animation with its controls in full screen;
     onKey may call toggleFullscreen(element) for the F key.
   Wrap a view's own code in a function, so that its names do not collide with those of the other scripts. */
'use strict';

// ---------- errors: shown in the box under the header, and counted by the self-test ----------
const ERRORS = [];
function showError(e) {
  const msg = String((e && e.message) || e);
  ERRORS.push(msg);
  const b = document.getElementById('err');
  if (b) { b.style.display = 'block'; b.textContent += `Error: ${msg}\n`; }
}
window.addEventListener('error', e => showError(`${e.message} (${String(e.filename || '').split('/').pop()}:${e.lineno})`));
window.addEventListener('unhandledrejection', e => showError(e.reason || 'a promise was rejected'));

// ---------- data parts: data.js is the core; each view loads its own part when it is first opened ----------
// Parts are loaded with a script tag rather than fetch(), which browsers block for pages opened from files.
window.DATA_PARTS = window.DATA_PARTS || {};
const PARTS = window.DATA_PARTS;
// Records that share their fields are written once ({__f: fields, __r: [values]}, pack() in build_data.py); unpack them
// into ordinary objects when a part is loaded.
function unpack(x) {
  if (Array.isArray(x)) return x.map(unpack);
  if (x && typeof x === 'object') {
    if (x.__f && x.__r) return x.__r.map(v => Object.fromEntries(x.__f.map((f, i) => [f, unpack(v[i])])));
    for (const k of Object.keys(x)) x[k] = unpack(x[k]);
  }
  return x;
}
const UNPACKED = new Set(), PENDING = {};
let BOOTED = false;
function adopt(name) {
  if (!UNPACKED.has(name)) { PARTS[name] = unpack(PARTS[name]); UNPACKED.add(name); if (BOOTED && PARTS[name] && PARTS[name].glossary) renderGlossary(); }
  return PARTS[name];
}
if (!PARTS.core) showError('data.js did not load; keep it next to index.html');
const D = adopt('core') || {};
// the file of a part, with its content version when build_data.py has stamped the page
function partFile(name) {
  const f = name === 'core' ? 'data.js' : name === 'studies' ? 'data/studies/index.js' : `data/${name}.js`;
  const v = (window.PART_VERSIONS || {})[f];
  return v ? `${f}?v=${v}` : f;
}
function loadPart(name) {
  if (PARTS[name]) return Promise.resolve(adopt(name));
  if (!PENDING[name]) PENDING[name] = new Promise((resolve, reject) => {
    const s = document.createElement('script'), file = partFile(name).split('?')[0];
    s.src = partFile(name);
    s.onload = () => PARTS[name] ? resolve(adopt(name)) : reject(new Error(`${file} did not set DATA_PARTS["${name}"]`));
    s.onerror = () => { delete PENDING[name]; s.remove(); reject(new Error(`Could not load ${file}; keep the data folder next to index.html`)); };
    document.head.append(s);
  });
  return PENDING[name];
}
const need = (...names) => Promise.all(names.flat().map(loadPart));
const loaded = name => UNPACKED.has(name);
const part = name => loaded(name) ? PARTS[name] : undefined;

// ---------- formatting, in American style: 12,126 and 0.95 ----------
const $ = id => document.getElementById(id);
const MINUS = '−';
const NF = {};
const nf = d => NF[d] || (NF[d] = new Intl.NumberFormat('en-US', { minimumFractionDigits: d, maximumFractionDigits: d }));
// a number with d decimals and thousands separators; a true minus sign, but none on a value that rounds to zero
function fmt(v, d = 0) {
  if (typeof v === 'string' && v.trim() !== '' && isFinite(+v)) v = +v;
  if (v == null || typeof v !== 'number' || !isFinite(v)) return '–';
  const s = nf(d).format(Math.abs(v));
  return (v < 0 && /[1-9]/.test(s) ? MINUS : '') + s;
}
const fx = (v, d = 2) => fmt(v, d);
const pct = (share, d = 1) => share == null || !isFinite(share) ? '–' : fmt(100 * share, d) + '%';
const sgn = (v, d = 0) => (v > 0 && /[1-9]/.test(fmt(v, d)) ? '+' : '') + fmt(v, d);
// a number as stored, with all its decimals (0.988 stays 0.988)
const num = v => fmt(v, (String(v).split('.')[1] || '').length);
const WORDS = ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine'];
// small whole numbers spelled out in running text, as American style has it; Spell for the start of a sentence
const spell = n => Number.isInteger(n) && n >= 0 && n < 10 ? WORDS[n] : fmt(n);
const Spell = n => { const s = spell(n); return s[0].toUpperCase() + s.slice(1); };
// 'a, b and c'
const andList = items => items.length < 2 ? items.join('') : `${items.slice(0, -1).join(', ')} and ${items[items.length - 1]}`;
const plural = (n, one, many) => `${fmt(n)} ${n === 1 ? one : many || one + 's'}`;
const MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'];
// '2026-10-11' -> 'October 11, 2026'
const fdate = iso => { const m = String(iso || '').match(/^(\d{4})-(\d{2})(?:-(\d{2}))?/); return m ? (m[3] ? `${MONTHS[+m[2] - 1]} ${+m[3]}, ${m[1]}` : `${MONTHS[+m[2] - 1]} ${m[1]}`) : '–'; };
const esc = s => String(s == null ? '' : s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);
// an en dash in displayed ranges ('5-49 years' is a data key, '5–49 years' its label)
const dash = s => typeof s === 'string' ? s.replace(/(\d)-(?=\d)/g, '$1–') : s;
const ordinal = n => n + (n % 100 >= 11 && n % 100 <= 13 ? 'th' : ['th', 'st', 'nd', 'rd'][n % 10] || 'th');
const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
const sum = a => a.reduce((t, v) => t + v, 0);
// labels of the data keys: lab('class', 'sepsis') -> 'Sepsis'; lab('tissue', 'pbmc') -> 'PBMC'
const lab = (kind, key) => ((D.labels || {})[kind] || {})[key] || String(key).replace(/_/g, ' ');
const labLower = (kind, key) => { const s = lab(kind, key); return /^[A-Z][a-z]/.test(s) ? s[0].toLowerCase() + s.slice(1) : s; };

// ---------- colors: read from the CSS tokens of the current theme ----------
// C.bg, C.ink, C.accent ... and the data colors C.cls[class], C.fam[technology], C.tier[anchored|unanchored].
const C = { cls: {}, fam: {}, tier: {} };
function hexRGB(h) { h = String(h).trim(); if (h.startsWith('rgb')) return h.match(/[\d.]+/g).slice(0, 3).map(Number); if (h.length === 4) h = '#' + [...h.slice(1)].map(c => c + c).join(''); return [1, 3, 5].map(k => parseInt(h.slice(k, k + 2), 16)); }
function toHex(rgb) { return '#' + rgb.map(v => Math.round(clamp(v, 0, 255)).toString(16).padStart(2, '0')).join(''); }
function mixHex(a, b, f) { const A = hexRGB(a), B = hexRGB(b); return toHex(A.map((v, i) => v + (B[i] - v) * f)); }
const tint = (hex, amount) => mixHex(hex, C.bg, amount);
const rgba = (hex, a) => `rgba(${hexRGB(hex).join(',')},${a})`;
function readTheme() {
  const cs = getComputedStyle(document.documentElement), v = n => cs.getPropertyValue(`--${n}`).trim();
  ['bg', 'surface', 'raised', 'ink', 'ink2', 'muted', 'line', 'grid', 'band', 'accent', 'mark', 'hilite'].forEach(n => { C[n] = v(n); });
  C.lineStrong = v('line-strong'); C.accentInk = v('accent-ink');
  const order = D.order || {};
  (order.class || []).forEach(k => { C.cls[k] = v(`cls-${k}`); });
  (order.family || []).forEach(k => { C.fam[k] = v(`fam-${k}`); });
  (order.tier || []).forEach(k => { C.tier[k] = v(`tier-${k}`); });
  C.dark = sum(hexRGB(C.bg)) < 300;
}

// ---------- SVG charts ----------
const NS = 'http://www.w3.org/2000/svg';
function S(tag, attrs, parent) {
  const n = document.createElementNS(NS, tag);
  for (const k in attrs || {}) { if (k === 'text') n.textContent = attrs[k]; else if (attrs[k] != null) n.setAttribute(k, attrs[k]); }
  if (parent) parent.appendChild(n);
  return n;
}
// an empty SVG of viewBox w x h in the element with this id (or the element itself); label is its aria-label
function svg(id, w, h, label) { const box = typeof id === 'string' ? $(id) : id; box.innerHTML = ''; const s = S('svg', { viewBox: `0 0 ${w} ${h}`, role: 'img', 'aria-label': label || '' }); box.appendChild(s); return s; }
// scales: lin(d0, d1, r0, r1) maps a value to a position (f.inv the reverse); lg the same on a log10 scale
const lin = (d0, d1, r0, r1) => { const f = v => r0 + (v - d0) / (d1 - d0 || 1) * (r1 - r0); f.inv = p => d0 + (p - r0) / (r1 - r0) * (d1 - d0); return f; };
const lg = (d0, d1, r0, r1) => { const a = Math.log10(d0), b = Math.log10(d1); const f = v => r0 + (Math.log10(Math.max(v, d0 * 1e-3)) - a) / (b - a) * (r1 - r0); f.inv = p => 10 ** (a + (p - r0) / (r1 - r0) * (b - a)); return f; };
// round tick values between lo and hi, about n of them
function ticks(lo, hi, n = 5) {
  const span = hi - lo; if (!(span > 0)) return [lo];
  const s0 = Math.pow(10, Math.floor(Math.log10(span / n))), e = span / n / s0;
  const step = s0 * (e >= 7.5 ? 10 : e >= 3.5 ? 5 : e >= 1.5 ? 2 : 1);
  const t = []; for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9 * step; v += step) t.push(+v.toFixed(10));
  return t;
}
function niceMax(v) { const t = ticks(0, v * 1.0001, 5); return t[t.length - 1] < v ? t[t.length - 1] + (t[1] - t[0]) : t[t.length - 1]; }
// a text in an SVG; a number keeps the decimals it has, so pass a formatted string where rounding is wanted
const T = (g, x, y, text, attrs = {}) => S('text', { x, y, 'font-size': 12, fill: C.ink2, text: typeof text === 'number' ? num(text) : dash(text), ...attrs }, g);
function gridY(g, y, tks, x0, x1, f = fmt) {
  tks.forEach(t => { S('line', { x1: x0, x2: x1, y1: y(t), y2: y(t), stroke: C.grid }, g); T(g, x0 - 6, y(t) + 4, f(t), { 'text-anchor': 'end', 'font-size': 11, fill: C.muted }); });
}
function gridX(g, x, tks, y0, y1, f = fmt) {
  tks.forEach(t => { S('line', { x1: x(t), x2: x(t), y1: y0, y2: y1, stroke: C.grid }, g); T(g, x(t), y1 + 15, f(t), { 'text-anchor': 'middle', 'font-size': 11, fill: C.muted }); });
}
// axis titles: below the x axis, or rotated along the y axis
function axisTitle(g, x, y, text, vertical = false) { return T(g, x, y, text, { 'text-anchor': 'middle', 'font-size': 12, fill: C.muted, transform: vertical ? `rotate(-90 ${x} ${y})` : null }); }
function rrect(g, x, y, w, h, fill, r = 3, extra = {}) { return S('rect', { x, y, width: Math.max(0, w), height: Math.max(0, h), rx: Math.min(r, Math.abs(w) / 2, Math.abs(h) / 2), fill, ...extra }, g); }
function rowBands(g, n, y0, rowH, x0, x1) { for (let i = 1; i < n; i += 2) S('rect', { x: x0, y: y0 + i * rowH, width: x1 - x0, height: rowH, fill: C.band }, g); }

// ---------- canvas ----------
// Size a canvas to its CSS width and the given height (or aspect ratio, width / height), sharp on high-density
// screens; returns {ctx, w, h} in CSS pixels.
function canvasFit(cv, opts = {}) {
  const w = Math.max(1, cv.clientWidth || cv.parentNode.clientWidth), h = Math.round(opts.height || w / (opts.aspect || 1.6)), r = window.devicePixelRatio || 1;
  if (cv.width !== Math.round(w * r) || cv.height !== Math.round(h * r)) { cv.width = Math.round(w * r); cv.height = Math.round(h * r); }
  cv.style.height = `${h}px`;
  const ctx = cv.getContext('2d'); ctx.setTransform(r, 0, 0, r, 0, 0);
  return { ctx, w, h, r };
}
// A dot drawn once and stamped many times (much faster than arc() for thousands of points): fill at the given opacity
// with a thin rim of the same hue, darker on a light page and lighter on a dark one. Draw it centered with dot.r.
function dotSprite(color, radius, alpha = 0.72) {
  const r = window.devicePixelRatio || 1, pad = 1, size = Math.ceil((radius + pad) * 2 * r);
  const cv = document.createElement('canvas'); cv.width = cv.height = size;
  const ctx = cv.getContext('2d'); ctx.scale(r, r);
  const c = (radius + pad);
  ctx.beginPath(); ctx.arc(c, c, radius, 0, 2 * Math.PI);
  ctx.fillStyle = rgba(color, alpha); ctx.fill();
  ctx.lineWidth = Math.max(0.5, radius / 4); ctx.strokeStyle = rgba(mixHex(color, C.dark ? '#ffffff' : '#000000', C.dark ? 0.25 : 0.3), Math.min(1, alpha + 0.15)); ctx.stroke();
  cv.r = c; cv.css = 2 * c;
  return cv;
}
function stamp(ctx, sprite, x, y) { ctx.drawImage(sprite, x - sprite.r, y - sprite.r, sprite.css, sprite.css); }

// ---------- keys and legends ----------
// keys([[label, color, kind]]) -> HTML; kind: '' (square), 'dot', 'line' or 'ring'
function keys(items) { return items.map(([label, col, kind]) => `<span class="key"><i class="${kind || ''}" style="${kind === 'line' || kind === 'ring' ? `border-color:${col}` : `background:${col}`}"></i>${dash(label)}</span>`).join(''); }
// Keys that pick out their category when hovered or focused: items = [{key, label, color, count}]; onHover(key or null)
function legend(el, items, onHover, kind = 'dot') {
  el = typeof el === 'string' ? $(el) : el;
  el.innerHTML = items.map(it => `<button type="button" class="key" data-k="${esc(it.key)}"><i class="${kind}" style="background:${it.color}"></i>${esc(it.label)}${it.count != null ? ` <b>${fmt(it.count)}</b>` : ''}</button>`).join('');
  if (onHover && !el.dataset.wired) {
    el.dataset.wired = '1';
    const on = e => { const b = e.target.closest('[data-k]'); el.__hover && el.__hover(b ? b.dataset.k : null); };
    el.addEventListener('mouseover', on); el.addEventListener('focusin', on);
    el.addEventListener('mouseleave', () => el.__hover && el.__hover(null)); el.addEventListener('focusout', () => el.__hover && el.__hover(null));
  }
  el.__hover = onHover;
}

// ---------- controls ----------
// segmented control: <div class="seg" id><button data-v>...; returns a function that reads the chosen value
function seg(id, onChange) {
  const box = typeof id === 'string' ? $(id) : id;
  box.addEventListener('click', e => { const b = e.target.closest('button[data-v]'); if (!b || b.classList.contains('on')) return; setSeg(box, b.dataset.v); onChange(b.dataset.v); });
  return () => (box.querySelector('button.on') || {}).dataset.v;
}
function setSeg(id, v) { (typeof id === 'string' ? $(id) : id).querySelectorAll('button[data-v]').forEach(x => { const on = x.dataset.v === String(v); x.classList.toggle('on', on); x.setAttribute('aria-pressed', on); }); }
// the buttons of a segmented control: [[value, label]]
const segHTML = (options, value) => options.map(([v, label]) => `<button type="button" data-v="${esc(v)}" class="${String(v) === String(value) ? 'on' : ''}" aria-pressed="${String(v) === String(value)}">${label}</button>`).join('');
const ICON_PLAY = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M7 5l12 7-12 7z"/></svg>';
const ICON_PAUSE = '<svg viewBox="0 0 24 24" aria-hidden="true"><rect x="6" y="5" width="4" height="14" rx="1"/><rect x="14" y="5" width="4" height="14" rx="1"/></svg>';
const ICON_REPLAY = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M12 5V2L7 6l5 4V7a5 5 0 1 1-5 5H5a7 7 0 1 0 7-7z"/></svg>';
const reducedMotion = () => matchMedia('(prefers-reduced-motion: reduce)').matches;

// ---------- full screen: an animation with its controls fills the screen; the button, F or Escape leaves ----------
const ICON_FULL = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M7 14H5v5h5v-2H7zm-2-4h2V7h3V5H5zm12 7h-3v2h5v-5h-2zM14 5v2h3v3h2V5z"/></svg>';
const ICON_FULL_EXIT = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M5 16h3v3h2v-5H5zm3-8H5v2h5V5H8zm6 11h2v-3h3v-2h-5zm2-11V5h-2v5h5V8z"/></svg>';
const fullscreenEl = () => document.fullscreenElement || document.webkitFullscreenElement || null;
const canFullscreen = () => !!(document.fullscreenEnabled || document.webkitFullscreenEnabled);
function leaveFullscreen() { if (fullscreenEl()) Promise.resolve((document.exitFullscreen || document.webkitExitFullscreen).call(document)).catch(() => {}); }
function toggleFullscreen(target) {
  if (fullscreenEl() === target) { leaveFullscreen(); return; }
  Promise.resolve((target.requestFullscreen || target.webkitRequestFullscreen).call(target)).catch(e => showError(`full screen was refused: ${e.message || e}`));
}
// A button that puts the element targetId in full screen and takes it out again; onChange(on) runs after every change.
// Where the browser cannot show an element in full screen (Safari on iPhone), the button is hidden.
function fullscreenButton(buttonId, targetId, onChange) {
  const b = $(buttonId), target = $(targetId);
  if (!canFullscreen()) { b.hidden = true; return; }
  const render = () => {
    const on = fullscreenEl() === target, label = on ? 'Leave full screen (F or Esc)' : 'Full screen (F)';
    b.innerHTML = on ? ICON_FULL_EXIT : ICON_FULL; b.setAttribute('aria-label', label); b.title = label; b.setAttribute('aria-pressed', on);
  };
  b.onclick = () => toggleFullscreen(target);
  ['fullscreenchange', 'webkitfullscreenchange'].forEach(ev => document.addEventListener(ev, () => {
    render(); untip();
    if (onChange) requestAnimationFrame(() => onChange(fullscreenEl() === target));
  }));
  render();
}

// ---------- tooltip ----------
const tipEl = () => $('tip');
// In full screen only the element on screen is drawn, so the tooltip moves into it while it lasts.
function tip(e, html) { const t = tipEl(), host = fullscreenEl() || document.body; if (t.parentNode !== host) host.appendChild(t); t.innerHTML = html; t.style.opacity = 1; moveTip(e); }
function moveTip(e) {
  const t = tipEl(), r = t.getBoundingClientRect(); let x = e.clientX + 14, y = e.clientY + 14;
  if (x + r.width > innerWidth - 8) x = e.clientX - r.width - 14; if (y + r.height > innerHeight - 8) y = e.clientY - r.height - 14;
  t.style.transform = `translate(${Math.max(4, x)}px,${Math.max(4, y)}px)`;
}
function untip() { tipEl().style.opacity = 0; }
// a tooltip on an element (html: text or a function that returns it)
function hover(node, html) { node.addEventListener('mousemove', e => tip(e, typeof html === 'function' ? html() : html)); node.addEventListener('mouseleave', untip); return node; }
function hitRow(g, x, y, w, h, html) { return hover(S('rect', { x, y, width: w, height: h, fill: 'transparent' }, g), html); }
// a colored dot for a tooltip line
const sw = col => `<i class="sw" style="background:${col}"></i>`;

// ---------- animation ----------
const ease = p => p < 0.5 ? 4 * p * p * p : 1 - Math.pow(-2 * p + 2, 3) / 2;
// count a number up or down in an element's text
function tween(el, to, f = fmt, ms = 450) {
  const from = el.dataset.v == null ? to : +el.dataset.v, t0 = performance.now(); el.dataset.v = to;
  if (reducedMotion()) ms = 1;
  const step = t => { const p = Math.min(1, (t - t0) / ms), e = 1 - Math.pow(1 - p, 3); el.textContent = f(from + (to - from) * e); if (p < 1) requestAnimationFrame(step); };
  requestAnimationFrame(step);
}
// animate numeric SVG attributes of many nodes at once: items = [[node, {attr: [from, to]}]]; returns a stop function
function animate(items, ms, done) {
  const t0 = performance.now(); let id;
  if (reducedMotion()) ms = 1;
  const step = t => {
    const p = Math.min(1, (t - t0) / ms), e = ease(p);
    items.forEach(([n, at]) => { for (const a in at) n.setAttribute(a, at[a][0] + (at[a][1] - at[a][0]) * e); });
    if (p < 1) id = requestAnimationFrame(step); else if (done) done();
  };
  id = requestAnimationFrame(step);
  return () => cancelAnimationFrame(id);
}

// ---------- statistics ----------
function rng(seed) { let a = seed >>> 0; return () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
function gauss(r) { let u = 0; while (u === 0) u = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * r()); }
const mean = a => a.length ? sum(a) / a.length : NaN;
const sd = a => { const m = mean(a); return Math.sqrt(sum(a.map(v => (v - m) ** 2)) / (a.length - 1)); };
const quant = (sorted, q) => { if (!sorted.length) return NaN; const p = (sorted.length - 1) * q, i = Math.floor(p); return sorted[i] + (sorted[Math.min(i + 1, sorted.length - 1)] - sorted[i]) * (p - i); };
const median = a => quant(a.slice().sort((x, y) => x - y), 0.5);
// AUROC: the probability that a positive sample scores above a negative one, ties counting one half (Mann-Whitney)
function auroc(pos, neg) {
  const all = pos.map(v => [v, 1]).concat(neg.map(v => [v, 0])).sort((a, b) => a[0] - b[0]);
  let rankSum = 0;
  for (let i = 0; i < all.length;) { let j = i; while (j < all.length && all[j][0] === all[i][0]) j++; const r = (i + j + 1) / 2; for (let k = i; k < j; k++) if (all[k][1]) rankSum += r; i = j; }
  return (rankSum - pos.length * (pos.length + 1) / 2) / (pos.length * neg.length);
}

// ---------- views, chapters and addresses ----------
// #start, #story, #resource/<chapter>, #studies/<study>/<chapter>, #methods/<chapter>, or the id of any element inside a
// view (#glossary). The back and forward buttons of the browser move between them.
const VIEWS = ['start', 'story', 'resource', 'studies', 'methods'];
const VIEW_NAME = { start: 'Start', story: 'How it works', resource: 'The resource', studies: 'Studies', methods: 'Method in detail' };
// where the last chapter of a book leads
const BOOK_NEXT = { resource: ['studies', 'Studies'], studies: ['methods', 'Method in detail'], methods: ['start', 'Back to the start'] };
const SPECS = {};
const ROUTE = { view: 'start', chapter: null, args: [] };
const CHAP = {};
function registerView(name, spec) {
  if (!VIEWS.includes(name)) throw new Error(`registerView: unknown view ${name}`);
  SPECS[name] = Object.assign({}, spec, { ready: false, loading: null });
}
const isBook = v => !!($(v) && $(v).querySelector('.book'));
const secs = v => [...$(v).querySelectorAll('.pages > section.sec')];
const groupOf = s => s.dataset.group || 'main';
const groupLabel = (v, g) => { const f = (SPECS[v] || {}).groupLabel; const l = f ? f(g) : undefined; return l === undefined ? VIEW_NAME[v] : l; };
const navName = (v, s) => { const f = (SPECS[v] || {}).navName; return (f && f(s)) || s.dataset.nav; };
function hashFor(v, sec) { const f = (SPECS[v] || {}).hashFor; return (sec && f && f(sec)) || (sec ? `${v}/${sec.id}` : v); }
const setHash = h => { if (location.hash.slice(1) !== h) history.pushState(null, '', '#' + h); };

// the rail, the chapter numbers and the line above each chapter's heading ("The resource · 2 of 8")
function refreshBook(v) {
  const list = secs(v); let last = null;
  list.forEach(s => {
    const g = groupOf(s), gl = list.filter(x => groupOf(x) === g), n = s.querySelector('.sec-h .n');
    if (n && !n.dataset.fixed) { if (n.textContent.trim()) n.dataset.fixed = '1'; else n.textContent = `${gl.indexOf(s) + 1}.`; }
    let pos = s.querySelector('.chap-pos');
    if (!pos) { s.insertAdjacentHTML('afterbegin', '<p class="eyebrow chap-pos"></p>'); pos = s.querySelector('.chap-pos'); }
    const name = groupLabel(v, g) || VIEW_NAME[v];
    pos.textContent = gl.length > 1 ? `${name} · ${gl.indexOf(s) + 1} of ${gl.length}` : name;
  });
  const groups = new Set(list.map(groupOf));
  $(v).querySelector('.rail ol').innerHTML = list.map((s, i) => {
    const g = groupOf(s), name = groupLabel(v, g);
    const head = groups.size > 1 && g !== last && name ? `<li class="grp">${name}</li>` : ''; last = g;
    const n = ((s.querySelector('.sec-h .n') || {}).textContent || '').replace(/\.$/, '');
    return head + `<li><button type="button" data-i="${i}" aria-current="${s.classList.contains('on')}"><span class="n">${n}</span><span class="t">${navName(v, s)}</span><span class="s">${s.dataset.sub || ''}</span></button></li>`;
  }).join('');
  const on = list.findIndex(s => s.classList.contains('on'));
  if (on >= 0) renderPager(v, on);
}
function renderPager(v, i) {
  const list = secs(v), prev = list[i - 1], next = list[i + 1], out = BOOK_NEXT[v] || ['start', 'Back to the start'];
  $(v).querySelector('.pager').innerHTML =
    (prev ? `<button type="button" class="prev" data-i="${i - 1}"><small>← Previous</small><b>${navName(v, prev)}</b></button>` : '') +
    (next ? `<button type="button" class="next" data-i="${i + 1}"><small>Next →</small><b>${navName(v, next)}</b></button>`
      : `<button type="button" class="next" data-go="${out[0]}"><small>Next →</small><b>${out[1]}</b></button>`);
}
function activateChapter(v, i) {
  const list = secs(v); i = clamp(i, 0, list.length - 1); CHAP[v] = i;
  list.forEach((s, k) => s.classList.toggle('on', k === i));
  const rail = $(v).querySelector('.rail'), btns = [...rail.querySelectorAll('button[data-i]')];
  btns.forEach(b => b.setAttribute('aria-current', +b.dataset.i === i));
  const cur = btns.find(b => +b.dataset.i === i);
  if (cur && rail.scrollWidth > rail.clientWidth + 4) rail.scrollLeft = cur.offsetLeft - rail.clientWidth / 2 + cur.clientWidth / 2;
  renderPager(v, i);
  return list[i];
}
function showViewFrame(v) {
  VIEWS.forEach(t => { const on = t === v, b = $(`tab-${t}`); $(t).hidden = !on; b.setAttribute('aria-selected', on); b.tabIndex = on ? 0 : -1; });
}
// load a view's part and run its init() once; a book gets its rail afterwards, since init() may add chapters
function ensureReady(v) {
  const spec = SPECS[v];
  if (!spec) { if (isBook(v) && !$(v).dataset.built) { refreshBook(v); $(v).dataset.built = '1'; } return Promise.resolve(); }
  if (spec.ready) return Promise.resolve();
  if (!spec.loading) spec.loading = (async () => {
    $(v).classList.add('busy');
    try {
      if (spec.part) await need(spec.part);
      if (spec.init) await spec.init();
      if (isBook(v)) refreshBook(v);
      spec.ready = true;
    } finally { $(v).classList.remove('busy'); }
  })().catch(e => { spec.loading = null; throw e; });
  return spec.loading;
}
let navToken = 0;
// go to an address; push=false when the address is already in the location bar (start, back and forward)
async function goTo(h, push = true) {
  const token = ++navToken;
  h = decodeURIComponent(String(h || '').replace(/^#/, ''));
  let [v, ...args] = h.split('/'), target = null;
  if (!VIEWS.includes(v)) {
    const el = h && $(h), main = el && el.closest('main.view');
    if (main) { v = main.id; target = el; } else v = 'start';
    args = [];
  }
  // a full-screen animation ends when the address leaves its view
  if (fullscreenEl() && !$(v).contains(fullscreenEl())) leaveFullscreen();
  showViewFrame(v);
  try { await ensureReady(v); } catch (e) { showError(e); return; }
  if (token !== navToken) return;
  const spec = SPECS[v] || {};
  let sec = null;
  if (isBook(v)) {
    let id = target && target.closest('section.sec') ? target.closest('section.sec').id : null;
    if (!id && spec.resolve) { try { id = await spec.resolve(args); } catch (e) { showError(e); } if (token !== navToken) return; }
    if (!id && args.length) id = args[0];
    const list = secs(v), i = list.findIndex(s => s.id === id);
    sec = activateChapter(v, i >= 0 ? i : (CHAP[v] || 0));
  }
  Object.assign(ROUTE, { view: v, chapter: sec ? sec.id : null, args });
  if (push) setHash(target && target !== sec ? h : hashFor(v, sec));
  if (spec.show) { try { await spec.show({ ...ROUTE }); } catch (e) { showError(e); } }
  if (token !== navToken) return;
  if (target && target !== sec && !target.matches('main.view')) target.scrollIntoView(); else window.scrollTo(0, 0);
}
const goChapter = (v, i) => goTo(hashFor(v, secs(v)[clamp(i, 0, secs(v).length - 1)]));
// redraw the view on screen (after a change of theme or of the window's width)
let themeEpoch = 0;
function redrawView() {
  const spec = SPECS[ROUTE.view] || {};
  try { if (spec.redraw) spec.redraw(); else if (spec.show && spec.ready) spec.show({ ...ROUTE }); } catch (e) { showError(e); }
}

// ---------- theme: system, light or dark; the choice is remembered on this computer ----------
const THEME_KEY = 'gts-theme';
const THEME_ICONS = {
  system: '<svg viewBox="0 0 24 24"><rect x="3" y="4" width="18" height="12" rx="2"/><path d="M8 20h8M12 16v4"/></svg>',
  light: '<svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/></svg>',
  dark: '<svg viewBox="0 0 24 24"><path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/></svg>'
};
const THEME_NEXT = { system: 'light', light: 'dark', dark: 'system' };
let themePref = 'system';
function applyTheme(pref, { persist = true, redraw = true } = {}) {
  themePref = pref;
  if (pref === 'system') delete document.documentElement.dataset.theme; else document.documentElement.dataset.theme = pref;
  if (persist) { try { localStorage.setItem(THEME_KEY, pref); } catch (e) { } }
  const b = $('theme'), name = { system: 'System', light: 'Light', dark: 'Dark' }[pref];
  b.innerHTML = `${THEME_ICONS[pref]}<span>${name}</span>`;
  b.setAttribute('aria-label', `Theme: ${pref}. Switch to ${THEME_NEXT[pref]}.`); b.title = `Theme: ${pref}`;
  readTheme(); themeEpoch++;
  if (redraw) redrawView();
}
function initTheme() {
  const q = new URLSearchParams(location.search).get('theme');
  if (q === 'light' || q === 'dark') { applyTheme(q, { persist: false, redraw: false }); return; }
  try { const t = localStorage.getItem(THEME_KEY); if (t in THEME_NEXT) themePref = t; } catch (e) { }
  applyTheme(themePref, { persist: false, redraw: false });
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { if (themePref === 'system') { readTheme(); themeEpoch++; redrawView(); } });
}

// ---------- glossary: the core terms and those of every part loaded so far, in alphabetical order ----------
function glossaryEntries() {
  const seen = new Set(), out = [];
  ['core', ...Object.keys(PARTS).filter(k => k !== 'core').sort()].forEach(k => {
    if (!UNPACKED.has(k)) return;
    ((PARTS[k] || {}).glossary || []).forEach(([term, def]) => { const key = term.toLowerCase(); if (!seen.has(key)) { seen.add(key); out.push([term, def, k]); } });
  });
  return out.sort((a, b) => a[0].localeCompare(b[0], 'en', { sensitivity: 'base' }));
}
function renderGlossary() {
  const html = glossaryEntries().map(([t, d]) => `<div class="gi"><dt>${t}</dt><dd>${d}</dd></div>`).join('');
  ['glossStart', 'glossPop'].forEach(id => { if ($(id)) $(id).innerHTML = html; });
}
function openGloss(open) { $('gloss').hidden = !open; $('glossBtn').setAttribute('aria-expanded', open); if (open) $('glossX').focus(); }

// ---------- interaction shared by all views ----------
function initNavigation() {
  VIEWS.forEach(t => { $(`tab-${t}`).onclick = () => goTo(t); });
  document.querySelector('.tabs').addEventListener('keydown', e => {
    if (!['ArrowLeft', 'ArrowRight'].includes(e.key)) return;
    const next = VIEWS[(VIEWS.indexOf(ROUTE.view) + (e.key === 'ArrowRight' ? 1 : VIEWS.length - 1)) % VIEWS.length];
    goTo(next); $(`tab-${next}`).focus();
  });
  VIEWS.filter(isBook).forEach(v => {
    $(v).querySelector('.rail').addEventListener('click', e => { const b = e.target.closest('button[data-i]'); if (b) goChapter(v, +b.dataset.i); });
    $(v).querySelector('.pager').addEventListener('click', e => { const b = e.target.closest('button[data-i]'); if (b) goChapter(v, +b.dataset.i); });
  });
  $('theme').onclick = () => applyTheme(THEME_NEXT[themePref]);
  $('glossBtn').onclick = () => openGloss($('gloss').hidden);
  $('glossX').onclick = () => openGloss(false);
  document.addEventListener('click', e => {
    if (!$('gloss').hidden && !e.target.closest('#gloss, #glossBtn')) openGloss(false);
    const g = e.target.closest('[data-go]'); if (g) { e.preventDefault(); goTo(g.dataset.go); return; }
    const a = e.target.closest('a[href^="#"]'); if (!a) return;
    const h = decodeURIComponent(a.getAttribute('href').slice(1));
    if (VIEWS.includes(h.split('/')[0]) || $(h)) { e.preventDefault(); goTo(h); }
  });
  window.addEventListener('popstate', () => goTo(location.hash, false));
  document.addEventListener('keydown', e => {
    if (e.key === 'Escape' && !$('gloss').hidden) { openGloss(false); $('glossBtn').focus(); return; }
    if (e.altKey || e.ctrlKey || e.metaKey || e.target.closest('input, select, textarea, .tabs, [role="slider"]')) return;
    const spec = SPECS[ROUTE.view] || {};
    if (spec.onKey && spec.onKey(e)) { e.preventDefault(); return; }
    if (isBook(ROUTE.view) && (e.key === 'ArrowRight' || e.key === 'ArrowLeft')) {
      const i = (CHAP[ROUTE.view] || 0) + (e.key === 'ArrowRight' ? 1 : -1);
      if (i >= 0 && i < secs(ROUTE.view).length) { e.preventDefault(); goChapter(ROUTE.view, i); }
    }
  });
  let width = innerWidth, timer = null;
  window.addEventListener('resize', () => { clearTimeout(timer); timer = setTimeout(() => { if (innerWidth !== width) { width = innerWidth; redrawView(); } }, 150); });
}

// ---------- self-test: open index.html?selftest=1 ----------
// It loads every part (and every study file listed in DATA_PARTS.studies.ids), opens every view and chapter in both
// themes, looking for JavaScript errors and for NaN, undefined or Infinity in the text, then runs every registered check.
// The result is shown in a panel and in the title: SELFTEST PASS n, or SELFTEST FAIL k/n with k checks failed of n.
const SELFTESTS = [];
function selftest(partName, name, fn) { SELFTESTS.push({ part: partName, name, fn }); }
const pause = () => new Promise(r => setTimeout(r, 0));
async function runSelftest() {
  const box = $('selftestBox'), results = [], themeBefore = themePref, start = location.hash;
  box.innerHTML = '<section class="selftest" id="selftest"><h2>Self-test running…</h2><p class="small" id="selftestNow"></p></section>';
  const now = t => { $('selftestNow').textContent = t; };
  const add = (p, name, ok, detail = '') => results.push({ part: p, name, ok: !!ok, detail: String(detail || '') });
  const errorsAtStart = ERRORS.length;
  // 1. every part, and every study file
  const parts = [...new Set(Object.values(SPECS).flatMap(s => [].concat(s.part || [])))];
  for (const p of parts) {
    now(`Loading ${partFile(p).split('?')[0]}`);
    try { await need(p); add('page', `loads ${partFile(p).split('?')[0]}`, true); } catch (e) { add('page', `loads ${partFile(p).split('?')[0]}`, false, e.message); }
  }
  const ids = ((part('studies') || {}).ids) || [];
  if (ids.length) {
    const failed = [];
    for (const id of ids) { try { await need(`studies/${id}`); } catch (e) { failed.push(id); } }
    add('page', `loads the file of every study (${fmt(ids.length)})`, !failed.length, failed.length ? `missing: ${failed.slice(0, 5).join(', ')}` : '');
  }
  // 2. every view and chapter, in both themes
  const bad = /\bNaN\b|\bundefined\b|\bInfinity\b|\[object /;
  for (const theme of ['light', 'dark']) {
    applyTheme(theme, { persist: false, redraw: false });
    for (const v of VIEWS) {
      const e0 = ERRORS.length, spec = SPECS[v] || {}, problems = [];
      now(`Drawing ${VIEW_NAME[v]} (${theme})`);
      await goTo(v, false);
      const routes = spec.routes ? await spec.routes() : isBook(v) ? secs(v).map(s => hashFor(v, s)) : [v];
      for (const r of routes) {
        await goTo(r, false); await pause();
        const shown = isBook(v) ? $(v).querySelector('section.sec.on') : $(v);
        const txt = shown ? shown.textContent : '';
        if (bad.test(txt)) problems.push(`${r}: "${txt.match(new RegExp(`.{0,40}(${bad.source}).{0,30}`))[0].trim()}"`);
      }
      const errs = ERRORS.slice(e0);
      add(v === 'start' ? 'core' : v, `every chapter draws in the ${theme} theme (${plural(routes.length, 'address', 'addresses')})`, !errs.length && !problems.length,
        [...errs.slice(0, 3), ...problems.slice(0, 3)].join('; '));
    }
  }
  applyTheme(themeBefore, { persist: false, redraw: false });
  await goTo(start || 'start', false);
  // 3. the registered checks
  for (const t of SELFTESTS) {
    now(`Checking ${t.part}: ${t.name}`);
    try {
      const r = await t.fn();
      if (r && typeof r === 'object') add(t.part, t.name, r.ok, r.detail); else add(t.part, t.name, r === true, r === true ? '' : 'returned false');
    } catch (e) { add(t.part, t.name, false, `${e.message} ${(e.stack || '').split('\n')[1] || ''}`); }
  }
  add('page', 'no JavaScript errors', ERRORS.length === 0, ERRORS.slice(0, 5).join('; ') + (errorsAtStart ? ` (${errorsAtStart} before the test)` : ''));
  // 4. the result
  const failed = results.filter(r => !r.ok), n = results.length;
  const head = failed.length ? `SELFTEST FAIL ${failed.length}/${n}` : `SELFTEST PASS ${n}`;
  document.title = head;
  box.innerHTML = `<section class="selftest ${failed.length ? 'fail' : 'pass'}" id="selftest"><h2>${head}</h2>
    <p class="small">${failed.length ? `${plural(failed.length, 'check')} of ${fmt(n)} failed; they are listed first.` : `All ${fmt(n)} checks passed.`}</p>
    <ol>${[...failed, ...results.filter(r => r.ok)].map(r => `<li class="${r.ok ? '' : 'f'}">${r.ok ? 'PASS' : 'FAIL'} · ${esc(r.part)} · ${esc(r.name)}${r.detail ? ` <span class="m">${esc(r.detail)}</span>` : ''}</li>`).join('')}</ol>
    <pre hidden id="selftest-result">${head}\n${results.map(r => `${r.ok ? 'PASS' : 'FAIL'}\t${r.part}\t${r.name}\t${r.detail}`).join('\n')}</pre></section>`;
  return { passed: n - failed.length, failed: failed.length, results };
}

// ---------- start: called by main.js once every script has registered its view ----------
function boot() {
  if ('scrollRestoration' in history) history.scrollRestoration = 'manual';
  BOOTED = true;
  initTheme();
  renderGlossary();
  initNavigation();
  goTo(location.hash, false);
  if (new URLSearchParams(location.search).get('selftest')) setTimeout(runSelftest, 300);
}
