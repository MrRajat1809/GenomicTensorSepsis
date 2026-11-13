/* Studies: the list of the included studies (sortable, searchable, beside a chart of all of them) and the excluded
   candidates with their reasons; then the chapters of the chosen study at addresses such as #studies/GSE65682/checks:
   overview, samples, platform and processing, genes, checks, correction and benchmark, each set beside the whole resource.
   Data: data/studies/index.js (the list, the candidates and the comparisons with the resource) and data/studies/<id>.js
   (one study), both written by studies() in build_data.py. Every sentence about a study is a template filled from these
   files, so that each study reads in the same way. */
'use strict';
(() => {
  const CHAPTERS = ['overview', 'samples', 'platform', 'genes', 'checks', 'correction', 'benchmark'];
  const FIRST = 'GSE65682';                   // the study shown before one is chosen: the guide's worked example
  const N = D.numbers;
  let P = null, LIST = [], BY = {}, study = null, find = null;
  const view = { sort: 'date', dir: 1, q: '', xq: '', xreason: null, hover: null, dots: {} };
  const cur = () => study ? part(`studies/${study}`) : undefined;

  // ======================================================================= words
  const REPO = { geo: 'GEO', arrayexpress: 'ArrayExpress', ena: 'ENA', ega: 'EGA', dbgap: 'dbGaP', node: 'NODE' };
  const famKey = i => D.order.family[i];
  const famText = i => lab('family', famKey(i));
  const famColor = i => C.fam[famKey(i)];
  const CLS_ONE = { sepsis: 'sepsis', healthy_control: 'healthy control', non_septic_control: 'non-septic control', other_condition: 'other condition', not_recorded: 'not recorded' };
  const clsCount = (n, k) => k === 'sepsis' || k === 'not_recorded' ? `${fmt(n)} ${CLS_ONE[k]}` : plural(n, CLS_ONE[k]);
  const clsList = counts => andList(D.order.class.map((k, i) => [k, counts[i]]).filter(([, n]) => n > 0).map(([k, n]) => clsCount(n, k)));
  // sample types in running text
  const TISSUE = { whole_blood: 'whole blood', blood_not_further_specified: 'blood not further specified', leukocytes: 'leukocytes', pbmc: 'PBMC', plasma: 'plasma', nasopharyngeal_swab: 'nasopharyngeal swabs', not_recorded: 'samples of unrecorded type' };
  const tisText = k => TISSUE[k] || labLower('tissue', k);
  const AGE_PL = { adult: 'adults', child: 'children', neonate: 'neonates' };
  const TASK = { sepsis_vs_healthy_control: 'Sepsis against healthy controls', sepsis_vs_non_septic_control: 'Sepsis against non-septic controls', mortality_in_sepsis: 'Death against survival in sepsis' };
  const TASK_SIDES = { sepsis_vs_healthy_control: ['sepsis', 'healthy controls'], sepsis_vs_non_septic_control: ['sepsis', 'non-septic controls'], mortality_in_sepsis: ['deaths', 'survivors'] };
  const SCHEME = { leave_one_group_out: 'Leave one group out', technology_transfer: 'Technology transfer', temporal: 'Temporal' };
  const GROUP = { bulk_blood: 'bulk-blood', pbmc: 'PBMC' };
  const METHOD = { none: 'No correction', plain_combat: 'Plain ComBat', flat_genes: 'Flat genes', matched_groups: 'Matched groups', healthy_anchored: 'Healthy controls' };
  const NO_VALUES = {
    absent_from_table: 'absent from the deposited table', no_own_data_file: 'without a data file of their own',
    not_processed: 'microRNA libraries, which hold no messenger RNA and were not processed', not_in_ncbi_counts: 'missing from the counts that NCBI generated', incomplete_array_file: 'with an incomplete array file',
  };
  const LEDGER = {
    repeat_deposit: 'repeat deposits whose first record is in the table', sample_type_not_corrected: 'of a sample type that is not corrected', newborn_study_not_corrected: 'of a newborn study, not corrected',
    no_class_to_line_up: 'without a class to line up', value_scale_not_comparable: 'on a value scale that cannot be aligned', study_has_too_few_samples: 'of a study platform with too few samples',
    no_check_applicable: 'of a study platform to which no mixing check applied',
  };
  const BENCH_REASON = {
    not_primary_record: 'are not the primary record of their sample', not_one_per_patient: 'are not the one sample marked for their patient',
    quality_flag: 'carry a quality flag', tissue_outside_task: 'are not of bulk blood', phase_outside_task: 'were drawn before the onset of sepsis or in recovery',
    class_outside_task: 'are of a class outside the task', label_not_recorded: 'have no recorded outcome',
  };
  const BENCH_SHORT = {
    not_primary_record: 'Not the primary record', not_one_per_patient: 'Not the sample of its patient', quality_flag: 'Quality flag',
    tissue_outside_task: 'Not bulk blood', phase_outside_task: 'Before onset or in recovery', class_outside_task: 'Class outside the task', label_not_recorded: 'No recorded outcome',
  };
  const ROLE = { train: 'trains', test: 'is tested', held_out_not_scored: 'is held out but not scored' };
  const SCALE = {
    log2: 'on the log₂ scale, where each step of one is a doubling', log2_like: 'log₂-like values normalized by the authors, kept as deposited',
    'log2 counts per million': 'log₂ counts per million, where each step of one is a doubling', linear: 'linear values (such as FPKM or normalized counts), kept as deposited',
    other: 'on a scale of the authors\' own, kept as deposited',
  };
  const repoName = d => REPO[d.src] || d.src;
  const studyLink = id => BY[id] ? `<a href="#studies/${id}">${id}</a>` : esc(id);
  const platLabel = p => p.p || 'no platform accession';
  const cleanTitle = t => String(t || '').replace(/^\[[^\]]*\]\s*/, '');
  const yr = d => d ? d.slice(0, 4) : '–';
  const howto = t => `<p class="howto"><span>How to read it</span><span>${t}</span></p>`;
  const sexColor = k => getComputedStyle($('studies')).getPropertyValue(`--sex-${k === 'female' || k === 'male' ? k : k === 'ambiguous' ? 'ambiguous' : 'not_recorded'}`).trim();
  const mixC = (a, f) => mixHex(a, C.bg, f);
  // a share as a percentage, with a decimal only where it would otherwise read 0% or 100%
  const pctS = v => v == null ? '–' : v > 0 && v < 0.01 ? pct(v, 1) : v < 1 && v > 0.99 ? pct(v, 1) : pct(v, 0);
  const pct1 = v => v != null && v >= 1 ? '100%' : pct(v, 1);
  const orList = items => items.length < 2 ? items.join('') : `${items.slice(0, -1).join(', ')} or ${items[items.length - 1]}`;
  const rankOf = (id, key) => 1 + LIST.filter(r => r[key] > BY[id][key]).length;

  // ======================================================================= the chart of all studies
  const yearPos = r => { const m = /^(\d{4})-(\d{2})-(\d{2})/.exec(r.date || ''); return m ? +m[1] + ((+m[2] - 1) * 30.44 + (+m[3] - 1)) / 365.25 : r.year + 0.5; };
  function niceLog(v, up) {
    const e = Math.floor(Math.log10(v)), steps = [1, 2, 3, 5, 10];
    if (up) { for (const s of steps) if (s * 10 ** e >= v) return s * 10 ** e; }
    for (const s of steps.slice().reverse()) if (s * 10 ** e <= v) return s * 10 ** e;
    return 10 ** e;
  }
  function studyTip(r) {
    return `<b>${r.id}</b> <span class="m">· ${repoName(r)} · ${yr(r.date)}</span><br>${esc(r.title)}<br>` +
      `${r.fam.map(i => `${sw(famColor(i))}${famText(i)}`).join(' ')}<br>${plural(r.records, 'sample record')}, ${plural(r.distinct, 'distinct sample')}` +
      `<br><span class="m">${clsList(r.cls)}</span>`;
  }
  // every study by release date (x) and number of records (y, logarithmic); hi is the study to pick out
  function drawMap(el, hi, compact) {
    const W = Math.max(300, el.clientWidth || 900), small = W < 600;
    const H = compact ? (small ? 230 : 250) : (small ? 270 : 330);
    const m = { l: small ? 46 : 58, r: 18, t: 18, b: 42 };
    const s = svg(el, W, H, `Scatter plot of the ${fmt(LIST.length)} included studies by release date and number of sample records, colored by technology${hi ? `; ${hi} is highlighted` : ''}.`);
    const x = lin(N.first_year, N.last_year + 1, m.l, W - m.r);
    const recs = LIST.map(r => r.records), lo = niceLog(Math.min(...recs), false), top = niceLog(Math.max(...recs), true);
    const y = lg(lo, top, H - m.b, m.t);
    const ytk = []; for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(top)); e++) [1, 3].forEach(k => { const v = k * 10 ** e; if (v >= lo && v <= top) ytk.push(v); });
    ytk.forEach(v => { S('line', { x1: m.l, x2: W - m.r, y1: y(v), y2: y(v), stroke: C.grid }, s); T(s, m.l - 6, y(v) + 4, fmt(v), { 'text-anchor': 'end', 'font-size': 11, fill: C.muted }); });
    const step = small ? 4 : 2;
    for (let v = Math.ceil(N.first_year / step) * step; v <= N.last_year; v += step) {
      S('line', { x1: x(v), x2: x(v), y1: m.t, y2: H - m.b, stroke: C.grid }, s);
      T(s, x(v), H - m.b + 15, String(v), { 'text-anchor': 'middle', 'font-size': 11, fill: C.muted });
    }
    axisTitle(s, (m.l + W - m.r) / 2, H - 6, 'Release date of the study');
    axisTitle(s, 13, (m.t + H - m.b) / 2, 'Sample records (log scale)', true);
    const r0 = small ? 4.2 : 5.4, dots = {};
    const order = LIST.slice().sort((a, b) => (a.id === hi) - (b.id === hi));
    order.forEach(r => {
      const isHi = r.id === hi, col = famColor(r.fam[0]);
      const c = S('circle', {
        cx: x(yearPos(r)), cy: y(r.records), r: isHi ? r0 + 2.5 : r0, fill: col, 'fill-opacity': hi && !isHi ? 0.45 : 0.82,
        stroke: isHi ? C.ink : mixHex(col, C.dark ? '#ffffff' : '#000000', 0.3), 'stroke-width': isHi ? 2 : 0.8, tabindex: 0, role: 'link',
        'aria-label': `${r.id}, ${yr(r.date)}, ${plural(r.records, 'record')}`, style: 'cursor:pointer',
      }, s);
      hover(c, () => studyTip(r));
      c.addEventListener('click', () => { untip(); openStudy(r.id); });
      c.addEventListener('keydown', e => { if (e.key === 'Enter') openStudy(r.id); });
      c.addEventListener('mouseenter', () => { if (!compact) hoverStudy(r.id, 'map'); });
      c.addEventListener('mouseleave', () => { if (!compact) hoverStudy(null, 'map'); });
      dots[r.id] = c;
    });
    if (hi && BY[hi]) {
      const r = BY[hi], cx = x(yearPos(r)), cy = y(r.records), right = cx < W - 110;
      T(s, cx + (right ? 11 : -11), cy - 9, hi, { 'text-anchor': right ? 'start' : 'end', 'font-size': 12.5, 'font-weight': 600, fill: C.ink });
    }
    return dots;
  }
  function famKeys(el) {
    const n = D.order.family.map((k, i) => LIST.filter(r => r.fam[0] === i).length);
    legend(el, D.order.family.map((k, i) => ({ key: k, label: lab('family', k), color: C.fam[k], count: n[i] })).filter(x => x.count > 0));
  }
  // a study picked out in the chart and the list together
  function hoverStudy(id, from) {
    const prev = view.hover; view.hover = id;
    if (prev && view.dots[prev]) { const r = view.dots[prev]; r.setAttribute('r', r.dataset.r0 || r.getAttribute('r')); r.setAttribute('stroke-width', 0.8); r.setAttribute('stroke', r.dataset.rim || r.getAttribute('stroke')); }
    if (id && view.dots[id]) {
      const r = view.dots[id];
      if (!r.dataset.r0) { r.dataset.r0 = r.getAttribute('r'); r.dataset.rim = r.getAttribute('stroke'); }
      r.setAttribute('r', +r.dataset.r0 + 2.5); r.setAttribute('stroke', C.ink); r.setAttribute('stroke-width', 2); r.parentNode.appendChild(r);
    }
    if (from === 'map') $('stList').querySelectorAll('tr.hi').forEach(t => t.classList.remove('hi'));
    if (from === 'map' && id) { const tr = $('stList').querySelector(`tr[data-id="${id}"]`); if (tr) tr.classList.add('hi'); }
  }

  // ======================================================================= choosing a study
  const COLS = [
    { k: 'id', label: 'Study', cls: 'acc', key: r => r.id, text: true },
    { k: 'title', label: 'Title', cls: 'c-title', key: r => r.title.toLowerCase(), text: true },
    { k: 'date', label: 'Released', cls: '', key: r => r.date || '', text: true },
    { k: 'fam', label: 'Technology', cls: 'c-tech', key: r => r.fam[0] },
    { k: 'records', label: 'Records', cls: 'r', key: r => r.records },
    { k: 'distinct', label: 'Distinct samples', cls: 'r', key: r => r.distinct },
    { k: 'cls', label: 'Classes', cls: 'c-cls', key: r => r.records ? r.cls[0] / r.records : 0 },
  ];
  const clsBar = (counts, wide) => {
    const t = sum(counts);
    return `<span class="st-cbar${wide ? ' wide' : ''}" data-cls="${counts.join(',')}">${t ? D.order.class.map((k, i) => counts[i] ? `<i style="width:${(100 * counts[i] / t).toFixed(2)}%;background:${C.cls[k]}"></i>` : '').join('') : ''}</span>`;
  };
  const clsTip = counts => { const t = sum(counts); return D.order.class.map((k, i) => counts[i] ? `${sw(C.cls[k])}${lab('class', k)}: <b>${fmt(counts[i])}</b> <span class="m">(${pctS(counts[i] / t)})</span>` : '').filter(Boolean).join('<br>'); };
  const matches = (text, q) => q.split(/\s+/).filter(Boolean).every(w => text.includes(w));
  const searchText = r => `${r.id} ${r.title} ${r.plats.join(' ')} ${r.year || ''} ${r.fam.map(famText).join(' ')} ${repoName(r)}`.toLowerCase();
  function renderList() {
    const c = COLS.find(x => x.k === view.sort) || COLS[2];
    const rows = LIST.filter(r => matches(searchText(r), view.q)).sort((a, b) => {
      const p = c.key(a), q = c.key(b);
      const o = c.text ? String(p).localeCompare(String(q), 'en', { numeric: true }) : p - q;
      return (o || a.id.localeCompare(b.id, 'en', { numeric: true })) * view.dir;
    });
    $('stList').innerHTML = `<caption class="small" style="text-align:left;caption-side:bottom;padding-top:8px">Records are the deposited sample records of a study; distinct samples count each array or library once, in the study that deposited it first. Click a column heading to sort.</caption>
      <thead><tr>${COLS.map(x => `<th class="${x.cls === 'r' ? 'r' : x.cls}" aria-sort="${x.k === view.sort ? (view.dir > 0 ? 'ascending' : 'descending') : 'none'}"><button type="button" data-k="${x.k}">${x.label}</button></th>`).join('')}</tr></thead>
      <tbody>${rows.length ? rows.map(r => `<tr data-id="${r.id}" class="${r.id === study ? 'sel' : ''}">
        <td class="acc"><a href="#studies/${r.id}">${r.id}</a></td>
        <td class="c-title"><span class="ttl" title="${esc(r.title)}">${esc(r.title)}</span></td>
        <td>${yr(r.date)}</td>
        <td class="c-tech"><span class="tech" title="${r.fam.map(famText).join(', ')}"><i style="background:${famColor(r.fam[0])}"></i><span class="tl">${famText(r.fam[0])}${r.fam.length > 1 ? ` <span class="small">+${r.fam.length - 1}</span>` : ''}</span></span></td>
        <td class="r">${fmt(r.records)}</td><td class="r">${fmt(r.distinct)}</td>
        <td class="c-cls">${clsBar(r.cls)}</td></tr>`).join('') : `<tr class="none"><td colspan="${COLS.length}">No included study matches “${esc(view.q)}”.</td></tr>`}</tbody>`;
    $('stCount').textContent = view.q ? `${fmt(rows.length)} of ${fmt(LIST.length)} studies match` : `${fmt(LIST.length)} studies`;
  }
  function renderExcluded() {
    const X = P.excluded, R = P.reasons, q = view.xq;
    const n = {}; X.forEach(x => { n[x.reason] = (n[x.reason] || 0) + 1; });
    const max = Math.max(...Object.values(n));
    $('stReasons').innerHTML = R.order.filter(k => n[k]).map(k => `<button type="button" data-r="${k}" aria-pressed="${view.xreason === k}" title="${esc(R.text[k])}"><span>${esc(R.label[k])}</span><span class="b${R.availability.includes(k) ? ' av' : ''}" style="width:${(100 * n[k] / max).toFixed(1)}%"></span><span class="v">${fmt(n[k])}</span></button>`).join('');
    const rank = k => R.order.indexOf(k);
    const rows = X.filter(x => (!view.xreason || x.reason === view.xreason) && matches(`${x.id} ${x.title} ${x.basis} ${R.label[x.reason]} ${x.year || ''} ${x.assay}`.toLowerCase(), q))
      .sort((a, b) => rank(a.reason) - rank(b.reason) || a.id.localeCompare(b.id, 'en', { numeric: true }));
    $('stExclList').innerHTML = `<thead><tr><th>Candidate</th><th class="c-title">Title</th><th>Year</th><th>Reason</th><th class="c-basis">Recorded basis</th></tr></thead>
      <tbody>${rows.length ? rows.map(x => `<tr data-x="${esc(x.id)}" class="${find && x.id === find ? 'sel' : ''}"><td class="acc">${x.url ? `<a href="${esc(x.url)}" target="_blank" rel="noopener">${esc(x.id)}</a>` : esc(x.id)}</td>
        <td class="c-title"><span class="ttl" title="${esc(x.title)}">${esc(x.title) || '<span class="small">no title</span>'}</span></td><td>${x.year || '–'}</td>
        <td>${esc(R.label[x.reason])}${x.overlaps ? `<br><span class="small">with ${x.overlaps.split(/[;,]\s*/).map(studyLink).join(', ')}</span>` : ''}</td>
        <td class="basis c-basis">${esc(x.basis)}${x.flag === 'uncertain' ? ' <span class="chip">uncertain</span>' : ''}</td></tr>`).join('') : `<tr class="none"><td colspan="5">No excluded candidate matches.</td></tr>`}</tbody>`;
    $('stExclCount').textContent = q || view.xreason ? `${fmt(rows.length)} of ${fmt(X.length)} candidates shown` : `${fmt(X.length)} candidates`;
  }
  function renderChooser() {
    const X = P.excluded, av = X.filter(x => P.reasons.availability.includes(x.reason)).length;
    const T_ = N.studies_by_technology;
    $('stChooseLead').textContent = `All ${fmt(LIST.length)} included studies, released from ${N.first_year} to ${N.last_year}: ${fmt(T_.array)} measured with microarrays and ${fmt(T_.rna_seq)} with RNA-seq. Choose one in the chart or the list to see its samples, its platform and processing, its genes, the checks, the correction and its place in the benchmark, each set beside the whole resource. Below the list are the ${fmt(X.length)} candidate studies that were excluded, each with its reason.`;
    $('stChooseHow').textContent = 'Each dot is one study, placed by the date of its release and the number of its sample records, on a logarithmic scale where each step of the grid is about three times the last. Its color is its technology (that of most of its records, where it used several). Hover over a dot or a row to see the study, and click it to open it.';
    view.dots = drawMap($('stMap'), null, false); view.hover = null;
    famKeys('stMapKeys');
    $('stSearch').value = view.q;
    renderList();
    $('stExclSum').innerHTML = `Excluded candidates <small>${fmt(X.length)} of the ${fmt(P.candidates)} candidate studies, each with the first reason that applied</small>`;
    $('stExclNote').textContent = `Every candidate was screened against the eligibility rules in a fixed order and excluded for the first reason that applied, so each is counted once: ${fmt(X.length - av)} for a scientific reason and ${fmt(av)} for availability (green bars), beside the ${fmt(LIST.length)} included studies. The studies excluded for availability met the scientific rules, but their data are under controlled access, not public, available only as raw reads or deposited without an expression file. Click a reason to list its studies; the basis is the note recorded for each decision.`;
    if (find) {
      const x = P.excluded.find(e => e.id === find), box = $('stFind') || (() => { const p = document.createElement('p'); p.id = 'stFind'; p.className = 'result'; $('stChooseLead').after(p); return p; })();
      box.innerHTML = x ? `<b>${esc(find)}</b> is not an included study. It was a candidate${x.title ? `, “${esc(x.title)}”` : ''}${x.year ? ` (${x.year})` : ''}, excluded for the first reason that applied: ${esc(P.reasons.text[x.reason].replace(/^[A-Z](?=[a-z])/, c => c.toLowerCase())).replace(/\.$/, '')}. The basis recorded for the decision reads “${esc(x.basis)}”.${x.url ? ` <a href="${esc(x.url)}" target="_blank" rel="noopener">${esc(find)} in its repository</a>.` : ''} It is shown in the list of excluded candidates below.`
        : `${esc(find)} is neither an included study nor one of the ${fmt(P.candidates)} candidate studies.`;
      if (x) { $('stExcl').open = true; view.xq = find.toLowerCase(); view.xreason = null; $('stExclSearch').value = find; }
    } else if ($('stFind')) $('stFind').remove();
    if ($('stExcl').open) renderExcluded();
  }
  function initChooser() {
    $('stList').addEventListener('click', e => {
      const b = e.target.closest('th button[data-k]');
      if (b) { const k = b.dataset.k; view.dir = view.sort === k ? -view.dir : (['records', 'distinct', 'cls'].includes(k) ? -1 : 1); view.sort = k; renderList(); return; }
      if (e.target.closest('a')) return;
      const tr = e.target.closest('tr[data-id]'); if (tr) openStudy(tr.dataset.id);
    });
    $('stList').addEventListener('mouseover', e => { const tr = e.target.closest('tr[data-id]'); const id = tr ? tr.dataset.id : null; if (id !== view.hover) hoverStudy(id, 'list'); });
    $('stList').addEventListener('mouseleave', () => hoverStudy(null, 'list'));
    $('stSearch').addEventListener('input', e => { view.q = e.target.value.trim().toLowerCase(); renderList(); });
    $('stExcl').addEventListener('toggle', () => { if ($('stExcl').open) renderExcluded(); });
    $('stExclSearch').addEventListener('input', e => { view.xq = e.target.value.trim().toLowerCase(); renderExcluded(); });
    $('stReasons').addEventListener('click', e => { const b = e.target.closest('button[data-r]'); if (!b) return; view.xreason = view.xreason === b.dataset.r ? null : b.dataset.r; renderExcluded(); });
    // the class bars of the list, and of any chapter, explain themselves on hover
    $('studies').addEventListener('mousemove', e => { const b = e.target.closest('.st-cbar[data-cls]'); if (b) tip(e, clsTip(b.dataset.cls.split(',').map(Number))); });
    $('studies').addEventListener('mouseout', e => { if (e.target.closest('.st-cbar[data-cls]')) untip(); });
    // choosing another study inside its chapters
    $('studies').addEventListener('change', e => { if (e.target.matches('select.st-pick')) openStudy(e.target.value, chapterKey()); });
    $('studies').addEventListener('click', e => {
      const b = e.target.closest('button[data-step]'); if (!b) return;
      const i = P.ids.indexOf(study) + +b.dataset.step; if (i >= 0 && i < P.ids.length) openStudy(P.ids[i], chapterKey());
    });
  }
  const chapterKey = () => (ROUTE.chapter || 'st-overview').replace(/^st-/, '').replace('choose', 'overview');
  const openStudy = (id, ch = 'overview') => goTo(`studies/${id}/${ch}`);

  // ======================================================================= small charts shared by the chapters
  // stacked horizontal bars: rows = [{label, sub, parts: [{key, label, n, color}]}]; every bar runs to 100%
  function stackBars(el, rows, label) {
    const W = Math.max(300, el.clientWidth || 800), small = W < 560, lw = small ? 0 : Math.min(250, W * 0.3);
    const bh = 20, gap = small ? 34 : 12, top = small ? 16 : 4, H = top + rows.length * (bh + gap);
    const s = svg(el, W, H, label), x0 = lw + (small ? 0 : 10), x1 = W - 4;
    rows.forEach((row, i) => {
      const y = top + i * (bh + gap), total = sum(row.parts.map(p => p.n));
      if (small) T(s, 0, y - 5, `${row.label}${row.sub ? ` · ${row.sub}` : ''}`, { 'font-size': 12, fill: C.ink2 });
      else { T(s, lw, y + 9, row.label, { 'text-anchor': 'end', 'font-size': 12.5, fill: C.ink2 }); if (row.sub) T(s, lw, y + 21, row.sub, { 'text-anchor': 'end', 'font-size': 11, fill: C.muted }); }
      if (!total) { rrect(s, x0, y, x1 - x0, bh, C.band, 3); T(s, x0 + 8, y + 14, 'none', { 'font-size': 11.5, fill: C.muted }); return; }
      let x = x0;
      row.parts.forEach(p => {
        if (!p.n) return;
        const w = (x1 - x0) * p.n / total;
        // a row of shares (rows.shares) shows percentages, a row of counts its counts
        hover(S('rect', { x, y, width: Math.max(0.5, w - 1), height: bh, fill: p.color }, s), row.shares ? `${sw(p.color)}${p.label}: <b>${pctS(p.n / total)}</b>` : `${sw(p.color)}${p.label}: <b>${fmt(p.n)}</b> <span class="m">(${pctS(p.n / total)} of ${fmt(total)})</span>`);
        const t = row.shares ? pctS(p.n / total) : fmt(p.n);
        if (w > t.length * 7 + 8) T(s, x + w / 2, y + 14, t, { 'text-anchor': 'middle', 'font-size': 11, fill: sum(hexRGB(p.color)) > 420 ? '#222' : '#fff', 'pointer-events': 'none' });
        x += w;
      });
    });
  }
  const clsParts = counts => D.order.class.map((k, i) => ({ key: k, label: lab('class', k), n: counts[i] || 0, color: C.cls[k] }));
  const techKeys = () => `<div class="keys">${keys(D.order.family.map(k => [lab('family', k), C.fam[k], 'dot']))}</div>`;
  const clsKeys = () => `<div class="keys">${keys(D.order.class.map(k => [lab('class', k), C.cls[k], 'dot']))}</div>`;
  // a strip of every table or stratum of the resource (values from 'get'), with the study's own picked out
  function stripChart(el, items, { get, lo, hi, axis, label, tipFn, isMine, color, fmtv = v => fx(v, 2), mark = null }) {
    const W = Math.max(300, el.clientWidth || 800), H = 92, m = { l: 14, r: 14, t: 22, b: 36 };
    const s = svg(el, W, H, label), x = lin(lo, hi, m.l, W - m.r);
    gridX(s, x, ticks(lo, hi, W < 560 ? 4 : 7), m.t - 6, H - m.b + 4, fmtv);
    axisTitle(s, (m.l + W - m.r) / 2, H - 4, axis);
    if (mark != null) S('line', { x1: x(mark), x2: x(mark), y1: m.t - 8, y2: H - m.b + 4, stroke: C.ink2, 'stroke-dasharray': '3 3' }, s);
    const mid = (m.t + H - m.b) / 2, rr = rng(7);
    const sorted = items.slice().sort((a, b) => isMine(a) - isMine(b));
    sorted.forEach(it => {
      const v = get(it); if (v == null || !isFinite(v)) return;
      const mine = isMine(it), cy = mine ? mid : mid + (rr() - 0.5) * 22;
      const c = S('circle', { cx: x(clamp(v, lo, hi)), cy, r: mine ? 7 : 4.2, fill: mine ? C.hilite : (color ? color(it) : C.lineStrong), 'fill-opacity': mine ? 0.95 : 0.6, stroke: mine ? C.ink : 'none', 'stroke-width': mine ? 1.5 : 0 }, s);
      hover(c, () => tipFn(it));
      if (mine) T(s, x(clamp(v, lo, hi)), m.t - 9, fmtv(v), { 'text-anchor': 'middle', 'font-size': 12, 'font-weight': 600, fill: C.ink });
    });
  }

  // ======================================================================= chapter 1: overview
  function tissuePhrase(t) { const k = Object.keys(t).filter(x => x !== 'not_recorded'); return k.length ? andList(k.map(tisText)) : 'samples of unrecorded type'; }
  function agePhrase(a) { const k = ['adult', 'child', 'neonate'].filter(x => a[x]); return k.length ? ` of ${andList(k.map(x => AGE_PL[x]))}` : ''; }
  function platformPhrase(plats) {
    const one = p => `${cleanTitle(p.title) || famText(p.fam)}${p.p ? ` (${p.p})` : ''}`;
    if (plats.length === 1) { const p = plats[0]; return famKey(p.fam) === 'rna_seq' ? `by RNA-seq on the ${one(p)}` : `on the ${one(p)}`; }
    return `on ${spell(plats.length)} platforms: ${andList(plats.map(one))}`;
  }
  const COVER = { sex: 'sex', age: 'age', patient: 'a patient identifier', time: 'a time point', outcome: 'an outcome', shock: 'septic shock' };
  function coverSentence(cov, records) {
    const present = Object.keys(COVER).filter(k => cov[k] > 0), absent = Object.keys(COVER).filter(k => !(cov[k] > 0));
    if (!present.length) return `Its records give none of sex, age, a patient identifier, a time point, an outcome or septic shock.`;
    const s = `Its records give ${andList(present.map(k => `${COVER[k]} (${cov[k] >= 1 ? 'all' : pctS(cov[k])})`))}`;
    return absent.length ? `${s}, but no ${orList(absent.map(k => COVER[k].replace(/^an? /, '')))}.` : `${s}.`;
  }
  function overviewHTML(d) {
    const L = BY[d.id], repeats = d.records - d.distinct, tables = d.platforms.filter(p => p.genes).length;
    const dist = repeats === 0 ? ', each of them a distinct sample'
      : d.distinct === 0 ? `, all of them repeat deposits of arrays or libraries first deposited in ${d.repeat_same ? 'other records' : 'other studies'}`
        : `: ${plural(d.distinct, 'distinct sample')} first deposited here and ${fmt(repeats)} repeat ${repeats === 1 ? 'deposit' : 'deposits'} of samples first deposited ${d.repeat_same && !d.repeat_other ? 'in another record of this study' : d.repeat_same ? 'earlier, here or in other studies' : 'in other studies'}`;
    const lead = `${d.id} is a ${repoName(d)} study released on ${fdate(d.date)}. It measured ${tissuePhrase(d.samples.tissue)}${agePhrase(d.samples.age_group)} ${platformPhrase(d.platforms)}, and holds ${plural(d.records, 'sample record')}${dist}${d.without ? `; ${plural(d.without, 'further record')} without expression values ${d.without === 1 ? 'is' : 'are'} listed apart with ${d.without === 1 ? 'its' : 'their'} reason` : ''}. By class, its records are ${clsList(d.cls)}. ${coverSentence(d.samples.coverage)}`;
    const rank = rankOf(d.id, 'records'), ranked = LIST.filter(r => r.records === d.records).length > 1 ? ' (shared)' : '';
    const links = [`<a href="${esc(d.url)}" target="_blank" rel="noopener">${d.id} in ${repoName(d)}</a>`]
      .concat(d.also_in.map(a => /^E-/.test(a) ? `<a href="https://www.ebi.ac.uk/biostudies/arrayexpress/studies/${esc(a)}" target="_blank" rel="noopener">copy ${esc(a)} in ArrayExpress</a>`
        : /^PRJ/.test(a) ? `<a href="https://www.ebi.ac.uk/ena/browser/view/${esc(a)}" target="_blank" rel="noopener">reads ${esc(a)} in ENA</a>` : esc(a)));
    const pub = d.pubmed.length ? `Articles in PubMed: ${d.pubmed.map(p => `<a href="https://pubmed.ncbi.nlm.nih.gov/${esc(p)}/" target="_blank" rel="noopener">${esc(p)}</a>`).join(', ')}.` : 'The record links no article in PubMed.';
    const rel = [d.superseries.length ? `It is a SuperSeries of ${andList(d.superseries.map(studyLink))}.` : '', d.subseries_of.length ? `It is part of the SuperSeries ${andList(d.subseries_of.map(esc))}.` : ''].join(' ');
    const others = [...new Set(d.group_members.map(g => g.id))].filter(id => id !== d.id);
    return `
      <p class="st-title">${esc(d.title)}</p>
      <p class="st-meta">${links.join(' · ')}</p>
      <p class="lead">${lead}</p>
      <div class="cards">
        <div class="card"><div class="lab">Sample records</div><div class="big">${fmt(d.records)}</div><div class="lab">${rank === 1 ? 'the largest' : `the ${ordinal(rank)} largest`}${ranked} of the ${fmt(LIST.length)} studies</div></div>
        <div class="card"><div class="lab">Distinct samples</div><div class="big">${fmt(d.distinct)}</div><div class="lab">first deposited here; ${pctS(d.distinct / N.distinct)} of the resource's ${fmt(N.distinct)}</div></div>
        <div class="card"><div class="lab">Expression tables</div><div class="big">${fmt(tables)}</div><div class="lab">of ${plural(d.platforms.length, 'study platform')}${d.without ? `; ${plural(d.without, 'record')} without values` : ''}</div></div>
        <div class="card"><div class="lab">Released</div><div class="big">${yr(d.date)}</div><div class="lab">${fdate(d.date)}; the resource spans ${N.first_year}–${N.last_year}</div></div>
      </div>
      <h3>Where it stands among the studies</h3>
      <div class="chart" id="stOvMap"></div>
      <div class="keys" id="stOvKeys"></div>
      ${howto(`Every included study as a dot, by release date and number of sample records (logarithmic scale), colored by technology; ${d.id} is ringed and labeled. Click any dot to open that study.`)}
      <h3>Platforms</h3>
      <div class="tablewrap"><table class="t"><thead><tr><th>Study platform</th><th>Technology</th><th class="r">Records</th><th class="r">Distinct samples</th><th class="r">Without values</th><th>Classes</th></tr></thead><tbody>
        ${d.platforms.map(p => `<tr><td><b>${esc(platLabel(p))}</b><br><span class="small">${esc(cleanTitle(p.title))}</span></td><td><span class="st-dot" style="background:${famColor(p.fam)}"></span>${famText(p.fam)}</td><td class="r">${fmt(p.records)}</td><td class="r">${fmt(p.distinct)}</td><td class="r">${fmt(p.without)}</td><td>${clsBar(p.cls)}</td></tr>`).join('')}
      </tbody></table></div>
      <p class="note">A study platform is the part of a study measured on one platform, the unit in which values were processed. ${d.platforms.length > 1 ? `${d.id} has ${spell(d.platforms.length)}.` : ''}</p>
      <h3>Where its samples also appear</h3>
      ${sharesHTML(d, others)}
      <h3>Why it was included</h3>
      <p class="body">Screening records one short basis for every decision. For ${d.id} it reads:</p>
      <blockquote class="st-quote">${esc(d.basis) || 'No basis recorded.'}</blockquote>
      ${d.flag === 'uncertain' ? '<p class="note">The decision is marked uncertain: the records were ambiguous, and a reader might decide otherwise. The mark does not change the decision.</p>' : ''}
      <p class="small">${pub} ${rel}</p>
      ${d.summary || d.design ? `<details class="st-own"><summary>The study's own summary, as deposited</summary>${d.summary ? `<p>${esc(d.summary)}${d.summary.length >= 3000 ? '…' : ''}</p>` : ''}${d.design ? `<p><b>Design.</b> ${esc(d.design)}${d.design.length >= 1500 ? '…' : ''}</p>` : ''}</details>` : ''}`;
  }
  function sharesHTML(d, others) {
    const total = sum(d.shares.map(s => s.samples));
    if (!d.shares.length) return `<p class="body">None of the arrays or libraries of ${d.id} is deposited in another included study${d.repeat_same ? `; ${plural(d.repeat_same, 'record')} repeat${d.repeat_same === 1 ? 's' : ''} a sample deposited earlier within the study itself` : ''}. ${others.length ? `It shares an identified patient with ${andList(others.map(studyLink))}.` : ''}</p>`;
    const first = s => s.first_here === s.samples ? 'here' : s.first_there === s.samples ? 'there' : andList([s.first_here ? `here (${fmt(s.first_here)})` : '', s.first_there ? `there (${fmt(s.first_there)})` : '', s.first_elsewhere ? `in a third study (${fmt(s.first_elsewhere)})` : ''].filter(Boolean));
    const groups = d.cohort_groups.length === 1 ? `cohort group ${d.cohort_groups[0]}` : `cohort groups ${andList(d.cohort_groups)}`;
    return `<p class="body">Its records share arrays or libraries with ${plural(d.shares.length, 'other included study', 'other included studies')}. Records of one physical sample were linked by evidence in their raw files, such as the scan time of an array or identical counts, and each record is kept. Together the studies that share samples form ${groups}; a sample counts once, in the study that deposited it first.</p>
      <div class="tablewrap"><table class="t"><thead><tr><th>Other study</th><th class="r">Shared samples</th><th>First deposited</th><th>Linked by</th></tr></thead><tbody>
      ${d.shares.map(s => `<tr><td>${studyLink(s.id)}</td><td class="r">${fmt(s.samples)}</td><td>${first(s)}</td><td class="small">${Object.entries(s.by).map(([k, n]) => `${esc(P.links[k] || k.replace(/_/g, ' '))}${Object.keys(s.by).length > 1 ? ` (${fmt(n)})` : ''}`).join('; ')}</td></tr>`).join('')}
      </tbody></table></div>
      <p class="note">Shared samples are counted once for every other study that holds them, so a sample deposited in three studies appears in two rows; ${fmt(total)} in all.</p>`;
  }
  function overviewDraw(d) { drawMap($('stOvMap'), d.id, true); famKeys('stOvKeys'); }

  // ======================================================================= chapter 2: samples
  const factList = (title, counts, labelOf) => `<dl><dt>${title}</dt>${Object.keys(counts).length ? Object.entries(counts).map(([k, n]) => `<dd><span>${esc(labelOf(k))}</span><b>${fmt(n)}</b></dd>`).join('') : '<dd><span>none</span></dd>'}</dl>`;
  const SEX = { female: 'Female', male: 'Male', other: 'Other', not_recorded: 'Not recorded' };
  const timeLabel = t => t.unit === 'day' ? `Day ${num(t.order)}` : t.unit === 'hour' ? `Hour ${num(t.order)}` : t.unit === 'ordinal' ? `${ordinal(t.order)} sample` : (t.words[0] || 'Not stated');
  const VAR = { sample_class: 'Class', tissue: 'Sample type', sex: 'Sex', age: 'Age', time: 'Time point', outcome: 'Outcome', septic_shock: 'Septic shock', patient: 'Patient' };
  function valueLabel(v, x) {
    if (x === 'not_recorded') return 'Not recorded';
    if (v === 'sample_class') return x === 'infection_sepsis_not_stated' ? 'Other condition (an infection; sepsis not stated)' : lab('class', x);
    if (v === 'tissue') return lab('tissue', x);
    if (v === 'age') return x !== '' && isFinite(+x) ? `${num(+x)} years` : lab('age', x);
    if (v === 'time') return `order ${x !== '' && isFinite(+x) ? num(+x) : x}`;
    if (v === 'outcome') return { survivor: 'Survivor', non_survivor: 'Non-survivor' }[x] || x;
    if (v === 'septic_shock') return { yes: 'Septic shock', no: 'No septic shock' }[x] || x;
    if (v === 'sex') return SEX[x] || x;
    return x;
  }
  function wordingHTML(d) {
    const W = d.samples.wording, rows = [];
    Object.keys(VAR).forEach(v => {
      const w = W[v]; if (!w) return;
      const items = w.items.filter(i => !i.more), more = w.items.find(i => i.more);
      if (!items.length) { rows.push(`<tr><td><b>${VAR[v]}</b></td><td colspan="5" class="small">Not recorded in any of the ${fmt(w.records)} records.</td></tr>`); return; }
      items.forEach((it, j) => rows.push(`<tr>${j === 0 ? `<td rowspan="${items.length + (more ? 1 : 0) + (w.not_recorded ? 1 : 0)}"><b>${VAR[v]}</b></td>` : ''}
        <td class="small">${esc(it.field)}</td><td class="w"><code>${esc(it.recorded) || '<span class="small">(empty)</span>'}</code></td><td>${esc(valueLabel(v, it.value))}</td><td class="r">${(v === 'age' || v === 'patient') && j > 0 && items[j - 1].field === it.field ? '' : fmt(it.n)}</td>
        <td class="ev">${esc(it.ex)} · ${esc(P.basis[it.basis] || it.basis)}${it.statement ? `: “${esc(it.statement)}”` : ''}${it.note ? `<br>${esc(it.note)}` : ''}</td></tr>`));
      if (more) rows.push(`<tr><td colspan="5" class="small">… and ${plural(more.more, 'other wording')} (${plural(more.n, 'record')})</td></tr>`);
      if (w.not_recorded) rows.push(`<tr><td colspan="5" class="small">Not recorded in ${plural(w.not_recorded, 'record')}</td></tr>`);
    });
    return `<div class="tablewrap"><table class="t st-words"><thead><tr><th>Variable</th><th>Field</th><th>Recorded wording</th><th>Curated value</th><th class="r">Records</th><th>Example and evidence</th></tr></thead><tbody>${rows.join('')}</tbody></table></div>`;
  }
  function samplesHTML(d) {
    const s = d.samples, repeats = d.records - d.distinct;
    const basis = s.basis, bsum = [['sample_record', 'on the sample record alone'], ['sample_record_and_study_description', 'on the record together with a quoted sentence of the study'], ['study_description', 'on the study description alone']]
      .filter(([k]) => basis[k]).map(([k, t]) => `${t} for ${plural(basis[k], 'record')}`);
    const sepsisN = s.sepsis, oc = s.outcome, sh = s.shock, ph = s.phase;
    const defs = Object.keys(s.outcome_definition);
    const outcome = sepsisN ? (oc.survivor || oc.non_survivor ? `Of the ${fmt(sepsisN)} sepsis records, ${fmt(oc.survivor || 0)} are of survivors and ${fmt(oc.non_survivor || 0)} of patients who died${defs.length ? ` (${andList(defs.map(x => x === 'not stated' ? 'definition not stated' : `${x} mortality`))})` : ''}${oc.not_recorded ? `; ${fmt(oc.not_recorded)} have no recorded outcome` : ''}.` : `None of the ${fmt(sepsisN)} sepsis records gives an outcome.`) : 'The study has no sepsis records, so outcome and septic shock do not apply.';
    const shock = sepsisN ? (sh.yes || sh.no ? `${fmt(sh.yes || 0)} sepsis records state septic shock and ${fmt(sh.no || 0)} state its absence${sh.not_recorded ? `; ${fmt(sh.not_recorded)} say nothing about it` : ''}.` : 'No record states whether the patient was in septic shock.') : '';
    const PH = [['during_sepsis', 'during sepsis'], ['before_onset', 'before its onset'], ['onset_not_recorded', 'at a time not known relative to the onset'], ['recovery_or_follow_up', 'in recovery or follow-up']];
    const phase = sepsisN && Object.keys(ph).some(k => k !== 'during_sepsis') ? `By phase, ${andList(PH.filter(([k]) => ph[k]).map(([k, t]) => `${fmt(ph[k])} were drawn ${t}`))}.` : '';
    const defn = s.sepsis3 ? `${plural(s.sepsis3, 'record')} rest on a study that states the Sepsis-3 definition.` : 'The study does not state in its records which definition of sepsis it used.';
    const flags = Object.entries(s.flags);
    const times = s.times;
    return `
      <p class="lead">The ${plural(d.records, 'sample record')} of ${d.id}, by what their records say. Each label was read from the wording of a record, or from a quoted sentence of the study's description where the record is silent; nothing was inferred from the expression values${repeats ? `. ${plural(repeats, 'record')} repeat${repeats === 1 ? 's' : ''} a sample first deposited elsewhere, so the counts of records and of distinct samples differ` : ''}.</p>
      <h3>Class, beside the whole resource</h3>
      <div class="chart" id="stSaCls"></div>${clsKeys()}
      ${howto('Each bar is divided by class in proportion to its samples; hover over a part to see its count. The last bar is every distinct sample of the resource, for comparison.')}
      <p class="body">The class rests ${andList(bsum)}${basis.not_recorded ? `, and is not recorded for ${plural(basis.not_recorded, 'record')}` : ''}.${s.infection_not_stated ? ` ${plural(s.infection_not_stated, 'record')} describe${s.infection_not_stated === 1 ? 's' : ''} an infection without stating sepsis; ${s.infection_not_stated === 1 ? 'it is an other condition' : 'they are other conditions'}, marked as such.` : ''}${s.class_differs ? ` ${plural(s.class_differs, 'record')} belong${s.class_differs === 1 ? 's' : ''} to samples whose records in different studies disagree on the class; each record keeps its own label.` : ''}</p>
      <h3>What the records state, beside the resource</h3>
      <div class="chart" id="stSaCov"></div>
      <div class="keys">${keys([[`${d.id}`, C.hilite], ['All records of the resource', C.lineStrong]])}</div>
      ${howto(`The share of ${d.id}'s records that give each variable (upper bar), and the share over all ${fmt(N.records)} records of the resource (lower bar). A variable that a record does not give is marked as not recorded, never guessed.`)}
      <h3>Sample type, age group and sex</h3>
      <div class="st-facts">${factList('Sample type', s.tissue, k => lab('tissue', k))}${factList('Age group', s.age_group, k => lab('age', k))}${factList('Sex', s.sex, k => SEX[k] || k)}</div>
      ${Object.keys(s.ages).length ? `<h3>Recorded age by class <span class="m">in years</span></h3><div class="chart" id="stSaAge"></div>
        ${howto(`Each box spans the middle half of the recorded ages of a class, the line inside it is the median and the whiskers reach the youngest and the oldest. The colored box is ${d.id}; the grey box under it is every distinct sample of the resource with a recorded age.`)}` : `<p class="body">The records give no age in years${Object.keys(s.age_group).some(k => k !== 'not_recorded') ? ', only an age group' : ''}.</p>`}
      <h3>Time points</h3>
      ${times.length ? `<div class="chart" id="stSaTime"></div>${howto('Records by sampling time, in the order the records give. Days and hours are counted as the study counts them, from admission, from the onset of sepsis or from enrolment; an order without a unit is the first, second and later sample of a patient.')}
        <p class="body">${plural(s.time_records, 'record')} give a time point${s.time_patients ? `, from ${plural(s.time_patients, 'identified patient')}` : ''}. The recorded wordings include ${andList([...new Set(times.flatMap(t => t.words))].slice(0, 4).map(w => `“${esc(w)}”`))}.</p>`
        : '<p class="body">No record gives a time point.</p>'}
      <h3>Outcome, septic shock and phase</h3>
      <p class="body">${outcome} ${shock} ${phase} ${defn}</p>
      <h3>Patients and quality flags</h3>
      <p class="body">${s.patients ? `Patient identifiers link ${plural(s.with_patient, 'record')} to ${plural(s.patients, 'patient')}; ` : 'The records name no patient, so each sample counts as its own patient; '}one sample is marked for each patient (${plural(s.one_per_patient, 'record')}${s.one_per_patient_primary !== s.one_per_patient ? `, ${fmt(s.one_per_patient_primary)} of them primary records` : ''}).
        ${s.flagged ? `${plural(s.flagged, 'record')} carr${s.flagged === 1 ? 'ies' : 'y'} a quality flag: ${andList(flags.map(([k, n]) => `${esc((P.flags[k] || k).toLowerCase().replace(/^no /, 'no '))} (${fmt(n)})`))}. Flagged records are kept and marked, so that each user can decide.` : 'No record carries a quality flag.'}</p>
      <h3>The wording behind the labels</h3>
      <p class="body">The most frequent recorded wordings of ${d.id} and the curated value each gave, with one sample that carries it and the evidence: the sample record, or a sentence of the study record quoted word for word.</p>
      ${wordingHTML(d)}`;
  }
  function samplesDraw(d) {
    const rows = [{ label: `${d.id}, records`, sub: plural(d.records, 'record'), parts: clsParts(d.cls) }];
    if (d.distinct !== d.records) rows.push({ label: `${d.id}, distinct samples`, sub: `${fmt(d.distinct)} first deposited here`, parts: clsParts(d.clsd) });
    rows.push({ label: 'The resource, distinct samples', sub: fmt(N.distinct), parts: clsParts(D.order.class.map(k => N.classes[k])) });
    stackBars($('stSaCls'), rows, `Stacked bars of the classes of ${d.id} beside the whole resource.`);
    coverageChart($('stSaCov'), d);
    if ($('stSaAge')) ageChart($('stSaAge'), d);
    if ($('stSaTime')) timeChart($('stSaTime'), d);
  }
  function coverageChart(el, d) {
    const cov = d.samples.coverage, R = P.resource.coverage, keysV = Object.keys(COVER);
    const W = Math.max(300, el.clientWidth || 800), small = W < 560, lw = small ? 92 : 150, rowH = 30, H = 10 + keysV.length * rowH + 26;
    const s = svg(el, W, H, `Bars of the share of records of ${d.id} that give each variable, beside the whole resource.`), x = lin(0, 1, lw + 8, W - 50);
    gridX(s, x, [0, 0.25, 0.5, 0.75, 1], 6, H - 26, v => pct(v, 0));
    keysV.forEach((k, i) => {
      const y = 10 + i * rowH, name = { sex: 'Sex', age: 'Age', patient: 'Patient', time: 'Time point', outcome: 'Outcome', shock: 'Septic shock' }[k];
      T(s, lw, y + 12, name, { 'text-anchor': 'end', 'font-size': 12.5, fill: C.ink2 });
      hover(rrect(s, x(0), y + 1, x(cov[k]) - x(0), 10, C.hilite, 2), `<b>${name}</b><br>${d.id}: ${pctS(cov[k])} of its records`);
      hover(rrect(s, x(0), y + 13, x(R[k]) - x(0), 7, C.lineStrong, 2), `<b>${name}</b><br>Resource: ${pctS(R[k])} of all records`);
      T(s, x(cov[k]) + 5, y + 10, pctS(cov[k]), { 'font-size': 11, fill: C.ink });
    });
  }
  function ageChart(el, d) {
    const A = d.samples.ages, R = P.resource.ages, ks = D.order.class.filter(k => A[k]);
    const W = Math.max(300, el.clientWidth || 800), small = W < 560, lw = small ? 92 : 140, rowH = 40, H = 8 + ks.length * rowH + 30;
    const s = svg(el, W, H, `Box plots of the recorded age of each class in ${d.id}, beside the whole resource.`);
    const hi = Math.max(100, ...ks.map(k => A[k][5]), ...ks.map(k => (R[k] || [0, 0, 0, 0, 0, 0])[5]));
    const x = lin(0, hi, lw + 8, W - 14);
    gridX(s, x, ticks(0, hi, small ? 4 : 6), 4, H - 30);
    axisTitle(s, (lw + W) / 2, H - 4, 'Age in years');
    const box = (q, y, h, col, alpha, label) => {
      const [n, mn, q1, md, q3, mx] = q;
      S('line', { x1: x(mn), x2: x(mx), y1: y + h / 2, y2: y + h / 2, stroke: col, 'stroke-width': 1.2 }, s);
      [mn, mx].forEach(v => S('line', { x1: x(v), x2: x(v), y1: y + 2, y2: y + h - 2, stroke: col, 'stroke-width': 1.2 }, s));
      hover(rrect(s, x(q1), y, Math.max(2, x(q3) - x(q1)), h, rgba(col, alpha), 2, { stroke: col, 'stroke-width': 1 }), `<b>${label}</b><br>${plural(n, 'sample')} with a recorded age<br>median ${num(md)}; middle half ${num(q1)}–${num(q3)}; range ${num(mn)}–${num(mx)} years`);
      S('line', { x1: x(md), x2: x(md), y1: y, y2: y + h, stroke: C.ink, 'stroke-width': 2 }, s);
    };
    ks.forEach((k, i) => {
      const y = 8 + i * rowH;
      T(s, lw, y + 13, lab('class', k), { 'text-anchor': 'end', 'font-size': 12.5, fill: C.ink2 });
      box(A[k], y + 2, 14, C.cls[k], 0.55, `${d.id}, ${labLower('class', k)}`);
      if (R[k]) box(R[k], y + 21, 9, C.lineStrong, 0.4, `The resource, ${labLower('class', k)}`);
    });
  }
  function timeChart(el, d) {
    const t = d.samples.times, W = Math.max(300, el.clientWidth || 800), H = 190, m = { l: 46, r: 10, t: 14, b: 50 };
    const s = svg(el, W, H, `Bar chart of the records of ${d.id} at each time point.`);
    const max = niceMax(Math.max(...t.map(x => x.n))), y = lin(0, max, H - m.b, m.t), bw = (W - m.l - m.r) / t.length;
    gridY(s, y, ticks(0, max, 4), m.l, W - m.r);
    axisTitle(s, 12, (m.t + H - m.b) / 2, 'Records', true);
    t.forEach((x, i) => {
      const cx = m.l + i * bw;
      hover(rrect(s, cx + bw * 0.15, y(x.n), bw * 0.7, y(0) - y(x.n), C.hilite, 2), `<b>${esc(timeLabel(x))}</b><br>${plural(x.n, 'record')}<br><span class="m">recorded as ${x.words.map(w => `“${esc(w)}”`).join(', ')}</span>`);
      if (bw > 26 || i % Math.ceil(30 / bw) === 0) T(s, cx + bw / 2, H - m.b + 14, timeLabel(x).replace(' sample', ''), { 'text-anchor': bw < 46 ? 'end' : 'middle', 'font-size': 11, fill: C.muted, transform: bw < 46 ? `rotate(-40 ${cx + bw / 2} ${H - m.b + 14})` : null });
    });
  }

  // ======================================================================= chapter 3: platform and processing
  function routeText(d, p) {
    const fam = famKey(p.fam);
    if (p.tier === 'not_processed') return `This study platform was not processed: ${esc(p.basis)}.`;
    if (p.tier === 'raw_array_files' && p.affy) {
      const a = p.affy;
      return `Read from the CEL file, the scanner's own file, of every array: ${plural(p.records, `${esc(a.chip)} array`)}${a.first ? ` scanned between ${fdate(a.first)} and ${fdate(a.last)}` : ''}. Robust multi-array average (RMA) turned the ${fmt(a.probes)} probes of each array into ${fmt(a.sets)} probe-set values: it corrects each array for its background light, quantile-normalizes the arrays of the study so that they share one distribution of intensities, and combines the probes of each set by a median polish, which a single bad probe cannot move. The design of the chip came from the Bioconductor package ${esc(a.package)}.`;
    }
    if (p.tier === 'raw_array_files') {
      const a = p.single || {};
      return `Read from the scanner file of every array${a.software ? ` (${esc(a.software)})` : ''}${a.first ? `, scanned between ${fdate(a.first)} and ${fdate(a.last)}` : ''}. Each spot was corrected for the background measured around it with the normexp model of limma, which keeps dim spots small but positive; the spots of one probe were averaged, and the arrays were quantile-normalized and put on the log₂ scale.`;
    }
    if (p.tier === 'unnormalized_intensities') {
      const c = (p.inten || {}).correction, file = p.source && p.source !== 'per_sample_files' ? ` (${esc(p.source.split(';')[0])}${p.source.includes(';') ? ' and others' : ''})` : '';
      const how = c === 'detection_p_values' ? 'The background of each array was inferred from the detection p-values of the export, which say how many of the chip\'s negative-control probes were brighter than each probe; with it, limma\'s neqc corrected the background as it would with the negative controls themselves.'
        : c === 'intensity_distribution' ? 'The export holds no detection p-values, so the background of each array was estimated from the distribution of its own intensities (the RMA-75 estimate of limma) and removed with the normexp model.'
          : 'The signals had already been background-corrected by the scanner software.';
      return `Read from the table of probe intensities that the authors exported before normalizing them${file}. ${how} The arrays were then quantile-normalized and put on the log₂ scale.`;
    }
    if (p.tier === 'counts') {
      const c = p.counts || {}, src = c.source === 'ncbi_counts' ? 'the table of counts that NCBI generates for GEO studies with one pipeline of its own' : c.source === 'per_sample_files' ? 'one count file per sample' : `the authors' table of counts${p.source && !['ncbi_counts', 'per_sample_files'].includes(p.source) ? ` (${esc(p.source.split(';')[0].split('/')[0])})` : ''}`;
      const lib = c.lib ? ` Library sizes, the total count of a sample, run from ${fmt(c.lib[0] / 1e6, 1)} to ${fmt(c.lib[2] / 1e6, 1)} million reads (median ${fmt(c.lib[1] / 1e6, 1)} million).` : '';
      return `RNA-seq read counts per gene, from ${src}.${lib} The counts are released as deposited${c.kind === 'decimals' ? ', including counts that are not whole numbers, as some counting programs estimate them' : ''}, and as log₂ counts per million after TMM normalization, which corrects each library size for the few genes that can take a large share of the reads.`;
    }
    if (p.tier === 'author_normalized') {
      const a = p.author || {};
      return a.values ? `Only values normalized by the authors were deposited (${esc(a.values)}). Another laboratory's pipeline cannot be undone from its output, so they are kept exactly as deposited, with their scale recorded.` : `Only values normalized by the authors were deposited, and their tables could not be read: ${esc(p.basis)}.`;
    }
    return esc(p.basis);
  }
  function checkText(p) {
    const c = p.check; if (!c) return '';
    if (c.status === 'compared') return `As a check of which column belongs to which record, the reprocessed values were compared with the authors' series matrix: ${fmt(c.own)} of ${plural(c.compared, 'sample')} were most similar to their own column${c.exception ? '. The rest are recorded as an exception that concerns the deposit, not the processing' : ''}.`;
    if (c.status === 'no_series_matrix_values') return 'GEO holds no series matrix with values for it, so the column-by-column check against the authors\' values could not be made.';
    if (c.status === 'same_source') return 'Its values were read from the series matrix itself, so the check against it does not apply.';
    if (c.status === 'too_few_shared_rows') return 'Its series matrix shares too few rows with the reprocessed table for the column-by-column check.';
    return '';
  }
  function platformHTML(d) {
    const many = d.platforms.length > 1;
    return `<p class="lead">How ${d.id} was read from the files its authors deposited. Each study platform was processed on its own, from the lowest data level that covers its samples, and never normalized together with another study, so that adding or removing a study never changes the numbers of another.</p>
      ${d.platforms.map((p, i) => {
        const nv = Object.entries(p.no_values);
        return `<div class="st-plat"><h3>${esc(platLabel(p))} <span class="m">${esc(cleanTitle(p.title))}</span> <span class="chip"><span class="st-dot" style="background:${famColor(p.fam)}"></span>${famText(p.fam)}</span></h3>
        <p class="body"><b>Route: ${esc({ raw_array_files: 'raw array files', unnormalized_intensities: 'unnormalized intensities', counts: 'read counts', author_normalized: 'values normalized by the authors', not_processed: 'not processed' }[p.tier] || p.tier)}.</b> ${routeText(d, p)}</p>
        ${p.scale ? `<p class="body">The released values are ${SCALE[p.scale] || esc(p.scale)}${p.features ? `, in ${plural(p.features, 'row')} of ${famKey(p.fam) === 'rna_seq' ? 'deposited gene identifiers' : p.tier === 'author_normalized' ? 'deposited identifiers' : 'probes or probe sets'}${p.id_example ? ` (such as ${esc(p.id_example)})` : ''}` : ''}.</p>` : ''}
        <p class="body">${nv.length ? `${p.records ? `${fmt(p.records)} of the ${plural(p.records + p.without, 'record')} have values. The other ${fmt(p.without)} are not in the sample table; they are listed in the ledger of records without values` : `None of the ${plural(p.without, 'record')} has values; they are not in the sample table, but listed in the ledger of records without values`}, each with its reason: ${andList(nv.map(([k, n]) => `${fmt(n)} ${NO_VALUES[k] || k.replace(/_/g, ' ')}`))}.` : p.records ? `All ${plural(p.records, 'record')} have values.` : 'It holds no records.'} ${checkText(p)}</p>
        ${p.flag === 'uncertain' ? '<p class="note">The choice of this route is marked uncertain in the records of the processing stage: the deposited files left room for another reading.</p>' : ''}
        ${p.scale ? `<div class="chart" id="stPlScale${i}"></div>` : ''}</div>`;
      }).join('')}
      ${d.platforms.some(p => p.scale) ? techKeys() + howto(`Every table of the resource on the same kind of scale, sorted by its median value: each vertical bar runs from the 5th to the 95th percentile of the released values of one table, the box spans the middle half and the dot is the median. ${d.id}${many ? '\'s tables are' : '\'s table is'} picked out. The tables differ in level and spread because they come from different technologies and pipelines; the correction deals with this for analyses across studies.`) : ''}`;
  }
  const LOG_SCALES = ['log2', 'log2_like', 'log2 counts per million'];
  function scaleChart(el, d, p) {
    const logGroup = LOG_SCALES.includes(p.scale);
    const all = P.resource.scales.filter(t => LOG_SCALES.includes(t.scale) === logGroup).sort((a, b) => a.q[2] - b.q[2]);
    const W = Math.max(300, el.clientWidth || 800), H = 260, m = { l: 46, r: 10, t: 22, b: 30 };
    const s = svg(el, W, H, `Ranges of the released values of every table on ${logGroup ? 'a logarithmic' : 'another'} scale; ${d.id} ${p.p} is picked out.`);
    const lo = Math.min(...all.map(t => t.q[0])), hi = Math.max(...all.map(t => t.q[4]));
    const y = lin(Math.floor(lo), Math.ceil(hi), H - m.b, m.t), bw = (W - m.l - m.r) / all.length;
    gridY(s, y, ticks(Math.floor(lo), Math.ceil(hi), 5), m.l, W - m.r);
    axisTitle(s, 12, (m.t + H - m.b) / 2, logGroup ? 'Value (log₂ units)' : 'Value', true);
    T(s, m.l, H - 8, `${fmt(all.length)} tables, sorted by median`, { 'font-size': 11, fill: C.muted });
    const mine = t => t.id === d.id && t.p === p.p;
    all.forEach((t, i) => {
      const cx = m.l + (i + 0.5) * bw, me = mine(t), fam = (P.resource.genes.find(g => g.id === t.id && g.p === t.p) || {}).fam;
      const col = me ? C.hilite : (fam != null ? famColor(fam) : C.lineStrong), a = me ? 1 : 0.5;
      S('line', { x1: cx, x2: cx, y1: y(t.q[0]), y2: y(t.q[4]), stroke: rgba(col, a), 'stroke-width': me ? 2 : 1 }, s);
      S('rect', { x: cx - Math.max(1, bw * 0.32), y: y(t.q[3]), width: Math.max(2, bw * 0.64), height: Math.max(1, y(t.q[1]) - y(t.q[3])), fill: rgba(col, me ? 0.9 : 0.35), stroke: me ? C.ink : 'none' }, s);
      S('circle', { cx, cy: y(t.q[2]), r: me ? 3.5 : 1.6, fill: me ? C.ink : col }, s);
      hitRow(s, cx - bw / 2, m.t, bw, H - m.t - m.b, `<b>${t.id}</b> <span class="m">${esc(t.p)}</span><br>${esc(t.scale)}<br>median ${fx(t.q[2], 2)}; middle half ${fx(t.q[1], 2)}–${fx(t.q[3], 2)}; 5th–95th percentile ${fx(t.q[0], 2)}–${fx(t.q[4], 2)}`);
      if (me) T(s, clamp(cx, m.l + 30, W - m.r - 30), m.t - 8, `${d.id}${d.platforms.length > 1 ? ` ${p.p}` : ''}`, { 'text-anchor': 'middle', 'font-size': 12, 'font-weight': 600, fill: C.ink });
    });
  }
  function platformDraw(d) { d.platforms.forEach((p, i) => { if ($(`stPlScale${i}`)) scaleChart($(`stPlScale${i}`), d, p); }); }

  // ======================================================================= chapter 4: genes
  function genesHTML(d) {
    const M = P.mapping, withGenes = d.platforms.filter(p => p.genes);
    const RULE = { highest_mean: 'the row with the highest mean expression was kept for each', sum_counts: 'their counts were summed', sum: 'their values were summed' };
    return `<p class="lead">The rows of ${d.id}'s ${withGenes.length > 1 ? 'tables' : 'table'} (probes, probe sets or the gene identifiers its authors deposited) were mapped to genes of the HGNC list, which gives every human gene one approved symbol and a fixed identifier. A row whose annotation names several genes was not assigned to any of them, and the gene-level table keeps one row per gene, so that every table of the resource can be read by the same gene keys.</p>
      ${d.platforms.map((p, i) => {
        if (!p.genes) return `<div class="st-plat"><h3>${esc(platLabel(p))}</h3><p class="body">No gene-level table: ${p.tier === 'not_processed' ? 'the study platform was not processed' : 'its tables could not be read'}.</p></div>`;
        const g = p.genes, mp = g.map || [];
        const parts = M.order.map((k, j) => [k, mp[j] || 0]).filter(([, n]) => n > 0);
        const MAP_TEXT = { one_gene: 'map to one gene', several_genes: 'to several genes and were left unassigned', no_gene_in_annotation: 'have no gene in the annotation', gene_not_in_hgnc: 'name a gene outside the HGNC list', accession_not_resolved: 'carry an accession that could not be resolved to a gene', not_in_annotation: 'are missing from the annotation' };
        return `<div class="st-plat"><h3>${esc(platLabel(p))} <span class="m">${esc(cleanTitle(p.title))}</span></h3>
          <div class="cards">
            <div class="card"><div class="lab">Rows of the table</div><div class="big">${fmt(g.features)}</div><div class="lab">annotated from ${esc(g.source || 'its platform table')}</div></div>
            <div class="card"><div class="lab">Rows with one gene</div><div class="big">${fmt(g.rows_used)}</div><div class="lab">${pctS(g.rows_used / g.features)} of the rows</div></div>
            <div class="card"><div class="lab">Genes in the table</div><div class="big">${fmt(g.genes)}</div><div class="lab">${fmt(g.coding)} of them protein-coding</div></div>
          </div>
          <div class="chart" id="stGeMap${i}"></div>
          <p class="body">Of the ${plural(g.features, 'row')}, ${andList(parts.map(([k, n]) => `${fmt(n)} ${MAP_TEXT[k] || esc((M.label[k] || k).toLowerCase())}`))}. ${g.several_rows ? `${plural(g.several_rows, 'gene is', 'genes are')} measured by more than one row; ${RULE[g.rule] || 'one value was kept for each'}${g.rule === 'highest_mean' ? ' gene' : ''}.` : 'No gene is measured by more than one row.'} The row-level table is released as well, with the gene of every row.</p></div>`;
      }).join('')}
      ${withGenes.length ? `<div class="keys">${keys(M.order.map(k => [M.label[k], M.color[k]]))}</div>${howto(`Each bar divides the rows of a table by what their annotation gives: one gene of the HGNC list, several genes, no gene, a gene that is not in the HGNC list, an accession that could not be resolved, or a row missing from the annotation. The lower bar is every table of the same technology, for comparison. Rows without a gene on arrays are mostly control probes and probes for unnamed transcripts; on RNA-seq, identifiers of genes outside the HGNC list.`)}
        <h3>Genes per table across the resource</h3><div class="chart" id="stGeAll"></div>${techKeys()}
        ${howto(`Every gene-level table of the resource as a dot, by the number of genes it holds, colored by technology; ${d.id} is picked out. RNA-seq tables hold more genes because sequencing reads every gene, while an array reads only the genes its probes were designed for.`)}` : ''}`;
  }
  function genesDraw(d) {
    const M = P.mapping;
    d.platforms.forEach((p, i) => {
      if (!p.genes || !$(`stGeMap${i}`)) return;
      const parts = arr => M.order.map((k, j) => ({ key: k, label: M.label[k], n: arr[j] || 0, color: M.color[k] }));
      const fam = famKey(p.fam);
      stackBars($(`stGeMap${i}`), [{ label: `${d.id}${d.platforms.length > 1 ? ` ${p.p}` : ''}`, sub: plural(p.genes.features, 'row'), parts: parts(p.genes.map || []) },
        { label: `All ${lab('family', fam)} tables`, sub: 'share of their rows', shares: true, parts: parts(M.by_family[fam] || []) }], `Stacked bars of the mapping of rows to genes in ${d.id}, beside all tables of its technology.`);
    });
    if ($('stGeAll')) genesAllChart($('stGeAll'), d);
  }
  function genesAllChart(el, d) {
    const all = P.resource.genes.slice().sort((a, b) => a.genes - b.genes), W = Math.max(300, el.clientWidth || 800), H = 200, m = { l: 56, r: 12, t: 22, b: 26 };
    const s = svg(el, W, H, `Dot plot of the number of genes in every table; ${d.id} is picked out.`);
    const max = niceMax(Math.max(...all.map(t => t.genes))), y = lin(0, max, H - m.b, m.t), bw = (W - m.l - m.r) / all.length;
    gridY(s, y, ticks(0, max, 4), m.l, W - m.r);
    axisTitle(s, 12, (m.t + H - m.b) / 2, 'Genes', true);
    T(s, m.l, H - 6, `${fmt(all.length)} tables, sorted by their number of genes`, { 'font-size': 11, fill: C.muted });
    all.forEach((t, i) => {
      const me = t.id === d.id, cx = m.l + (i + 0.5) * bw;
      const c = S('circle', { cx, cy: y(t.genes), r: me ? 6 : 3.4, fill: me ? C.hilite : famColor(t.fam), 'fill-opacity': me ? 1 : 0.75, stroke: me ? C.ink : 'none', 'stroke-width': 1.5 }, s);
      hover(c, `<b>${t.id}</b> <span class="m">${esc(t.p)}</span><br>${sw(famColor(t.fam))}${famText(t.fam)}<br>${fmt(t.genes)} genes, ${fmt(t.coding)} protein-coding`);
      if (me) T(s, clamp(cx, m.l + 30, W - m.r - 30), y(t.genes) - 11, fmt(t.genes), { 'text-anchor': 'middle', 'font-size': 12, 'font-weight': 600, fill: C.ink });
    });
  }

  // ======================================================================= chapter 5: checks
  const sexRules = () => P.resource.sex_rules;
  const TISSUE_SHORT = k => lab('tissue', k);
  function sexHTML(d) {
    const many = d.checks.sex.length > 1;
    if (!d.checks.sex.length) return `<h3>Sex read from the genes</h3><p class="body">Not tested: ${d.id} has no table with values.</p>`;
    return d.checks.sex.map((t, i) => {
      const head = `<h3>Sex read from the genes${many ? ` <span class="m">${esc(platLabel(t))}</span>` : ''}</h3>`;
      const R = sexRules();
      const STATUS = {
        too_few_samples: `Not tested: the table holds ${plural(t.samples, 'sample')}, fewer than the ${fmt(R.minimum_samples)} the check needs.`,
        genes_missing: 'Not tested: the table does not hold the sex genes that the check reads.',
        scale_not_comparable: 'Not tested: the values of this table are on a scale that cannot be read as log₂ values.',
        one_group: 'The split found no second group: the sex genes of these samples do not fall into two groups far enough apart, so no sample was called. A study of one sex looks like this.',
        no_tissue_shows_two_sexes: 'The table holds several sample types and was split within each; none showed two groups far enough apart, so no sample was called.',
      };
      if (STATUS[t.status]) return head + `<p class="body">${STATUS[t.status]}</p>`;
      const compared = t.concordant + t.discordant;
      const seps = t.splits.map(s => s.separation), sizes = t.splits.flatMap(s => [s.lower, s.higher]);
      const splits = t.splits.length > 2 ? [`within each of its ${spell(t.splits.length)} sample types, into groups of ${fmt(Math.min(...sizes))} to ${fmt(Math.max(...sizes))} samples that lie ${fx(Math.min(...seps), 1)} to ${fx(Math.max(...seps), 1)} pooled standard deviations apart`]
        : t.splits.map(s => `${s.tissue ? `${TISSUE_SHORT(s.tissue)}: ` : ''}${fmt(s.lower)} and ${fmt(s.higher)} samples, ${fx(s.separation, 1)} pooled standard deviations apart`);
      const disc = (t.points || []).filter(p => p.a === 'discordant');
      return head + `<p class="body">For every sample the mean of the Y-chromosome genes${t.xist ? ' minus XIST' : ' (the table has no XIST)'} was computed, and the samples were split into two groups by this number alone, without the recorded sex (a two-means split): ${andList(splits)}${t.splits.length ? ` (a split stands only from ${fmt(R.minimum_separation)})` : ''}. ${fmt(t.female)} samples read as female, ${fmt(t.male)} as male${t.ambiguous ? ` and ${fmt(t.ambiguous)} lie too close to the midpoint to call` : ''}.
        ${t.recorded && compared ? `Of the ${plural(compared, 'called sample')} with a recorded sex, ${fmt(t.concordant)} agree with the record and ${fmt(t.discordant)} do not: ${pct1(t.conc)}, against ${pct1(N.sex_concordance)} across the resource.` : 'The records give no sex, so the reading could not be compared with them; it is released with every sample all the same.'}
        ${t.unreliable ? ` With fewer than ${pct(R.unreliable_below, 0)} in agreement, the recorded sex of this table is marked as unreliable: a disagreement of this size points to the deposit, such as a reversed or misattached field, rather than to single samples.` : ''}</p>
        <div class="st-grid2">${t.plane ? `<div><h4>The sex plane</h4><div class="chart" id="stChPlane${i}"></div></div>` : ''}<div><h4>The split</h4><div class="chart" id="stChSplit${i}"></div></div></div>
        <div class="keys">${keys([['Recorded female', sexColor('female'), 'dot'], ['Recorded male', sexColor('male'), 'dot'], ['Sex not recorded', sexColor('not_recorded'), 'dot'], ['Disagrees with the record', C.ink, 'ring']])}</div>
        ${howto(`${t.plane ? 'In the plane, each dot is a sample placed by its XIST (across) and its Y-chromosome genes (up), each scaled from the mean of the samples read as male (0) to that of the samples read as female or male (1) in this table: women sit at the lower right, men at the upper left. ' : ''}In the split, each dot is a sample placed by its Y-genes-minus-XIST score; the dashed line is the midpoint between the two group means, and the shaded band around it is the zone where a sample is too close to call. Colors show the recorded sex; a dark ring marks a sample whose reading disagrees with its record.`)}
        ${disc.length ? `<p class="body">${plural(disc.length, 'sample disagrees', 'samples disagree')} with ${disc.length === 1 ? 'its' : 'their'} record, each listed for inspection; a disagreement can be a mistake in the record, a sample swapped in the laboratory, or biology, such as the loss of the Y chromosome from the blood cells of some older men:</p>
          <div class="st-chips">${disc.slice(0, 40).map(p => `<span class="chip"><b>${esc(p.s)}</b> recorded ${esc(p.r)}, reads ${esc(p.e)}</span>`).join('')}${disc.length > 40 ? `<span class="chip">and ${fmt(disc.length - 40)} more</span>` : ''}</div>` : ''}`;
    }).join('');
  }
  function planeChart(el, t, id) {
    const pts = t.points.filter(p => p.x != null && p.y != null), W = Math.max(280, el.clientWidth || 500), H = Math.round(Math.min(W, 420) * 0.9), m = { l: 44, r: 10, t: 10, b: 40 };
    const s = svg(el, W, H, `Scatter plot of the samples of ${id} by XIST and the Y-chromosome genes, colored by recorded sex.`);
    const ext = a => [Math.min(-0.3, ...a), Math.max(1.3, ...a)];
    const [xa, xb] = ext(pts.map(p => clamp(p.x, -1, 2))), [ya, yb] = ext(pts.map(p => clamp(p.y, -1, 2)));
    const x = lin(xa, xb, m.l, W - m.r), y = lin(ya, yb, H - m.b, m.t);
    gridX(s, x, [0, 0.5, 1], m.t, H - m.b, v => fx(v, 1)); gridY(s, y, [0, 0.5, 1], m.l, W - m.r, v => fx(v, 1));
    axisTitle(s, (m.l + W - m.r) / 2, H - 6, 'XIST, male (0) to female (1) mean');
    axisTitle(s, 12, (m.t + H - m.b) / 2, 'Y genes, female (0) to male (1)', true);
    const r = pts.length > 400 ? 2.6 : 3.4;
    pts.slice().sort((a, b) => (a.a === 'discordant') - (b.a === 'discordant')).forEach(p => {
      const rec = p.r === 'female' || p.r === 'male' ? p.r : 'not_recorded', bad = p.a === 'discordant';
      hover(S('circle', { cx: x(clamp(p.x, xa, xb)), cy: y(clamp(p.y, ya, yb)), r: bad ? r + 1.4 : r, fill: sexColor(rec), 'fill-opacity': 0.7, stroke: bad ? C.ink : 'none', 'stroke-width': bad ? 1.6 : 0 }, s),
        `<b>${esc(p.s)}</b><br>recorded ${esc(SEX[p.r] ? SEX[p.r].toLowerCase() : p.r)}, reads ${esc(p.e)}<br><span class="m">score ${fx(p.d, 2)}</span>`);
    });
  }
  function splitChart(el, t, id) {
    const R = sexRules(), splits = t.splits.length ? t.splits : [{ tissue: '', midpoint: null, gap: null }];
    const W = Math.max(280, el.clientWidth || 500), rowH = splits.length > 2 ? 60 : 74, m = { l: 12, r: 12, t: 8, b: 36 }, H = m.t + splits.length * (rowH + (splits[0].tissue ? 12 : 0)) + m.b;
    const s = svg(el, W, H, `Strip of the sex scores of the samples of ${id}, with the two-means split.`);
    const pts = t.points || [], ds = pts.map(p => p.d);
    const lo = Math.min(...ds), hi = Math.max(...ds), pad = (hi - lo) * 0.05 || 1;
    const x = lin(lo - pad, hi + pad, m.l, W - m.r);
    gridX(s, x, ticks(lo - pad, hi + pad, W < 420 ? 4 : 6), m.t, H - m.b, v => fx(v, 0));
    axisTitle(s, W / 2, H - 4, t.xist ? 'Y-chromosome genes minus XIST (log₂ units)' : 'Mean of the Y-chromosome genes (log₂ units)');
    const rr = rng(11);
    splits.forEach((sp, i) => {
      const y0 = m.t + i * (rowH + (sp.tissue ? 12 : 0)) + (sp.tissue ? 12 : 0), mid = y0 + rowH / 2 + 4;
      if (sp.tissue) T(s, m.l, y0 + 2, TISSUE_SHORT(sp.tissue), { 'font-size': 11.5, 'font-weight': 600, fill: C.ink2 });
      if (sp.midpoint != null) {
        const z = R.ambiguous_zone * sp.gap;
        S('rect', { x: x(sp.midpoint - z), y: y0 + 14, width: x(sp.midpoint + z) - x(sp.midpoint - z), height: rowH - 18, fill: C.band }, s);
        S('line', { x1: x(sp.midpoint), x2: x(sp.midpoint), y1: y0 + 14, y2: y0 + rowH - 4, stroke: C.ink2, 'stroke-dasharray': '4 3' }, s);
        [[sp.midpoint - sp.gap / 2, 'female mean'], [sp.midpoint + sp.gap / 2, 'male mean']].forEach(([v, l]) => {
          S('line', { x1: x(v), x2: x(v), y1: y0 + 14, y2: y0 + 20, stroke: C.ink, 'stroke-width': 1.5 }, s);
          if (W > 420 && i === 0 && !sp.tissue) T(s, x(v), y0 + 11, l, { 'text-anchor': 'middle', 'font-size': 10.5, fill: C.muted });
        });
      }
      pts.filter(p => !sp.tissue || p.t === sp.tissue).sort((a, b) => (a.a === 'discordant') - (b.a === 'discordant')).forEach(p => {
        const rec = p.r === 'female' || p.r === 'male' ? p.r : 'not_recorded', bad = p.a === 'discordant';
        hover(S('circle', { cx: x(p.d), cy: mid + (rr() - 0.5) * (rowH - 34), r: bad ? 4 : 2.8, fill: sexColor(rec), 'fill-opacity': 0.7, stroke: bad ? C.ink : 'none', 'stroke-width': bad ? 1.5 : 0 }, s),
          `<b>${esc(p.s)}</b><br>score ${fx(p.d, 2)}; reads ${esc(p.e)}<br>recorded ${esc(SEX[p.r] ? SEX[p.r].toLowerCase() : p.r)}`);
      });
    });
  }
  function repeatsHTML(d) {
    const R = d.checks.repeats;
    if (!R.length) return `<h3>Two deposits of one sample</h3><p class="body">No record of ${d.id} is linked to a second deposit of the same array or library in a table of another study, so this check does not apply.</p>`;
    const pairs = sum(R.map(r => r.pairs));
    return `<h3>Two deposits of one sample</h3>
      <p class="body">${plural(pairs, 'record')} of ${d.id} ${pairs === 1 ? 'is' : 'are'} linked to a second deposit of the same array or library. The two deposits were processed apart, each with the other samples of its own study, so their values need not be equal, but the order of the genes from low to high should be nearly the same. Spearman's correlation compares that order: 1 means the same order.</p>
      <div class="tablewrap"><table class="t"><thead><tr><th>Second deposit in</th><th class="r">Pairs</th><th class="r">Lowest Spearman</th><th class="r">Median</th><th class="r">Each other's closest match</th></tr></thead><tbody>
      ${R.map(r => `<tr><td>${r.id === d.id ? 'this study' : studyLink(r.id)}</td><td class="r">${fmt(r.pairs)}</td><td class="r">${fx(r.lowest, 3)}</td><td class="r">${fx(r.median, 3)}</td><td class="r">${fmt(r.best)} of ${fmt(r.pairs)}</td></tr>`).join('')}
      </tbody></table></div>`;
  }
  function rocPoints(pos, neg) {
    const all = pos.map(v => [v, 1]).concat(neg.map(v => [v, 0])).sort((a, b) => b[0] - a[0]), out = [[0, 0]];
    let tp = 0, fp = 0;
    for (let i = 0; i < all.length;) { let j = i; while (j < all.length && all[j][0] === all[i][0]) { if (all[j][1]) tp++; else fp++; j++; } out.push([fp / neg.length, tp / pos.length]); i = j; }
    return out;
  }
  function scoreHTML(d) {
    const min = P.resource.score_minimum, strata = d.checks.scores, ns = d.checks.not_scored, all = P.resource.scores;
    const intro = `<h3>The sepsis score</h3><p class="body">The Sepsis MetaScore combines eleven genes, six that rise and five that fall in sepsis, by the formula of its authors. Within one study it should be higher in sepsis than in healthy controls. It was computed for every table and sample type that holds at least ${fmt(min)} sepsis samples drawn during sepsis and ${fmt(min)} healthy controls. The comparison tests whether known biology shows in the data, not the clinical value of the score.</p>`;
    const med = median(all.map(x => x.auroc));
    const body = strata.map((t, i) => `<h4>${lab('tissue', t.tissue)}${d.platforms.length > 1 ? ` · ${esc(t.p)}` : ''}</h4>
      <p class="body">Over ${clsCount(t.sepsis, 'sepsis')} samples and ${plural(t.healthy, 'healthy control')}, the area under the ROC curve (AUROC) is <b>${fx(t.auroc, 3)}</b>${t.low != null ? ` (95% interval ${fx(t.low, 3)}–${fx(t.high, 3)}, from a bootstrap over patients)` : ''}: the chance that a random sepsis sample scores above a random healthy control. ${t.sepsis_pp !== t.sepsis || t.healthy_pp !== t.healthy ? `Over one sample per patient (${fmt(t.sepsis_pp)} and ${fmt(t.healthy_pp)}) it is ${fx(t.auroc_pp, 3)}. ` : ''}${t.left_out ? `${plural(t.left_out, 'sepsis sample')} drawn before the onset or at an unknown time ${t.left_out === 1 ? 'was' : 'were'} left out. ` : ''}Across the ${fmt(all.length)} tables and sample types tested, the median AUROC is ${fx(med, 3)}. ${t.genes ? `${fmt(t.as_published)} of the ${plural(t.genes, 'gene')} of the score${' '}and HLA-DRA that the table holds changed in their published direction.` : ''}${t.role ? ` ${d.id} was one of the ${t.role} cohorts of the score's authors, so this result is not independent of the score.` : ''}</p>
      <div class="st-grid2"><div><h4>Score by class</h4><div class="chart" id="stChScore${i}"></div></div><div><h4>ROC curve</h4><div class="chart" id="stChRoc${i}"></div></div></div>`).join('');
    const where = n => `the ${labLower('tissue', n.tissue)}${d.platforms.length > 1 ? ` of ${esc(n.p || 'the table')}` : ''}`;
    const few = ns.filter(n => n.comparable), other = ns.filter(n => !n.comparable);
    const notes = (few.length === 1 ? `<p class="body">Not computed for ${where(few[0])}: it holds ${clsCount(few[0].sepsis, 'sepsis')} samples drawn during sepsis and ${plural(few[0].healthy, 'healthy control')}, and the score needs at least ${fmt(min)} of each.</p>`
      : few.length ? `<p class="body">Not computed for ${andList(few.map(n => `${where(n)} (${fmt(n.sepsis)} and ${fmt(n.healthy)})`))}: the numbers are the sepsis samples drawn during sepsis and the healthy controls of each, and the score needs at least ${fmt(min)} of each.</p>` : '')
      + other.map(n => `<p class="body">Not computed for ${where(n)}: its values (${esc(n.scale || 'of another scale')}) cannot be read on the log₂ scale.</p>`).join('');
    const none = !strata.length && !ns.length ? `<p class="body">Not computed: ${d.id} has no table with values.</p>` : '';
    return intro + body + notes + none + (strata.length ? `${howto('In the score panel each dot is a sample, sepsis above and healthy controls below; the more the two rows are apart, the better the score separates them. The ROC curve runs over every threshold of the score: across, the share of healthy controls above it; up, the share of sepsis samples above it. A curve along the diagonal is chance, one through the upper left corner perfect separation.')}
      <h4>Among all tables tested</h4><div class="chart" id="stChAll"></div>${howto(`Each dot is the AUROC of one table and sample type of the resource; ${d.id} is picked out. Dots to the right separate sepsis from healthy controls better.`)}` : '');
  }
  function scoreStrip(el, t, id) {
    const W = Math.max(280, el.clientWidth || 500), H = 200, m = { l: 74, r: 12, t: 12, b: 36 };
    const s = svg(el, W, H, `Strip plot of the sepsis score of the sepsis samples and healthy controls of ${id}.`);
    const all = t.pos.concat(t.neg), lo = Math.min(...all), hi = Math.max(...all), pad = (hi - lo) * 0.05 || 1;
    const x = lin(lo - pad, hi + pad, m.l, W - m.r);
    gridX(s, x, ticks(lo - pad, hi + pad, W < 420 ? 4 : 6), m.t, H - m.b, v => fx(v, 1));
    axisTitle(s, (m.l + W) / 2, H - 4, 'Sepsis MetaScore');
    const rr = rng(5), rowH = (H - m.t - m.b) / 2;
    [['sepsis', t.pos, t.pos_s], ['healthy_control', t.neg, t.neg_s]].forEach(([k, vals, ids], j) => {
      const cy = m.t + rowH * (j + 0.5);
      T(s, m.l - 8, cy + 4, lab('class', k), { 'text-anchor': 'end', 'font-size': 12, fill: C.ink2 });
      const r = vals.length > 300 ? 2.2 : 3;
      vals.forEach((v, q) => hover(S('circle', { cx: x(v), cy: cy + (rr() - 0.5) * rowH * 0.7, r, fill: C.cls[k], 'fill-opacity': 0.6 }, s), `<b>${esc(ids[q])}</b><br>${lab('class', k)}; score ${fx(v, 2)}`));
      const md = median(vals);
      S('line', { x1: x(md), x2: x(md), y1: cy - rowH * 0.42, y2: cy + rowH * 0.42, stroke: C.ink, 'stroke-width': 2 }, s);
    });
  }
  function rocChart(el, t, id) {
    const W = Math.max(240, Math.min(el.clientWidth || 400, 330)), H = W, m = { l: 44, r: 12, t: 12, b: 40 };
    const s = svg(el, W, H, `ROC curve of the sepsis score in ${id}.`); s.style.maxWidth = `${W}px`;
    const x = lin(0, 1, m.l, W - m.r), y = lin(0, 1, H - m.b, m.t);
    gridX(s, x, [0, 0.5, 1], m.t, H - m.b, v => fx(v, 1)); gridY(s, y, [0, 0.5, 1], m.l, W - m.r, v => fx(v, 1));
    axisTitle(s, (m.l + W - m.r) / 2, H - 6, 'Healthy controls above the threshold');
    axisTitle(s, 12, (m.t + H - m.b) / 2, 'Sepsis samples above', true);
    S('line', { x1: x(0), y1: y(0), x2: x(1), y2: y(1), stroke: C.lineStrong, 'stroke-dasharray': '4 4' }, s);
    const pts = rocPoints(t.pos, t.neg);
    S('path', { d: 'M' + pts.map(([a, b], i) => `${i ? 'L' : ''}${x(a).toFixed(1)},${y(b).toFixed(1)}`).join(' ') + ` L${x(1)},${y(1)} L${x(1)},${y(0)} Z`, fill: rgba(C.cls.sepsis, 0.12), stroke: 'none' }, s);
    S('path', { d: 'M' + pts.map(([a, b], i) => `${i ? 'L' : ''}${x(a).toFixed(1)},${y(b).toFixed(1)}`).join(' '), fill: 'none', stroke: C.cls.sepsis, 'stroke-width': 2 }, s);
    const a = auroc(t.pos, t.neg);
    T(s, x(0.97), y(0.06), `AUROC ${fx(a, 3)}`, { 'text-anchor': 'end', 'font-size': 13, 'font-weight': 600, fill: C.ink });
  }
  function checksHTML(d) {
    return `<p class="lead">Three checks compare ${d.id}'s data with facts that the curation did not decide: the sex that the genes of each sample show, the agreement of two deposits of one sample, and the behavior of a published sepsis score. Each works inside the study's own table: nothing is pooled with other studies, no value is changed and a disagreement is reported, never removed.</p>
      ${sexHTML(d)}${repeatsHTML(d)}${scoreHTML(d)}
      ${d.checks.markers.length ? `<h3>Marker genes of sorted cells</h3><p class="body">In a study of sorted cells, the marker genes of each cell type should be higher in the fraction they mark than in the other fractions of the same study. The AUROC says how well they separate them.</p>
        <div class="tablewrap"><table class="t"><thead><tr><th>Cell type</th><th class="r">Markers read</th><th class="r">Samples</th><th class="r">Other samples</th><th class="r">AUROC</th></tr></thead><tbody>
        ${d.checks.markers.map(k => `<tr><td>${lab('tissue', k.cell)}</td><td class="r">${fmt(k.markers)}</td><td class="r">${fmt(k.samples)}</td><td class="r">${fmt(k.others)}</td><td class="r">${fx(k.auroc, 3)}</td></tr>`).join('')}</tbody></table></div>` : ''}`;
  }
  function checksDraw(d) {
    d.checks.sex.forEach((t, i) => { if ($(`stChPlane${i}`)) planeChart($(`stChPlane${i}`), t, d.id); if ($(`stChSplit${i}`)) splitChart($(`stChSplit${i}`), t, d.id); });
    d.checks.scores.forEach((t, i) => { scoreStrip($(`stChScore${i}`), t, d.id); rocChart($(`stChRoc${i}`), t, d.id); });
    if ($('stChAll')) {
      const all = P.resource.scores, lo = Math.min(0.5, ...all.map(x => x.auroc));
      stripChart($('stChAll'), all, {
        get: x => x.auroc, lo: Math.floor(lo * 10) / 10, hi: 1, axis: 'AUROC of the sepsis score, sepsis against healthy controls', label: `Strip of the AUROC of every table tested; ${d.id} is picked out.`,
        isMine: x => x.id === d.id, color: () => C.cls.sepsis, fmtv: v => fx(v, 2),
        tipFn: x => `<b>${x.id}</b> <span class="m">${esc(x.p)} · ${labLower('tissue', x.tissue)}</span><br>AUROC ${fx(x.auroc, 3)}${x.low != null ? ` (${fx(x.low, 3)}–${fx(x.high, 3)})` : ''}<br>${clsCount(x.sepsis, 'sepsis')}, ${plural(x.healthy, 'healthy control')}`,
      });
    }
  }

  // ======================================================================= chapter 6: correction
  function unitText(d, u) {
    const minH = P.resource.minimum_healthy, g = GROUP[u.group] || u.group;
    if (u.role === 'check_only') return `Its ${plural(u.records, 'record')} are all repeat deposits of samples whose first record is kept in another study, so this study platform was corrected only to compare the two deposits of each sample. Its records are not in the corrected table.`;
    if (u.status === 'left_out' || u.mstatus === 'left_out') {
      const why = u.mreason || u.reason;
      const W = {
        newborn_study_not_corrected: 'Not corrected: most of its samples are of newborns. No study with healthy controls is of newborns, so a correction of newborns could not be tested with hidden controls, and lining newborn sepsis up with older patients would remove what is specific to the newborn response.',
        no_class_to_line_up: 'Not corrected: its samples are all of other conditions, so matched groups have no class to line up with the reference studies.',
        study_has_too_few_samples: `Not corrected: fewer than ${fmt(P.resource.minimum_unit)} of its samples could estimate its effects.`,
        no_check_applicable: 'Corrected by matched groups, then left out: no study-mixing check applied to it, because it has too few samples of any class shared with the reference studies, so the correction could not be checked.',
      };
      return W[why] || `Left out of the corrected ${g} table${why ? ` (${esc(why.replace(/_/g, ' '))})` : ''}.`;
    }
    if (u.method === 'healthy_anchored') return `Corrected on its ${plural(u.anchors, 'healthy control')}, as COCONUT does: the level and spread of each gene in the study were estimated from the controls alone, shrunk toward those of the other studies (empirical Bayes), and removed from every sample. The patients took no part in the estimate, so their differences from the controls are kept. ${plural(u.estimation, 'sample')} estimated the study's effects, over ${fmt(u.genes)} genes.`;
    if (u.method === 'matched_groups') return `With ${u.anchors ? `only ${plural(u.anchors, 'healthy control')}` : 'no healthy controls'}, fewer than the ${fmt(minH)} needed to correct on them, it was corrected by matched groups: the mean and spread of each of its classes were set to those of the same class and age group in the studies corrected on healthy controls, from ${plural(u.estimation, 'sample')} over ${fmt(u.genes)} genes. This removes real differences in severity or timing along with the study effect, so part of the separation of classes across such studies comes from the correction itself.`;
    return esc(u.status);
  }
  const CHECK = { study_mixing: 'Study mixing', sex_genes: 'Sex genes', sepsis_score: 'Sepsis score' };
  const checkRule = c => c.check === 'study_mixing' ? `rise ≤ ${fx(c.limit, 2)}` : `≥ ${fx(c.limit, 2)}`;
  function checksTable(u) {
    const rows = u.checks.filter(c => c.applicable);
    if (!rows.length) return '<p class="note">No check applied to this study platform.</p>';
    return `<div class="tablewrap"><table class="t"><thead><tr><th>Check</th><th>Class</th><th class="r">Samples</th><th class="r">Before</th><th class="r">After</th><th>Rule</th><th>Result</th></tr></thead><tbody>
      ${rows.map(c => `<tr><td>${CHECK[c.check] || c.check}</td><td>${c.cls ? lab('class', c.cls) : '–'}</td><td class="r">${fmt(c.samples)}</td><td class="r">${fx(c.before, 3)}</td><td class="r">${fx(c.after, 3)}</td><td class="small" style="white-space:nowrap">${checkRule(c)}</td>
        <td>${c.passed ? '<span class="st-pass">passed</span>' : c.gate ? '<span class="st-fail">failed</span>' : '<span class="st-fail">beyond the limit</span>'}${c.gate ? '' : ' <span class="small">(reported only)</span>'}</td></tr>`).join('')}</tbody></table></div>
      <p class="note">Study mixing is the excess share of a sample's nearest neighbors of the same class that come from its own study; it should fall. The sex genes check the agreement of one split of the Y genes against XIST over all studies with the sex read within this study; the sepsis score checks that its sepsis samples still score above the healthy controls of the other studies (AUROC). For studies corrected on healthy controls, mixing is reported but cannot remove the study, since these studies are the reference.</p>`;
  }
  function correctionHTML(d) {
    const C2 = d.correction, inT = C2.in_tables, units = C2.units;
    const inAny = Object.values(inT).some(n => n > 0);
    const led = C2.ledger;
    const ledgerText = andList(led.map(l => `${fmt(l.n)} ${l.outcome === 'in_corrected_table' ? `in the corrected ${GROUP[l.group] || l.group} table` : l.outcome === 'corrected_for_checks_only' ? 'corrected only to compare repeat deposits' : `left out, ${LEDGER[l.reason] || esc(l.reason.replace(/_/g, ' '))}`}`));
    const lead = inAny ? `${andList(Object.entries(inT).filter(([, n]) => n).map(([g, n]) => `${plural(n, 'record')} of ${d.id} ${n === 1 ? 'is' : 'are'} in the corrected ${GROUP[g]} table`))}. The correction removes the differences between studies so that samples of many studies share one scale; the per-study table of ${d.id} is kept as processed, and remains the reference for any result that depends on the correction.`
      : `None of the records of ${d.id} is in a corrected table. Its per-study table is kept as processed and can be analyzed on its own.`;
    const held = units.filter(u => u.held.length);
    return `<p class="lead">${lead}</p>
      <h3>Every record of the study in the correction</h3>
      <div class="chart" id="stCoLedger"></div>
      <p class="body">${led.length === 1 ? `All ${plural(d.records, 'record')} of ${d.id} are ${ledgerText.replace(/^[\d,]+ /, '')}` : `Of its ${plural(d.records, 'record')}, ${ledgerText}`}.</p>
      ${units.map((u, i) => `<div class="st-plat"><h3>${esc(platLabel(u))} <span class="m">${GROUP[u.group] ? `${GROUP[u.group]} group` : ''} · ${u.method ? (u.method === 'healthy_anchored' ? 'corrected on healthy controls' : 'matched groups') : 'not corrected'}</span></h3>
        <p class="body">${unitText(d, u)}</p>${u.status === 'in_unit' && u.mstatus === 'corrected' ? checksTable(u) : u.checks.some(c => c.applicable) ? checksTable(u) : ''}</div>`).join('')}
      ${!units.length ? `<p class="body">No study platform of ${d.id} entered the correction: ${led.length ? `its records are ${andList(led.map(l => `${LEDGER[l.reason] || esc(l.reason.replace(/_/g, ' '))} (${fmt(l.n)})`))}` : 'it has no records with values'}.</p>` : ''}
      ${Object.entries(inT).filter(([, n]) => n).map(([g]) => `<h3>Its samples in the corrected ${GROUP[g]} table, before and after</h3>
        <div class="st-grid2"><div class="st-canvas"><canvas id="stCoV0-${g}" role="img"></canvas><span class="st-cap">Before correction</span></div><div class="st-canvas"><canvas id="stCoV1-${g}" role="img"></canvas><span class="st-cap">After correction</span></div></div>
        ${clsKeys()}
        ${howto(`Each dot is a sample of the corrected ${GROUP[g]} table, placed by the first two principal components of its most variable genes; ${d.id}'s samples are colored by class and the others are grey. Before the correction a study usually sits apart, as an island of its own; after it, its samples should lie among those of other studies of the same class. Hover over a colored dot to see its sample.`)}`).join('')}
      ${held.length ? `<h3>The held-out test</h3>${held.map((u, i) => heldHTML(d, u, i)).join('')}` : ''}`;
  }
  function heldHTML(d, u, i) {
    const by = Object.fromEntries(u.held.map(h => [h.method, h])), none = by.none, mg = by.matched_groups, R = (P.resource.held_out[u.group] || {});
    const units = (R.none || {}).units;
    return `<p class="body">${d.id}${d.platforms.length > 1 ? ` (${esc(u.p)})` : ''} was one of the ${units ? fmt(units) : ''} ${GROUP[u.group]} study platforms with healthy controls held out in turn: its ${plural(u.held[0].hidden, 'healthy control')} were hidden, the study was corrected from its other ${plural(u.held[0].estimation, 'sample')} by each method, and the hidden controls should then land among the healthy controls of the other studies. The location error is their remaining distance, in log₂ units averaged over genes.${mg && none ? ` Matched groups ${mg.location < none.location ? 'brought them closer than' : 'did not bring them closer than'} no correction (${fx(mg.location, 2)} against ${fx(none.location, 2)}); the sampling floor, the error expected from sampling alone, is ${fx(u.held[0].floor, 2)}.` : ''}</p>
      <div class="chart" id="stCoHeld${i}"></div>
      ${howto(`Bars are ${d.id}'s location error under each method; the ring on each bar is the median over all held-out study platforms, and the dashed line is the sampling floor of ${d.id}. Shorter bars are better.`)}`;
  }
  function heldChart(el, d, u) {
    const order = P.methods.filter(m => u.held.some(h => h.method === m)), R = P.resource.held_out[u.group] || {};
    const W = Math.max(300, el.clientWidth || 700), small = W < 560, lw = small ? 100 : 140, rowH = 30, H = 10 + order.length * rowH + 30;
    const s = svg(el, W, H, `Bars of the location error of ${d.id} under each correction method in the held-out test.`);
    const max = niceMax(Math.max(...u.held.map(h => h.location), ...order.map(m => (R[m] || {}).location || 0)));
    const x = lin(0, max, lw + 8, W - 50);
    gridX(s, x, ticks(0, max, 5), 6, H - 30, v => fx(v, 1));
    axisTitle(s, (lw + W) / 2, H - 4, small ? 'Location error (log₂ units)' : 'Location error of the hidden healthy controls (log₂ units)');
    const floor = u.held[0].floor;
    order.forEach((mth, i) => {
      const h = u.held.find(z => z.method === mth), y = 10 + i * rowH, med = (R[mth] || {}).location;
      T(s, lw, y + 13, METHOD[mth] || mth, { 'text-anchor': 'end', 'font-size': 12.5, fill: C.ink2 });
      hover(rrect(s, x(0), y + 3, x(h.location) - x(0), 14, mth === 'matched_groups' ? C.hilite : C.lineStrong, 2), `<b>${METHOD[mth]}</b><br>${d.id}: ${fx(h.location, 3)}<br>median of all study platforms: ${fx(med, 3)}<br><span class="m">dispersion left among the training studies ${fx(h.dispersion, 3)}</span>`);
      T(s, x(h.location) + 6, y + 14, fx(h.location, 2), { 'font-size': 11, fill: C.ink });
      if (med != null) S('circle', { cx: x(med), cy: y + 10, r: 4.5, fill: 'none', stroke: C.ink, 'stroke-width': 1.5 }, s);
    });
    S('line', { x1: x(floor), x2: x(floor), y1: 4, y2: H - 30, stroke: C.ink2, 'stroke-dasharray': '3 3' }, s);
  }
  // the whole corrected table before and after, with the study's samples in color (bulk blood: the start page's points)
  function viewOf(g) { return g === 'bulk_blood' ? D.hero : P.pbmc_view; }
  function drawViews(d, g) {
    const V = viewOf(g), si = V.studies.indexOf(d.id); if (si < 0) return;
    const n = V.x0.length, k = V.scale, idx = [];
    for (let i = 0; i < n; i++) if (V.study[i] === si) idx.push(i);
    const ext = (a, b) => { let lo = Infinity, hi = -Infinity; for (let i = 0; i < n; i++) { lo = Math.min(lo, a[i], b[i]); hi = Math.max(hi, a[i], b[i]); } const p = (hi - lo) * 0.04; return [(lo - p) / k, (hi + p) / k]; };
    const XD = ext(V.x0, V.x1), YD = ext(V.y0, V.y1);
    [0, 1].forEach(st => {
      const cv = $(`stCoV${st}-${g}`); if (!cv) return;
      const w = cv.clientWidth || 400, m = { l: 36, r: 10, t: 34, b: 34 };
      const h = Math.round(clamp((w - m.l - m.r) * (YD[1] - YD[0]) / (XD[1] - XD[0]) + m.t + m.b, 220, 520));
      const { ctx } = canvasFit(cv, { height: h });
      const kk = Math.min((w - m.l - m.r) / (XD[1] - XD[0]), (h - m.t - m.b) / (YD[1] - YD[0]));
      const cx = m.l + (w - m.l - m.r) / 2, cy = m.t + (h - m.t - m.b) / 2, xm = (XD[0] + XD[1]) / 2, ym = (YD[0] + YD[1]) / 2;
      const sx = v => cx + (v / k - xm) * kk, sy = v => cy - (v / k - ym) * kk;
      const X = st ? V.x1 : V.x0, Y = st ? V.y1 : V.y0;
      ctx.clearRect(0, 0, w, h);
      ctx.font = `11px ${getComputedStyle(document.body).fontFamily}`; ctx.fillStyle = C.muted; ctx.textAlign = 'center';
      const vv = (V.variance || {})[st ? 'after' : 'before'];
      ctx.fillText(`Component 1${vv ? ` (${pct(vv[0], 0)} of the variation)` : ''}`, cx, h - 10);
      ctx.save(); ctx.translate(13, cy); ctx.rotate(-Math.PI / 2); ctx.fillText(`Component 2${vv ? ` (${pct(vv[1], 0)})` : ''}`, 0, 0); ctx.restore();
      const grey = dotSprite(C.lineStrong, 1.6, 0.4);
      for (let i = 0; i < n; i++) if (V.study[i] !== si) stamp(ctx, grey, sx(X[i]), sy(Y[i]));
      const sprites = D.order.class.map(c => dotSprite(C.cls[c], 3.2, 0.85));
      idx.forEach(i => stamp(ctx, sprites[+V.cls[i]], sx(X[i]), sy(Y[i])));
      cv.setAttribute('aria-label', `Scatter plot of the ${fmt(n)} samples of the corrected ${GROUP[g]} table ${st ? 'after' : 'before'} the correction, with the ${fmt(idx.length)} samples of ${d.id} colored by class.`);
      cv.onmousemove = e => {
        const r = cv.getBoundingClientRect(), px = e.clientX - r.left, py = e.clientY - r.top; let best = -1, bd = 64;
        idx.forEach(i => { const dd = (sx(X[i]) - px) ** 2 + (sy(Y[i]) - py) ** 2; if (dd < bd) { bd = dd; best = i; } });
        if (best >= 0) { const c = D.order.class[+V.cls[best]]; tip(e, `<b>${esc(V.sample[best])}</b> <span class="m">· ${d.id}</span><br>${sw(C.cls[c])}${lab('class', c)}`); } else untip();
      };
      cv.onmouseleave = untip;
    });
  }
  function correctionDraw(d) {
    const led = d.correction.ledger, OUT = { in_corrected_table: 0, corrected_for_checks_only: 1, left_out: 2 };
    const parts = led.slice().sort((a, b) => OUT[a.outcome] - OUT[b.outcome] || b.n - a.n).map((l, i) => ({
      key: `${l.outcome}${l.reason}`, n: l.n,
      label: l.outcome === 'in_corrected_table' ? `In the corrected ${GROUP[l.group] || l.group} table` : l.outcome === 'corrected_for_checks_only' ? 'Corrected only to compare repeat deposits' : `Left out: ${LEDGER[l.reason] || l.reason.replace(/_/g, ' ')}`,
      color: l.outcome === 'in_corrected_table' ? (l.method === 'healthy_anchored' ? C.tier.anchored : C.tier.unanchored) : l.outcome === 'corrected_for_checks_only' ? mixC(C.tier.anchored, 0.55) : mixHex(C.lineStrong, C.bg, i % 2 ? 0.35 : 0),
    }));
    stackBars($('stCoLedger'), [{ label: `${d.id}`, sub: plural(d.records, 'record'), parts }], `Stacked bar of the outcome of every record of ${d.id} in the correction.`);
    Object.entries(d.correction.in_tables).forEach(([g, nn]) => { if (nn) drawViews(d, g); });
    d.correction.units.filter(u => u.held.length).forEach((u, i) => { if ($(`stCoHeld${i}`)) heldChart($(`stCoHeld${i}`), d, u); });
  }

  // ======================================================================= chapter 7: benchmark
  function rolesText(b, r) {
    const R = P.resource, min = R.bench_minimum, tech = (b.groups.find(g => g.group === r.group) || {}).tech;
    const first = (b.groups.find(g => g.group === r.group) || {}).first;
    if (r.scheme === 'leave_one_group_out') return r.test ? `held out and tested once, in its own fold, and trains in the other ${plural(r.train, 'fold')}` : r.held_out_not_scored ? `held out once but not scored, and trains in the other ${plural(r.train, 'fold')}` : `trains in all ${plural(r.train, 'fold')}; with fewer than ${fmt(min)} samples of a class, it is never held out`;
    if (r.scheme === 'technology_transfer') { const role = r.test ? ROLE.test : r.held_out_not_scored ? ROLE.held_out_not_scored : ROLE.train; return `as ${tech === 'rna_seq' ? 'an RNA-seq' : 'an array'} group it trains when ${tech === 'rna_seq' ? 'RNA-seq trains' : 'arrays train'} and ${role} when ${tech === 'rna_seq' ? 'arrays train' : 'RNA-seq trains'}`; }
    if (r.scheme === 'temporal') { const role = r.test ? ROLE.test : r.held_out_not_scored ? ROLE.held_out_not_scored : ROLE.train; return `first released on ${fdate(first)}, ${first < R.bench_split_date ? 'before' : 'after'} ${fdate(R.bench_split_date)}, so it ${role}`; }
    return '';
  }
  function benchHTML(d) {
    const R = P.resource, others = d.psg_studies;
    return `<p class="lead">The benchmark defines ${spell(R.tasks.length)} tasks on the per-study tables. Every split puts whole patient-safe groups on one side, so that no deposit of a sample and no identified patient lies in both training and test. ${d.id} belongs to patient-safe group${d.psg.length > 1 ? 's' : ''} ${andList(d.psg)}${others.length ? `, together with ${andList(others.map(studyLink))}` : ', on its own'}.</p>
      ${d.bench.map((b, i) => {
        const [pl, nl] = TASK_SIDES[b.task], rs = Object.entries(b.reasons);
        return `<h3>${TASK[b.task] || b.task}</h3>
          <div class="chart" id="stBe${i}"></div>
          <p class="body">${b.eligible ? `Of its ${plural(b.records, 'record')}, ${fmt(b.eligible)} ${b.eligible === 1 ? 'is a sample' : 'are samples'} of the task: ${fmt(b.pos)} ${pl} and ${fmt(b.neg)} ${nl}.` : `None of its ${plural(b.records, 'record')} is a sample of the task.`} ${rs.length ? `The others are left out where they ${andList(rs.map(([k, n]) => `${BENCH_REASON[k] || k.replace(/_/g, ' ')} (${fmt(n)})`))}, each for the first reason that applies.` : ''}
          ${b.groups.map(g => `Group ${g.group} gives the task ${fmt(g.pos)} ${pl} and ${fmt(g.neg)} ${nl}.`).join(' ')}</p>
          ${b.roles.length ? `<div class="chart" id="stBeF${i}"></div><ul class="find">${b.roles.map(r => `<li><b>${SCHEME[r.scheme]}${d.psg.length > 1 ? ` (${r.group})` : ''}:</b> ${rolesText(b, r)}.</li>`).join('')}</ul>`
            : `<p class="body">${b.groups.length ? '' : `Group${d.psg.length > 1 ? 's' : ''} ${andList(d.psg)} ${d.psg.length > 1 ? 'give' : 'gives'} this task no samples, so ${d.psg.length > 1 ? 'they take' : 'it takes'} no part in its folds.`}</p>`}`;
      }).join('')}
      ${howto(`The bar divides ${d.id}'s records for each task: samples of the positive and negative class first, then the records left out, by the first reason that applies. The squares are the folds of each split scheme: filled green where the group is tested, outlined where it is held out but not scored, grey where it trains. A model is scored within each held-out group, never over held-out groups pooled.`)}
      <p class="note">The corrected tables are not inputs for the sepsis tasks: their correction used the healthy controls of every study, so a model tested on a held-out study would have seen its controls through the correction.</p>`;
  }
  function benchBar(el, d, b) {
    const [pl, nl] = TASK_SIDES[b.task], greys = Object.keys(b.reasons);
    const parts = [{ key: 'pos', label: `Sample of the task: ${pl}`, n: b.pos, color: b.task === 'mortality_in_sepsis' ? mixHex(C.cls.sepsis, '#000000', 0.45) : C.cls.sepsis },
      { key: 'neg', label: `Sample of the task: ${nl}`, n: b.neg, color: b.task === 'mortality_in_sepsis' ? mixC(C.cls.sepsis, 0.45) : b.task === 'sepsis_vs_healthy_control' ? C.cls.healthy_control : C.cls.non_septic_control }]
      .concat(greys.map((k, i) => ({ key: k, label: `Left out: ${BENCH_SHORT[k] || k}`, n: b.reasons[k], color: mixHex(C.lineStrong, C.bg, (i % 3) * 0.3) })));
    stackBars(el, [{ label: d.id, sub: plural(b.records, 'record'), parts }], `Stacked bar of the records of ${d.id} in the task ${TASK[b.task]}.`);
  }
  function foldStrip(el, d, b) {
    const R = P.resource, folds = R.folds[b.task] || {}, rows = b.roles;
    const W = Math.max(300, el.clientWidth || 700), small = W < 560, lw = small ? 0 : 170, sq = 13, gap = 3, rowH = small ? 40 : 26, H = 6 + rows.length * rowH;
    const s = svg(el, W, H, `The folds of each split scheme of ${TASK[b.task]}, with the role of the group of ${d.id}.`);
    rows.forEach((r, i) => {
      const y = 6 + i * rowH + (small ? 14 : 0), total = folds[r.scheme] || (r.train + r.test + r.held_out_not_scored);
      T(s, small ? 0 : lw, small ? y - 4 : y + 10, `${SCHEME[r.scheme]}${d.psg.length > 1 ? ` (${r.group})` : ''} · ${plural(total, 'fold')}`, { 'text-anchor': small ? 'start' : 'end', 'font-size': 12, fill: C.ink2 });
      const per = Math.max(1, Math.floor((W - lw - 14) / (sq + gap)));
      const cells = [...Array(r.test).fill('test'), ...Array(r.held_out_not_scored).fill('held_out_not_scored'), ...Array(Math.max(0, total - r.test - r.held_out_not_scored)).fill('train')];
      cells.slice(0, per).forEach((c, j) => {
        const x = lw + (small ? 0 : 10) + j * (sq + gap);
        hover(rrect(s, x, y, sq, sq, c === 'test' ? C.hilite : c === 'train' ? C.lineStrong : 'none', 2, c === 'held_out_not_scored' ? { stroke: C.hilite, 'stroke-width': 1.5 } : {}),
          `<b>${SCHEME[r.scheme]}</b><br>${c === 'train' ? 'a fold in which the group trains' : c === 'test' ? `the fold in which the group is tested${r.test_folds[0] ? ` (${esc(r.test_folds[0])})` : ''}` : 'the fold in which the group is held out but not scored'}`);
      });
    });
  }
  function benchDraw(d) { d.bench.forEach((b, i) => { benchBar($(`stBe${i}`), d, b); if ($(`stBeF${i}`)) foldStrip($(`stBeF${i}`), d, b); }); }

  // ======================================================================= the chapters
  const RENDER = {
    overview: [overviewHTML, overviewDraw], samples: [samplesHTML, samplesDraw], platform: [platformHTML, platformDraw], genes: [genesHTML, genesDraw],
    checks: [checksHTML, checksDraw], correction: [correctionHTML, correctionDraw], benchmark: [benchHTML, benchDraw],
  };
  function switcher() {
    const i = P.ids.indexOf(study);
    return `<div class="st-switch"><span>Study</span><select class="st-pick" aria-label="Choose another study">${P.ids.map(id => `<option value="${id}"${id === study ? ' selected' : ''}>${id}</option>`).join('')}</select>
      <button type="button" class="btn ghost" data-step="-1"${i <= 0 ? ' disabled' : ''} aria-label="Previous study">← ${i > 0 ? P.ids[i - 1] : ''}</button><button type="button" class="btn ghost" data-step="1"${i >= P.ids.length - 1 ? ' disabled' : ''} aria-label="Next study">${i < P.ids.length - 1 ? P.ids[i + 1] : ''} →</button>
      <a href="#studies">All studies</a></div>`;
  }
  function renderChapter(key) {
    const sec = $(`st-${key}`), body = sec.querySelector(".st-body"), d = cur();
    sec.querySelector('.st-ctx').textContent = d ? `${d.id} · ${famText((BY[d.id] || { fam: [0] }).fam[0])} · ${yr(d.date)}` : '';
    if (!d) { body.innerHTML = '<p class="st-empty">Choose a study first.</p>'; return; }
    const [html, draw] = RENDER[key];
    body.innerHTML = switcher() + html(d);
    draw(d);
  }

  registerView('studies', {
    part: 'studies',
    async init() {
      P = part('studies'); LIST = P.list; BY = Object.fromEntries(LIST.map(r => [r.id, r]));
      study = P.ids.includes(FIRST) ? FIRST : P.ids[0];
      await need(`studies/${study}`);
      initChooser();
    },
    // #studies -> the list; #studies/<id> or #studies/<id>/<chapter> -> a chapter of that study (section st-<chapter>);
    // #studies/<candidate> -> the list, with that excluded candidate shown
    async resolve(args) {
      find = null;
      const a = args[0];
      if (!a || a === 'choose' || a === 'st-choose') return 'st-choose';
      if (P.ids.includes(a)) { study = a; await need(`studies/${a}`); return `st-${CHAPTERS.includes(args[1]) ? args[1] : 'overview'}`; }
      find = a; return 'st-choose';
    },
    hashFor(sec) { return sec.dataset.group === 'profile' && study ? `studies/${study}/${sec.id.slice(3)}` : 'studies'; },
    groupLabel(g) { return g === 'profile' ? (study || 'No study chosen') : g === 'choose' ? '' : undefined; },
    routes() { return ['studies', ...['GSE65682', 'GSE66099', 'GSE60424', 'GSE296830', 'E-MTAB-17576'].filter(id => P.ids.includes(id)).flatMap(id => CHAPTERS.map(c => `studies/${id}/${c}`))]; },
    show(route) {
      refreshBook('studies');
      if (route.chapter === 'st-choose') renderChooser(); else renderChapter(route.chapter.slice(3));
    },
  });

  // ======================================================================= checks (?selftest=1)
  const studyFile = id => part(`studies/${id}`);
  selftest('studies', 'the list holds every included study once, and with the excluded candidates every candidate', () => {
    const ids = new Set(P.list.map(r => r.id)), x = new Set(P.excluded.map(r => r.id));
    const ok = ids.size === N.studies && P.ids.length === N.studies && P.ids.every(id => ids.has(id)) && P.candidates === N.candidates && ids.size + x.size === N.candidates && [...x].every(id => !ids.has(id));
    return { ok, detail: `${fmt(ids.size)} included and ${fmt(x.size)} excluded of ${fmt(N.candidates)} candidates` };
  });
  selftest('studies', 'in every study the class counts add up to its records and to its distinct samples', () => {
    const bad = P.ids.filter(id => { const r = BY[id], f = studyFile(id); return sum(r.cls) !== r.records || sum(r.clsd) !== r.distinct || sum(f.cls) !== f.records || sum(f.clsd) !== f.distinct || f.records !== r.records || f.distinct !== r.distinct; });
    return { ok: !bad.length, detail: bad.slice(0, 5).join(', ') };
  });
  selftest('studies', 'the records, distinct samples and classes of the studies add up to the totals of the resource', () => {
    const rec = sum(LIST.map(r => r.records)), dis = sum(LIST.map(r => r.distinct)), cls = D.order.class.map((k, i) => sum(LIST.map(r => r.clsd[i])));
    const gone = sum(LIST.map(r => r.without));
    const ok = rec === N.records && dis === N.distinct && D.order.class.every((k, i) => cls[i] === N.classes[k]) && gone === N.records_without_values
      && rec + gone === N.included_records;
    return { ok, detail: `records ${fmt(rec)} against ${fmt(N.records)}; distinct ${fmt(dis)} against ${fmt(N.distinct)}; without values ${fmt(gone)} against ${fmt(N.records_without_values)}` };
  });
  selftest('studies', 'the study platforms of every study add up to its records with values, and to its records without values', () => {
    const bad = P.ids.filter(id => { const f = studyFile(id), r = BY[id]; return sum(f.platforms.map(p => p.records)) !== f.records || sum(f.platforms.map(p => p.without)) !== f.without || f.without !== r.without; });
    return { ok: !bad.length, detail: bad.slice(0, 5).join(', ') };
  });
  selftest('studies', 'the AUROC recomputed in the page from every study\'s scores equals the released one', () => {
    let n = 0, worst = 0, where = '';
    P.ids.forEach(id => studyFile(id).checks.scores.forEach(t => {
      n++; const a = auroc(t.pos, t.neg), e = Math.abs(a - t.auroc);
      if (t.pos.length !== t.sepsis || t.neg.length !== t.healthy) { worst = Infinity; where = `${id} ${t.tissue}: ${t.pos.length}/${t.neg.length} scores`; }
      else if (e > worst) { worst = e; where = `${id} ${t.tissue}`; }
    }));
    return { ok: n > 0 && worst <= 1e-3, detail: `${fmt(n)} strata; largest difference ${worst === Infinity ? 'count mismatch' : worst.toExponential(1)} (${where})` };
  });
  selftest('studies', 'the area under the ROC curve drawn in the page equals the AUROC (GSE65682)', () => {
    const t = (studyFile('GSE65682') || studyFile(P.ids[0])).checks.scores[0]; if (!t) return { ok: false, detail: 'no stratum' };
    const pts = rocPoints(t.pos, t.neg); let area = 0; for (let i = 1; i < pts.length; i++) area += (pts[i][0] - pts[i - 1][0]) * (pts[i][1] + pts[i - 1][1]) / 2;
    return { ok: Math.abs(area - auroc(t.pos, t.neg)) < 1e-9, detail: `trapezoid ${fx(area, 6)}, rank ${fx(auroc(t.pos, t.neg), 6)}` };
  });
  selftest('studies', 'in every sex check the concordant and discordant samples match the samples drawn', () => {
    const bad = [];
    P.ids.forEach(id => studyFile(id).checks.sex.forEach(t => {
      const pts = t.points || [], c = pts.filter(p => p.a === 'concordant').length, x = pts.filter(p => p.a === 'discordant').length;
      if (c !== t.concordant || x !== t.discordant || (t.status.startsWith('tested') && t.female + t.male + t.ambiguous !== pts.filter(p => ['female', 'male', 'ambiguous'].includes(p.e)).length)) bad.push(`${id} ${t.p}`);
    }));
    return { ok: !bad.length, detail: bad.slice(0, 5).join(', ') };
  });
  selftest('studies', 'in every task the samples and the records left out add up to the study\'s records', () => {
    const bad = [];
    P.ids.forEach(id => { const f = studyFile(id); f.bench.forEach(b => { if (b.pos + b.neg !== b.eligible || b.eligible + sum(Object.values(b.reasons)) !== b.records || b.records !== f.records) bad.push(`${id} ${b.task}`); }); });
    return { ok: !bad.length, detail: bad.slice(0, 5).join(', ') };
  });
  selftest('studies', 'the outcomes of every study in the correction add up to its records, and its samples in the corrected tables match the points drawn', () => {
    const bad = [];
    P.ids.forEach(id => {
      const f = studyFile(id), led = sum(f.correction.ledger.map(l => l.n)), inLed = sum(f.correction.ledger.filter(l => l.outcome === 'in_corrected_table').map(l => l.n));
      const pts = g => { const V = viewOf(g), si = V.studies.indexOf(id); return si < 0 ? 0 : V.study.filter(s => s === si).length; };
      const inT = sum(Object.values(f.correction.in_tables));
      if (led !== f.records || inLed !== inT || Object.entries(f.correction.in_tables).some(([g, n]) => pts(g) !== n)) bad.push(id);
    });
    return { ok: !bad.length, detail: bad.slice(0, 5).join(', ') };
  });
  selftest('studies', 'every chapter of every study draws without NaN, undefined or a script error', async () => {
    const before = location.hash, e0 = ERRORS.length, bad = /\bNaN\b|\bundefined\b|\bInfinity\b|\[object /, problems = [];
    for (const id of P.ids) for (const c of CHAPTERS) {
      await goTo(`studies/${id}/${c}`, false);
      const txt = $(`st-${c}`).textContent;
      if (bad.test(txt)) problems.push(`${id}/${c}: "${txt.match(new RegExp(`.{0,40}(${bad.source}).{0,30}`))[0].trim()}"`);
    }
    await goTo(before || 'studies', false);
    const errs = ERRORS.slice(e0);
    return { ok: !problems.length && !errs.length, detail: `${fmt(P.ids.length * CHAPTERS.length)} chapters; ${[...problems.slice(0, 3), ...errs.slice(0, 2)].join('; ')}` };
  });
  selftest('studies', 'the search of the list finds a study by accession, title word and platform', () => {
    const q = view.q, r = BY[P.ids[0]], words = [r.id.toLowerCase(), r.title.toLowerCase().split(/\s+/).find(w => w.length > 4) || r.id.toLowerCase(), (r.plats[0] || r.id).toLowerCase()];
    const ok = words.every(w => LIST.filter(x => matches(searchText(x), w)).some(x => x.id === r.id));
    view.q = q;
    return { ok, detail: words.join(', ') };
  });
})();
