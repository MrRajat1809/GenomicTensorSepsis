/* Method in detail: one chapter per step of the pipeline, each with its Equation, an interactive Simulation, a Worked
   example on real samples and the Result of the step for the whole resource. Every number comes from data/methods.js
   (methods() in build_data.py) or from the core data; the constants of the invented examples are marked as invented in
   the text. Every worked example is recomputed here from its exported inputs and checked against its exported result
   with ?selftest=1. Uses common.js. */
'use strict';
(() => {
  let M = null;                       // the part, once loaded
  const ST = {};                      // the state of the simulations
  // ---------- small builders ----------
  const V = s => `<span class="v">${s}</span>`;
  const NB = s => `<span class="n">${s}</span>`;
  const frac = (a, b) => `<span class="frac"><span>${a}</span><span>${b}</span></span>`;
  const how = t => `<p class="howto"><span>How to read it</span><span>${t}</span></p>`;
  const H3 = (kind, label, title) => `<h3><span class="tag t-${kind}">${label}</span>${title ? ` ${title}` : ''}</h3>`;
  const EQ = (rows, note) => `<div class="eq">${rows.map(r => `<span class="row">${r}</span>`).join('')}${note ? `<span class="c">${note}</span>` : ''}</div>`;
  const RES = html => `<div class="result">${html}</div>`;
  const code = s => `<code>${esc(s)}</code>`;
  const sub = s => `<sub>${s}</sub>`;
  const neg = v => v < 0 ? `${MINUS}${String(v).replace(/^-/, '')}` : v;
  // the drawing width of a chart: its box on screen, so that text keeps its size on a phone
  const wOf = id => clamp(($(id) && $(id).clientWidth) || 860, 300, 980);
  const small = w => w < 560;
  const close = (a, b, tol) => Math.abs(a - b) <= tol;
  const ok = (cond, detail) => ({ ok: !!cond, detail: detail || '' });
  const dayOf = iso => fdate(String(iso).slice(0, 10));
  const TIME = iso => String(iso).slice(11, 19);
  const clsDot = k => `<i class="mdot" style="background:${C.cls[k]}"></i>`;
  const studyLink = s => `<a href="#studies/${esc(s)}">${esc(s)}</a>`;
  // a chart box: the SVG and its aria label
  function chart(id, w, h, label) { const s = svg(id, w, h, label); return s; }
  // number inputs of the invented examples
  const numIn = (id, v, attrs = '') => `<input class="mnum" type="number" id="${id}" value="${v}" step="any" ${attrs} aria-label="${id}">`;

  // ======================================================================= 0. overview
  function buildOver() {
    const f = M.overview.flow, N = D.numbers, s = M.search, T7 = M.correction.tables, b = M.bench;
    const open = s.by_source.filter(r => ['geo', 'arrayexpress', 'ena'].includes(r.source)), cat = s.by_source.filter(r => !open.includes(r));
    const SRC = { geo: 'GEO', arrayexpress: 'ArrayExpress', ena: 'ENA', ega: 'EGA', dbgap: 'dbGaP', node: 'NODE' };
    const box = (href, n, title, line) => `<div class="fbox core"><a href="#methods/${href}"><b>${n}</b> ${title}<span class="fsub">${line}</span></a></div>`;
    const folds = sum(b.tasks.map(t => sum(Object.values(t.folds))));
    $('mbOver').innerHTML = `
      ${H3('sim', 'Diagram', 'From the repositories to the released tables')}
      <div class="flowgrid mflow">
        <div class="flowcol"><h4>Public repositories</h4>
          ${open.map(r => `<div class="fbox"><b>${SRC[r.source]}</b><span class="fsub">${fmt(r.records)} records found, open data</span></div>`).join('')}
          <div class="fbox"><b>${andList(cat.map(r => SRC[r.source]))}</b><span class="fsub">${fmt(sum(cat.map(r => r.records)))} records, catalogued only: controlled or unstated access</span></div>
        </div>
        <div class="flowarrow" aria-hidden="true">→</div>
        <div class="flowcol"><h4>Steps, one chapter each</h4>
          ${box('m-search', 1, 'Search and screening', `${fmt(f.records_found)} records → ${fmt(f.candidate_studies)} candidates → ${fmt(f.included_studies)} included studies`)}
          ${box('m-arrays', 2, 'Arrays to numbers', `${fmt(M.arrays.affymetrix.arrays)} Affymetrix arrays read from their scanner files, the other arrays from intensities`)}
          ${box('m-rnaseq', 3, 'RNA-seq to numbers', `${fmt(M.counts.samples)} libraries of ${fmt(M.counts.study_platforms)} study platforms as log₂ counts per million`)}
          ${box('m-genes', 4, 'Probes to genes', `${fmt(M.genes.features.total)} rows of ${fmt(M.genes.tables.count)} matrices related to HGNC genes`)}
          ${box('m-labels', 5, 'Labels from records', `${fmt(f.sample_records_with_values)} records with expression values labeled by ${fmt(M.labels.rule_sets)} sets of rules`)}
          ${box('m-repeats', 6, 'Repeat deposits and patients', `${fmt(f.sample_records_with_values)} records → ${fmt(f.distinct_samples)} distinct samples → ${fmt(f.one_sample_per_patient)} one per patient`)}
          ${box('m-sex', 7, 'Sex from expression', `${fmt(M.sex.resource.samples.compared)} samples compared with their recorded sex`)}
          ${box('m-auroc', 8, 'The sepsis score and the AUROC', `${fmt(M.score.resource.strata)} study strata scored`)}
          ${box('m-correct', 9, 'The correction', `${fmt(T7.bulk_blood.samples)} bulk-blood and ${fmt(T7.pbmc.samples)} PBMC samples on one scale`)}
          ${box('m-heldout', 10, 'The held-out test', `${fmt(M.heldout.medians.bulk_blood.units + M.heldout.medians.pbmc.units)} units corrected with their healthy controls hidden`)}
          ${box('m-bench', 11, 'Benchmark splits and leakage', `${spell(b.tasks.length)} tasks, ${fmt(folds)} folds of patient-safe groups`)}
        </div>
        <div class="flowarrow" aria-hidden="true">→</div>
        <div class="flowcol"><h4>Released tables</h4>
          <div class="fbox"><b>Per-study tables</b><span class="fsub">${fmt(N.tables)} matrices as processed, each with its gene-level table</span></div>
          <div class="fbox"><b>The sample table</b><span class="fsub">${fmt(N.records)} records: curated labels with their sources, repeat deposits, patients and flags</span></div>
          <div class="fbox"><b>The checks</b><span class="fsub">sex read from the genes, the pairs of repeat deposits, the sepsis score of every stratum</span></div>
          <div class="fbox"><b>Corrected tables</b><span class="fsub">bulk blood, ${fmt(T7.bulk_blood.studies)} studies; PBMC, ${fmt(T7.pbmc.studies)} studies</span></div>
          <div class="fbox"><b>Benchmark</b><span class="fsub">every record of every task, and the folds</span></div>
        </div>
      </div>
      ${how('Read from left to right. Each box in the middle opens its chapter; the numbers are those of the whole resource. The per-study tables are never changed after processing: labels, links and checks are added beside them, and the corrected values go into tables of their own.')}`;
  }

  // ======================================================================= 1. search and screening
  const SRCN = { geo: 'GEO', arrayexpress: 'ArrayExpress', ena: 'ENA', ega: 'EGA', dbgap: 'dbGaP', node: 'NODE' };
  function buildSearch() {
    const s = M.search, sci = sum(s.reasons.filter(r => r.group === 'scientific').map(r => r.studies)), av = sum(s.reasons.filter(r => r.group === 'availability').map(r => r.studies));
    ST.terms = s.terms.map(() => true);
    ST.ladder = 'GSE103119'; ST.invented = s.reasons.map(r => r.key === 'single_cell' || r.key === 'sorted_cells');
    const tot = k => sum(s.by_source.map(r => r[k]));
    $('mbSearch').innerHTML = `
      ${H3('eq', 'Equation')}
      ${EQ([`${V('C')} = Σ${sub('r')} ${V('E')}${sub('r')} + ${V('I')} &nbsp;&nbsp; ${NB(fmt(s.candidates_total))} = ${NB(fmt(sci))} + ${NB(fmt(av))} + ${NB(fmt(s.included))}`,
        `recall = ${frac('known studies found', 'known studies')} = ${frac(NB(fmt(s.recall.found)), NB(fmt(s.recall.checked)))} = ${NB(fx(s.recall.found / s.recall.checked, 2))}`],
        `C: candidate studies, each study counted once across the repositories. E<sub>r</sub>: studies excluded for reason r, the first reason in a fixed order that applies to the study, so that each study counts once; the ${fmt(sci)} are excluded for a scientific reason and the ${fmt(av)} for availability, which comes last. I: included studies. The known studies are the cohorts of the multi-cohort analyses of sepsis transcriptomes by Sweeney and colleagues whose own records name sepsis; ${spell(s.recall.outside_scope)} further cohorts of those analyses never name sepsis and are outside the rules.`)}
      ${H3('sim', 'Simulation a', 'Switch search terms off')}
      <p class="note mtop0">Every term was searched separately in every repository, on ${dayOf(s.searched)}. Click a term to leave it out of the search and see which candidates and included studies it alone reached.</p>
      <div class="controls"><button type="button" class="btn ghost" id="msAll">All ${fmt(s.terms.length)} terms</button><button type="button" class="btn ghost" id="msTwo">Only sepsis and septic</button></div>
      <div class="chart" id="msTerms"></div>
      <div class="readout" id="msTermsOut" aria-live="polite"></div>
      ${how('One row per search term. The grey bar counts the candidate studies the term reached, the green part the included studies among them; a study reached by several terms appears in each of their rows. A study counts as reached while at least one of the terms switched on reached it. Hover over a row for the hits of the term in each repository.')}
      <h3>Counting each study once</h3>
      <div class="tablewrap"><table class="t" id="msDup"><thead><tr><th>Repository</th><th class="r">Records found</th><th class="r">Set aside</th><th class="r">Copies of a study already found</th><th class="r">Candidate studies</th></tr></thead><tbody>
        ${s.by_source.map(r => `<tr><td>${SRCN[r.source]}</td><td class="r">${fmt(r.records)}</td><td class="r">${fmt(r.no_runs + r.not_expression)}</td><td class="r">${fmt(r.duplicates)}</td><td class="r">${fmt(r.candidates)}</td></tr>`).join('')}
        <tr class="total"><td>All</td><td class="r">${fmt(tot('records'))}</td><td class="r">${fmt(tot('no_runs') + tot('not_expression'))}</td><td class="r">${fmt(tot('duplicates'))}</td><td class="r">${fmt(tot('candidates'))}</td></tr></tbody></table></div>
      <p class="note">The repositories are taken in a fixed order (${andList(s.sources.map(k => SRCN[k]))}), and a record is a copy when a repository earlier in the order holds the same study by a shared accession or a declared link. ${s.duplicates.map(d => `${fmt(d.records)} ${SRCN[d.source]} ${d.records === 1 ? 'record is a copy' : 'records are copies'} of ${SRCN[d.of]} ${d.of === 'geo' ? 'series' : 'studies'}`).join('; ')}. Set aside at once: ${fmt(tot('no_runs'))} ENA projects without a human RNA-seq run, and ${fmt(tot('not_expression'))} ArrayExpress studies that are not expression studies by their recorded type.</p>
      ${H3('sim', 'Simulation b', 'The first reason that applies')}
      <p class="note mtop0">Each candidate is taken through the reasons in order, and the first that applies is its reason. Choose a real study, or describe an invented one by ticking what is true of it.</p>
      <div class="controls"><div class="seg" id="msWhich">${segHTML([['GSE103119', 'GSE103119 (excluded)'], ['GSE65682', 'GSE65682 (included)'], ['invented', 'An invented study']], ST.ladder)}</div></div>
      <div class="mladder" id="msLadder"></div>
      ${H3('ex', 'Worked example', 'GSE103119 and GSE65682')}
      <div class="grid2"><div class="readout" id="msEx1"></div><div class="readout" id="msEx2"></div></div>
      ${H3('res', 'Result', 'From candidates to included studies')}
      <div class="chart" id="msFunnel"></div>
      ${how('One row per reason, in the order in which they are checked. The grey bar is the number of candidates still in after the reason; the colored end is the number the reason removed (dark: a scientific reason, orange: availability). The last row is the included studies. Hover over a row for the reason in full.')}
      <div id="msRes"></div>`;
    $('msAll').onclick = () => { ST.terms = s.terms.map(() => true); drawSearch(); };
    $('msTwo').onclick = () => { ST.terms = s.terms.map(t => t === 'sepsis' || t === 'septic'); drawSearch(); };
    seg('msWhich', v => { ST.ladder = v; drawLadder(); });
    $('msLadder').addEventListener('change', e => { const b = e.target.closest('input[data-r]'); if (b) { ST.invented[+b.dataset.r] = b.checked; drawLadder(); } });
    const ex = s.excluded, inc = s.included_example;
    $('msEx1').innerHTML = `<h4 class="exhead">GSE103119 <span class="m">excluded</span></h4>A GEO series released in ${ex.year}, with ${fmt(ex.records)} sample records on ${esc(ex.platform)} (${esc(ex.platform_title)}). The search reached it by ${spell(ex.terms.length)} ${ex.terms.length === 1 ? 'term' : 'terms'} only, <i>${andList(ex.terms.map(esc))}</i>; none of the terms sepsis, septic or septicemia did. Its records give a condition for every sample: ${andList(ex.conditions.map(([k, v]) => `${fmt(v)} “${esc(k)}”`))}, all whole blood. It passes every reason up to the one that asks for a blood sample from a patient the study diagnoses with sepsis, and stops there. The recorded basis: “${esc(ex.basis)}”.`;
    $('msEx2').innerHTML = `<h4 class="exhead">GSE65682 <span class="m">included</span></h4>A series of the MARS consortium released in ${inc.year}, with ${fmt(inc.records)} sample records on ${esc(inc.platform_title)}. It matched ${spell(inc.terms.length)} of the ${spell(s.terms.length)} terms: ${andList(inc.terms.map(t => `<i>${esc(t)}</i>`))}. ArrayExpress holds a copy, ${esc(inc.copy)}, marked as a duplicate, so the study counts once, under GEO. ${inc.recall ? 'It is one of the known cohorts of the recall check, and the search found it.' : ''} It passes every reason; the recorded basis: “${esc(inc.basis)}”. It becomes one study platform, ${esc(inc.study_platform)}, whose records curation later sorts into ${fmt(inc.classes.sepsis)} sepsis samples, ${fmt(inc.classes.healthy_control)} healthy controls, ${fmt(inc.classes.non_septic_control)} non-septic controls and ${fmt(inc.classes.other_condition)} other conditions; ${fmt(inc.outcome_recorded)} give an outcome. Inclusion is a decision about the study: its records of other conditions are kept with the rest.`;
    const byS = {}; s.candidates.forEach(c => { if (c.reason === 0) byS[c.source] = (byS[c.source] || 0) + 1; });
    $('msRes').innerHTML = RES(`<b>${fmt(s.included)}</b> of the ${fmt(s.candidates_total)} candidate studies were included: ${andList(Object.entries(byS).map(([k, v]) => `${fmt(v)} from ${SRCN[s.sources[k]]}`))}. Together they hold <b>${fmt(s.included_records)}</b> sample records. For every candidate, excluded or not, the study record and its sample records were stored (${fmt(s.harvested_records)} sample records), so that each decision rests on recorded evidence. The search found all <b>${fmt(s.recall.found)}</b> of the ${fmt(s.recall.checked)} known cohorts.`);
  }
  // candidates reached by the terms switched on: a candidate without a recorded term (an ENA match inside a word) is always reached
  const termMask = () => ST.terms.reduce((m, on, i) => on ? m | (1 << i) : m, 0);
  const reached = (c, mask) => c.terms === 0 || (c.terms & mask) !== 0;
  function drawSearch() {
    const s = M.search, mask = termMask(), w = wOf('msTerms'), sm = small(w);
    const rows = s.terms.map((t, i) => { const own = s.candidates.filter(c => c.terms & (1 << i)); return { t, i, all: own.length, inc: own.filter(c => c.reason === 0).length }; });
    const rowH = 26, top = 20, left = sm ? 150 : 210, right = sm ? 58 : 130, H = top + rowH * rows.length + 30;
    const g = chart('msTerms', w, H, 'Bar chart: candidate and included studies reached by each search term; click a row to switch the term on or off.');
    const xmax = niceMax(Math.max(...rows.map(r => r.all))), x = lin(0, xmax, left, w - right);
    gridX(g, x, ticks(0, xmax, sm ? 3 : 5), top - 6, top + rowH * rows.length);
    axisTitle(g, (left + w - right) / 2, H - 2, 'Candidate studies reached by the term');
    rows.forEach((r, k) => {
      const y = top + k * rowH, on = ST.terms[r.i];
      const grp = S('g', { class: 'mrow', role: 'button', tabindex: 0, 'aria-pressed': on, 'aria-label': `${r.t}: ${on ? 'on' : 'off'}` }, g);
      S('rect', { x: 0, y, width: w, height: rowH, fill: 'transparent' }, grp);
      S('rect', { x: 4, y: y + 7, width: 12, height: 12, rx: 3, fill: on ? C.ink : 'none', stroke: C.lineStrong }, grp);
      if (on) S('path', { d: `M${7},${y + 13} l3,3 l6,-7`, stroke: C.bg, 'stroke-width': 1.8, fill: 'none' }, grp);
      T(grp, 22, y + 17, sm && r.t.length > 19 ? `${r.t.slice(0, 18)}…` : r.t, { 'font-size': sm ? 11.5 : 12.5, fill: on ? C.ink : C.muted, 'text-decoration': on ? null : 'line-through' });
      rrect(grp, x(0), y + 5, x(r.all) - x(0), rowH - 10, on ? C.lineStrong : C.grid, 3);
      rrect(grp, x(0), y + 5, x(r.inc) - x(0), rowH - 10, on ? C.accent : tint(C.accent, 0.6), 3);
      T(grp, x(r.all) + 5, y + 17, sm ? `${fmt(r.all)} · ${fmt(r.inc)}` : `${fmt(r.all)} · ${fmt(r.inc)} included`, { 'font-size': 11.5, fill: C.muted });
      hover(grp, () => `<b>${esc(r.t)}</b><br>${s.sources.map((k, j) => `${SRCN[k]}: ${fmt(s.hits[j][r.i])} hits`).join('<br>')}<br><span class="m">${fmt(r.all)} candidate studies, ${fmt(r.inc)} included</span>`);
      const toggle = () => { ST.terms[r.i] = !ST.terms[r.i]; drawSearch(); };
      grp.addEventListener('click', toggle);
      grp.addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); toggle(); } });
    });
    const cand = s.candidates.filter(c => reached(c, mask)), inc = cand.filter(c => c.reason === 0);
    const lost = s.candidates.filter(c => c.reason === 0 && !reached(c, mask)).map(c => c.accession);
    const n = ST.terms.filter(Boolean).length;
    $('msTermsOut').innerHTML = `With ${plural(n, 'term')} the searches reach <b>${fmt(cand.length)}</b> of the ${fmt(s.candidates_total)} candidate studies and <b>${fmt(inc.length)}</b> of the ${fmt(s.included)} included studies.` +
      (lost.length ? ` Lost: ${lost.slice(0, 12).map(studyLink).join(', ')}${lost.length > 12 ? ` and ${fmt(lost.length - 12)} more` : ''}.` : ' No included study is lost.');
  }
  function drawLadder() {
    const s = M.search;
    let stop = -1, header;
    if (ST.ladder === 'invented') { stop = ST.invented.findIndex(Boolean); header = 'An invented study: tick what is true of it'; }
    else { const c = s.candidates.find(x => x.accession === ST.ladder); stop = c.reason - 1; header = ST.ladder; }
    const inv = ST.ladder === 'invented';
    $('msLadder').innerHTML = `<p class="mlad-h">${esc(header)} <span class="small">✓ does not apply · ✗ the first reason that applies</span></p><ol>${s.reasons.map((r, i) => {
      const state = stop < 0 || i < stop ? 'pass' : i === stop ? 'stop' : 'after';
      const mark = state === 'pass' ? '✓' : state === 'stop' ? '✗' : '·';
      return `<li class="${state} ${r.group}"><span class="mk" aria-hidden="true">${mark}</span>${inv ? `<label><input type="checkbox" data-r="${i}" ${ST.invented[i] ? 'checked' : ''}> ${esc(r.label)}</label>` : `<span>${esc(r.label)}</span>`}${state === 'stop' ? ' <b>← its reason</b>' : ''}</li>`;
    }).join('')}<li class="${stop < 0 ? 'stop incl' : 'after'}"><span class="mk" aria-hidden="true">${stop < 0 ? '✓' : '·'}</span><span>Included</span></li></ol>
    <p class="note">${stop < 0 ? 'No reason applies, so the study is included.' : `The study is counted under “${esc(s.reasons[stop].label)}” only, even if a later reason would also apply.`}${inv && ST.invented.filter(Boolean).length > 1 ? ' A study of single sorted cells, for example, counts as single-cell, not as sorted cells.' : ''}</p>`;
  }
  function drawFunnel() {
    const s = M.search, w = wOf('msFunnel'), sm = small(w), steps = s.funnel.slice(1), byKey = Object.fromEntries(s.reasons.map(r => [r.key, r]));
    const rowH = 23, top = 24, left = sm ? 160 : 300, right = sm ? 40 : 60, H = top + rowH * (steps.length + 2) + 28;
    const g = chart('msFunnel', w, H, `Funnel: ${fmt(s.candidates_total)} candidate studies reduced reason by reason to ${fmt(s.included)} included studies.`);
    const x = lin(0, s.candidates_total, left, w - right);
    gridX(g, x, ticks(0, s.candidates_total, sm ? 3 : 6), top - 6, top + rowH * (steps.length + 2) - 4);
    axisTitle(g, (left + w - right) / 2, H - 2, 'Candidate studies still in');
    const row = (k, label, remain, removed, color, html) => {
      const y = top + k * rowH;
      T(g, left - 8, y + 15, label, { 'text-anchor': 'end', 'font-size': 12, fill: C.ink2 });
      rrect(g, x(0), y + 4, x(remain) - x(0), rowH - 8, k === steps.length + 1 ? C.accent : C.lineStrong, 2);
      if (removed) rrect(g, x(remain), y + 4, Math.max(2, x(remain + removed) - x(remain)), rowH - 8, color, 2);
      T(g, x(remain + (removed || 0)) + 5, y + 15, removed ? `${MINUS}${fmt(removed)}` : fmt(remain), { 'font-size': 11.5, fill: C.muted });
      hitRow(g, 0, y, w, rowH, html);
    };
    row(0, 'Candidate studies', s.candidates_total, 0, null, `<b>${fmt(s.candidates_total)} candidate studies</b>`);
    steps.forEach((st, k) => { const r = byKey[st.step]; row(k + 1, sm ? r.label.replace(/^(.{22}).+$/, '$1…') : r.label, st.remaining, st.removed, r.group === 'scientific' ? C.ink2 : C.tier.unanchored,
      `<b>${esc(r.label)}</b><br>${esc(r.description)}<br><span class="m">${fmt(st.removed)} ${st.removed === 1 ? 'study' : 'studies'}${st.records != null ? `, ${fmt(st.records)} sample records` : ''}; ${fmt(st.remaining)} still in</span>`); });
    row(steps.length + 1, 'Included studies', s.included, 0, null, `<b>${fmt(s.included)} included studies</b>`);
  }

  // ======================================================================= 2. arrays to numbers
  // the standard normal: density, and the share below z from the complementary error function (fractional error below 1.2e-7)
  function erfc(x) {
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  }
  const phi = z => Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI), Phi = z => 0.5 * erfc(-z / Math.SQRT2);
  // normexp: the expected signal of an observed intensity x, given the background mean mu, its spread sigma and the mean signal alpha
  function normexp(x, mu, sigma, alpha) { const m = x - mu - sigma * sigma / alpha; return { m, s: m + sigma * phi(m / sigma) / Phi(m / sigma) }; }
  // quantile normalization of the columns of a table (rows are probes): each value takes the mean of the values of its rank;
  // tied values share the mean of the targets of their ranks
  function quantileNormalize(cols) {
    const n = cols[0].length, sorted = cols.map(c => c.slice().sort((a, b) => a - b));
    const target = sorted[0].map((_, k) => mean(sorted.map(s => s[k])));
    return { sorted, target, out: cols.map(c => c.map(v => { const lo = c.filter(u => u < v).length, ties = c.filter(u => u === v).length; return mean(target.slice(lo, lo + ties)); })), n };
  }
  // Tukey's median polish (rows are probes, columns are arrays), as R's medpolish: returns the common level, the probe and array
  // effects, the remainders and the number of rounds; RMA reports common level + array effect for every array
  function medianPolish(x, maxit = 10, eps = 0.01) {
    const nr = x.length, nc = x[0].length, z = x.map(r => r.slice()), r = Array(nr).fill(0), c = Array(nc).fill(0);
    let t = 0, old = 0, rounds = 0;
    for (let it = 0; it < maxit; it++) {
      rounds++;
      for (let i = 0; i < nr; i++) { const d = median(z[i]); for (let j = 0; j < nc; j++) z[i][j] -= d; r[i] += d; }
      let d = median(c); for (let j = 0; j < nc; j++) c[j] -= d; t += d;
      for (let j = 0; j < nc; j++) { const cd = median(z.map(row => row[j])); for (let i = 0; i < nr; i++) z[i][j] -= cd; c[j] += cd; }
      d = median(r); for (let i = 0; i < nr; i++) r[i] -= d; t += d;
      const now = sum(z.flat().map(Math.abs));
      if (now === 0 || Math.abs(now - old) < eps * now) break;
      old = now;
    }
    return { t, r, c, z, rounds, values: c.map(v => t + v) };
  }
  const QN0 = [[1, 4, 3], [6, 2, 5], [3, 8, 10], [12, 16, 20]];   // invented: four probes (rows) on three arrays A, B, C
  const MP0 = [[6, 8, 7], [9, 11, 10], [5, 7, 12]];               // invented: three probes of one probe set on arrays A, B, C
  function buildArrays() {
    const a = M.arrays, r = a.rma, bg = a.background;
    ST.bgx = bg.probes.find(p => p.probe === 'median').x;
    const inTable = (prefix, m) => `<table class="t mgrid"><thead><tr><th>Probe</th><th class="r">A</th><th class="r">B</th><th class="r">C</th></tr></thead><tbody>${m.map((row, i) => `<tr><td>${prefix}${sub(i + 1)}</td>${row.map((v, j) => `<td class="r">${numIn(`${prefix}${i}${j}`, v, 'min="0"')}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
    $('mbArrays').innerHTML = `
      ${H3('eq', 'Equation', 'The three steps of RMA')}
      ${EQ([`1. background: &nbsp;${V('ŝ')} = ${V('m')} + ${V('σ')} ${frac(`φ(${V('m')}/${V('σ')})`, `Φ(${V('m')}/${V('σ')})`)}, &nbsp; ${V('m')} = ${V('x')} − ${V('μ')} − ${frac(`${V('σ')}²`, V('α'))}`,
        `2. quantile normalization: &nbsp;${V('x')}′${sub('pj')} = ${frac('1', V('J'))} Σ${sub('k')} ${V('x')}${sub(`(r${sub('pj')}) k`)}`,
        `3. median polish: &nbsp;log₂ ${V('y')}${sub('pj')} = ${V('μ')}′ + ${V('a')}${sub('p')} + ${V('b')}${sub('j')} + ${V('e')}${sub('pj')}, &nbsp; value of array ${V('j')} = ${V('μ')}′ + ${V('b')}${sub('j')}`],
        'x: the intensity of a probe as scanned. μ and σ: the mean and spread of the background, the light of a spot to which nothing has bound. α: the mean signal. φ is the density of the standard normal distribution and Φ the share of it below a value. In step 2, r<sub>pj</sub> is the rank of probe p in array j, x<sub>(r)k</sub> the value of rank r in array k, and J the number of arrays: every value is replaced by the mean of the values of the same rank. In step 3, y<sub>pj</sub> are the normalized intensities of the probes of one probe set; the common level μ′, the probe effects a<sub>p</sub> and the array effects b<sub>j</sub> are fitted with medians, and each array gets μ′ + b<sub>j</sub>.')}
      ${H3('sim', 'Simulation a', `Background correction with the model of ${bg.array}`)}
      <p class="note mtop0">limma's normexp, used for Illumina and Agilent arrays, rests on the same model as the first step of RMA. The background below is the one estimated for ${bg.array} of ${bg.study} from its detection p-values. Move the intensity.</p>
      <div class="controls"><label>Observed intensity ${V('x')} = <b id="maBgV"></b> <input type="range" id="maBgX" min="60" max="300" step="0.1" value="${ST.bgx}"></label></div>
      <div class="chart" id="maBgChart"></div>
      <div class="eq" id="maBgEq"></div>
      ${how('The curve is the corrected value ŝ for every observed intensity; the dashed line is plain subtraction of the background mean, x − μ. Far above the background the two agree. Near and below it, plain subtraction goes to zero and below, where a value has no logarithm, while ŝ stays small but positive. The marked points are the median probe of the array and an invented dim probe.')}
      ${H3('sim', 'Simulation b', 'Quantile normalization, on three invented arrays')}
      <p class="note mtop0">An invented example of four probes on three arrays A, B and C; change any intensity.</p>
      <div class="grid2"><div>${inTable('p', QN0)}</div><div class="chart" id="maQnChart"></div></div>
      <div class="readout" id="maQnOut"></div>
      ${how('The table holds the intensities; the chart draws each array before (open dots) and after (filled dots) normalization, on the log₂ scale, with lines that join a probe to its new value. After normalization every array holds the same values, and within each array the probes keep their order.')}
      ${H3('sim', 'Simulation c', 'Median polish, on one invented probe set')}
      <p class="note mtop0">An invented probe set of three probes, already normalized and on the log₂ scale. The probe q<sub>3</sub> misbehaves on array C; change it, or any other value.</p>
      <div class="grid2"><div>${inTable('q', MP0)}</div><div class="readout" id="maMpOut"></div></div>
      ${H3('ex', 'Worked example', `a real array of ${r.study}`)}
      <div class="readout" id="maEx"></div>
      <div class="chart" id="maDots"></div>
      ${how(`One dot per array of ${r.study}. Top row: the median of its raw probe intensities, on the log₂ scale; bottom row: the median of its ${fmt(r.probe_sets)} values after RMA, on the same scale. The lines join the two values of one array. ${r.array} is the larger dot. Hover over a dot for its array.`)}
      <h4 class="exhead">The background of ${bg.array} <span class="m">${bg.study}, column ${esc(bg.label)} of the deposited table</span></h4>
      <div class="tablewrap" id="maBgTable"></div>
      ${H3('res', 'Result', 'How each study platform was read')}
      <div class="chart" id="maRoutes"></div>
      <div id="maRes"></div>`;
    $('maBgX').addEventListener('input', e => { ST.bgx = +e.target.value; drawBackground(); });
    ['p', 'q'].forEach(pre => $('mbArrays').querySelector(`#${pre}00`).closest('table').addEventListener('input', () => pre === 'p' ? drawQN() : drawMP()));
  }
  const readGrid = (pre, rows) => Array.from({ length: rows }, (_, i) => [0, 1, 2].map(j => { const v = parseFloat($(`${pre}${i}${j}`).value); return isFinite(v) ? v : 0; }));
  function drawBackground() {
    const bg = M.arrays.background, { mu, sigma, alpha, offset } = bg, w = wOf('maBgChart'), h = small(w) ? 230 : 260, m = { l: 52, r: 16, t: 14, b: 40 };
    $('maBgV').textContent = fx(ST.bgx, 1);
    const g = chart('maBgChart', w, h, 'Curve of the background-corrected signal against the observed intensity, with plain subtraction for comparison.');
    const X = lin(60, 300, m.l, w - m.r), ymax = 300 - mu, Y = lin(-60, ymax, h - m.b, m.t);
    gridY(g, Y, ticks(-60, ymax, 5), m.l, w - m.r); gridX(g, X, ticks(60, 300, small(w) ? 4 : 7), h - m.b, h - m.b);
    S('line', { x1: m.l, x2: w - m.r, y1: Y(0), y2: Y(0), stroke: C.lineStrong }, g);
    axisTitle(g, (m.l + w - m.r) / 2, h - 6, 'Observed intensity x'); axisTitle(g, 14, (m.t + h - m.b) / 2, 'Corrected value', true);
    const xs = []; for (let x = 60; x <= 300; x += 1) xs.push(x);
    S('path', { d: xs.map((x, i) => `${i ? 'L' : 'M'}${X(x)},${Y(x - mu)}`).join(''), stroke: C.muted, 'stroke-dasharray': '5 4', fill: 'none', 'stroke-width': 1.5 }, g);
    S('path', { d: xs.map((x, i) => `${i ? 'L' : 'M'}${X(x)},${Y(normexp(x, mu, sigma, alpha).s)}`).join(''), stroke: C.accentInk, fill: 'none', 'stroke-width': 2.4 }, g);
    S('line', { x1: X(mu), x2: X(mu), y1: m.t, y2: h - m.b, stroke: C.lineStrong, 'stroke-dasharray': '2 3' }, g);
    T(g, X(mu) + 4, m.t + 10, `background mean μ = ${fx(mu, 1)}`, { 'font-size': 11, fill: C.muted });
    T(g, X(290), Y(290 - mu) + 16, 'x − μ', { 'font-size': 11, fill: C.muted, 'text-anchor': 'end' });
    T(g, X(240), Y(normexp(240, mu, sigma, alpha).s) - 10, 'ŝ (normexp)', { 'font-size': 11.5, fill: C.accentInk, 'font-weight': 600, 'text-anchor': 'end' });
    bg.probes.filter(p => p.x <= 300).forEach(p => {
      const e = normexp(p.x, mu, sigma, alpha);
      hover(S('circle', { cx: X(p.x), cy: Y(e.s), r: 4.5, fill: C.bg, stroke: C.ink, 'stroke-width': 1.5 }, g), `<b>${p.probe === 'median' ? 'Median probe of the array' : 'Invented dim probe'}</b><br>x = ${fx(p.x, 1)}, ŝ = ${fx(e.s, 1)}`);
      T(g, X(p.x) - 6, Y(e.s) - 10, p.probe === 'median' ? 'median probe' : 'invented dim probe', { 'font-size': 11, fill: C.ink2, 'text-anchor': 'end' });
    });
    const cur = normexp(ST.bgx, mu, sigma, alpha);
    S('circle', { cx: X(ST.bgx), cy: Y(cur.s), r: 6, fill: C.accentInk, stroke: C.bg, 'stroke-width': 2 }, g);
    $('maBgEq').innerHTML = `${V('m')} = ${NB(fx(ST.bgx, 1))} − ${NB(fx(mu, 1))} − ${frac(`${NB(fx(sigma, 1))}²`, NB(fx(alpha, 1)))} = ${NB(fx(cur.m, 1))} &nbsp;→&nbsp; ${V('ŝ')} = ${NB(fx(cur.s, 1))} &nbsp; (plain subtraction: ${NB(fx(ST.bgx - mu, 1))}) &nbsp;→&nbsp; log₂(${V('ŝ')} + ${offset}) = ${NB(fx(Math.log2(cur.s + offset), 2))}` +
      `<span class="c">μ = ${fx(mu, 1)}, σ = ${fx(sigma, 1)} and α = ${fx(alpha, 1)} for ${bg.array}. The offset of ${fmt(offset)} keeps values near the background from spreading over many log₂ units.</span>`;
  }
  function drawQN() {
    const rows = readGrid('p', 4), cols = [0, 1, 2].map(j => rows.map(r => r[j])), q = quantileNormalize(cols), w = wOf('maQnChart'), h = 230, m = { l: 34, r: 8, t: 22, b: 26 };
    const g = chart('maQnChart', w, h, 'Dot chart of the invented arrays before and after quantile normalization.');
    const lv = v => Math.log2(Math.max(v, 0.25)), all = [...rows.flat(), ...q.target].map(lv), Y = lin(Math.min(...all) - 0.3, Math.max(...all) + 0.3, h - m.b, m.t);
    gridY(g, Y, ticks(Y.inv(h - m.b), Y.inv(m.t), 4), m.l, w - m.r, v => fx(v, 0));
    T(g, 4, m.t - 8, 'log₂', { 'font-size': 11, fill: C.muted });
    const colW = (w - m.l - m.r) / 3;
    ['A', 'B', 'C'].forEach((name, j) => {
      const x0 = m.l + colW * j + colW * 0.3, x1 = m.l + colW * j + colW * 0.7;
      T(g, (x0 + x1) / 2, h - 8, name, { 'text-anchor': 'middle', 'font-size': 12, fill: C.ink2, 'font-weight': 600 });
      cols[j].forEach((v, p) => {
        S('line', { x1: x0, y1: Y(lv(v)), x2: x1, y2: Y(lv(q.out[j][p])), stroke: C.lineStrong }, g);
        S('circle', { cx: x0, cy: Y(lv(v)), r: 4, fill: C.bg, stroke: C.ink2, 'stroke-width': 1.4 }, g);
        hover(S('circle', { cx: x1, cy: Y(lv(q.out[j][p])), r: 4, fill: C.accentInk }, g), `<b>p${p + 1}, array ${name}</b><br>${num(v)} → ${fx(q.out[j][p], 2)}`);
      });
    });
    const tbl = `<table class="t mgrid"><thead><tr><th>Rank</th><th class="r">A</th><th class="r">B</th><th class="r">C</th><th class="r">Mean: the target</th></tr></thead><tbody>${q.target.map((t, k) => `<tr><td>${k + 1}</td>${q.sorted.map(s => `<td class="r">${num(s[k])}</td>`).join('')}<td class="r"><b>${fx(t, 2)}</b></td></tr>`).join('')}</tbody></table>`;
    const out = `<table class="t mgrid"><thead><tr><th>Probe</th><th class="r">A</th><th class="r">B</th><th class="r">C</th></tr></thead><tbody>${rows.map((_, p) => `<tr><td>p${sub(p + 1)}</td>${[0, 1, 2].map(j => `<td class="r">${q.out[j][p] > 0 ? fx(Math.log2(q.out[j][p]), 2) : '–'}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
    $('maQnOut').innerHTML = `<div class="grid2"><div><b>Sorted, and the mean of each rank</b>${tbl}</div><div><b>Normalized, on the log₂ scale</b>${out}</div></div>`;
  }
  function drawMP() {
    const x = readGrid('q', 3), p = medianPolish(x), means = [0, 1, 2].map(j => mean(x.map(r => r[j])));
    const order = v => [0, 1, 2].slice().sort((a, b) => v[b] - v[a]).map(j => 'ABC'[j]).join(' > ');
    $('maMpOut').innerHTML = `<b>Median polish</b> (${spell(p.rounds)} ${p.rounds === 1 ? 'round' : 'rounds'}): common level ${V('μ')}′ = <b>${fx(p.t, 2)}</b>, probe levels ${p.r.map((v, i) => `q${sub(i + 1)} ${fx(p.t + v, 2)}`).join(', ')}; array effects ${p.c.map((v, j) => `${'ABC'[j]} ${sgn(v, 2)}`).join(', ')}.<br>
      Values reported by RMA: ${p.values.map((v, j) => `${'ABC'[j]} <b>${fx(v, 2)}</b>`).join(', ')} (order ${order(p.values)}).<br>
      Plain means of the columns: ${means.map((v, j) => `${'ABC'[j]} ${fx(v, 2)}`).join(', ')} (order ${order(means)}).<br>
      <span class="small">A median is not moved by one extreme value, so a single bad probe changes the plain means but not what RMA reports.</span>`;
  }
  function buildArraysText() {
    const r = M.arrays.rma, raw = r.each.map(e => e.raw), rmaV = r.each.map(e => e.rma), lo = Math.min(...raw), hi = Math.max(...raw);
    const mlo = Math.min(...rmaV), mhi = Math.max(...rmaV);
    $('maEx').innerHTML = `${r.study}, a study of children with septic shock and of healthy controls, holds ${fmt(r.arrays)} ${esc(r.chip)} arrays scanned between ${r.first_year} and ${r.last_year}. Its design package, ${code(r.package)}, lists ${fmt(r.pm_probes)} probes in ${fmt(r.probe_sets)} probe sets, about ${fmt(Math.round(r.pm_probes / r.probe_sets))} probes per set. The first array to be scanned, ${r.array}, is the file ${code(r.cel)}, scanned on ${dayOf(r.scan_day)}. The median intensity of its probes is ${fmt(r.raw_median)}, or ${fx(Math.log2(r.raw_median), 2)} on the log₂ scale.
      Across the ${fmt(r.arrays)} arrays the median raw intensity runs from ${fmt(lo)} to ${fmt(hi)}, a ${fx(hi / lo, 1)}-fold range, from ${fx(Math.log2(lo), 2)} to ${fx(Math.log2(hi), 2)} on the log₂ scale. After RMA the ${fmt(r.pm_probes)} intensities of ${r.array} have become ${fmt(r.probe_sets)} values with a median of <b>${fx(r.array_median, 2)}</b>, and the medians of all arrays lie between ${fx(mlo, 2)} and ${fx(mhi, 2)}. They are not identical, because quantile normalization gives the probes one distribution, while the summaries of the probe sets need not share one.`;
    const bg = M.arrays.background;
    $('maBgTable').innerHTML = `<table class="t"><thead><tr><th>Probe</th><th class="r">${V('x')}</th><th class="r">${V('x')} − ${V('μ')}</th><th class="r">${V('m')}</th><th class="r">${V('ŝ')}</th><th class="r">log₂(${V('ŝ')} + ${bg.offset})</th></tr></thead><tbody>${bg.probes.map(p => {
      const e = normexp(p.x, bg.mu, bg.sigma, bg.alpha);
      return `<tr><td>${{ median: 'Median probe of the array', bright: 'Bright probe (top 1%)', dim: 'Invented dim probe' }[p.probe]}</td><td class="r">${fx(p.x, 1)}</td><td class="r">${fx(p.x - bg.mu, 1)}</td><td class="r">${fx(e.m, 1)}</td><td class="r">${fx(e.s, 1)}</td><td class="r">${fx(Math.log2(e.s + bg.offset), 2)}</td></tr>`;
    }).join('')}</tbody></table><p class="note">limma estimated μ = ${fx(bg.mu, 1)}, σ = ${fx(bg.sigma, 1)} and α = ${fx(bg.alpha, 1)}, so σ²/α = ${fx(bg.sigma ** 2 / bg.alpha, 2)}. Only ${pct(bg.detected, 1)} of the probes of this array have a detection p-value below ${num(bg.detected_below)}: on a blood array most probes lie near the background. For the bright probe the correction subtracts the background and the small term σ²/α, nothing more; the median probe keeps more than plain subtraction leaves; the invented dim probe, below the background mean, gets a small positive value.</p>`;
    const rt = M.arrays.routes, LAB = { raw_array_files: 'Raw array files', unnormalized_intensities: 'Unnormalized intensities', counts: 'Read counts', author_normalized: 'Values normalized by the authors', not_processed: 'Not processed' };
    ST.routes = Object.keys(LAB).filter(k => k in rt).map(k => [LAB[k], rt[k], k]);
    $('maRes').innerHTML = RES(`Each study platform is processed on its own, from the lowest data level that covers its samples and from one source: <b>${fmt(rt.raw_array_files)}</b> study platforms from the raw files of every array (RMA for the ${fmt(M.arrays.affymetrix.study_platforms)} Affymetrix study platforms, with ${fmt(M.arrays.affymetrix.arrays)} arrays, and normexp for the scanner files of other arrays), <b>${fmt(rt.unnormalized_intensities)}</b> from intensities exported before normalization, <b>${fmt(rt.counts)}</b> from RNA-seq read counts (next chapter) and <b>${fmt(rt.author_normalized)}</b> from values normalized by their authors, kept as deposited. ${Spell(rt.not_processed)} holds microRNA only and was not processed.`);
  }
  function drawArrayDots() {
    const r = M.arrays.rma, w = wOf('maDots'), h = 190, m = { l: small(w) ? 12 : 170, r: 16, t: 16, b: 38 };
    const g = chart('maDots', w, h, `Dot chart of the median intensity of every array of ${r.study} before and after RMA.`);
    const v1 = r.each.map(e => Math.log2(e.raw)), v2 = r.each.map(e => e.rma), lo = Math.floor(Math.min(...v1, ...v2) * 2) / 2, hi = Math.ceil(Math.max(...v1, ...v2) * 2) / 2;
    const X = lin(lo, hi, m.l, w - m.r), y1 = m.t + 30, y2 = h - m.b - 22, rnd0 = rng(26440);
    gridX(g, X, ticks(lo, hi, 6), m.t, h - m.b, v => fx(v, 1));
    axisTitle(g, (m.l + w - m.r) / 2, h - 4, 'Median, log₂ scale');
    if (!small(w)) { T(g, m.l - 10, y1 + 4, 'Raw intensities', { 'text-anchor': 'end', 'font-size': 12, fill: C.ink2 }); T(g, m.l - 10, y2 + 4, 'After RMA', { 'text-anchor': 'end', 'font-size': 12, fill: C.ink2 }); }
    else { T(g, m.l, y1 - 18, 'Raw intensities', { 'font-size': 11, fill: C.muted }); T(g, m.l, y2 + 26, 'After RMA', { 'font-size': 11, fill: C.muted }); }
    const jit = r.each.map(() => [(rnd0() - 0.5) * 18, (rnd0() - 0.5) * 18]);
    r.each.forEach((e, i) => S('line', { x1: X(v1[i]), y1: y1 + jit[i][0], x2: X(v2[i]), y2: y2 + jit[i][1], stroke: C.grid }, g));
    r.each.forEach((e, i) => {
      const me = e.accession === r.array, html = `<b>${esc(e.accession)}</b><br>raw median ${fmt(e.raw, e.raw % 1 ? 1 : 0)} (log₂ ${fx(v1[i], 2)})<br>median after RMA ${fx(e.rma, 2)}`;
      hover(S('circle', { cx: X(v1[i]), cy: y1 + jit[i][0], r: me ? 6 : 3, fill: me ? C.cls.sepsis : C.ink2, 'fill-opacity': me ? 1 : 0.55 }, g), html);
      hover(S('circle', { cx: X(v2[i]), cy: y2 + jit[i][1], r: me ? 6 : 3, fill: me ? C.cls.sepsis : C.accentInk, 'fill-opacity': me ? 1 : 0.6 }, g), html);
    });
    const k = r.each.findIndex(e => e.accession === r.array);
    T(g, X(v1[k]) + 9, y1 + jit[k][0] - 7, r.array, { 'font-size': 11.5, fill: C.ink, 'font-weight': 600 });
    // the routes, as bars
    const rows = ST.routes, w2 = wOf('maRoutes'), left = small(w2) ? 150 : 230, rowH = 26, H2 = rowH * rows.length + 34;
    const g2 = chart('maRoutes', w2, H2, 'Bar chart of the study platforms by processing route.'), xm = niceMax(Math.max(...rows.map(x => x[1]))), X2 = lin(0, xm, left, w2 - 40);
    gridX(g2, X2, ticks(0, xm, 4), 4, rowH * rows.length + 4);
    axisTitle(g2, (left + w2 - 40) / 2, H2 - 2, 'Study platforms');
    rows.forEach(([label, n], i) => { const y = 6 + i * rowH; T(g2, left - 8, y + 14, label, { 'text-anchor': 'end', 'font-size': 12, fill: C.ink2 }); rrect(g2, X2(0), y + 3, X2(n) - X2(0), rowH - 8, C.accent, 3); T(g2, X2(n) + 5, y + 14, fmt(n), { 'font-size': 11.5, fill: C.muted }); hitRow(g2, 0, y, w2, rowH, `<b>${esc(label)}</b><br>${fmt(n)} study platforms`); });
  }

  // ======================================================================= 3. RNA-seq to numbers
  function logCPM(y, L, f, meanN, prior) { const N = f * L, c = prior * N / meanN; return { N, c, ratio: (y + c) / (N + 2 * c) * 1e6, value: Math.log2((y + c) / (N + 2 * c) * 1e6), zero: Math.log2(c / (N + 2 * c) * 1e6) }; }
  const zeroOf = (meanMillions, prior) => Math.log2(prior / (meanMillions * 1e6 + 2 * prior) * 1e6);
  function buildCounts() {
    const c = M.counts, k = c.cpm, e = c.expressed;
    ST.cpm = { y: k.count, L: k.library, f: k.factor, prior: true };
    $('mbCounts').innerHTML = `
      ${H3('eq', 'Equation')}
      ${EQ([`log₂ CPM${sub('gj')} = log₂ ( ${frac(`${V('y')}${sub('gj')} + ${V('c')}${sub('j')}`, `${V('N')}${sub('j')} + 2${V('c')}${sub('j')}`)} × 10⁶ ), &nbsp; ${V('N')}${sub('j')} = ${V('f')}${sub('j')} ${V('L')}${sub('j')}, &nbsp; ${V('c')}${sub('j')} = ${V('c')} ${frac(`${V('N')}${sub('j')}`, `${V('N̄')}`)}`,
        `a gene without reads: &nbsp; ${frac(`${V('c')}${sub('j')}`, `${V('N')}${sub('j')} + 2${V('c')}${sub('j')}`)} = ${frac(V('c'), `${V('N̄')} + 2${V('c')}`)}`],
        `y<sub>gj</sub>: the count of reads of gene g in sample j. L<sub>j</sub>: the library size, the total count of the sample. f<sub>j</sub>: its TMM factor, which corrects the library size where a few genes take a large share of the reads; the factors of the samples normalized together multiply to one. N<sub>j</sub>: the effective library size, and N̄ its mean over the samples normalized together. c = ${fmt(k.prior)}: the prior count, scaled to each library as c<sub>j</sub>, so that it is the same small amount relative to the size of every library.`)}
      ${H3('sim', 'Simulation', `Change the count and the library, starting from ${k.gene} in ${k.sample}`)}
      <div class="controls">
        <label>Reads of the gene ${V('y')} = <b id="mcYv"></b> <input type="range" id="mcY" min="0" max="600" step="1"></label>
        <label>Library size ${V('L')} = <b id="mcLv"></b> <input type="range" id="mcL" min="600" max="800" step="1"></label>
        <label>TMM factor ${V('f')} = <b id="mcFv"></b> <input type="range" id="mcF" min="0.5" max="2" step="0.001"></label>
        <label><input type="checkbox" id="mcP" checked> add the prior count</label>
        <button type="button" class="btn ghost" id="mcReset">Back to ${k.sample}</button>
      </div>
      <div class="chart" id="mcChart"></div>
      <div class="eq" id="mcEq"></div>
      ${how(`The curve is the log₂ count per million of a gene in a library of this size, for every number of reads (the axis is logarithmic, with zero at the left edge); the dot is the gene. Without the prior count (dashed) a gene with no reads has no value and a gene with a few reads swings widely. The dotted line is the value of zero reads, the same for every library of a normalization group. N̄ stays at the mean of ${k.study}, ${fx(k.mean_effective / 1e6, 2)} million.`)}
      ${H3('ex', 'Worked example', `${k.gene} in ${k.sample} of ${k.study}`)}
      <ol class="msteps" id="mcEx"></ol>
      ${H3('eq', 'Equation', 'The expressed-gene flag')}
      ${EQ([`#{ ${V('j')} : ${frac(`${V('y')}${sub('gj')}`, `${V('N')}${sub('j')}`)} × 10⁶ ≥ ${frac(V('k'), V('Ñ'))} × 10⁶ } ≥ ${V('n')}* &nbsp; and &nbsp; Σ${sub('j')} ${V('y')}${sub('gj')} ≥ ${V('t')}`],
        `The rule of filterByExpr in edgeR, with its default limits. Ñ: the median effective library size of the study; k = ${fmt(e.minimum_count)}, so the limit is the count per million that ${fmt(e.minimum_count)} reads give in a library of median size; t = ${fmt(e.minimum_total)} reads in all. n*: the number of samples in the smallest group, or, when that group holds more than ${fmt(e.large_n)}, ${fmt(e.large_n)} plus ${pct(e.proportion, 0)} of the excess. The flag removes nothing; it tells which values rest on reads.`)}
      ${H3('ex', 'Worked example', `the flag in ${k.study}`)}
      <div class="readout" id="mcFlag"></div>
      <div class="chart" id="mcFlagChart"></div>
      ${H3('res', 'Result', 'The value of a gene without reads, in every normalization group')}
      <div class="chart" id="mcZero"></div>
      ${how('One dot per normalization group of the count tables: the mean effective library size of its samples and the log₂ count per million that a gene with no reads receives there. Every dot lies on the curve c/(N̄ + 2c) × 10⁶, so the value depends only on how deeply the group was sequenced. A value near the bottom of a table may rest on no reads at all.')}
      <div id="mcRes"></div>`;
    const setInputs = () => { $('mcY').value = Math.round(100 * Math.log10(ST.cpm.y + 1)); $('mcL').value = Math.round(100 * Math.log10(ST.cpm.L)); $('mcF').value = ST.cpm.f; $('mcP').checked = ST.cpm.prior; };
    setInputs();
    $('mcY').addEventListener('input', e => { ST.cpm.y = Math.round(10 ** (+e.target.value / 100) - 1); drawCPM(); });
    $('mcL').addEventListener('input', e => { ST.cpm.L = Math.round(10 ** (+e.target.value / 100)); drawCPM(); });
    $('mcF').addEventListener('input', e => { ST.cpm.f = +e.target.value; drawCPM(); });
    $('mcP').addEventListener('change', e => { ST.cpm.prior = e.target.checked; drawCPM(); });
    $('mcReset').onclick = () => { ST.cpm = { y: k.count, L: k.library, f: k.factor, prior: true }; setInputs(); drawCPM(); };
    const r = logCPM(k.count, k.library, k.factor, k.mean_effective, k.prior);
    $('mcEx').innerHTML = `
      <li>The library holds ${V('L')} = ${fmt(k.library)} reads, and ${k.gene} (row ${esc(k.feature)}, its Entrez gene) has ${V('y')} = ${fmt(k.count)} of them. Without any correction that is ${fmt(k.count)}/${fx(k.library / 1e6, 2)} = ${fx(k.count / k.library * 1e6, 1)} counts per million.</li>
      <li>The TMM factor of the library is ${V('f')} = ${fx(k.factor, 3)}, so it counts as a library of ${V('N')} = ${fx(k.factor, 3)} × ${fx(k.library / 1e6, 2)} = ${fx(r.N / 1e6, 2)} million reads.</li>
      <li>The mean effective library size of the ${fmt(k.libraries)} libraries of ${k.study}, normalized as one group, is ${V('N̄')} = ${fx(k.mean_effective / 1e6, 2)} million, so the prior count of this library is ${V('c')}${sub('j')} = ${fmt(k.prior)} × ${fx(r.N / 1e6, 2)}/${fx(k.mean_effective / 1e6, 2)} = ${fx(r.c, 2)}.</li>
      <li>The equation gives (${fmt(k.count)} + ${fx(r.c, 2)})/(${fx(r.N / 1e6, 2)} million + 2 × ${fx(r.c, 2)}) × 10⁶ = ${fx(r.ratio, 1)}, and log₂ ${fx(r.ratio, 1)} = <b>${fx(r.value, 2)}</b>, the value in the table (${fx(k.value, 2)}).</li>
      <li>A gene with no reads in this library gets ${fx(r.c, 2)}/(${fx(r.N / 1e6, 2)} million + 2 × ${fx(r.c, 2)}) × 10⁶ = ${fx(Math.exp(r.zero * Math.LN2), 3)}, or <b>${fx(r.zero, 2)}</b> on the log₂ scale, as in every other library of the study.</li>`;
    const need = e.large_n + (k.libraries - e.large_n) * e.proportion, cut = e.minimum_count / e.median_effective * 1e6, G = Object.fromEntries(e.checked.map(x => [x.gene, x]));
    ST.flag = { need, cut };
    $('mcFlag').innerHTML = `${k.study} is one group of ${fmt(k.libraries)} libraries, more than ${fmt(e.large_n)}, so ${V('n')}* = ${fmt(e.large_n)} + ${num(e.proportion)} × (${fmt(k.libraries)} − ${fmt(e.large_n)}) = <b>${fx(need, 1)}</b>: a gene must reach the limit in at least that many libraries. The median effective library size is ${V('Ñ')} = ${fx(e.median_effective / 1e6, 2)} million, so the limit is ${fmt(e.minimum_count)}/${fx(e.median_effective / 1e6, 2)} = <b>${fx(cut, 2)}</b> counts per million. Of the ${fmt(e.genes)} genes of the table, <b>${fmt(e.flagged)}</b> are flagged as expressed.
      ${k.gene} reaches the limit in all ${fmt(G[k.gene].samples)} libraries and is flagged. XIST, expressed from the inactive X chromosome of women, reaches it in ${fmt(G.XIST.samples)}; the study has ${fmt(e.women)} women, too few to reach ${fx(need, 1)}, so XIST is not flagged although it is plainly expressed in them. IL10 reaches it in ${fmt(G.IL10.samples)} libraries and narrowly misses the flag.`;
    const src = c.sources;
    $('mcRes').innerHTML = RES(`<b>${fmt(c.study_platforms)}</b> study platforms provide counts, <b>${fmt(c.samples)}</b> libraries in all: ${fmt(src.ncbi_counts)} from the counts that NCBI computes for GEO with one pipeline of its own, ${fmt(src.table)} from the authors' count tables and ${fmt(src.per_sample_files)} from one file per sample. Across the ${fmt(c.zeros.length)} normalization groups a gene without reads receives between ${fx(Math.min(...c.zeros.map(z => z.zero)), 2)} and ${fx(Math.max(...c.zeros.map(z => z.zero)), 2)}. No gene was filtered out.`);
  }
  function drawCPM() {
    const k = M.counts.cpm, s = ST.cpm, prior = s.prior ? k.prior : 0, w = wOf('mcChart'), h = small(w) ? 240 : 270, m = { l: 52, r: 18, t: 16, b: 42 };
    $('mcYv').textContent = fmt(s.y); $('mcLv').textContent = `${fx(s.L / 1e6, 1)} million`; $('mcFv').textContent = fx(s.f, 3);
    const g = chart('mcChart', w, h, 'Curve of the log2 count per million against the number of reads of a gene, for the chosen library.');
    const X = lin(0, 6, m.l, w - m.r), Y = lin(-8, 18, h - m.b, m.t);
    gridY(g, Y, [-8, -4, 0, 4, 8, 12, 16], m.l, w - m.r);
    [0, 1, 10, 100, 1000, 1e4, 1e5, 1e6].forEach(t => { const xx = X(Math.log10(t + 1)); S('line', { x1: xx, x2: xx, y1: m.t, y2: h - m.b, stroke: C.grid }, g); if (!small(w) || [0, 100, 1e4, 1e6].includes(t)) T(g, xx, h - m.b + 15, fmt(t), { 'text-anchor': 'middle', 'font-size': 11, fill: C.muted }); });
    axisTitle(g, (m.l + w - m.r) / 2, h - 6, 'Reads of the gene, y'); axisTitle(g, 14, (m.t + h - m.b) / 2, 'log₂ counts per million', true);
    const ys = []; for (let v = 0; v <= 6.0001; v += 0.02) ys.push(10 ** v - 1);
    const withP = logCPM(0, s.L, s.f, k.mean_effective, k.prior), N = withP.N;
    S('path', { d: ys.filter(y => y > 0).map((y, i) => `${i ? 'L' : 'M'}${X(Math.log10(y + 1))},${Y(Math.max(-8, Math.log2(y / N * 1e6)))}`).join(''), stroke: C.muted, 'stroke-dasharray': '5 4', fill: 'none', 'stroke-width': 1.5 }, g);
    S('path', { d: ys.map((y, i) => `${i ? 'L' : 'M'}${X(Math.log10(y + 1))},${Y(logCPM(y, s.L, s.f, k.mean_effective, k.prior).value)}`).join(''), stroke: C.fam.rna_seq, fill: 'none', 'stroke-width': 2.4 }, g);
    if (s.prior) { S('line', { x1: m.l, x2: w - m.r, y1: Y(withP.zero), y2: Y(withP.zero), stroke: C.fam.rna_seq, 'stroke-dasharray': '1 3' }, g); T(g, w - m.r, Y(withP.zero) - 5, `zero reads: ${fx(withP.zero, 2)}`, { 'text-anchor': 'end', 'font-size': 11, fill: C.ink2 }); }
    const r = logCPM(s.y, s.L, s.f, k.mean_effective, prior), val = prior ? r.value : (s.y > 0 ? Math.log2(s.y / r.N * 1e6) : -Infinity);
    if (isFinite(val)) S('circle', { cx: X(Math.log10(s.y + 1)), cy: Y(clamp(val, -8, 18)), r: 6, fill: C.cls.sepsis, stroke: C.bg, 'stroke-width': 2 }, g);
    $('mcEq').innerHTML = `${V('N')} = ${fx(s.f, 3)} × ${fx(s.L / 1e6, 2)} million = ${fx(r.N / 1e6, 2)} million; &nbsp; ${V('c')}${sub('j')} = ${prior ? `${fmt(k.prior)} × ${fx(r.N / 1e6, 2)}/${fx(k.mean_effective / 1e6, 2)} = ${fx(r.c, 2)}` : '0 (no prior count)'}<br>log₂ CPM = log₂ ( ${frac(`${fmt(s.y)} + ${fx(r.c, 2)}`, `${fx(r.N / 1e6, 2)} million + ${fx(2 * r.c, 2)}`)} × 10⁶ ) = ${NB(isFinite(val) ? fx(val, 2) : 'no value: the logarithm of zero does not exist')}` +
      `<span class="c">${prior ? 'With the prior count every gene has a value, and a gene with zero or one read does not swing widely.' : 'Without the prior count, a gene with no reads has no value, and a gene with one read more or less moves by a whole log₂ unit when counts are small.'}</span>`;
  }
  function drawFlagChart() {
    const e = M.counts.expressed, k = M.counts.cpm, w = wOf('mcFlagChart'), left = 90, rowH = 30, top = 18, H = rowH * e.checked.length + 40 + top, g = chart('mcFlagChart', w, H, 'Bar chart of the libraries in which three genes reach the limit of the expressed flag.');
    const X = lin(0, k.libraries, left, w - 70);
    gridX(g, X, ticks(0, k.libraries, 5), top + 6, top + rowH * e.checked.length + 6);
    axisTitle(g, (left + w - 70) / 2, H - 2, `Libraries of ${k.study} at or above ${fx(ST.flag.cut, 2)} counts per million`);
    e.checked.forEach((x, i) => { const y = top + 8 + i * rowH; T(g, left - 8, y + 15, x.gene, { 'text-anchor': 'end', 'font-size': 12.5, fill: C.ink2, 'font-weight': 600 }); rrect(g, X(0), y + 4, X(x.samples) - X(0), rowH - 10, x.flagged ? C.fam.rna_seq : C.lineStrong, 3); T(g, X(x.samples) + 5, y + 16, `${fmt(x.samples)} · ${x.flagged ? 'flagged' : 'not flagged'}`, { 'font-size': 11.5, fill: C.muted }); hitRow(g, 0, y, w, rowH, `<b>${x.gene}</b><br>${fmt(x.samples)} libraries reach the limit; ${fmt(x.total)} reads in all`); });
    S('line', { x1: X(ST.flag.need), x2: X(ST.flag.need), y1: top - 2, y2: top + rowH * e.checked.length + 8, stroke: C.cls.sepsis, 'stroke-width': 1.6, 'stroke-dasharray': '4 3' }, g);
    T(g, X(ST.flag.need), top - 6, `n* = ${fx(ST.flag.need, 1)}`, { 'text-anchor': 'middle', 'font-size': 11, fill: C.cls.sepsis });
  }
  function drawZeros() {
    const c = M.counts, z = c.zeros, w = wOf('mcZero'), h = small(w) ? 240 : 270, m = { l: 52, r: 18, t: 16, b: 42 };
    const g = chart('mcZero', w, h, 'Dot chart of the value of a zero count against the mean library size of every normalization group, on the curve that sets it.');
    const xs = z.map(r => r.mean_millions), lo = Math.min(...xs) / 1.5, hi = Math.max(...xs) * 1.5, X = lg(lo, hi, m.l, w - m.r);
    const vs = z.map(r => r.zero), Y = lin(Math.floor(Math.min(...vs)) - 0.5, Math.ceil(Math.max(...vs)) + 0.5, h - m.b, m.t);
    gridY(g, Y, ticks(Y.inv(h - m.b), Y.inv(m.t), 5), m.l, w - m.r);
    [0.1, 0.3, 1, 3, 10, 30, 100].filter(t => t >= lo && t <= hi).forEach(t => { S('line', { x1: X(t), x2: X(t), y1: m.t, y2: h - m.b, stroke: C.grid }, g); T(g, X(t), h - m.b + 15, num(t), { 'text-anchor': 'middle', 'font-size': 11, fill: C.muted }); });
    axisTitle(g, (m.l + w - m.r) / 2, h - 6, 'Mean effective library size, millions of reads'); axisTitle(g, 14, (m.t + h - m.b) / 2, 'log₂ CPM of zero reads', true);
    const pts = []; for (let v = Math.log10(lo); v <= Math.log10(hi); v += 0.02) pts.push(10 ** v);
    S('path', { d: pts.map((t, i) => `${i ? 'L' : 'M'}${X(t)},${Y(zeroOf(t, c.cpm.prior))}`).join(''), stroke: C.lineStrong, fill: 'none', 'stroke-width': 1.5 }, g);
    const ends = [z.reduce((a, b) => b.zero < a.zero ? b : a), z.reduce((a, b) => b.zero > a.zero ? b : a)];
    z.forEach(r => hover(S('circle', { cx: X(r.mean_millions), cy: Y(r.zero), r: r.study === c.cpm.study ? 6 : 4, fill: r.study === c.cpm.study ? C.cls.sepsis : C.fam.rna_seq, 'fill-opacity': 0.8, stroke: C.bg }, g),
      `<b>${esc(r.study)}</b> <span class="m">${esc(r.platform)}${r.group ? ` · ${esc(r.group)}` : ''}</span><br>${fmt(r.libraries)} libraries, mean ${fx(r.mean_millions, 2)} million reads<br>zero reads: ${fx(r.zero, 2)}`));
    ends.forEach((r, i) => T(g, X(r.mean_millions) + (i ? 8 : 0), Y(r.zero) + (i ? 4 : 20), r.study, { 'font-size': 11, fill: C.ink2, 'text-anchor': i ? 'start' : 'end' }));
  }

  // ======================================================================= 4. probes to genes
  const STATUS = ['one_gene', 'several_genes', 'no_gene_in_annotation', 'gene_not_in_hgnc', 'accession_not_resolved'];
  const statusLabel = k => (M.vis.mapping_labels || {})[k] || String(k).replace(/_/g, ' ');
  // the figure grey of rows without a gene would vanish on the paper background; the page's grey of 'not recorded' stands in
  const statusColor = k => { const c = (M.vis.mapping_colors || {})[k] || C.muted; return /^#e0e0e0$/i.test(c) ? C.cls.not_recorded : c; };
  function buildGenes() {
    const g = M.genes, h = g.hla, kept = h.probe_sets.find(p => p.feature === h.kept), other = h.probe_sets.find(p => p.feature !== h.kept);
    ST.rule = 'highest';
    $('mbGenes').innerHTML = `
      ${H3('eq', 'Equation')}
      ${EQ([`arrays and log-scale values: &nbsp; ${V('x')}${sub('gj')} = ${V('x')}${sub('r*j')}, &nbsp; ${V('r')}* = the row of gene ${V('g')} with the highest mean ${frac('1', V('n'))} Σ${sub('j')} ${V('x')}${sub('rj')}`,
        `RNA-seq counts: &nbsp; ${V('y')}${sub('gj')} = Σ${sub('r')} ${V('y')}${sub('rj')} (identical rows counted once), then log₂ CPM again`],
        'g: an HGNC gene; r: a row (feature) of the matrix that maps to that gene alone; j: a sample, and n the samples of the study. On arrays one measured row is kept unchanged, the one with the highest mean in the study, so that every value of the gene-level table was measured and is not an average of probes that may behave differently. Values deposited on a linear scale are added up, like counts. A row that maps to several genes is kept in the row-level table but left out of the gene-level table.')}
      ${H3('sim', 'Simulation', `The two probe sets of ${h.gene} in ${h.study}`)}
      <div class="controls"><span>Gene-level value from</span><div class="seg" id="mgRule">${segHTML([['highest', 'the row with the higher mean (used)'], ['average', 'the average of the two rows'], ['other', 'the other row']], ST.rule)}</div></div>
      <div class="chart" id="mgChart"></div>
      <div class="readout" id="mgOut"></div>
      ${how(`One point per array of ${h.study}, in the order of scanning, for each probe set of ${h.gene}; the thick line is the gene-level value under the chosen rule, and the vertical line marks ${h.array}. Hover for the values of an array.`)}
      ${H3('ex', 'Worked example', `${h.gene} in ${h.study}, and two rows with several genes`)}
      <div class="readout" id="mgEx"></div>
      <div class="chart" id="mgStatus"></div>
      ${how('The rows of one matrix by what their annotation says they measure. Only the rows that map to one gene enter the gene-level table; every row keeps its status in the map of its matrix.')}
      <div class="readout" id="mgCopy"></div>
      ${H3('res', 'Result', 'Every row of every matrix')}
      <div class="chart" id="mgAll"></div>
      <div id="mgRes"></div>`;
    seg('mgRule', v => { ST.rule = v; drawGenes(); });
    const st = g.statuses, totalRows = sum(Object.values(st));
    $('mgEx').innerHTML = `${h.gene}, part of the HLA-DR molecule of monocytes and other immune cells, is measured by ${spell(h.probe_sets.length)} probe sets of the ${h.study} chip, and both map to it alone, through their Entrez gene. Over the ${fmt(h.arrays.length)} arrays, ${code(kept.feature)} has a mean of ${fx(mean(kept.values), 2)} and ${code(other.feature)} a mean of ${fx(mean(other.values), 2)}; the first is kept. In ${h.array}, the gene-level value of ${h.gene} is therefore <b>${fx(kept.values[h.arrays.indexOf(h.array)], 2)}</b>, the value of ${code(kept.feature)}, and not ${fx(other.values[h.arrays.indexOf(h.array)], 2)}, nor an average of the two. The gene-level table of ${h.study} holds ${fmt(g.genes)} genes; ${fmt(g.genes_with_several_rows)} of them had more than one row to choose from.
      <br><br>Some rows match genes whose sequences are almost the same, and their intensity is the sum of what they bind. ${g.several.map(x => `${code(x.feature)} maps to ${andList(x.genes)}`).join('; ')}. Choosing one of the genes would credit it with RNA that may belong to the others, so such rows stay out of the gene-level table. Of the ${fmt(totalRows)} rows of ${h.study}, ${fmt(st.one_gene)} map to one gene, ${fmt(st.several_genes)} to several, ${fmt(st.no_gene_in_annotation)} have no gene in the annotation and ${fmt(st.gene_not_in_hgnc)} lead to genes outside the HGNC set.`;
    const cp = g.copy;
    $('mgCopy').innerHTML = `<h4 class="exhead">A count that is written twice <span class="m">${cp.study}</span></h4>In the count table of ${cp.study}, ${cp.gene}, which lies in the region shared by the X and Y chromosomes, has two rows, ${code(cp.rows[0])} and its Y copy ${code(cp.rows[1])}, with the same count in every sample: in ${cp.sample} both read ${fmt(cp.count)}. They are one measurement written twice, so the gene-level count is <b>${fmt(cp.count)}</b>, not ${fmt(2 * cp.count)}. This study has ${fmt(cp.identical)} such copies.`;
    const f = g.features;
    $('mgRes').innerHTML = RES(`The ${fmt(f.total)} rows of the ${fmt(g.tables.count)} matrices: <b>${fmt(f.one_gene)}</b> (${pct(f.one_gene / f.total, 0)}) map to one gene, ${fmt(f.several_genes)} to several, ${fmt(f.no_gene_in_annotation)} have no gene in their annotation, ${fmt(f.gene_not_in_hgnc)} lead to a gene the HGNC set does not list and ${fmt(f.accession_not_resolved)} carry only accessions that lead to no gene. A symbol is matched first as an approved symbol, then as a previous symbol and then as an alias: ${fmt(g.matched.previous_symbol)} rows were matched through a previous symbol and ${fmt(g.matched.alias_symbol)} through an alias (the author table of GSE310929, for example, writes ${andList(Object.entries(g.renamed).map(([o, n]) => `${o} for ${n}`))}). The gene-level tables hold between ${fmt(g.tables.fewest)} and ${fmt(g.tables.most)} genes, with a median of ${fmt(g.tables.median)}.`);
  }
  function stackBar(id, parts, total, label) {
    const w = wOf(id), H = 40, g = chart(id, w, H, label), X = lin(0, total, 0, w);
    let x0 = 0;
    parts.forEach(([k, n]) => { if (!n) return; const x1 = x0 + n; hover(rrect(g, X(x0), 8, Math.max(1, X(x1) - X(x0) - 1), 26, statusColor(k), 2), `<b>${esc(statusLabel(k))}</b><br>${fmt(n)} rows (${pct(n / total, 1)})`); x0 = x1; });
    const key = parts.filter(([, n]) => n);
    $(id).insertAdjacentHTML('beforeend', `<div class="keys">${keys(key.map(([k, n]) => [`${statusLabel(k)} ${fmt(n)}`, statusColor(k)]))}</div>`);
  }
  function drawGenes() {
    const h = M.genes.hla, w = wOf('mgChart'), hh = small(w) ? 230 : 260, m = { l: 46, r: 14, t: 14, b: 36 }, a = h.probe_sets[0], b = h.probe_sets[1];
    const gene = a.values.map((v, j) => ST.rule === 'highest' ? v : ST.rule === 'other' ? b.values[j] : (v + b.values[j]) / 2);
    const g = chart('mgChart', w, hh, `Line chart of the two probe sets of ${h.gene} across the arrays of ${h.study} and the gene-level value.`);
    const all = [...a.values, ...b.values], X = lin(0, h.arrays.length - 1, m.l, w - m.r), Y = lin(Math.floor(Math.min(...all)), Math.ceil(Math.max(...all)), hh - m.b, m.t);
    gridY(g, Y, ticks(Y.inv(hh - m.b), Y.inv(m.t), 5), m.l, w - m.r);
    axisTitle(g, (m.l + w - m.r) / 2, hh - 8, `Arrays of ${h.study}, in the order of scanning`); axisTitle(g, 12, (m.t + hh - m.b) / 2, 'log₂ value', true);
    const k = h.arrays.indexOf(h.array);
    S('line', { x1: X(k), x2: X(k), y1: m.t, y2: hh - m.b, stroke: C.lineStrong, 'stroke-dasharray': '3 3' }, g);
    T(g, X(k) + 4, m.t + 10, h.array, { 'font-size': 11, fill: C.muted });
    [[a, C.fam.affymetrix], [b, C.fam.agilent]].forEach(([p, col]) => p.values.forEach((v, j) => S('circle', { cx: X(j), cy: Y(v), r: 2.4, fill: col, 'fill-opacity': 0.75 }, g)));
    S('path', { d: gene.map((v, j) => `${j ? 'L' : 'M'}${X(j)},${Y(v)}`).join(''), stroke: C.ink, 'stroke-width': 1.8, fill: 'none', 'stroke-opacity': 0.85 }, g);
    h.arrays.forEach((acc, j) => hitRow(g, X(j) - (X(1) - X(0)) / 2, m.t, X(1) - X(0), hh - m.b - m.t, `<b>${esc(acc)}</b><br>${sw(C.fam.affymetrix)}${a.feature}: ${fx(a.values[j], 2)}<br>${sw(C.fam.agilent)}${b.feature}: ${fx(b.values[j], 2)}<br>gene-level value: ${fx(gene[j], 2)}`));
    $('mgChart').insertAdjacentHTML('beforeend', `<div class="keys">${keys([[`${a.feature} (mean ${fx(mean(a.values), 2)})`, C.fam.affymetrix, 'dot'], [`${b.feature} (mean ${fx(mean(b.values), 2)})`, C.fam.agilent, 'dot'], ['gene-level value', C.ink, 'line']])}</div>`);
    const kept = ST.rule === 'highest';
    $('mgOut').innerHTML = `Gene-level value of ${h.gene} in ${h.array}: <b>${fx(gene[k], 2)}</b>; its mean over the arrays: ${fx(mean(gene), 2)}. ${kept ? `This is the rule of the resource: ${code(a.feature)}, measured as it is.` : ST.rule === 'average' ? 'An average of two probe sets is a value that no probe measured; where the probes behave differently, as these do in some arrays, it blurs both.' : `The other probe set reads ${h.gene} lower in most arrays.`}`;
    const st = M.genes.statuses;
    stackBar('mgStatus', STATUS.filter(s => s in st).map(s => [s, st[s]]), sum(Object.values(st)), `Stacked bar of the rows of ${h.study} by mapping status.`);
    const f = M.genes.features;
    stackBar('mgAll', STATUS.map(s => [s, f[s] || 0]), f.total, 'Stacked bar of every row of every matrix by mapping status.');
  }

  // ======================================================================= 5. labels from records
  const FIELD_KEY = { title: 'title', 'abdominal sepsis and controls': 'abdominal', 'endotype class': 'endotype', 'pneumonia diagnoses': 'pneumonia' };
  // the class that a list of rules gives a record, and the number of the rule that decided (0: none)
  function applyRules(rules, rec) {
    for (let i = 0; i < rules.length; i++) {
      const r = rules[i], word = rec[FIELD_KEY[r.from]];
      if (word == null || word === '') continue;
      if (r.patterns) { const hit = r.patterns.find(([p]) => new RegExp(p, 'i').test(word)); if (hit) return [hit[1], r.n || i + 1]; continue; }
      const v = r.values[String(word).trim().toLowerCase()];
      if (v == null || v === 'next') continue;
      return [v, r.n || i + 1];
    }
    return ['not_recorded', 0];
  }
  const VALUE_LABEL = v => v === 'infection_sepsis_not_stated' ? 'other condition, infection without stated sepsis' : v === 'next' ? 'next rule' : labLower('class', v);
  function buildLabels() {
    const L = M.labels, rec = L.record, rules = L.rules.map((r, i) => Object.assign({ n: i + 1 }, r));
    ST.lab = { title: rec.title, abdominal: rec['abdominal sepsis and controls'], endotype: rec['endotype class'], pneumonia: rec['pneumonia diagnoses'], swap: false };
    const opts = (field, list) => `<select id="ml_${field}">${list.map(v => `<option value="${esc(v)}">${esc(v)}</option>`).join('')}</select>`;
    const words = r => r.patterns ? r.patterns.map(([p, v]) => `begins with ${code(p.replace(/^\^/, ''))} → ${esc(VALUE_LABEL(v))}`).join('; ') : Object.entries(r.values).map(([k, v]) => `${code(k)} → ${esc(VALUE_LABEL(v))}`).join('; ');
    const titles = [...new Set(L.records.map(r => r.title.replace(/\s*\[.*$/, '')))].sort();
    const choice = { abdominal: [...new Set(L.records.map(r => r.abdominal))].sort(), endotype: [...new Set(L.records.map(r => r.endotype))].sort(), pneumonia: [...new Set(L.records.map(r => r.pneumonia))].sort() };
    $('mbLabels').innerHTML = `
      ${H3('eq', 'Equation')}
      ${EQ([`label(${V('i')}) = ${V('v')}${sub(`k*`)}(${V('i')}), &nbsp; ${V('k')}* = the first ${V('k')} for which rule ${V('k')} gives record ${V('i')} a value other than “next”; &nbsp; no such ${V('k')}: not recorded`],
        'A rule reads one field of the record and gives a value from a fixed list, by a list of wordings or by a pattern; a wording a list does not name stops the whole run, unless the rule says what any other wording gives. A rule may also pass a record on with the value “next”. Where the record does not say enough but the study does, a rule may quote a sentence of the study record, which must be found there word for word. Nothing is inferred from the expression values.')}
      ${H3('sim', 'Simulation', `The class rules of ${L.study}, on a record you can change`)}
      <div class="mrules" id="mlRules"></div>
      <div class="controls">
        <label>title begins with ${opts('title', titles)}</label>
        <label>abdominal sepsis and controls ${opts('abdominal', choice.abdominal)}</label>
        <label>endotype class ${opts('endotype', choice.endotype)}</label>
        <label>pneumonia diagnoses ${opts('pneumonia', choice.pneumonia)}</label>
      </div>
      <div class="controls"><label><input type="checkbox" id="mlSwap"> try the pneumonia rule before the endotype rule</label><button type="button" class="btn ghost" id="mlReset">Back to ${L.sample}</button></div>
      <div class="readout" id="mlOut" aria-live="polite"></div>
      ${how(`The rules are tried from the top; the first that gives a value decides (green), the ones above it give nothing or pass the record on, and the ones below are not reached. The record starts as ${L.sample}. Changing the order of the rules applies it to all ${fmt(L.records.length)} records of the study, to show that the order is a decision.`)}
      ${H3('ex', 'Worked example', `the record ${L.sample}`)}
      <div class="readout" id="mlEx"></div>
      ${H3('res', 'Result')}
      <div class="tablewrap" id="mlTable"></div>
      <div class="chart" id="mlBasis"></div>
      <div id="mlRes"></div>`;
    ST.labRules = rules; ST.labWords = words;
    ['title', 'abdominal', 'endotype', 'pneumonia'].forEach(f => { const el = $(`ml_${f}`); el.value = f === 'title' ? ST.lab.title.replace(/\s*\[.*$/, '') : ST.lab[f]; el.onchange = () => { ST.lab[f] = el.value; drawLabels(); }; });
    $('mlSwap').onchange = e => { ST.lab.swap = e.target.checked; drawLabels(); };
    $('mlReset').onclick = () => { ST.lab = { title: rec.title, abdominal: rec['abdominal sepsis and controls'], endotype: rec['endotype class'], pneumonia: rec['pneumonia diagnoses'], swap: false }; $('mlSwap').checked = false; ['title', 'abdominal', 'endotype', 'pneumonia'].forEach(f => { $(`ml_${f}`).value = f === 'title' ? ST.lab.title.replace(/\s*\[.*$/, '') : ST.lab[f]; }); drawLabels(); };
    const lb = L.labels, days = lb.age * 365.25, adult = L.age_groups.adult_from_years;
    $('mlEx').innerHTML = `${L.study} holds ${fmt(L.records.length)} whole-blood arrays of the MARS consortium, from healthy subjects and from patients admitted to intensive care, some of them with sepsis. The record ${L.sample} holds, among other fields, the title ${code(rec.title)}, ${code(`gender: ${rec.gender}`)}, ${code(`age: ${rec.age}`)}, ${code(`abdominal sepsis and controls: ${rec['abdominal sepsis and controls']}`)}, ${code(`endotype class: ${rec['endotype class']}`)}, ${code(`pneumonia diagnoses: ${rec['pneumonia diagnoses']}`)} and ${code(`${lb.outcome_rule.from}: ${rec[lb.outcome_rule.from]}`)}.
      <ol class="msteps"><li><b>First rule.</b> The title does not begin with ${code(rules[0].patterns[0][0].replace(/^\^/, ''))}, so the rule gives nothing.</li>
      <li><b>Second rule.</b> The field holds ${code(rec['abdominal sepsis and controls'])}, which passes the record on.</li>
      <li><b>Third rule.</b> The endotype class ${code(rec['endotype class'])} gives <b>sepsis</b>. The word sepsis is not in the record; the Mars endotypes were defined by the investigators among patients with sepsis, so the rule quotes the study summary, “${esc(rules[2].study_states)}”, and the label rests on the record and the study description.</li>
      <li>The fourth and fifth rules are not reached. Had they been, ${code(rec['pneumonia diagnoses'])} (community-acquired pneumonia) would have made the record an other condition, an infection without stated sepsis.</li></ol>
      The other variables come from other rules. The wording ${code(rec.gender)} becomes ${esc(lb.sex)} through a list of wordings that all studies share. The age is read as ${fmt(lb.age)} years, and the age group is decided on the exact number of days, ${fmt(lb.age)} × 365.25 = ${fx(days, 1)}, at least ${fmt(adult)} years, so ${esc(lb.age_group)}. The outcome field gives ${labLower('outcome', lb.outcome).replace(/_/g, ' ')}, with the outcome definition ${esc(lb.outcome_definition)}. The record names no patient, time point or septic shock, so these are not recorded.`;
    $('mlRes').innerHTML = RES(`Each of the <b>${fmt(L.rule_sets)}</b> study platforms has its own rules, written after reading the complete list of its recorded values. Over all ${fmt(sum(Object.values(L.basis)))} records, the class rests on the sample record alone for <b>${fmt(L.basis.sample_record)}</b>, on the record together with a quoted study sentence for <b>${fmt(L.basis.sample_record_and_study_description)}</b> and on the study description alone for <b>${fmt(L.basis.study_description)}</b>${L.basis.not_recorded ? `; ${spell(L.basis.not_recorded)} could not be decided` : ''}. The rules quote ${fmt(L.study_sentences)} study sentences. Infections without stated sepsis, kept as other conditions, account for ${fmt(L.infection_not_stated)} records.`);
  }
  function labelRecord() { return { title: ST.lab.title, abdominal: ST.lab.abdominal, endotype: ST.lab.endotype, pneumonia: ST.lab.pneumonia }; }
  const ruleOrder = swap => { const r = ST.labRules.slice(); if (swap) [r[2], r[3]] = [r[3], r[2]]; return r; };
  function drawLabels() {
    const L = M.labels, rules = ruleOrder(ST.lab.swap), [value, n] = applyRules(rules, labelRecord()), at = rules.findIndex(r => r.n === n);
    $('mlRules').innerHTML = rules.map((r, i) => {
      const word = labelRecord()[FIELD_KEY[r.from]], state = at < 0 ? 'none' : i < at ? 'pass' : i === at ? 'stop' : 'after';
      return `<div class="mrule ${state}"><span class="rn">${ordinal(r.n)} rule</span><span class="rf">${esc(r.from)}: ${code(word)}</span><span class="rw">${ST.labWords(r)}</span>${r.study_states ? `<span class="rq">quotes the study: “${esc(r.study_states)}”</span>` : ''}</div>`;
    }).join('');
    const records = L.records.map(r => ({ title: r.title, abdominal: r.abdominal, endotype: r.endotype, pneumonia: r.pneumonia, stored: r.value }));
    const now = records.map(r => applyRules(rules, r)[0]), changed = now.filter((v, i) => v !== records[i].stored).length;
    const sepsis = now.filter(v => v === 'sepsis').length;
    $('mlOut').innerHTML = `This record is <b>${esc(VALUE_LABEL(value))}</b>${n ? `, decided by the ${ordinal(n)} rule` : ', since no rule decides'}. ${ST.lab.swap ? `With the pneumonia rule before the endotype rule, <b>${fmt(changed)}</b> of the ${fmt(records.length)} records of ${L.study} change class, and ${fmt(sepsis)} are sepsis instead of ${fmt(records.filter(r => r.stored === 'sepsis').length)}.` : `In the order of the rule set, the rules give all ${fmt(records.length)} records of ${L.study} the class they have in the sample table.`}`;
    // the records decided by each rule, from the real records
    const tally = {}; L.records.forEach(r => { const k = `${r.rule}|${r.value}`; tally[k] = (tally[k] || 0) + 1; });
    $('mlTable').innerHTML = `<table class="t"><thead><tr><th>Rule</th><th>Field</th><th>Value</th><th class="r">Records of ${L.study} decided</th></tr></thead><tbody>${Object.entries(tally).sort((a, b) => a[0].localeCompare(b[0], 'en', { numeric: true })).map(([k, c]) => { const [r, v] = k.split('|'); return `<tr><td>${+r ? ordinal(+r) : '–'}</td><td>${+r ? esc(ST.labRules[+r - 1].from) : ''}</td><td>${esc(VALUE_LABEL(v))}</td><td class="r">${fmt(c)}</td></tr>`; }).join('')}<tr class="total"><td colspan="3">All records</td><td class="r">${fmt(L.records.length)}</td></tr></tbody></table>`;
    const B = L.basis, parts = [['sample_record', 'The sample record'], ['sample_record_and_study_description', 'The record and a study sentence'], ['study_description', 'The study description'], ['not_recorded', 'Not decided']].filter(([k]) => B[k]);
    const w = wOf('mlBasis'), total = sum(parts.map(([k]) => B[k])), g = chart('mlBasis', w, 40, 'Stacked bar of all records by what their class rests on.'), X = lin(0, total, 0, w), cols = [C.accentInk, C.accent, tint(C.accent, 0.55), C.lineStrong];
    let x0 = 0;
    parts.forEach(([k, label], i) => { hover(rrect(g, X(x0), 8, Math.max(1, X(x0 + B[k]) - X(x0) - 1), 26, cols[i], 2), `<b>${label}</b><br>${fmt(B[k])} records (${pct(B[k] / total, 1)})`); x0 += B[k]; });
    $('mlBasis').insertAdjacentHTML('beforeend', `<div class="keys">${keys(parts.map(([k, label], i) => [`${label} ${fmt(B[k])}`, cols[i]]))}</div>`);
  }

  // ======================================================================= 6. repeat deposits and patients
  // the sample marked for a patient: without a flag, earliest in time, of whole blood, lowest accession
  function markOne(samples) {
    const key = s => [s.flag ? 1 : 0, s.time_order == null ? Infinity : s.time_order, s.tissue === 'whole_blood' ? 0 : 1, s.record];
    return samples.slice().sort((a, b) => { const A = key(a), B = key(b); for (let i = 0; i < A.length; i++) if (A[i] !== B[i]) return A[i] < B[i] ? -1 : 1; return 0; })[0].record;
  }
  function buildRepeats() {
    const R = M.repeats, P = R.patient, T0 = R.totals;
    ST.ns = P.samples.map(s => Object.assign({}, s));
    $('mbRepeats').innerHTML = `
      ${H3('eq', 'Equation')}
      ${EQ([`same sample: &nbsp; (chip type, scan time, median raw intensity)${sub('a')} = (…)${sub('b')}, or (slide, medians), or (barcode, position), or the digest of the counts`,
        `primary record = the record of the study released first (on the same day: the lower accession); &nbsp; sample identifier = its study : accession`,
        `one per patient = the first of the patient's samples by: without a flag › earliest time order › whole blood › lowest study and accession`],
        'The raw files hold facts about the physical array or library that do not change when it is deposited again: an Affymetrix CEL file records the chip type and the moment of the scan to the second, Agilent and CodeLink files the slide, Illumina exports the barcode of the chip and the position on it; for RNA-seq a short fingerprint of the whole table of counts is computed. Study platforms that share a sample form a cohort group; cohort groups that share an identified patient form one patient-safe group.')}
      ${H3('sim', 'Simulation', `One sample per patient: ${P.patient} of ${P.study}`)}
      <p class="note mtop0">Patient ${esc(P.patient)} was sampled on ${spell(P.samples.length)} days. Mark a sample as flagged, and see which one the preference order marks.</p>
      <div class="mdays" id="mrDays"></div>
      <div class="readout" id="mrDaysOut" aria-live="polite"></div>
      ${H3('ex', 'Worked example', `one array in ${spell(R.deposits.length)} studies: ${R.sample_id}`)}
      <div class="readout" id="mrEx"></div>
      <div class="tablewrap" id="mrTable"></div>
      <div class="chart" id="mrTime"></div>
      ${how('One mark per deposit of the array, at the release date of its study; every file gives the same chip type, scan time and median intensity, so the six records are one array. The first released is its primary record.')}
      <div class="readout" id="mrEx2"></div>
      <div class="mgroups" id="mrSafe"></div>
      ${H3('res', 'Result')}
      <div class="chart" id="mrEvidence"></div>
      ${how('The linked pairs of records by the evidence that links them, over the whole resource. Pairs linked by expression identity were near-identical profiles in studies known to share samples.')}
      <div id="mrRes"></div>`;
    $('mrDays').addEventListener('click', e => { const b = e.target.closest('button[data-i]'); if (!b) return; const s = ST.ns[+b.dataset.i]; s.flag = !s.flag; drawDays(); });
    const d = R.deposits, first = d.find(x => x.primary), same = d.every(x => x.scan_time === first.scan_time && x.median === first.median && x.chip === first.chip);
    const gse66099 = d.find(x => x.study === 'GSE66099'), healthy = d.filter(x => x.class === 'healthy_control').length;
    const spellings = [...new Set(d.map(x => x.patient_id))].sort((a, b) => b.includes('-') - a.includes('-'));
    $('mrEx').innerHTML = `The Affymetrix array ${first.record} of ${first.study} holds the blood of a child recruited as a control, patient ${esc(spellings[0])}. Its CEL file, of chip type ${code(first.chip)}, was scanned on ${dayOf(first.scan_time)} at ${TIME(first.scan_time)} and has a median raw intensity of ${fmt(first.median)}. CEL files in ${spell(d.length - 1)} later series give ${same ? 'the same chip type, scan time and median' : 'other values'}, so the ${spell(d.length)} records are one array${gse66099 ? `; the file in GSE66099 even carries the first accession in its name, ${code(gse66099.file)}` : ''}. ${first.study} was released first, on ${dayOf(first.released)}, so its record is the primary record and the sample identifier is <b>${esc(R.sample_id)}</b>.
      <br><br>The records do not agree on the class. ${first.study} calls the child a control and describes its controls as “${esc(first.statement)}”; since nothing there calls them healthy, its rule makes the record a non-septic control. The later deposits call the array normal, or quote a sentence about healthy controls, so their ${spell(healthy)} records are healthy controls. Each record keeps its own label and the disagreement is marked (${fmt(T0.class_differs)} distinct samples have records that differ in class). The patient is written ${andList(spellings.map(s => code(s)))} in different studies; because the records are one sample, the spellings are one patient.`;
    $('mrTable').innerHTML = `<table class="t"><thead><tr><th>Study</th><th>Released</th><th>Record</th><th>CEL file</th><th>Recorded wording</th><th>Curated class</th></tr></thead><tbody>${d.map(x => `<tr class="${x.primary ? 'sel' : ''}"><td>${studyLink(x.study)}${x.primary ? ' <span class="chip">primary</span>' : ''}</td><td>${fdate(x.released)}</td><td>${esc(x.record)}</td><td><code>${esc(x.file)}</code></td><td><code>${esc(x.wording)}</code></td><td>${clsDot(x.class)}${esc(labLower('class', x.class))}</td></tr>`).join('')}</tbody></table>`;
    const gr = R.group_rows, gRec = sum(gr.map(x => x.records)), gDist = sum(gr.map(x => x.primary));
    $('mrEx2').innerHTML = `<h4 class="exhead">The patients of ${P.study} <span class="m">and one patient in two cohort groups</span></h4>The array belongs to the cohort group ${R.group}: ${spell(gr.length)} study platforms (${andList(gr.map(x => x.study))}) with ${fmt(gRec)} records but only ${fmt(gDist)} distinct samples, ${fx(gRec / gDist, 1)} records per sample.
      ${P.study} holds ${fmt(P.records)} whole-blood arrays of ${fmt(P.patients)} patients, ${fmt(P.sepsis)} with sepsis and ${fmt(P.healthy)} healthy controls; the rules read the patient from the field ${code('group id')} and the day from ${code('group day')}. Patient ${esc(P.patient)}, who did not survive, was sampled from ${code(P.samples[0].time_point)} to ${code(P.samples[P.samples.length - 1].time_point)}; none of the ${spell(P.samples.length)} arrays is flagged, so the earliest, ${P.samples[0].record}, is marked. The healthy controls were sampled as ${andList(P.healthy_points.map(code))}. In all, ${fmt(P.with_most)} patients have ${spell(P.most)} samples each, and the ${fmt(P.records)} records come down to ${fmt(P.patients)} marked samples, one per patient.
      <br><br>A cohort group holds every deposit of a sample, but not always every sample of a patient. GSE13015 measured patient I001 on two versions of the Illumina BeadChip, and each array was deposited again in another study, so the patient's two arrays fall into two cohort groups. The records of the later deposits name no patient, but every record of a sample receives the patient that any of its records names. The two groups become one patient-safe group, ${esc([...new Set(R.safe.rows.map(r => r.group))].sort().join('+'))}; ${fmt(R.safe.patients)} patients join them in this way.`;
    const byG = {}; R.safe.rows.forEach(r => (byG[r.group] = byG[r.group] || []).push(r));
    $('mrSafe').innerHTML = Object.entries(byG).map(([g, rows]) => `<div class="mgroup"><b>Cohort group ${esc(g)}</b>${rows.map(r => `<div class="mrec ${r.primary ? 'p' : ''}">${esc(r.study)} · ${esc(r.platform)} · ${esc(r.record)}<span>${r.patient_id ? `patient ${esc(r.patient_id)}` : 'names no patient'}${r.primary ? ' · primary record' : ` · deposit of ${esc(r.sample_id)}`}</span></div>`).join('')}</div>`).join('<div class="mjoin" aria-hidden="true">+</div>') + `<p class="note mfull">Patient I001 joins the two groups into the patient-safe group ${Object.keys(byG).join('+')}.</p>`;
    const S0 = T0.samples, co = T0.cohorts;
    $('mrRes').innerHTML = RES(`The ${fmt(D.numbers.records)} records are <b>${fmt(S0.distinct)}</b> distinct samples: ${fmt(S0.deposited_more_than_once)} samples were deposited more than once, accounting for ${fmt(S0.repeated_records)} repeat records, and no sample in more than ${spell(S0.most_records_of_one_sample)} records. The ${fmt(co.study_platforms)} study platforms form <b>${fmt(co.cohort_groups)}</b> cohort groups, and patients join some of them into <b>${fmt(co.patient_safe_groups)}</b> patient-safe groups. ${fmt(T0.identified)} patients are identified, with ${fmt(T0.with_patient)} distinct samples between them, ${fmt(T0.in_several_studies)} of them in more than one study; another ${fmt(T0.without_patient)} distinct samples name no patient and count as patients of their own. One sample is marked for each: <b>${fmt(T0.one_per_patient)}</b>, an upper bound on the number of patients.`);
  }
  function drawDays() {
    const marked = markOne(ST.ns);
    $('mrDays').innerHTML = ST.ns.map((s, i) => `<div class="mday ${s.record === marked ? 'on' : ''}"><b>${esc(s.time_point)}</b><span>${esc(s.record)}</span><span class="small">day ${fmt(s.time_order)}${s.record === marked ? ' · marked' : ''}</span>
      <button type="button" class="btn ghost" data-i="${i}" data-k="flag" aria-pressed="${s.flag}">${s.flag ? 'Flagged' : 'No flag'}</button></div>`).join('');
    const m = ST.ns.find(s => s.record === marked), changed = ST.ns.some((s, i) => s.flag !== M.repeats.patient.samples[i].flag);
    $('mrDaysOut').innerHTML = `Marked: <b>${esc(m.record)}</b> (${esc(m.time_point)}). ${changed ? (m.flag ? 'Every sample is now flagged, so the earliest of them is marked.' : 'A flagged sample loses its place to the earliest sample without a flag.') : `As in the sample table: none of the ${spell(ST.ns.length)} is flagged, so the earliest is marked.`}`;
  }
  function drawRepeatCharts() {
    const d = M.repeats.deposits, w = wOf('mrTime'), h = 120, m = { l: 20, r: 20, t: 40, b: 30 };
    const yr = s => { const t = new Date(s); return t.getUTCFullYear() + (t - Date.UTC(t.getUTCFullYear(), 0, 1)) / 31557600000; };
    const scan = yr(d[0].scan_time), lo = Math.floor(Math.min(scan, ...d.map(x => yr(x.released)))), hi = Math.ceil(Math.max(...d.map(x => yr(x.released))));
    const g = chart('mrTime', w, h, 'Timeline of the scan and the six deposits of one array.'), X = lin(lo, hi, m.l, w - m.r);
    gridX(g, X, ticks(lo, hi, small(w) ? 4 : 8), m.t - 10, h - m.b, v => String(v));
    S('line', { x1: X(lo), x2: X(hi), y1: h - m.b, y2: h - m.b, stroke: C.lineStrong }, g);
    S('path', { d: `M${X(scan)},${h - m.b - 4} l5,-9 l-10,0z`, fill: C.ink }, g);
    T(g, X(scan), m.t - 22, 'scanned', { 'text-anchor': 'middle', 'font-size': 11, fill: C.muted });
    d.forEach((x, i) => {
      const cx = X(yr(x.released)), cy = h - m.b - 16 - (i % 2) * 14;
      hover(S('circle', { cx, cy, r: x.primary ? 7 : 5, fill: C.cls[x.class], stroke: C.bg, 'stroke-width': 1.5 }, g), `<b>${esc(x.study)}</b> ${esc(x.record)}<br>released ${fdate(x.released)}<br>${esc(labLower('class', x.class))}${x.primary ? ' · primary record' : ''}`);
      if (!small(w) || x.primary) T(g, cx, cy - 10, x.study, { 'text-anchor': 'middle', 'font-size': 10.5, fill: C.ink2 });
    });
    const ev = M.repeats.evidence, LAB = { chip_type_scan_time_and_median_intensity: 'Chip type, scan time and median', beadchip_barcode_and_position: 'BeadChip barcode and position', slide_identifier_and_median_intensities: 'Slide and median intensities', identical_counts: 'Identical counts', expression_identity: 'Expression identity', record_names_a_re_analysis: 'A record names a re-analysis', 'expression_identity; identical_counts': 'Expression identity and identical counts' };
    const w2 = wOf('mrEvidence'), left = small(w2) ? 170 : 260, rowH = 25, H2 = rowH * ev.length + 34, g2 = chart('mrEvidence', w2, H2, 'Bar chart of the linked pairs of records by evidence.'), xm = niceMax(Math.max(...ev.map(e => e[1]))), X2 = lin(0, xm, left, w2 - 50);
    gridX(g2, X2, ticks(0, xm, 4), 4, rowH * ev.length + 4); axisTitle(g2, (left + w2 - 50) / 2, H2 - 2, 'Linked pairs of records');
    ev.forEach(([k, n], i) => { const y = 6 + i * rowH; T(g2, left - 8, y + 14, LAB[k] || k.replace(/_/g, ' '), { 'text-anchor': 'end', 'font-size': 12, fill: C.ink2 }); rrect(g2, X2(0), y + 3, X2(n) - X2(0), rowH - 8, C.fam.affymetrix, 3); T(g2, X2(n) + 5, y + 14, fmt(n), { 'font-size': 11.5, fill: C.muted }); hitRow(g2, 0, y, w2, rowH, `<b>${esc(LAB[k] || k)}</b><br>${fmt(n)} pairs`); });
  }

  // ======================================================================= 7. sex from expression
  // recorded sex: magenta for women and amber for men, colors that no class or technology of the page uses
  const sexColors = () => C.dark ? { f: '#d77fbd', m: '#e3a84a', n: C.muted } : { f: '#b04f95', m: '#c98516', n: C.muted };
  // two-means in one dimension: centers at the lowest and the highest value, each sample to the nearer center, centers to
  // the group means, until no sample changes group; returns the steps for the animation and the final split
  function twoMeans(d, start) {
    let lo = start ? start[0] : Math.min(...d), hi = start ? start[1] : Math.max(...d), assign = null;
    const steps = [];
    for (let it = 0; it < 100; it++) {
      const next = d.map(v => Math.abs(v - lo) <= Math.abs(v - hi) ? 0 : 1);
      steps.push({ lo, hi, cut: (lo + hi) / 2 });
      if (assign && next.every((a, i) => a === assign[i])) break;
      assign = next;
      lo = mean(d.filter((_, i) => !assign[i])); hi = mean(d.filter((_, i) => assign[i]));
    }
    return { steps, assign, lo, hi };
  }
  // the groups on either side of a split t, and what the checks of the stage ask of them
  function splitStats(d, y, x, t, zone) {
    const L = [], H = []; d.forEach((v, i) => (v < t ? L : H).push(i));
    if (!L.length || !H.length) return null;
    const mL = mean(L.map(i => d[i])), mH = mean(H.map(i => d[i]));
    const ss = (idx, mu, arr) => sum(idx.map(i => (arr[i] - mu) ** 2));
    const pooled = a => Math.sqrt((ss(L, mean(L.map(i => a[i])), a) + ss(H, mean(H.map(i => a[i])), a)) / (L.length + H.length - 2));
    const sp = pooled(d), gap = mH - mL, mid = (mL + mH) / 2;
    const ySep = (mean(H.map(i => y[i])) - mean(L.map(i => y[i]))) / pooled(y), xSep = (mean(L.map(i => x[i])) - mean(H.map(i => x[i]))) / pooled(x);
    const call = d.map(v => Math.abs(v - mid) < zone * gap ? 'a' : v >= mid ? 'm' : 'f');
    return { L, H, mL, mH, sp, gap, mid, sep: gap / sp, ySep, xSep, call };
  }
  function buildSex() {
    const X = M.sex, t = X.table, p = X.policy, e = Object.fromEntries(X.each.map(r => [r.accession, r]));
    const run = twoMeans(X.each.map(r => r.d));
    ST.sex = { t: run.steps[run.steps.length - 1].cut, anim: null };
    const [mA, fA] = X.examples.map(a => e[a]);
    $('mbSex').innerHTML = `
      ${H3('eq', 'Equation')}
      ${EQ([`${V('d')}${sub('i')} = ${frac(`Σ${sub(`g ∈ Y`)} ${V('x')}${sub('gi')}`, `|${V('Y')}|`)} − ${V('x')}${sub(`XIST, i`)}`,
        `${V('s')}${sub('p')} = √ ${frac(`Σ${sub(`i∈L`)}(${V('d')}${sub('i')} − ${V('d̄')}${sub('L')})² + Σ${sub(`i∈H`)}(${V('d')}${sub('i')} − ${V('d̄')}${sub('H')})²`, `${V('n')}${sub('L')} + ${V('n')}${sub('H')} − 2`)}, &nbsp; separation = ${frac(`${V('d̄')}${sub('H')} − ${V('d̄')}${sub('L')}`, `${V('s')}${sub('p')}`)}`],
        `x<sub>gi</sub>: the log₂ expression of gene g in sample i. Y: the genes of the Y chromosome that the table holds among the ${spell(X.y_genes.length)} ${andList(X.y_genes)}, which only men express; XIST is expressed strongly by women and hardly by men. L and H: the lower and higher group of a two-means split of d, made without the recorded sex. A split stands only where the separation is at least ${fmt(p.minimum_separation)}, the smaller group holds at least ${fmt(p.minimum_group)} samples and ${pct(p.minimum_group_share, 0)} of the table, and the higher group has the Y genes higher and XIST lower, each by at least ${fmt(p.component_separation)} pooled standard deviations. A sample closer to the midpoint than ${num(p.ambiguous_zone)} of the gap is called ambiguous.`)}
      ${H3('sim', 'Simulation', `Move the split on the samples of ${X.study}`)}
      <div class="controls"><button type="button" class="tbtn play" id="msxPlay" aria-label="Run two-means from the two extremes"></button><span>Run two-means from the lowest and highest value, or drag the split</span>
        <label>split ${V('d')} = <b id="msxTv"></b> <input type="range" id="msxT" min="-12" max="8" step="0.01" aria-label="The split of d"></label></div>
      <div class="grid2 msx"><div class="chart" id="msxPlane"></div><div class="chart" id="msxStrip"></div></div>
      <div class="keys" id="msxKeys"></div>
      <div class="readout" id="msxOut" aria-live="polite"></div>
      ${how(`The first chart places every sample of ${X.study} by its XIST value and the mean of its Y genes, colored by its recorded sex, which takes no part in the split; a split of d is a diagonal line on this plane. The second chart lays the samples out along d, with the split (the line you can drag) and the two group means. The readout recomputes the groups, their separation and the agreement with the recorded sex for the current split.`)}
      ${H3('ex', 'Worked example', `sex in ${X.study}`)}
      <div class="readout" id="msxEx"></div>
      <div class="tablewrap" id="msxTable"></div>
      ${H3('res', 'Result')}
      <div id="msxRes"></div>`;
    const ex = splitStats(X.each.map(r => r.d), X.each.map(r => r.y), X.each.map(r => r.xist), run.steps[run.steps.length - 1].cut, p.ambiguous_zone);
    ST.sexRun = run; ST.sexFinal = ex;
    $('msxEx').innerHTML = `${X.study} holds ${fmt(t.samples)} arrays, all with a recorded sex, and its table has all ${spell(X.y_genes.length)} Y genes and XIST. For ${mA.accession}, recorded male, ${V('d')} = ${fx(mA.y, 2)} − ${fx(mA.xist, 2)} = <b>${fx(mA.d, 1)}</b>; for ${fA.accession}, the sepsis record of the chapter on labels, recorded female, ${V('d')} = ${fx(fA.y, 2)} − ${fx(fA.xist, 2)} = <b>${fx(fA.d, 1)}</b>.
      Two-means over all ${fmt(t.samples)} values stops after ${spell(run.steps.length - 1)} ${run.steps.length === 2 ? 'round' : 'rounds'} with ${fmt(ex.L.length)} samples in the lower group and ${fmt(ex.H.length)} in the higher, with means ${fx(ex.mL, 2)} and ${fx(ex.mH, 2)}. The gap is ${fx(ex.mH, 2)} − (${fx(ex.mL, 2)}) = ${fx(ex.gap, 2)} and the pooled standard deviation ${fx(ex.sp, 3)}, so the separation is ${fx(ex.gap, 2)}/${fx(ex.sp, 3)} = <b>${fx(ex.sep, 1)}</b>, far above ${fmt(p.minimum_separation)}; the Y genes separate the groups by ${fx(ex.ySep, 1)} and XIST by ${fx(ex.xSep, 1)} pooled standard deviations. The midpoint is ${fx(ex.mid, 2)}, and the ambiguous zone reaches ${fx(p.ambiguous_zone * ex.gap, 2)} on either side, from ${fx(ex.mid - p.ambiguous_zone * ex.gap, 1)} to ${fx(ex.mid + p.ambiguous_zone * ex.gap, 1)}. So ${fmt(ex.call.filter(c => c === 'f').length)} samples are called female, ${fmt(ex.call.filter(c => c === 'm').length)} male and ${fmt(ex.call.filter(c => c === 'a').length)} ambiguous.`;
    const rec = X.each.map(r => r.recorded), cnt = (r, c) => X.each.filter((x, i) => rec[i] === r && ex.call[i] === c).length;
    const agree = cnt('f', 'f') + cnt('m', 'm'), dis = cnt('f', 'm') + cnt('m', 'f');
    $('msxTable').innerHTML = `<table class="t mshort"><thead><tr><th>Recorded sex</th><th class="r">Read female</th><th class="r">Read male</th><th class="r">Ambiguous</th></tr></thead><tbody>${[['f', 'Female'], ['m', 'Male']].map(([k, l]) => `<tr><td>${l}</td><td class="r">${fmt(cnt(k, 'f'))}</td><td class="r">${fmt(cnt(k, 'm'))}</td><td class="r">${fmt(cnt(k, 'a'))}</td></tr>`).join('')}</tbody></table>
      <p class="note">Of the ${fmt(agree + dis)} called samples, <b>${fmt(agree)}</b> agree with the recorded sex and ${fmt(dis)} do not: an agreement of ${pct(agree / (agree + dis), 1)}. A single disagreement can be a mistake in the record, a sample swapped in the laboratory, or biology, such as the loss of the Y chromosome from the blood cells of some older men; the check lists every one for inspection.</p>`;
    const R = X.resource, sm = R.samples;
    $('msxRes').innerHTML = RES(`Of the ${fmt(R.matrices)} gene-level tables, <b>${fmt(R.two_sexes)}</b> showed both sexes. Over the ${fmt(sm.compared)} called samples that have a recorded sex, the two agreed for <b>${pct(sm.concordance, 1)}</b>. Where at least ${fmt(X.policy.minimum_compared)} samples of a table are compared and fewer than ${pct(X.policy.unreliable_below, 0)} agree, the recorded sex of the whole table is marked unreliable: ${andList(R.unreliable.map(u => esc(u.split('|')[0])))}. Outside them the agreement was <b>${pct(sm.concordance_outside_unreliable_matrices, 1)}</b>.`);
    $('msxT').addEventListener('input', ev => { stopSex(); ST.sex.t = +ev.target.value; drawSex(); });
    dragAlong('msxStrip', ST.sex, v => { ST.sex.t = v; drawSex(); }, stopSex);
    $('msxPlay').onclick = () => { if (ST.sex.anim) stopSex(); else playSex(); };
  }
  function stopSex() { if (ST.sex.anim) clearInterval(ST.sex.anim); ST.sex.anim = null; if ($('msxPlay')) $('msxPlay').innerHTML = ICON_PLAY; }
  function playSex() {
    const steps = ST.sexRun.steps; let k = 0;
    $('msxPlay').innerHTML = ICON_PAUSE;
    const step = () => { ST.sex.t = steps[k].cut; ST.sex.centers = [steps[k].lo, steps[k].hi]; drawSex(); k++; if (k >= steps.length) { stopSex(); ST.sex.centers = null; drawSex(); } };
    step(); if (!reducedMotion()) ST.sex.anim = setInterval(step, 700); else { ST.sex.t = steps[steps.length - 1].cut; ST.sex.centers = null; drawSex(); }
  }
  function drawSex() {
    const X = M.sex, rows = X.each, d = rows.map(r => r.d), y = rows.map(r => r.y), x = rows.map(r => r.xist), t = ST.sex.t, p = X.policy;
    $('msxTv').textContent = fx(t, 2); $('msxT').value = t;
    const st = splitStats(d, y, x, t, p.ambiguous_zone), colR = sexColors();
    // the plane
    const w = wOf('msxPlane'), h = Math.min(320, Math.max(240, w * 0.8)), m = { l: 44, r: 10, t: 10, b: 38 };
    const g = chart('msxPlane', w, h, `Scatter plot of XIST against the mean of the Y genes for the samples of ${X.study}, with the split line.`);
    const xs = lin(Math.floor(Math.min(...x)), Math.ceil(Math.max(...x)), m.l, w - m.r), ys = lin(Math.floor(Math.min(...y)), Math.ceil(Math.max(...y)), h - m.b, m.t);
    gridX(g, xs, ticks(xs.inv(m.l), xs.inv(w - m.r), 5), m.t, h - m.b); gridY(g, ys, ticks(ys.inv(h - m.b), ys.inv(m.t), 5), m.l, w - m.r);
    axisTitle(g, (m.l + w - m.r) / 2, h - 6, 'XIST, log₂'); axisTitle(g, 12, (m.t + h - m.b) / 2, 'Mean of the Y genes, log₂', true);
    const clip = `msxc${themeEpoch}`; const defs = S('defs', {}, g); S('rect', { x: m.l, y: m.t, width: w - m.l - m.r, height: h - m.t - m.b }, S('clipPath', { id: clip }, defs));
    const gg = S('g', { 'clip-path': `url(#${clip})` }, g);
    // the split d = t is the line y = x + t
    const x0 = xs.inv(m.l), x1 = xs.inv(w - m.r);
    S('path', { d: `M${xs(x0)},${ys(x0 + t)} L${xs(x1)},${ys(x1 + t)}`, stroke: C.ink, 'stroke-width': 1.8, 'stroke-dasharray': '6 4' }, gg);
    rows.forEach((r, i) => hover(S('circle', { cx: xs(r.xist), cy: ys(r.y), r: 3, fill: colR[r.recorded], 'fill-opacity': 0.7, stroke: st && st.call[i] === 'a' ? C.ink : 'none' }, gg), `<b>${esc(r.accession)}</b><br>Y genes ${fx(r.y, 2)}, XIST ${fx(r.xist, 2)}, d = ${fx(r.d, 2)}<br>recorded ${r.recorded === 'f' ? 'female' : r.recorded === 'm' ? 'male' : 'not recorded'}${st ? ` · read ${st.call[i] === 'f' ? 'female' : st.call[i] === 'm' ? 'male' : 'ambiguous'}` : ''}`));
    // the strip along d
    const w2 = wOf('msxStrip'), h2 = h, mm = { l: 14, r: 14, t: 16, b: 38 };
    const g2 = chart('msxStrip', w2, h2, 'Dot strip of the score d, with the split you can drag.');
    const dx = lin(-12, 8, mm.l, w2 - mm.r), r0 = rng(65682), band = (h2 - mm.t - mm.b) / 2;
    gridX(g2, dx, ticks(-12, 8, 5), mm.t, h2 - mm.b); axisTitle(g2, (mm.l + w2 - mm.r) / 2, h2 - 6, 'd = mean of the Y genes − XIST');
    rows.forEach((r, i) => S('circle', { cx: dx(r.d), cy: mm.t + band * (r.recorded === 'm' ? 1.5 : 0.5) + (r0() - 0.5) * band * 0.8, r: 2.6, fill: colR[r.recorded], 'fill-opacity': 0.7 }, g2));
    T(g2, mm.l + 2, mm.t + 10, 'recorded female', { 'font-size': 11, fill: C.muted }); T(g2, mm.l + 2, mm.t + band + 10, 'recorded male', { 'font-size': 11, fill: C.muted });
    if (st) {
      rrect(g2, dx(st.mid - p.ambiguous_zone * st.gap), mm.t, dx(st.mid + p.ambiguous_zone * st.gap) - dx(st.mid - p.ambiguous_zone * st.gap), h2 - mm.t - mm.b, rgba(C.ink, 0.06), 0);
      [st.mL, st.mH].forEach(v => S('line', { x1: dx(v), x2: dx(v), y1: mm.t, y2: h2 - mm.b, stroke: C.accentInk, 'stroke-width': 1.5 }, g2));
      T(g2, dx(st.mL), h2 - mm.b - 4, 'mean L', { 'text-anchor': 'middle', 'font-size': 10.5, fill: C.accentInk }); T(g2, dx(st.mH), h2 - mm.b - 4, 'mean H', { 'text-anchor': 'middle', 'font-size': 10.5, fill: C.accentInk });
    }
    if (ST.sex.centers) ST.sex.centers.forEach(v => S('path', { d: `M${dx(v)},${mm.t - 2} l-5,-8 l10,0z`, fill: C.cls.sepsis }, g2));
    const handle = S('g', { class: 'mdrag', 'aria-hidden': 'true' }, g2);
    S('line', { x1: dx(t), x2: dx(t), y1: mm.t - 4, y2: h2 - mm.b, stroke: C.ink, 'stroke-width': 2.2 }, handle);
    S('circle', { cx: dx(t), cy: mm.t - 4, r: 6, fill: C.ink }, handle);
    ST.sex.drag = { scale: dx, width: w2, lo: -12, hi: 8 };
    $('msxKeys').innerHTML = keys([['recorded female', colR.f, 'dot'], ['recorded male', colR.m, 'dot'], ['ambiguous call', C.ink, 'ring'], ['split', C.ink, 'line']]);
    if (!st) { $('msxOut').innerHTML = 'Every sample lies on one side of the split: there are no two groups.'; return; }
    const rec = rows.map(r => r.recorded), called = st.call.map((c, i) => c !== 'a' && rec[i] !== 'n' ? (c === rec[i]) : null).filter(v => v != null);
    const small_ = Math.min(st.L.length, st.H.length), stands = st.sep >= p.minimum_separation && small_ >= p.minimum_group && small_ >= p.minimum_group_share * d.length && st.ySep >= p.component_separation && st.xSep >= p.component_separation;
    $('msxOut').innerHTML = `Lower group <b>${fmt(st.L.length)}</b> (mean ${fx(st.mL, 2)}), higher group <b>${fmt(st.H.length)}</b> (mean ${fx(st.mH, 2)}); pooled standard deviation ${fx(st.sp, 3)}, separation <b>${fx(st.sep, 1)}</b>; Y genes ${fx(st.ySep, 1)} and XIST ${fx(st.xSep, 1)} pooled standard deviations apart. ${stands ? 'The split meets every condition.' : 'This split would <b>not</b> stand: it fails a condition, so no sample would be called.'} Called samples that agree with their recorded sex: <b>${pct(mean(called.map(Number)), 1)}</b> of ${fmt(called.length)}.`;
  }
  // dragging a vertical line along the horizontal axis of a chart: the listeners sit on the chart's box, which outlives the
  // SVG that every move redraws; state.drag = {scale, width, lo, hi} is set by the drawing, onValue(v) moves the line
  function dragAlong(id, state, onValue, before) {
    const box = $(id); let on = false;
    const set = ev => { const d = state.drag, s = box.querySelector('svg'); if (!d || !s) return; const r = s.getBoundingClientRect(); onValue(+clamp(d.scale.inv((ev.clientX - r.left) * d.width / r.width), d.lo, d.hi).toFixed(2)); };
    box.addEventListener('pointerdown', ev => { if (!state.drag || ev.button > 0) return; on = true; if (before) before(); box.setPointerCapture(ev.pointerId); set(ev); });
    box.addEventListener('pointermove', ev => { if (on) set(ev); });
    const end = () => { on = false; };
    box.addEventListener('pointerup', end); box.addEventListener('pointercancel', end);
  }

  // ======================================================================= 8. the sepsis score and the AUROC
  const geo = a => Math.exp(mean(a.map(Math.log)));
  function rocCurve(pos, neg) {
    const all = [...new Set([...pos, ...neg])].sort((a, b) => b - a), pts = [[0, 0]];
    all.forEach(t => pts.push([neg.filter(v => v >= t).length / neg.length, pos.filter(v => v >= t).length / pos.length]));
    return pts;
  }
  function buildScore() {
    const Sc = M.score, strata = Sc.strata, home = strata.find(s => s.study === Sc.study && s.platform === Sc.platform);
    ST.auc = { k: strata.indexOf(home), t: null };
    const label = s => `${s.study}${strata.filter(x => x.study === s.study).length > 1 ? ` · ${s.platform} · ${lab('tissue', s.tissue)}` : ''}`;
    $('mbScore').innerHTML = `
      ${H3('eq', 'Equation')}
      ${EQ([`${V('s')}${sub('i')} = ( ∏${sub(`g∈U`)} ${V('x')}′${sub('gi')} )${'<sup>1/|U|</sup>'} − ( ∏${sub(`g∈D`)} ${V('x')}′${sub('gi')} )${'<sup>1/|D|</sup>'}`,
        `AUROC = ${frac('1', `${V('n')}${sub('S')} ${V('n')}${sub('H')}`)} Σ${sub(`i∈S`)} Σ${sub(`j∈H`)} ( [${V('s')}${sub('i')} > ${V('s')}${sub('j')}] + ½ [${V('s')}${sub('i')} = ${V('s')}${sub('j')}] )`],
        `U and D: the ${spell(Sc.up.length)} genes reported higher in sepsis and the ${spell(Sc.down.length)} reported lower by Sweeney and colleagues in 2015 that the table holds (at least ${spell(Sc.minimum_genes)} of each). x′<sub>gi</sub>: the log₂ value of gene g in sample i, raised by one constant for the whole table where its lowest value is below one, so that every value is positive. Each term is a geometric mean. S: the sepsis samples drawn during sepsis; H: the healthy controls; n<sub>S</sub> and n<sub>H</sub> their numbers; [·] is 1 where the statement inside is true and 0 otherwise. 0.5 is chance; 1 means that every sepsis sample scores above every control. The direction of the score is fixed in advance.`)}
      ${H3('sim', 'Simulation', 'Drag the threshold and watch the ROC curve')}
      <div class="controls"><label>Study <select id="maucS">${strata.map((s, i) => `<option value="${i}">${esc(label(s))}</option>`).join('')}</select></label>
        <label>threshold = <b id="maucTv"></b> <input type="range" id="maucT" step="0.01" aria-label="Threshold of the score"></label></div>
      <div class="grid2 mauc"><div class="chart" id="maucDots"></div><div class="chart" id="maucRoc"></div></div>
      <div class="readout" id="maucOut" aria-live="polite"></div>
      ${how('The first chart shows the score of every sepsis sample and healthy control of the study; every sample at or above the threshold is called sepsis (filled). The second is the ROC curve: the share of sepsis samples called (sensitivity) against the share of controls called (false-positive rate), for every threshold; the dot is the current one. The area under the curve equals the share of sepsis-control pairs in which the sepsis sample scores higher.')}
      ${H3('ex', 'Worked example', `the score of ${Sc.sample}, and the AUROC of ${Sc.study} counted from its pairs`)}
      <div class="grid2"><div class="tablewrap" id="maucGenes"></div><div class="readout" id="maucEx"></div></div>
      ${H3('res', 'Result')}
      <div class="chart" id="maucAll"></div>
      ${how('One row per study stratum (study platform and sample type) with at least five sepsis samples and five healthy controls: its AUROC with its bootstrap interval, in which patients are drawn at random with replacement 2,000 times. Strata on which the score was derived or validated are marked. Click a row to open it in the simulation.'.replace('2,000', fmt(Sc.bootstrap.replicates)))}
      <div id="maucRes"></div>`;
    $('maucS').value = ST.auc.k;
    $('maucS').onchange = e => { ST.auc.k = +e.target.value; ST.auc.t = null; drawScore(); };
    $('maucT').addEventListener('input', e => { ST.auc.t = +e.target.value; drawScore(); });
    dragAlong('maucDots', ST.auc, v => { ST.auc.t = v; drawScore(); });
    const up = Sc.up.map(g => g.value), dn = Sc.down.map(g => g.value), gu = geo(up), gd = geo(dn);
    const rowsUp = Sc.up, rowsDn = Sc.down, n = Math.max(rowsUp.length, rowsDn.length);
    $('maucGenes').innerHTML = `<table class="t"><thead><tr><th>Up gene</th><th class="r">Value</th><th>Down gene</th><th class="r">Value</th></tr></thead><tbody>${Array.from({ length: n }, (_, i) => `<tr><td>${rowsUp[i] ? esc(rowsUp[i].symbol) + (rowsUp[i].published !== rowsUp[i].symbol ? ` <span class="ui">(${esc(rowsUp[i].published)})</span>` : '') : ''}</td><td class="r">${rowsUp[i] ? fx(rowsUp[i].value, 2) : ''}</td><td>${rowsDn[i] ? esc(rowsDn[i].symbol) + (rowsDn[i].published !== rowsDn[i].symbol ? ` <span class="ui">(${esc(rowsDn[i].published)})</span>` : '') : ''}</td><td class="r">${rowsDn[i] ? fx(rowsDn[i].value, 2) : ''}</td></tr>`).join('')}<tr class="total"><td>Geometric mean</td><td class="r">${fx(gu, 2)}</td><td>Geometric mean</td><td class="r">${fx(gd, 2)}</td></tr></tbody></table>`;
    const P = home.sepsis, Hh = home.healthy, pairs = P.length * Hh.length; let above = 0, ties = 0;
    P.forEach(a => Hh.forEach(b => { if (a > b) above++; else if (a === b) ties++; }));
    ST.aucEx = { gu, gd, score: gu - gd, above, ties, pairs };
    $('maucEx').innerHTML = `The lowest value of the ${Sc.study} table is not below one, so nothing is raised. For the sepsis record ${Sc.sample}, the up genes give (${up.map(v => fx(v, 2)).join(' × ')})<sup>1/${up.length}</sup> = ${fx(gu, 2)}, the down genes (${dn.map(v => fx(v, 2)).join(' × ')})<sup>1/${dn.length}</sup> = ${fx(gd, 2)}, and the score is ${fx(gu, 2)} − ${fx(gd, 2)} = <b>${fx(gu - gd, 2)}</b>. The old symbols C9orf95 and KIAA1370 of the publication are found under their current names.
      <br><br>In ${Sc.study}, the ${fmt(P.length)} sepsis samples and ${fmt(Hh.length)} healthy controls make ${fmt(P.length)} × ${fmt(Hh.length)} = ${fmt(pairs)} pairs. In ${fmt(above)} the sepsis sample has the higher score and in ${fmt(pairs - above - ties)} the control${ties ? `; ${fmt(ties)} are ties` : ''}, so the AUROC is ${fmt(above)}/${fmt(pairs)} = <b>${fx((above + ties / 2) / pairs, 3)}</b>. The mean score is ${fx(mean(P), 2)} for sepsis and ${fx(mean(Hh), 2)} for the controls. Drawing patients again gives the interval ${fx(home.low, 3)} to ${fx(home.high, 3)}.`;
    const R = Sc.resource, sc = R.score, smallest = strata.reduce((a, b) => b.auroc < a.auroc ? b : a);
    $('maucRes').innerHTML = RES(`<b>${fmt(R.strata)}</b> strata held at least ${spell(Sc.minimum_group)} sepsis samples drawn during sepsis and as many healthy controls. The AUROC of the score ranged from <b>${fx(sc.auroc_lowest, 2)}</b> to <b>${fx(Math.max(...strata.map(s => s.auroc)), 2)}</b>, with a median of <b>${fx(sc.auroc_median, 3)}</b>; in the ${fmt(R.independent.strata)} strata independent of its derivation the median was ${fx(R.independent.auroc_median, 3)}. The lowest, ${smallest.study}, has ${fmt(smallest.sepsis.length)} sepsis samples and ${fmt(smallest.healthy.length)} controls, and its interval runs from ${fx(smallest.low, 2)} to ${fx(smallest.high, 2)}: small groups give wide intervals. Gene by gene, the difference between sepsis and healthy controls had its published direction in ${pct(R.genes.share_as_published, 1)} of ${fmt(R.genes.comparisons)} comparisons. Sepsis against healthy volunteers is the easiest comparison in this field: a high AUROC shows that the data and the labels carry the expected biology, not that the score tells sepsis from other illness.`);
  }
  function drawScore() {
    const Sc = M.score, s = Sc.strata[ST.auc.k], pos = s.sepsis, neg = s.healthy, all = [...pos, ...neg], lo = Math.min(...all), hi = Math.max(...all);
    if (ST.auc.t == null) ST.auc.t = +((mean(pos) + mean(neg)) / 2).toFixed(2);
    const tIn = $('maucT'); tIn.min = (lo - 0.05).toFixed(2); tIn.max = (hi + 0.05).toFixed(2); tIn.value = ST.auc.t;
    const t = ST.auc.t; $('maucTv').textContent = fx(t, 2);
    const w = wOf('maucDots'), h = 260, m = { l: 14, r: 14, t: 24, b: 40 }, g = chart('maucDots', w, h, `Dot strips of the scores of sepsis samples and healthy controls of ${s.study}, with the threshold.`);
    const X = lin(lo - 0.08 * (hi - lo), hi + 0.08 * (hi - lo), m.l, w - m.r), band = (h - m.t - m.b) / 2, r0 = rng(ST.auc.k + 7);
    gridX(g, X, ticks(X.inv(m.l), X.inv(w - m.r), 5), m.t, h - m.b, v => fx(v, 1));
    axisTitle(g, (m.l + w - m.r) / 2, h - 6, 'Sepsis MetaScore');
    [[pos, C.cls.sepsis, 'Sepsis'], [neg, C.cls.healthy_control, 'Healthy control']].forEach(([vals, col, name], k) => {
      T(g, m.l + 2, m.t + band * k + 11, `${name} (${fmt(vals.length)})`, { 'font-size': 11.5, fill: C.ink2 });
      vals.forEach(v => S('circle', { cx: X(v), cy: m.t + band * (k + 0.55) + (r0() - 0.5) * band * 0.7, r: vals.length > 300 ? 2.2 : 3, fill: v >= t ? col : 'none', stroke: col, 'stroke-width': 1, 'fill-opacity': 0.75 }, g));
    });
    const thr = S('g', { class: 'mdrag', 'aria-hidden': 'true' }, g);
    S('line', { x1: X(t), x2: X(t), y1: m.t - 6, y2: h - m.b, stroke: C.ink, 'stroke-width': 2 }, thr);
    S('circle', { cx: X(t), cy: m.t - 6, r: 6, fill: C.ink }, thr);
    T(thr, X(t) + (X(t) > w * 0.7 ? -9 : 9), m.t - 2, `threshold ${fx(t, 2)}`, { 'text-anchor': X(t) > w * 0.7 ? 'end' : 'start', 'font-size': 11, fill: C.ink });
    ST.auc.drag = { scale: X, width: w, lo: lo - 0.05, hi: hi + 0.05 };
    // ROC
    const w2 = wOf('maucRoc'), s2 = Math.min(w2, 300), m2 = { l: 46, r: 12, t: 12, b: 38 }, g2 = chart('maucRoc', w2, s2, `ROC curve of the score in ${s.study}, with the point of the current threshold.`);
    const side = Math.min(w2 - m2.l - m2.r, s2 - m2.t - m2.b), X2 = lin(0, 1, m2.l, m2.l + side), Y2 = lin(0, 1, m2.t + side, m2.t);
    gridX(g2, X2, [0, 0.25, 0.5, 0.75, 1], m2.t, m2.t + side, v => fx(v, 2)); gridY(g2, Y2, [0, 0.25, 0.5, 0.75, 1], m2.l, m2.l + side, v => fx(v, 2));
    axisTitle(g2, m2.l + side / 2, m2.t + side + 32, 'False-positive rate'); axisTitle(g2, 12, m2.t + side / 2, 'Sensitivity', true);
    S('line', { x1: X2(0), y1: Y2(0), x2: X2(1), y2: Y2(1), stroke: C.lineStrong, 'stroke-dasharray': '4 4' }, g2);
    const roc = rocCurve(pos, neg), area = auroc(pos, neg);
    S('path', { d: `M${X2(0)},${Y2(0)} ` + roc.map(([fp, tp]) => `L${X2(fp)},${Y2(tp)}`).join(' ') + ` L${X2(1)},${Y2(1)} L${X2(1)},${Y2(0)} Z`, fill: rgba(C.cls.sepsis, 0.12) }, g2);
    S('path', { d: roc.map(([fp, tp], i) => `${i ? 'L' : 'M'}${X2(fp)},${Y2(tp)}`).join(' ') + ` L${X2(1)},${Y2(1)}`, stroke: C.cls.sepsis, 'stroke-width': 2, fill: 'none' }, g2);
    const tp = pos.filter(v => v >= t).length / pos.length, fp = neg.filter(v => v >= t).length / neg.length;
    S('circle', { cx: X2(fp), cy: Y2(tp), r: 6, fill: C.ink, stroke: C.bg, 'stroke-width': 2 }, g2);
    T(g2, X2(0.55), Y2(0.18), `area ${fx(area, 3)}`, { 'font-size': 12.5, fill: C.ink, 'font-weight': 600 });
    $('maucOut').innerHTML = `At ${fx(t, 2)}: ${fmt(pos.filter(v => v >= t).length)} of ${fmt(pos.length)} sepsis samples are called sepsis (sensitivity <b>${pct(tp, 1)}</b>), and ${fmt(neg.filter(v => v >= t).length)} of ${fmt(neg.length)} controls (false-positive rate <b>${pct(fp, 1)}</b>, specificity ${pct(1 - fp, 1)}). The area under the curve does not depend on the threshold: <b>${fx(area, 3)}</b> (stored ${fx(s.auroc, 3)}; interval ${s.low != null ? `${fx(s.low, 3)} to ${fx(s.high, 3)}` : 'not computed, too few patients'}).`;
    // every stratum
    const st = Sc.strata.slice().sort((a, b) => b.auroc - a.auroc), w3 = wOf('maucAll'), rowH = 15, left = small(w3) ? 92 : 140, H3_ = rowH * st.length + 44;
    const g3 = chart('maucAll', w3, H3_, 'Dot chart of the AUROC of the score in every stratum, with its interval.');
    const lo3 = Math.min(0.5, Math.floor(Math.min(...st.map(x => x.low != null ? x.low : x.auroc)) * 10) / 10), X3 = lin(lo3, 1, left, w3 - 20);
    gridX(g3, X3, ticks(lo3, 1, 5), 6, rowH * st.length + 8, v => fx(v, 1));
    axisTitle(g3, (left + w3 - 20) / 2, H3_ - 4, 'AUROC of the Sepsis MetaScore, sepsis against healthy controls');
    st.forEach((x, i) => {
      const y = 12 + i * rowH, me = Sc.strata.indexOf(x) === ST.auc.k;
      T(g3, left - 8, y + 4, x.study, { 'text-anchor': 'end', 'font-size': 11, fill: me ? C.ink : C.ink2, 'font-weight': me ? 600 : 400 });
      if (x.low != null) S('line', { x1: X3(x.low), x2: X3(x.high), y1: y, y2: y, stroke: C.lineStrong, 'stroke-width': 2 }, g3);
      S('circle', { cx: X3(x.auroc), cy: y, r: me ? 5 : 3.5, fill: x.derivation ? C.bg : C.cls.sepsis, stroke: C.cls.sepsis, 'stroke-width': 1.5 }, g3);
      const hit = hitRow(g3, 0, y - rowH / 2, w3, rowH, `<b>${esc(x.study)}</b> <span class="m">${esc(x.platform)} · ${esc(lab('tissue', x.tissue))}</span><br>AUROC ${fx(x.auroc, 3)}${x.low != null ? ` (${fx(x.low, 2)} to ${fx(x.high, 2)})` : ''}<br>${fmt(x.sepsis.length)} sepsis, ${fmt(x.healthy.length)} healthy controls${x.derivation ? `<br>used in the ${esc(x.derivation)} of the score` : ''}`);
      hit.style.cursor = 'pointer'; hit.addEventListener('click', () => { ST.auc.k = Sc.strata.indexOf(x); ST.auc.t = null; $('maucS').value = ST.auc.k; drawScore(); });
    });
    $('maucAll').insertAdjacentHTML('beforeend', `<div class="keys">${keys([['independent of the score', C.cls.sepsis, 'dot'], ['used to derive or validate the score', C.cls.sepsis, 'ring'], ['bootstrap interval', C.lineStrong, 'line']])}</div>`);
  }

  // ======================================================================= 9. the correction
  // the shrunk variance and location of empirical Bayes (ComBat's estimators)
  const shrinkVar = (own, n, a, b) => (0.5 * (n - 1) * own + b) / (0.5 * (n - 1) + a - 1);
  const shrinkWeight = (tau2, v, n) => tau2 / (tau2 + v / n);
  // the reference distribution at a quantile position (evenly spaced positions, linear between them)
  function refAt(ref, pos) {
    const k = (pos - ref.first) / (ref.last - ref.first) * (ref.value.length - 1);
    if (k <= 0) return ref.value[0]; if (k >= ref.value.length - 1) return ref.value[ref.value.length - 1];
    const i = Math.floor(k); return ref.value[i] + (ref.value[i + 1] - ref.value[i]) * (k - i);
  }
  const lgamma = z => { const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]; let x = z, y = z, tmp = x + 5.5; tmp -= (x + 0.5) * Math.log(tmp); let ser = 1.000000000190015; for (const v of c) ser += v / ++y; return -tmp + Math.log(2.5066282746310005 * ser / x); };
  const invGammaPdf = (v, a, b) => v <= 0 ? 0 : Math.exp(a * Math.log(b) - lgamma(a) - (a + 1) * Math.log(v) - b / v);
  function buildCorrection() {
    const K = M.correction, A = K.anchored, mt = K.matched, ref = K.reference, [es, ec] = A.examples, lim = K.limits;
    ST.n = A.n;
    const z = x => (x.aligned - A.sex_effect * x.female - A.anchor_mean) / A.anchor_sd;
    const corr = x => A.anchor_mean + A.anchor_sd * (z(x) - A.location) / A.scale + A.sex_effect * x.female;
    $('mbCorrect').innerHTML = `
      ${H3('eq', 'Equation')}
      ${EQ([`1. alignment: &nbsp; ${V('p')} = ${frac(`${V('r')} − 0.5`, V('n'))}, &nbsp; ${V('x̃')} = ${V('Q')}${sub('ref')}(${V('p')})`,
        `2. anchor units: &nbsp; ${V('z')}${sub('gi')} = ${frac(`${V('x̃')}${sub('gi')} − ${V('β')}${sub('g')} ${V('f')}${sub('i')} − ${V('m')}${sub('g')}`, `${V('s')}${sub('g')}`)}`,
        `3. shrinkage: &nbsp; (${V('σ')}*)² = ${frac(`½(${V('n')} − 1) ${V('σ̂')}² + ${V('b')}`, `½(${V('n')} − 1) + ${V('a')} − 1`)}, &nbsp; ${V('δ')}* = ${V('w')} ${V('δ̂')} + (1 − ${V('w')}) ${V('δ̄')}, &nbsp; ${V('w')} = ${frac(`${V('τ')}²`, `${V('τ')}² + (${V('σ')}*)²/${V('n')}`)}`,
        `4. corrected value: &nbsp; ${V('x')}*${sub('gi')} = ${V('m')}${sub('g')} + ${V('s')}${sub('g')} ${frac(`${V('z')}${sub('gi')} − ${V('δ')}*${sub('gu')}`, `${V('σ')}*${sub('gu')}`)} + ${V('β')}${sub('g')} ${V('f')}${sub('i')}`,
        `matched groups: &nbsp; ${V('σ̂')}²${sub('gu')} = ${frac(`Σ${sub('c')} ${V('SS')}${sub('gc')} / (${V('s')}${sub('g')}² ${V('ω')}${sub('gc')})`, `Σ${sub('c')} (${V('n')}${sub('c')} − 1)`)}, &nbsp; ${V('δ̂')}${sub('gu')} = ${frac(`Σ${sub('c')} ${V('n')}${sub('c')} (${V('z̄')}${sub('gc')} − ${V('σ')}*${sub('gu')} ${V('ρ')}${sub('gc')})`, `Σ${sub('c')} ${V('n')}${sub('c')}`)}`],
        `1. r: the rank of a gene among the n genes of the sample; p its quantile position; Q<sub>ref</sub> the reference distribution of the group (bulk blood or PBMC), the average over the units of their average distribution, read at ${fmt(lim.grid_points)} positions. 2. β<sub>g</sub>: the sex effect of gene g, the slope of its aligned values on a female indicator f<sub>i</sub> (1 for a woman, 0 for a man), set aside and added back at the end. m<sub>g</sub> and s<sub>g</sub>: the anchor mean and standard deviation, from the healthy controls of the units that have at least ${spell(lim.minimum_healthy_controls)}, the spread taken around each control's own unit mean. 3. δ̂ and σ̂²: the mean and variance of a unit's healthy controls in anchor units, its location and scale before shrinkage; a and b: the shape and rate of the inverse-gamma prior, and δ̄ and τ² the mean and variance of the locations, all fitted over the genes of the unit; n: its healthy controls. Matched groups, for units without enough healthy controls: SS<sub>gc</sub> is the sum of squares of class c around its mean, z̄<sub>gc</sub> its mean in anchor units, and ρ<sub>gc</sub> and ω<sub>gc</sub> the mean and variance of the same class and age group in the units corrected on healthy controls.`)}
      ${H3('sim', 'Simulation', `How many healthy controls? The shrinkage of ${A.gene} in ${A.study}`)}
      <div class="controls"><label>healthy controls ${V('n')} = <b id="mkNv"></b> <input type="range" id="mkN" min="3" max="120" step="1" value="${A.n}"></label><button type="button" class="btn ghost" id="mkReset">Back to the ${fmt(A.n)} of ${A.study}</button></div>
      <div class="grid2"><div class="chart" id="mkCurve"></div><div class="chart" id="mkPrior"></div></div>
      <div class="eq" id="mkEq"></div>
      ${how(`The first chart draws the shrunk scale σ* and the weight w that the gene's own location keeps, for every number of controls, with the gene's own estimates and the prior of ${A.study} held as they are; the dot is the current n. The second shows the variances of the ${fmt(A.genes)} genes of the unit, from which the prior (the curve) is fitted, with the gene's own variance and its shrunk value. With few controls the prior pulls the own estimate strongly; with many, the gene keeps its own.`)}
      ${H3('ex', 'Worked example a', `aligning ${A.gene} in ${es.sample}`)}
      <div class="grid2"><div class="readout" id="mkAlign"></div><div class="chart" id="mkRef"></div></div>
      ${H3('ex', 'Worked example b', `the location and scale of ${A.gene} in ${A.study}`)}
      <div class="readout" id="mkEB"></div>
      ${H3('ex', 'Worked example c', `correcting ${A.gene} in two samples of ${A.study}`)}
      <div class="readout" id="mkCorr"></div>
      ${H3('ex', 'Worked example d', `${A.gene} in ${mt.study}, by matched groups`)}
      <div class="readout" id="mkMatched"></div>
      <h3>What each method assumes</h3>
      <ul class="find" id="mkAssume"></ul>
      ${H3('res', 'Result')}
      <div id="mkRes"></div>`;
    $('mkN').addEventListener('input', e => { ST.n = +e.target.value; drawShrink(); });
    $('mkReset').onclick = () => { ST.n = A.n; $('mkN').value = A.n; drawShrink(); };
    $('mkAlign').innerHTML = `${es.sample} is from a woman with sepsis. In the gene-level table of ${A.study} its ${A.gene} value is ${fx(es.raw, 3)} on a log₂ scale. Among the ${fmt(es.genes)} genes of the bulk-blood group that the sample holds, ${A.gene} has rank ${fmt(es.rank)}, so its quantile position is (${fmt(es.rank)} − 0.5)/${fmt(es.genes)} = ${fx(es.position, 3)}. The reference distribution has the value <b>${fx(refAt(ref, es.position), 3)}</b> there, and this becomes the aligned value. The number changed; the gene's place among the genes of the sample did not.`;
    const v = shrinkVar(A.own_variance, A.n, A.prior_shape, A.prior_rate), w = shrinkWeight(A.prior_spread, v, A.n), d = w * A.own_location + (1 - w) * A.prior_location, half = 0.5 * (A.n - 1);
    $('mkEB').innerHTML = `In anchor units, the ${fmt(A.n)} healthy controls of ${A.study} have a mean ${A.gene} of ${V('δ̂')} = ${fx(A.own_location, 3)} and a variance of ${fx(A.own_variance, 3)} (standard deviation ${fx(Math.sqrt(A.own_variance), 3)}): the unit reads the gene well above the anchor, and its controls spread less than the anchor's.
      <br><b>Scale.</b> The prior fitted over the genes of the unit has shape ${V('a')} = ${fx(A.prior_shape, 3)} and rate ${V('b')} = ${fx(A.prior_rate, 3)}. With ½(${V('n')} − 1) = ${fx(half, 1)}: (${V('σ')}*)² = (${fx(half, 1)} × ${fx(A.own_variance, 3)} + ${fx(A.prior_rate, 3)}) / (${fx(half, 1)} + ${fx(A.prior_shape - 1, 3)}) = ${fx(half * A.own_variance + A.prior_rate, 3)}/${fx(half + A.prior_shape - 1, 3)} = ${fx(v, 3)}, so ${V('σ')}* = <b>${fx(Math.sqrt(v), 3)}</b>. With only ${fmt(A.n)} controls the prior pulls the low own estimate a long way toward the typical scale of the unit.
      <br><b>Location.</b> Over the genes the locations average ${V('δ̄')} = ${fx(A.prior_location, 3)} with variance ${V('τ')}² = ${fx(A.prior_spread, 2)}. The uncertainty of the gene's own mean is ${fx(v, 3)}/${fmt(A.n)} = ${fx(v / A.n, 3)}, so ${V('w')} = ${fx(A.prior_spread, 2)}/(${fx(A.prior_spread, 2)} + ${fx(v / A.n, 3)}) = ${fx(w, 3)} and ${V('δ')}* = ${fx(w, 3)} × ${fx(A.own_location, 3)} + (1 − ${fx(w, 3)}) × ${fx(A.prior_location, 3)} = <b>${fx(d, 3)}</b>. Locations differ so much from gene to gene that the prior is wide, and the gene keeps nearly all of its own estimate.`;
    const zs = z(es), zc = z(ec);
    $('mkCorr').innerHTML = `For ${A.gene} the anchor, built from ${fmt(A.anchor_units)} units, has ${V('m')} = ${fx(A.anchor_mean, 3)} and ${V('s')} = ${fx(A.anchor_sd, 3)}; the sex effect is ${V('β')} = ${fx(A.sex_effect, 3)}, and both samples are from women, so ${V('f')} = 1.
      <br><b>The sepsis sample</b> ${es.sample}, aligned value ${fx(es.aligned, 3)}: ${V('z')} = (${fx(es.aligned, 3)} − ${fx(A.sex_effect, 3)} − ${fx(A.anchor_mean, 3)})/${fx(A.anchor_sd, 3)} = ${fx(zs, 2)}. Removing the unit's effects, (${fx(zs, 2)} − ${fx(A.location, 3)})/${fx(A.scale, 3)} = ${fx((zs - A.location) / A.scale, 2)}, and back on the log₂ scale: ${fx(A.anchor_mean, 3)} + ${fx(A.anchor_sd, 3)} × ${fx((zs - A.location) / A.scale, 2)} + ${fx(A.sex_effect, 3)} = <b>${fx(corr(es), 2)}</b>, the value of ${es.sample} in the corrected bulk-blood table (${fx(es.corrected, 2)}).
      <br><b>The healthy control</b> ${ec.sample} has the value ${fx(ec.raw, 3)}, aligned ${fx(ec.aligned, 3)}, ${V('z')} = ${fx(zc, 2)} and corrected value <b>${fx(corr(ec), 2)}</b>. The control now lies close to the anchor mean ${fx(A.anchor_mean, 2)}, and the patient ${fx((zs - A.location) / A.scale, 2)} anchor standard deviations above it, as ${A.gene} is expected to be raised in sepsis. Every sample of the unit gets the same shift and scale for a gene, so the order of its samples is kept.`;
    const zbar = (mt.mean - A.anchor_mean) / A.anchor_sd, vraw = mt.spread / mt.omega, after = (zbar - mt.location) / mt.scale;
    $('mkMatched').innerHTML = `${mt.study}, a study of the GAinS consortium, holds ${fmt(mt.samples)} leukocyte samples of adults with sepsis and no controls, so the sums have one class. Their mean ${A.gene}, aligned and without the sex effect, is ${fx(mt.mean, 3)}, in anchor units ${V('z̄')} = (${fx(mt.mean, 3)} − ${fx(A.anchor_mean, 3)})/${fx(A.anchor_sd, 3)} = ${fx(zbar, 3)}. The reference for adults with sepsis, from ${fmt(mt.reference_units)} units corrected on healthy controls, has mean ${V('ρ')} = ${fx(mt.rho, 3)} and variance ${V('ω')} = ${fx(mt.omega, 3)} in anchor units: in sepsis ${A.gene} is raised and varies far more than in health. The patients of ${mt.study} have a variance of ${fx(mt.spread, 3)}, so ${V('σ̂')}² = ${fx(mt.spread, 3)}/${fx(mt.omega, 3)} = ${fx(vraw, 3)}, and shrinkage gives ${V('σ')}* = ${fx(mt.scale, 3)}. Then ${V('δ̂')} = ${fx(zbar, 3)} − ${fx(mt.scale, 3)} × ${fx(mt.rho, 3)} = ${fx(zbar - mt.scale * mt.rho, 3)}, and shrinkage gives ${V('δ')}* = ${fx(mt.location, 3)}. After correction the patients' mean is (${fx(zbar, 3)} − ${fx(mt.location, 3)})/${fx(mt.scale, 3)} = <b>${fx(after, 2)}</b> anchor units, or ${fx(A.anchor_mean + A.anchor_sd * after, 2)} on the log₂ scale: the reference mean. Had these patients truly had higher or lower ${A.gene} than the reference patients, because they were sicker or sampled later, that difference would now be gone.`;
    $('mkAssume').innerHTML = `<li><b>Correction on healthy controls</b> assumes that healthy people are alike in every study, so that what differs between the controls of two studies is technical. Real differences in the age, ancestry or sampling of the controls are removed with it; the patients, who take no part in the estimate, keep theirs.</li>
      <li><b>Matched groups</b> assume that the patients of each class in a unit look, on average, like the patients of that class and age group in the units corrected on healthy controls. Differences in severity, timing, pathogen or the share of an endotype are removed with the study effect, and part of the separation between classes across these units is produced by the correction itself. The sample table records which method corrected each unit, so that an analysis can leave these units out.</li>
      <li><b>Both</b> keep the order of the samples within a unit, and estimate the scale within classes, because patients with sepsis vary much more than healthy people. Newborn studies were not corrected: no unit corrected on healthy controls is of newborns.</li>`;
    const U = K.units, T7 = K.tables, mu = K.matched_units, kb = K.kept.bulk_blood, kp = K.kept.pbmc;
    $('mkRes').innerHTML = RES(`In bulk blood the correction worked on <b>${fmt(U.bulk_blood.units)}</b> study units with ${fmt(U.bulk_blood.estimation_samples)} estimation samples and ${fmt(U.bulk_blood.genes)} genes: ${fmt(U.bulk_blood.anchored_units)} with at least ${spell(lim.minimum_healthy_controls)} healthy controls and ${fmt(U.bulk_blood.newborn_units)} of newborns, which were left out. Of the ${fmt(mu.unanchored_units)} units without enough healthy controls, ${fmt(mu.corrected_by_matched_groups)} were corrected by matched groups, and ${fmt(kb.matched_groups)} kept after the checks. The corrected bulk-blood table holds <b>${fmt(T7.bulk_blood.samples)}</b> samples of ${fmt(T7.bulk_blood.studies)} studies (${fmt(T7.bulk_blood.by_method.healthy_anchored)} corrected on healthy controls, ${fmt(T7.bulk_blood.by_method.matched_groups)} by matched groups), and the PBMC table <b>${fmt(T7.pbmc.samples)}</b> samples of ${fmt(T7.pbmc.studies)} studies and ${fmt(T7.pbmc.genes)} genes (${spell(kp.healthy_anchored)} units corrected on healthy controls and ${spell(kp.matched_groups)} by matched groups). The largest sex effect is that of ${esc(K.sex_effects.bulk_blood[0][0])}, ${fx(K.sex_effects.bulk_blood[0][1], 2)} in bulk blood (lower in women); ${fmt(K.sex_over_one.bulk_blood)} genes have an effect above one log₂ unit, all on the Y chromosome. For ${A.gene} it is ${fx(A.sex_effect, 3)}.`);
  }
  function drawShrink() {
    const A = M.correction.anchored, n = ST.n; $('mkNv').textContent = fmt(n);
    const v = shrinkVar(A.own_variance, n, A.prior_shape, A.prior_rate), w = shrinkWeight(A.prior_spread, v, n), d = w * A.own_location + (1 - w) * A.prior_location;
    const es = A.examples[0], zs = (es.aligned - A.sex_effect * es.female - A.anchor_mean) / A.anchor_sd, xs = A.anchor_mean + A.anchor_sd * (zs - d) / Math.sqrt(v) + A.sex_effect * es.female;
    // the curves of sigma* and w over n
    const wd = wOf('mkCurve'), h = 250, m = { l: 44, r: 44, t: 14, b: 38 }, g = chart('mkCurve', wd, h, 'Curves of the shrunk scale and the weight of the own location against the number of healthy controls.');
    const X = lin(3, 120, m.l, wd - m.r), Y = lin(0, 1.2, h - m.b, m.t), ns = []; for (let k = 3; k <= 120; k++) ns.push(k);
    gridY(g, Y, [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2], m.l, wd - m.r, t => fx(t, 1)); gridX(g, X, [10, 20, 40, 60, 80, 100, 120], h - m.b, h - m.b);
    axisTitle(g, (m.l + wd - m.r) / 2, h - 6, 'Healthy controls in the unit, n');
    S('line', { x1: m.l, x2: wd - m.r, y1: Y(Math.sqrt(A.own_variance)), y2: Y(Math.sqrt(A.own_variance)), stroke: C.tier.anchored, 'stroke-dasharray': '4 4' }, g);
    T(g, wd - m.r, Y(Math.sqrt(A.own_variance)) + 14, `own σ̂ = ${fx(Math.sqrt(A.own_variance), 2)}`, { 'text-anchor': 'end', 'font-size': 11, fill: C.tier.anchored });
    S('path', { d: ns.map((k, i) => `${i ? 'L' : 'M'}${X(k)},${Y(Math.sqrt(shrinkVar(A.own_variance, k, A.prior_shape, A.prior_rate)))}`).join(''), stroke: C.tier.anchored, 'stroke-width': 2.4, fill: 'none' }, g);
    S('path', { d: ns.map((k, i) => `${i ? 'L' : 'M'}${X(k)},${Y(shrinkWeight(A.prior_spread, shrinkVar(A.own_variance, k, A.prior_shape, A.prior_rate), k))}`).join(''), stroke: C.tier.unanchored, 'stroke-width': 2.4, fill: 'none' }, g);
    S('line', { x1: X(A.n), x2: X(A.n), y1: m.t, y2: h - m.b, stroke: C.lineStrong, 'stroke-dasharray': '2 3' }, g);
    T(g, X(A.n) + 4, h - m.b - 6, `${A.study}: ${fmt(A.n)}`, { 'font-size': 11, fill: C.muted });
    S('circle', { cx: X(n), cy: Y(Math.sqrt(v)), r: 5.5, fill: C.tier.anchored, stroke: C.bg, 'stroke-width': 2 }, g);
    S('circle', { cx: X(n), cy: Y(w), r: 5.5, fill: C.tier.unanchored, stroke: C.bg, 'stroke-width': 2 }, g);
    T(g, X(110), Y(shrinkWeight(A.prior_spread, shrinkVar(A.own_variance, 110, A.prior_shape, A.prior_rate), 110)) - 8, 'weight w', { 'text-anchor': 'end', 'font-size': 11.5, fill: C.tier.unanchored, 'font-weight': 600 });
    T(g, X(110), Y(Math.sqrt(shrinkVar(A.own_variance, 110, A.prior_shape, A.prior_rate))) - 8, 'shrunk scale σ*', { 'text-anchor': 'end', 'font-size': 11.5, fill: C.tier.anchored, 'font-weight': 600 });
    // the variances of the genes of the unit and the inverse-gamma prior fitted to them
    const hv = A.variance_histogram, w2 = wOf('mkPrior'), m2 = { l: 14, r: 14, t: 14, b: 38 }, g2 = chart('mkPrior', w2, h, `Histogram of the variances of the healthy controls of ${A.study} over its genes, with the fitted prior.`);
    const bins = hv.counts.length, bw = (hv.hi - hv.lo) / bins, total = sum(hv.counts), X2 = lin(hv.lo, 2.5, m2.l, w2 - m2.r);
    const dens = hv.counts.map(c => c / total / bw), ymax = Math.max(...dens, invGammaPdf(A.prior_rate / (A.prior_shape + 1), A.prior_shape, A.prior_rate)) * 1.1, Y2 = lin(0, ymax, h - m2.b, m2.t);
    gridX(g2, X2, [0, 0.5, 1, 1.5, 2, 2.5], m2.t, h - m2.b, t => fx(t, 1));
    axisTitle(g2, (m2.l + w2 - m2.r) / 2, h - 6, 'Variance of the controls in anchor units, per gene');
    hv.counts.forEach((c, i) => { const x0 = hv.lo + i * bw; if (x0 >= 2.5) return; rrect(g2, X2(x0), Y2(dens[i]), Math.max(1, X2(x0 + bw) - X2(x0) - 1), Y2(0) - Y2(dens[i]), rgba(C.lineStrong, 0.7), 1); });
    const vs = []; for (let x = 0.01; x <= 2.5; x += 0.01) vs.push(x);
    S('path', { d: vs.map((x, i) => `${i ? 'L' : 'M'}${X2(x)},${Y2(invGammaPdf(x, A.prior_shape, A.prior_rate))}`).join(''), stroke: C.ink2, 'stroke-width': 1.8, fill: 'none' }, g2);
    T(g2, X2(1.6), Y2(invGammaPdf(1.6, A.prior_shape, A.prior_rate)) - 8, 'prior', { 'font-size': 11.5, fill: C.ink2 });
    [[A.own_variance, `own ${fx(A.own_variance, 2)}`, C.muted], [v, `shrunk ${fx(v, 2)}`, C.tier.anchored]].forEach(([x, label, col], k) => {
      S('line', { x1: X2(x), x2: X2(x), y1: m2.t + 14, y2: h - m2.b, stroke: col, 'stroke-width': 2 }, g2);
      T(g2, X2(x) + (k ? 5 : -5), m2.t + 10 + k * 14, label, { 'text-anchor': k ? 'start' : 'end', 'font-size': 11, fill: col });
    });
    $('mkEq').innerHTML = `${V('n')} = ${fmt(n)}: &nbsp; (${V('σ')}*)² = ${frac(`${fx(0.5 * (n - 1), 1)} × ${fx(A.own_variance, 3)} + ${fx(A.prior_rate, 3)}`, `${fx(0.5 * (n - 1), 1)} + ${fx(A.prior_shape - 1, 3)}`)} = ${NB(fx(v, 3))}, &nbsp; ${V('w')} = ${frac(fx(A.prior_spread, 2), `${fx(A.prior_spread, 2)} + ${fx(v, 3)}/${fmt(n)}`)} = ${NB(fx(w, 3))}, &nbsp; ${V('δ')}* = ${NB(fx(d, 3))}
      <span class="c">With these effects the sepsis sample ${es.sample} would get a corrected ${A.gene} of ${fx(xs, 2)}${n === A.n ? ', its value in the corrected table' : ` instead of ${fx(es.corrected, 2)}`}. The own estimates (δ̂ = ${fx(A.own_location, 3)}, σ̂² = ${fx(A.own_variance, 3)}) and the prior are held fixed; in the real unit they would move a little with other controls.</span>`;
    // the reference distribution with the sample's position
    const ref = M.correction.reference, w3 = wOf('mkRef'), h3 = 220, m3 = { l: 44, r: 12, t: 12, b: 38 }, g3 = chart('mkRef', w3, h3, 'Curve of the reference distribution of bulk blood, with the position of the gene in the sample.');
    const X3 = lin(0, 1, m3.l, w3 - m3.r), lo3 = Math.floor(ref.value[0]), hi3 = Math.ceil(ref.value[ref.value.length - 1]), Y3 = lin(lo3, hi3, h3 - m3.b, m3.t);
    gridX(g3, X3, [0, 0.25, 0.5, 0.75, 1], m3.t, h3 - m3.b, t => fx(t, 2)); gridY(g3, Y3, ticks(lo3, hi3, 5), m3.l, w3 - m3.r);
    axisTitle(g3, (m3.l + w3 - m3.r) / 2, h3 - 6, 'Quantile position'); axisTitle(g3, 12, (m3.t + h3 - m3.b) / 2, 'Reference value, log₂', true);
    const step = Math.max(1, Math.floor(ref.value.length / 400));
    S('path', { d: ref.value.filter((_, i) => i % step === 0 || i === ref.value.length - 1).map((v2, i, arr) => { const k = Math.min(i * step, ref.value.length - 1); const p = ref.first + (ref.last - ref.first) * k / (ref.value.length - 1); return `${i ? 'L' : 'M'}${X3(p)},${Y3(v2)}`; }).join(''), stroke: C.ink2, 'stroke-width': 2, fill: 'none' }, g3);
    const pa = es.position, va = refAt(ref, pa);
    S('path', { d: `M${X3(pa)},${h3 - m3.b} L${X3(pa)},${Y3(va)} L${m3.l},${Y3(va)}`, stroke: C.cls.sepsis, 'stroke-dasharray': '4 3', fill: 'none', 'stroke-width': 1.5 }, g3);
    S('circle', { cx: X3(pa), cy: Y3(va), r: 5.5, fill: C.cls.sepsis, stroke: C.bg, 'stroke-width': 2 }, g3);
    T(g3, X3(pa) - 8, Y3(va) - 8, `${A.gene}: p = ${fx(pa, 3)} → ${fx(va, 2)}`, { 'text-anchor': 'end', 'font-size': 11.5, fill: C.ink });
  }

  // ======================================================================= 10. the held-out test
  const METHOD_LABEL = { none: 'No correction', plain_combat: 'Plain ComBat', flat_genes: 'Flat genes', matched_groups: 'Matched groups' };
  const methodColor = k => ({ none: C.muted, plain_combat: C.cls.other_condition, flat_genes: C.fam.illumina_beadchip, matched_groups: C.tier.unanchored })[k];
  // the sampling floor of n hidden controls over G genes: half the root mean square distance between the means of two
  // random halves, in within-study standard deviations (invented controls drawn from a standard normal distribution)
  function samplingFloor(n, G, r) {
    let ss = 0;
    for (let g = 0; g < G; g++) { const v = Array.from({ length: n }, () => gauss(r)); const idx = v.map((_, i) => [r(), i]).sort((a, b) => a[0] - b[0]).map(x => x[1]); const h = Math.floor(n / 2); const a = mean(idx.slice(0, h).map(i => v[i])), b = mean(idx.slice(h).map(i => v[i])); ss += (a - b) ** 2; }
    return Math.sqrt(ss / G) / 2;
  }
  function buildHeldout() {
    const Hd = M.heldout, ex = Hd.units.find(u => u.study === Hd.example), md = Hd.medians;
    ST.ho = { n: ex.hidden, seed: 1, draws: [], method: 'matched_groups' };
    const mi = k => Hd.methods.indexOf(k);
    $('mbHeldout').innerHTML = `
      ${H3('eq', 'Equation')}
      ${EQ([`${V('d')}${sub('g')} = ${frac(`${V('h̄')}${sub('g')} − ${V('p̄')}${sub('g')}`, `${V('σ')}${sub('g')}`)}, &nbsp; location error = √( ${frac('1', V('G'))} Σ${sub('g')} ${V('d')}${sub('g')}² )`,
        `sampling floor = ½ √( ${frac('1', V('G'))} Σ${sub('g')} ( ${V('h̄')}${sub('g')}${'<sup>(1)</sup>'} − ${V('h̄')}${sub('g')}${'<sup>(2)</sup>'} )² / ${V('σ')}${sub('g')}² )`],
        'h̄<sub>g</sub>: the mean of the hidden healthy controls of the held-out unit for gene g, after a method has corrected the unit from its other samples; p̄<sub>g</sub>: the average of the control means of the other units, after the same method; σ<sub>g</sub>: the pooled within-unit standard deviation of those controls; G: the genes. The pool dispersion is the same distance for each unit of the training pool from the other units of the pool. The floor splits the hidden controls into two random halves, (1) and (2): it is the error that chance alone gives a group of that size.')}
      ${H3('sim', 'Simulation', 'The sampling floor, on invented controls')}
      <p class="note mtop0">Invented healthy controls drawn from one standard normal distribution for each of 200 genes, so that no method error exists; only chance remains. Choose how many controls are hidden.</p>
      <div class="controls"><label>hidden controls = <b id="mhNv"></b> <input type="range" id="mhN" min="4" max="120" step="1" value="${ex.hidden}"></label><button type="button" class="btn" id="mhDraw">Draw again</button></div>
      <div class="chart" id="mhSim"></div>
      <div class="readout" id="mhSimOut"></div>
      ${how('Each dot is the floor computed from one random draw of the invented controls; the curve is 1/√n, the standard error of a mean of n values. The more controls are hidden, the closer a perfect method could place their mean. The real floors of the held-out units are drawn as rings at their own numbers of controls.')}
      ${H3('ex', 'Worked example', `${ex.study} held out`)}
      <div class="readout" id="mhEx"></div>
      ${H3('res', 'Result', 'Every unit held out, by method')}
      <div class="controls"><span>Pick out</span><div class="seg" id="mhM">${segHTML(Hd.methods.map(k => [k, METHOD_LABEL[k]]), ST.ho.method)}</div></div>
      <div class="chart" id="mhUnits"></div>
      <div class="keys" id="mhKeys"></div>
      ${how('One row per held-out bulk-blood unit, sorted by its error without correction: the location error of its hidden controls after each method (dots) and its sampling floor (the bar). The chosen method is drawn larger. Lower is better; nothing can be expected below the floor.')}
      <div class="tablewrap" id="mhTable"></div>
      <div id="mhRes"></div>`;
    $('mhN').addEventListener('input', e => { ST.ho.n = +e.target.value; ST.ho.draws = []; drawFloor(); });
    $('mhDraw').onclick = () => { ST.ho.seed++; ST.ho.draws = []; drawFloor(); };
    seg('mhM', v => { ST.ho.method = v; drawUnits(); });
    const e = ex.error, ds = ex.dispersion;
    $('mhEx').innerHTML = `With its ${fmt(ex.hidden)} healthy controls hidden, ${ex.study} has ${fmt(ex.visible)} sepsis samples to be corrected from. Uncorrected, the hidden controls lie <b>${fx(e[mi('none')], 2)}</b> standard deviations from the other units' controls. Matched groups, using only the patients, bring them to <b>${fx(e[mi('matched_groups')], 2)}</b>; flat genes to ${fx(e[mi('flat_genes')], 2)}; plain ComBat to ${fx(e[mi('plain_combat')], 2)}. The pool dispersion is ${fx(ds[mi('matched_groups')], 2)} for matched groups and ${fx(ds[mi('flat_genes')], 2)} for flat genes, but ${fx(ds[mi('plain_combat')], 2)} for plain ComBat: by setting each training unit's overall mean to one value, whatever its mix of patients and controls, it leaves the controls of those units spread apart. The sampling floor of the ${fmt(ex.hidden)} hidden controls is ${fx(ex.floor, 2)}, so no method could be expected to come much closer than that.`;
    const row = (g, k) => md[g] && md[g][k] ? [fx(md[g][k].error, 2), fx(md[g][k].dispersion, 2)] : ['–', '–'];
    $('mhTable').innerHTML = `<table class="t"><thead><tr><th>Method, medians over the units</th><th class="r">Bulk blood: location error</th><th class="r">pool dispersion</th><th class="r">PBMC: location error</th><th class="r">pool dispersion</th></tr></thead><tbody>${Hd.methods.map(k => `<tr><td>${sw(methodColor(k))}${METHOD_LABEL[k]}</td><td class="r">${row('bulk_blood', k)[0]}</td><td class="r">${row('bulk_blood', k)[1]}</td><td class="r">${row('pbmc', k)[0]}</td><td class="r">${row('pbmc', k)[1]}</td></tr>`).join('')}<tr class="total"><td>Sampling floor</td><td class="r">${fx(md.bulk_blood.floor, 2)}</td><td></td><td class="r">${fx(md.pbmc.floor, 2)}</td><td></td></tr></tbody></table>
      <p class="note">${fmt(md.bulk_blood.units)} bulk-blood units of ${fmt(md.bulk_blood.studies)} studies and ${fmt(md.pbmc.units)} PBMC units; matched groups applied to ${fmt(md.bulk_blood.matched_groups.units)} of the bulk-blood units, the others having too few samples of a matched class. The flat-gene method predicts a unit's effects from the ${fmt(Hd.flat_genes)} genes that sepsis changes least; since sepsis changes nearly every gene of whole blood, it was kept as a comparison only.</p>`;
    const pc = Hd.paired;
    $('mhRes').innerHTML = RES(`Before correction the healthy controls of one study lie a median of <b>${fx(md.bulk_blood.none.dispersion, 2)}</b> within-study standard deviations from those of the other studies: two groups of healthy people look further apart than most healthy people within one study. Plain ComBat has the lowest median location error, ${fx(md.bulk_blood.plain_combat.error, 2)}, against ${fx(md.bulk_blood.matched_groups.error, 2)} for matched groups, but leaves the pool spread by ${fx(md.bulk_blood.plain_combat.dispersion, 2)}.${pc ? ` Unit by unit the two are close: the median paired difference, ComBat minus matched groups, is ${fx(pc.median_difference, 2)} (${pct(Hd.level, 0)} interval ${fx(pc.difference_low, 2)} to ${fx(pc.difference_high, 2)}; Wilcoxon signed-rank p = ${fx(pc.wilcoxon_p, 3)}), and ComBat is closer in ${fmt(pc.method_closer)} of ${fmt(pc.units)} units.` : ''} The error of matched groups is several times the floor of ${fx(md.bulk_blood.floor, 2)}: a residual difference between studies remains after correction, and analyses should allow for it. The alignment and the sex effect were estimated once with the held-out unit included, so the test is slightly optimistic for every method alike.`);
  }
  function drawFloor() {
    const Hd = M.heldout, w = wOf('mhSim'), h = 230, m = { l: 44, r: 14, t: 12, b: 38 }, n = ST.ho.n; $('mhNv').textContent = fmt(n);
    if (!ST.ho.draws.length) { const r = rng(1000 * ST.ho.seed + n); ST.ho.draws = Array.from({ length: 12 }, () => samplingFloor(n, 200, r)); }
    const g = chart('mhSim', w, h, 'Dots of simulated sampling floors against the number of hidden controls, with the curve 1 over the square root of n.');
    const X = lin(0, 120, m.l, w - m.r), Y = lin(0, 0.6, h - m.b, m.t);
    gridX(g, X, [0, 20, 40, 60, 80, 100, 120], m.t, h - m.b); gridY(g, Y, [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6], m.l, w - m.r, t => fx(t, 1));
    axisTitle(g, (m.l + w - m.r) / 2, h - 6, 'Hidden healthy controls, n'); axisTitle(g, 12, (m.t + h - m.b) / 2, 'Sampling floor', true);
    const ns = []; for (let k = 3; k <= 120; k++) ns.push(k);
    S('path', { d: ns.map((k, i) => `${i ? 'L' : 'M'}${X(k)},${Y(Math.min(0.6, 1 / Math.sqrt(k)))}`).join(''), stroke: C.lineStrong, 'stroke-width': 1.8, fill: 'none' }, g);
    Hd.units.filter(u => u.group === 'bulk_blood' && u.hidden <= 120).forEach(u => hover(S('circle', { cx: X(u.hidden), cy: Y(Math.min(0.6, u.floor)), r: 4, fill: 'none', stroke: C.tier.anchored, 'stroke-width': 1.5 }, g), `<b>${esc(u.study)}</b> <span class="m">${esc(u.platform)}</span><br>${fmt(u.hidden)} hidden controls, floor ${fx(u.floor, 2)}`));
    ST.ho.draws.forEach((f, i) => S('circle', { cx: X(n) + (i - 5.5) * 1.6, cy: Y(Math.min(0.6, f)), r: 3.5, fill: C.ink, 'fill-opacity': 0.7 }, g));
    $('mhSimOut').innerHTML = `With ${fmt(n)} hidden controls the floor comes out between <b>${fx(Math.min(...ST.ho.draws), 2)}</b> and <b>${fx(Math.max(...ST.ho.draws), 2)}</b> in ${fmt(ST.ho.draws.length)} draws, close to 1/√${fmt(n)} = ${fx(1 / Math.sqrt(n), 2)}.`;
  }
  function drawUnits() {
    const Hd = M.heldout, units = Hd.units.filter(u => u.group === 'bulk_blood').sort((a, b) => b.error[0] - a.error[0]), w = wOf('mhUnits'), rowH = 17, left = small(w) ? 92 : 150, H = rowH * units.length + 46;
    const g = chart('mhUnits', w, H, 'Dot chart of the location error of every held-out unit after each method, with its sampling floor.');
    const xmax = niceMax(Math.max(...units.flatMap(u => u.error.filter(v => v != null)))), X = lin(0, xmax, left, w - 16);
    gridX(g, X, ticks(0, xmax, 5), 6, rowH * units.length + 8, v => fx(v, 1));
    axisTitle(g, (left + w - 16) / 2, H - 4, 'Location error of the hidden controls, within-study standard deviations');
    units.forEach((u, i) => {
      const y = 14 + i * rowH;
      T(g, left - 8, y + 4, u.study, { 'text-anchor': 'end', 'font-size': 11, fill: C.ink2 });
      rrect(g, X(0), y - 3, X(u.floor) - X(0), 6, C.lineStrong, 2);
      Hd.methods.forEach((k, j) => { if (u.error[j] == null) return; const on = k === ST.ho.method; S('circle', { cx: X(u.error[j]), cy: y, r: on ? 5 : 3, fill: methodColor(k), 'fill-opacity': on ? 1 : 0.55, stroke: on ? C.bg : 'none' }, g); });
      hitRow(g, 0, y - rowH / 2, w, rowH, `<b>${esc(u.study)}</b> <span class="m">${esc(u.platform)} · ${esc(lab('family', u.family))}</span><br>${fmt(u.hidden)} hidden controls, ${fmt(u.visible)} other estimation samples<br>${Hd.methods.map((k, j) => u.error[j] == null ? '' : `${sw(methodColor(k))}${METHOD_LABEL[k]}: ${fx(u.error[j], 2)}`).filter(Boolean).join('<br>')}<br>sampling floor ${fx(u.floor, 2)}`);
    });
    $('mhKeys').innerHTML = keys([...Hd.methods.map(k => [METHOD_LABEL[k], methodColor(k), 'dot']), ['sampling floor', C.lineStrong]]);
  }

  // ======================================================================= 11. benchmark splits and leakage
  const TASK_LABEL = { sepsis_vs_healthy_control: 'Sepsis against healthy controls', sepsis_vs_non_septic_control: 'Sepsis against non-septic controls', mortality_in_sepsis: 'Death against survival in sepsis' };
  const SCHEME_LABEL = { leave_one_group_out: 'Leave one group out', technology_transfer: 'Technology transfer', temporal: 'Temporal' };
  const REASON_LABEL = { not_primary_record: 'not the primary record of its sample', not_one_per_patient: 'not the one sample marked for its patient', quality_flag: 'a quality flag', tissue_outside_task: 'not whole blood, blood not further specified or leukocytes', phase_outside_task: 'drawn before onset, or in recovery or follow-up', class_outside_task: 'a class outside the task', label_not_recorded: 'no recorded outcome' };
  // the invented example: two held-out studies whose classifier knows the study and, with 'signal' > 0, some sepsis
  function pooledExample(sa, sb, signal) {
    const A = { pos: sa, neg: 10 - sa, base: 0.8 }, B = { pos: sb, neg: 10 - sb, base: 0.2 };
    const scores = st => ({ pos: Array(st.pos).fill(st.base + signal), neg: Array(st.neg).fill(st.base) });
    const a = scores(A), b = scores(B);
    const within = s => s.pos.length && s.neg.length ? auroc(s.pos, s.neg) : null;
    return { a, b, A, B, auA: within(a), auB: within(b), pooled: auroc([...a.pos, ...b.pos], [...a.neg, ...b.neg]) };
  }
  function buildBench() {
    const Bn = M.bench, f = Bn.fold, t1 = Bn.tasks.find(t => t.task === f.task), R = M.repeats, P = R.patient;
    ST.pool = { a: 9, b: 1, signal: 0 };
    $('mbBench').innerHTML = `
      <ul class="find">
        <li><b>The same patient on both sides.</b> Patient ${esc(P.patient)} of ${P.study} was sampled on ${spell(P.samples.length)} days. Split at random, a classifier could learn this patient's own profile from one day and recognize it on another.</li>
        <li><b>Two deposits of one sample.</b> The array of ${esc(R.sample_id)} is deposited in ${spell(R.deposits.length)} studies. Split study by study, GSE66099 could train and GSE4607 test, and the classifier would be scored on an array it was trained on.</li>
        <li><b>Studies confounded with the class.</b> In the task of sepsis against healthy controls, ${fmt(Bn.sepsis_only_groups)} of the ${fmt(t1.groups)} patient-safe groups hold sepsis samples and no healthy controls, and ${spell(Bn.control_only_groups)} hold controls and no sepsis: a classifier can learn the study instead of the illness. Holding out whole groups removes that reward.</li>
        <li><b>Pooled scoring.</b> Even with whole groups held out, scoring them together can reward recognizing the study, as the simulation shows.</li>
      </ul>
      ${H3('eq', 'Equation')}
      ${EQ([`score of a fold = { AUROC${sub('k')} : held-out group ${V('k')} with ≥ ${fmt(Bn.policy.minimum_per_class)} samples of each class }, &nbsp; never one AUROC over the pooled groups`],
        `Every fold assigns whole patient-safe groups to one side, so that no deposit of a sample and no identified patient lies on both sides. Leave one group out holds out each group with at least ${spell(Bn.policy.minimum_per_class)} samples of each class in turn; technology transfer trains on the array groups and tests the RNA-seq groups, and the reverse; the temporal split trains on the groups first released before ${fdate(Bn.policy.temporal_split_date)} and tests the later ones.`)}
      ${H3('sim', 'Simulation', 'Pooled against per-study scoring, on two invented studies')}
      <p class="note mtop0">An invented example. Two held-out studies of ten samples each. The classifier scores every sample of study A 0.8 and of study B 0.2, plus a bonus for sepsis that you can set; at zero it has learned only which study a sample comes from.</p>
      <div class="controls">
        <label>sepsis samples in study A <b id="mbAv"></b> <input type="range" id="mbA" min="1" max="9" step="1" value="9"></label>
        <label>in study B <b id="mbBv"></b> <input type="range" id="mbB" min="1" max="9" step="1" value="1"></label>
        <label>what it knows of sepsis <b id="mbSv"></b> <input type="range" id="mbS" min="0" max="0.3" step="0.05" value="0"></label>
      </div>
      <div class="grid2"><div class="chart" id="mbPool"></div><div class="readout" id="mbPoolOut" aria-live="polite"></div></div>
      ${how('Each square is a sample, filled for sepsis and open for a control, placed at its score. Within each study the AUROC asks only whether its own sepsis samples score above its own controls; the pooled AUROC also compares sepsis samples of one study with controls of the other, and so rewards the difference between the studies.')}
      ${H3('ex', 'Worked example', `the fold that holds out ${f.group}`)}
      <div class="readout" id="mbEx"></div>
      ${H3('res', 'Result', 'The three tasks')}
      <div class="tablewrap" id="mbTasks"></div>
      <div class="chart" id="mbReasons"></div>
      ${how(`Every record is listed for every task, either as a sample of the task or with the first reason that excludes it; the bars are the records that each reason removes from the task of sepsis against healthy controls, in the order in which the reasons are checked.`)}
      <div id="mbRes"></div>`;
    ['A', 'B', 'S'].forEach(k => $(`mb${k}`).addEventListener('input', () => { ST.pool = { a: +$('mbA').value, b: +$('mbB').value, signal: +$('mbS').value }; drawPool(); }));
    const gr = R.group_rows, plats = f.platforms;
    $('mbEx').innerHTML = `In leave one group out for sepsis against healthy controls, one fold holds out the patient-safe group ${f.group}, the cohort group of the array of patient ${esc((R.deposits.find(x => x.patient_id.includes('-')) || R.deposits[0]).patient_id)}: ${spell(gr.length)} study platforms (${andList(gr.map(x => x.study))}) with ${fmt(sum(gr.map(x => x.records)))} records but ${fmt(sum(gr.map(x => x.primary)))} distinct samples. After the exclusions, ${f.group} gives the task <b>${fmt(f.test_positives + f.test_negatives)}</b> samples, ${fmt(f.test_positives)} with sepsis and ${fmt(f.test_negatives)} healthy controls, from ${spell(plats.length)} of its study platforms (${andList(plats.map(([s, , n]) => `${s} ${fmt(n)}`))}); ${andList(gr.filter(x => !plats.some(p => p[0] === x.study)).map(x => x.study))} give none, because all their records are repeat deposits whose primary records lie in the other studies of the group. The ${fmt(f.train_groups)} other groups train, with ${fmt(f.train_positives)} sepsis samples and ${fmt(f.train_negatives)} healthy controls. As a check, ${fmt(f.train_positives)} + ${fmt(f.test_positives)} = <b>${fmt(f.train_positives + f.test_positives)}</b> and ${fmt(f.train_negatives)} + ${fmt(f.test_negatives)} = <b>${fmt(f.train_negatives + f.test_negatives)}</b>, the totals of the task.
      <br><br>Had the split been made by study, the ${fmt(gr.find(x => x.study === 'GSE66099').records)} records of GSE66099 could have trained while the same arrays were tested through their primary records. In the temporal scheme ${f.group} is dated by its first release, ${fdate(f.first_release)}, so it trains. When RNA-seq trains and arrays are tested, ${f.group} and ${Bn.contrast.group}, the group of ${andList(Bn.contrast.studies)}, are both test groups; ${f.group} holds ${fmt(f.test_positives)} sepsis samples to ${fmt(f.test_negatives)} controls and ${Bn.contrast.group} ${fmt(Bn.contrast.positives)} to ${fmt(Bn.contrast.negatives)}. Pooled, a classifier that only recognized the two groups would look as if it recognized sepsis.`;
    $('mbTasks').innerHTML = `<table class="t"><thead><tr><th>Task</th><th class="r">Samples</th><th class="r">Positive</th><th class="r">Negative</th><th class="r">Patient-safe groups</th>${Bn.schemes.map(s => `<th class="r">Folds: ${SCHEME_LABEL[s].toLowerCase()}</th>`).join('')}</tr></thead><tbody>${Bn.tasks.map(t => `<tr><td>${TASK_LABEL[t.task] || esc(t.task)}</td><td class="r">${fmt(t.eligible)}</td><td class="r">${fmt(t.positives)}</td><td class="r">${fmt(t.negatives)}</td><td class="r">${fmt(t.groups)}</td>${Bn.schemes.map(s => `<td class="r">${fmt(t.folds[s] || 0)}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
    const folds = sum(Bn.tasks.map(t => sum(Object.values(t.folds))));
    $('mbRes').innerHTML = RES(`The ${spell(Bn.tasks.length)} tasks have <b>${fmt(folds)}</b> folds in all. The positive class is sepsis in the first two tasks and death in the third, whose labels keep each study's own definition of the outcome. The benchmark uses the per-study tables, not the corrected ones: the correction used the healthy controls of every study, and in matched-group units its sepsis samples too, so a model tested on corrected values would be given part of the answer. Any scaling across studies has to be fitted on the training groups of a fold, or applied to each sample on its own. In the fold that trains on arrays and tests on RNA-seq for the first task there are ${fmt(Bn.transfer.test_groups)} scored test groups, and so ${fmt(Bn.transfer.test_groups)} scores; ${fmt(Bn.transfer.not_scored)} further groups are held out but not scored.`);
  }
  function drawPool() {
    const p = ST.pool, r = pooledExample(p.a, p.b, p.signal);
    $('mbAv').textContent = fmt(p.a); $('mbBv').textContent = fmt(p.b); $('mbSv').textContent = fx(p.signal, 2);
    const w = wOf('mbPool'), h = 200, m = { l: 70, r: 14, t: 14, b: 36 }, g = chart('mbPool', w, h, 'Scores of the samples of two invented studies.');
    const X = lin(0, 1.2, m.l, w - m.r);
    gridX(g, X, [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2], m.t, h - m.b, v => fx(v, 1)); axisTitle(g, (m.l + w - m.r) / 2, h - 6, 'Score given by the classifier');
    [['Study A', r.a, 0], ['Study B', r.b, 1]].forEach(([name, s, k]) => {
      const y0 = m.t + 20 + k * 70; T(g, m.l - 10, y0 + 12, name, { 'text-anchor': 'end', 'font-size': 12, fill: C.ink2, 'font-weight': 600 });
      [[s.pos, 1], [s.neg, 0]].forEach(([vals, isP]) => vals.forEach((v, j) => { const col = isP ? C.cls.sepsis : C.cls.healthy_control; rrect(g, X(v) - 27 + (j % 5) * 11, y0 - 4 + Math.floor(j / 5) * 11 + (isP ? 0 : 26), 9, 9, isP ? col : 'none', 2, { stroke: col, 'stroke-width': 1.4 }); }));
    });
    const f = v => v == null ? 'not defined (one class only)' : fx(v, 2);
    $('mbPoolOut').innerHTML = `<b>Within study A:</b> AUROC ${f(r.auA)}<br><b>Within study B:</b> AUROC ${f(r.auB)}<br><b>Pooled:</b> AUROC <b>${fx(r.pooled, 2)}</b><br><span class="small">${p.signal === 0 ? 'The classifier knows nothing of sepsis, so within each study it is no better than chance; pooled, it looks good only because the studies differ in their share of sepsis.' : 'With some knowledge of sepsis the within-study AUROC rises; the pooled AUROC mixes it with the difference between the studies.'}</span>`;
  }
  function drawReasons() {
    const Bn = M.bench, t1 = Bn.tasks[0], rows = Bn.reasons.map(k => [k, t1.excluded[k] || 0]).filter(([, n]) => n > 0), w = wOf('mbReasons'), left = small(w) ? 170 : 330, rowH = 25, H = rowH * (rows.length + 1) + 34;
    const g = chart('mbReasons', w, H, 'Bar chart of the records excluded from the first task by each reason, and the samples of the task.'), xm = niceMax(Math.max(...rows.map(r => r[1]), t1.eligible)), X = lin(0, xm, left, w - 50);
    gridX(g, X, ticks(0, xm, 4), 4, rowH * (rows.length + 1) + 4); axisTitle(g, (left + w - 50) / 2, H - 2, 'Records');
    [...rows, ['eligible', t1.eligible]].forEach(([k, n], i) => { const y = 6 + i * rowH, lbl = k === 'eligible' ? 'Samples of the task' : `Excluded: ${REASON_LABEL[k] || k}`; T(g, left - 8, y + 14, small(w) ? lbl.replace(/^(.{26}).+$/, '$1…') : lbl, { 'text-anchor': 'end', 'font-size': 12, fill: k === 'eligible' ? C.ink : C.ink2, 'font-weight': k === 'eligible' ? 600 : 400 }); rrect(g, X(0), y + 3, X(n) - X(0), rowH - 8, k === 'eligible' ? C.accent : C.lineStrong, 3); T(g, X(n) + 5, y + 14, fmt(n), { 'font-size': 11.5, fill: C.muted }); hitRow(g, 0, y, w, rowH, `<b>${esc(lbl)}</b><br>${fmt(n)} records`); });
  }

  // ======================================================================= the view
  const DRAW = {
    'm-search': () => { drawSearch(); drawLadder(); drawFunnel(); },
    'm-arrays': () => { drawBackground(); drawQN(); drawMP(); drawArrayDots(); },
    'm-rnaseq': () => { drawCPM(); drawFlagChart(); drawZeros(); },
    'm-genes': drawGenes, 'm-labels': drawLabels,
    'm-repeats': () => { drawDays(); drawRepeatCharts(); },
    'm-sex': drawSex, 'm-auroc': drawScore, 'm-correct': drawShrink,
    'm-heldout': () => { drawFloor(); drawUnits(); },
    'm-bench': () => { drawPool(); drawReasons(); },
  };
  const draw = id => { if (ST.sex) stopSex(); const f = DRAW[id]; if (f) f(); };
  registerView('methods', {
    part: 'methods',
    init() {
      M = part('methods');
      [buildOver, buildSearch, buildArrays, buildArraysText, buildCounts, buildGenes, buildLabels, buildRepeats, buildSex, buildScore, buildCorrection, buildHeldout, buildBench].forEach(f => f());
    },
    show(route) { draw(route.chapter); },
    redraw() { draw(ROUTE.chapter); },
  });

  // ======================================================================= checks (?selftest=1)
  // every worked example is recomputed here from its exported inputs and compared with its exported result
  const T_ = (name, fn) => selftest('methods', name, fn);
  T_('every chapter has its Equation, Simulation, Worked example and Result', () => {
    const missing = secs('methods').filter(s => s.id !== 'm-over').map(s => [s.id, ['t-eq', 't-sim', 't-ex', 't-res'].filter(k => !s.querySelector(`h3 .tag.${k}`))]).filter(([, m]) => m.length);
    return ok(!missing.length, missing.map(([id, m]) => `${id}: ${m.join(', ')}`).join('; '));
  });
  T_('overview: the flow of the figure gives the main numbers of the core', () => {
    const f = M.overview.flow, N = D.numbers;
    return ok(f.included_studies === N.studies && f.candidate_studies === N.candidates && f.included_sample_records === N.included_records
      && f.records_without_values === N.records_without_values && f.sample_records_with_values === N.records && f.distinct_samples === N.distinct
      && f.one_sample_per_patient === N.one_per_patient,
      `${f.included_studies}/${N.studies}, ${f.sample_records_with_values}/${N.records}, ${f.distinct_samples}/${N.distinct}`);
  });
  T_('search: the excluded and included studies add up to the candidates, and each reason counts its candidates', () => {
    const s = M.search, byR = s.reasons.map((r, i) => s.candidates.filter(c => c.reason === i + 1).length);
    const bad = s.reasons.filter((r, i) => byR[i] !== r.studies).map(r => r.key);
    return ok(sum(s.reasons.map(r => r.studies)) + s.included === s.candidates_total && s.candidates.length === s.candidates_total && s.candidates.filter(c => c.reason === 0).length === s.included && !bad.length, bad.join(', '));
  });
  T_('search: every term on reaches every candidate; the funnel removes each reason in order down to the included studies', () => {
    const s = M.search, all = (1 << s.terms.length) - 1, f = s.funnel;
    const steps = f.slice(1).every((x, i) => x.remaining === f[i].remaining - x.removed && x.step === s.reasons[i].key && x.removed === s.reasons[i].studies);
    return ok(s.candidates.every(c => reached(c, all)) && steps && f[f.length - 1].remaining === s.included && f[0].remaining === s.candidates_total);
  });
  T_('search: GSE103119 stops at its recorded reason and GSE65682 passes every reason', () => {
    const s = M.search, c1 = s.candidates.find(c => c.accession === 'GSE103119'), c2 = s.candidates.find(c => c.accession === 'GSE65682');
    return ok(s.reasons[c1.reason - 1].key === s.excluded.reason && c2.reason === 0 && s.excluded.conditions.reduce((t, [, n]) => t + n, 0) === s.excluded.records);
  });
  T_('arrays: quantile normalization of the invented arrays gives the table of the guide', () => {
    const q = quantileNormalize([0, 1, 2].map(j => QN0.map(r => r[j]))), out = QN0.map((_, p) => [0, 1, 2].map(j => Math.log2(q.out[j][p])));
    const want = [[1, 2, 1], [3, 1, 2], [2, 3, 3], [4, 4, 4]];
    return ok(out.every((r, p) => r.every((v, j) => close(v, want[p][j], 1e-9))), JSON.stringify(out));
  });
  T_('arrays: median polish of the invented probe set gives 6, 8 and 7, with probe levels 7, 10 and 6', () => {
    const p = medianPolish(MP0);
    return ok(p.values.every((v, j) => close(v, [6, 8, 7][j], 1e-9)) && p.r.every((v, i) => close(p.t + v, [7, 10, 6][i], 1e-9)), p.values.join(', '));
  });
  T_('arrays: normexp recomputed from the background of GSM1317896 gives its three probes', () => {
    const b = M.arrays.background, bad = b.probes.filter(p => { const e = normexp(p.x, b.mu, b.sigma, b.alpha); return !close(e.s, p.signal, 0.01) || !close(Math.log2(e.s + b.offset), p.log, 0.002) || !close(e.m, p.m, 0.01); });
    return ok(!bad.length, bad.map(p => p.probe).join(', '));
  });
  T_('arrays: GSM648590 is the first array of GSE26440 scanned, with the medians of the worked example', () => {
    const r = M.arrays.rma, e = r.each[0];
    return ok(e.accession === r.array && e.raw === r.raw_median && close(e.rma, r.array_median, 1e-4) && r.each.length === r.arrays, `${e.accession}`);
  });
  T_('RNA-seq: log₂ CPM of HLA-DRA in GSM4683511 recomputed from its count, library, TMM factor and prior equals the table', () => {
    const k = M.counts.cpm, r = logCPM(k.count, k.library, k.factor, k.mean_effective, k.prior);
    return ok(close(r.value, k.value, 2e-3) && close(r.zero, k.zero, 2e-3), `${fx(r.value, 4)} against ${fx(k.value, 4)}; zero ${fx(r.zero, 4)} against ${fx(k.zero, 4)}`);
  });
  T_('RNA-seq: the value of a zero count in every normalization group is c/(N̄ + 2c) × 10⁶', () => {
    const c = M.counts, bad = c.zeros.filter(z => !close(zeroOf(z.mean_millions, c.cpm.prior), z.zero, 2e-3));
    return ok(!bad.length, bad.slice(0, 3).map(z => z.study).join(', '));
  });
  T_('RNA-seq: the expressed flag of GSE154918 follows n* and the limit', () => {
    const e = M.counts.expressed, need = e.large_n + (M.counts.cpm.libraries - e.large_n) * e.proportion;
    return ok(e.checked.every(x => x.flagged === (x.samples >= need && x.total >= e.minimum_total)), `n* = ${fx(need, 1)}`);
  });
  T_('genes: HLA-DRA keeps the probe set with the higher mean, and its gene-level value is that row', () => {
    const h = M.genes.hla, [a, b] = h.probe_sets, i = h.arrays.indexOf(h.array);
    return ok(a.feature === h.kept && mean(a.values) > mean(b.values) && close(mean(a.values), a.mean, 1e-3) && close(mean(b.values), b.mean, 1e-3) && close(a.values[i], h.value, 1e-4) && h.arrays.length === M.arrays.rma.arrays);
  });
  T_('genes: the rows of GSE26440 and of the resource add up by status', () => {
    const g = M.genes, f = g.features;
    return ok(sum(Object.values(g.statuses)) === M.arrays.rma.probe_sets && STATUS.reduce((t, k) => t + (f[k] || 0), 0) + (f.not_in_annotation || 0) === f.total);
  });
  T_('labels: the rules of GSE65682, applied here, give every one of its records its stored class and rule', () => {
    const L = M.labels, rules = L.rules.map((r, i) => Object.assign({ n: i + 1 }, r)), bad = L.records.filter(r => { const [v, n] = applyRules(rules, r); return v !== r.value || n !== r.rule; });
    const ex = L.records.find(r => r.accession === L.sample);
    return ok(!bad.length && ex.value === 'sepsis' && ex.rule === 3, bad.slice(0, 3).map(r => r.accession).join(', '));
  });
  T_('repeats: the six deposits share chip, scan time and median, and the primary record is the first released', () => {
    const d = M.repeats.deposits, first = d.slice().sort((a, b) => a.released < b.released ? -1 : a.released > b.released ? 1 : a.record < b.record ? -1 : 1)[0];
    return ok(d.every(x => x.scan_time === d[0].scan_time && x.median === d[0].median && x.chip === d[0].chip) && first.primary && d.filter(x => x.primary).length === 1 && M.repeats.sample_id === `${first.study}:${first.record}`);
  });
  T_('repeats: the preference order marks the stored sample of NS_19', () => {
    const s = M.repeats.patient.samples;
    return ok(markOne(s) === s.find(x => x.marked).record);
  });
  T_('sex: two-means recomputed on the scores of GSE65682 gives its stored split, calls and agreement', () => {
    const X = M.sex, t = X.table, d = X.each.map(r => r.d), run = twoMeans(d), st = splitStats(d, X.each.map(r => r.y), X.each.map(r => r.xist), run.steps[run.steps.length - 1].cut, X.policy.ambiguous_zone);
    const rec = X.each.map(r => r.recorded), agree = st.call.filter((c, i) => c !== 'a' && c === rec[i]).length, dis = st.call.filter((c, i) => c !== 'a' && rec[i] !== 'n' && c !== rec[i]).length;
    const calls = st.call.every((c, i) => c === X.each[i].call);
    return ok(st.L.length === t.lower_group && st.H.length === t.higher_group && close(st.gap, t.gap, 1e-3) && close(st.sep, t.separation, 0.01) && close(st.mid, t.midpoint, 1e-3) && calls && agree === t.concordant && dis === t.discordant,
      `${st.L.length}/${t.lower_group}, separation ${fx(st.sep, 3)}/${fx(t.separation, 3)}, ${agree}/${t.concordant}`);
  });
  T_('sex: d is the mean of the Y genes minus XIST for the two worked samples', () => {
    const e = M.sex.each.filter(r => M.sex.examples.includes(r.accession));
    return ok(e.length === 2 && e.every(r => close(r.y - r.xist, r.d, 2e-4)));
  });
  T_('score: the geometric means of the eleven genes of GSM1602802 give its stored score', () => {
    const S0 = M.score, v = geo(S0.up.map(g => g.value)) - geo(S0.down.map(g => g.value));
    return ok(close(v, S0.score, 1e-3), `${fx(v, 4)} against ${fx(S0.score, 4)}`);
  });
  T_('score: the AUROC counted from the pairs equals the stored AUROC in every stratum, and the area under the ROC curve', () => {
    const bad = M.score.strata.filter(s => { const a = auroc(s.sepsis, s.healthy); const roc = [...rocCurve(s.sepsis, s.healthy), [1, 1]]; let area = 0; for (let i = 1; i < roc.length; i++) area += (roc[i][0] - roc[i - 1][0]) * (roc[i][1] + roc[i - 1][1]) / 2; return !close(a, s.auroc, 1e-6) || !close(area, a, 1e-9); });
    return ok(!bad.length, bad.slice(0, 3).map(s => s.study).join(', '));
  });
  T_('correction: CEACAM1 in GSE54514 recomputed from the reference, the anchor, the shrinkage and the effects gives the corrected table', () => {
    const A = M.correction.anchored, ref = M.correction.reference, out = [];
    A.examples.forEach(x => {
      if (!close(refAt(ref, x.position), x.aligned, 2e-3)) out.push(`${x.sample} aligned`);
      const z = (x.aligned - A.sex_effect * x.female - A.anchor_mean) / A.anchor_sd, v = A.anchor_mean + A.anchor_sd * (z - A.location) / A.scale + A.sex_effect * x.female;
      if (!close(v, x.corrected, 2e-3)) out.push(`${x.sample} corrected ${fx(v, 4)} against ${fx(x.corrected, 4)}`);
      if (!close((x.rank - 0.5) / x.genes, x.position, 1e-6)) out.push(`${x.sample} position`);
    });
    const v = shrinkVar(A.own_variance, A.n, A.prior_shape, A.prior_rate), w = shrinkWeight(A.prior_spread, v, A.n);
    if (!close(Math.sqrt(v), A.scale, 1e-4)) out.push(`scale ${fx(Math.sqrt(v), 5)}`);
    if (!close(w * A.own_location + (1 - w) * A.prior_location, A.location, 1e-4)) out.push('location');
    return ok(!out.length, out.join('; '));
  });
  T_('correction: matched groups of CEACAM1 in E-MTAB-4421 line the patients up with the reference mean', () => {
    const A = M.correction.anchored, m = M.correction.matched, z = (m.mean - A.anchor_mean) / A.anchor_sd, after = (z - m.location) / m.scale;
    return ok(close(m.spread / m.omega, m.variance_raw, 1e-4) && close(z - m.scale * m.rho, m.location_raw, 1e-4) && close(after, m.rho, 0.01), `after ${fx(after, 3)}, reference ${fx(m.rho, 3)}`);
  });
  T_('held-out test: the medians over the units, recomputed, equal the stored medians', () => {
    const Hd = M.heldout, out = [];
    Object.entries(Hd.medians).forEach(([g, md]) => {
      const units = Hd.units.filter(u => u.group === g);
      Hd.methods.forEach((k, j) => { if (!md[k]) return; const v = median(units.map(u => u.error[j]).filter(x => x != null)), dd = median(units.map(u => u.dispersion[j]).filter(x => x != null)); if (!close(v, md[k].error, 1e-4) || !close(dd, md[k].dispersion, 1e-4)) out.push(`${g} ${k}`); });
      if (!close(median(units.map(u => u.floor)), md.floor, 1e-4) || units.length !== md.units) out.push(`${g} floor`);
    });
    return ok(!out.length, out.join('; '));
  });
  T_('benchmark: the fold of G024 and its training groups add up to the task, and the invented pooled example gives 0.90', () => {
    const f = M.bench.fold, t = M.bench.tasks.find(x => x.task === f.task), p = pooledExample(9, 1, 0);
    const ex = M.bench.tasks.every(x => sum(Object.values(x.excluded)) + x.eligible === x.records);
    return ok(f.train_positives + f.test_positives === t.positives && f.train_negatives + f.test_negatives === t.negatives && close(p.pooled, 0.9, 1e-12) && p.auA === 0.5 && p.auB === 0.5 && ex
      && sum(f.platforms.map(x => x[2])) === f.test_positives + f.test_negatives, `pooled ${p.pooled}`);
  });
})();
