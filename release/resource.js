/* The resource: what GenomicTensorSepsis holds, chapter by chapter, as a book with a rail of chapters. Each chapter has a
   lead paragraph, one or two charts drawn by hand with a "How to read it" note and tooltips, and a short note of what
   the charts show. The chapters are written here (chapterHTML) and drawn when shown (DRAW). Their data are in
   data/resource.js (resource() in build_data.py): the points of the paper's figures and the pipeline tables; the
   before-and-after view of the correction reuses the start page's points (data.js, hero). */
'use strict';
(() => {
  let R = null;
  const N = D.numbers;
  // choices made inside the chapters
  const V = {
    year: 'platforms', rec: 'age', genes: 'tables', map: 'share', rep: 'cc', score: 'roc', pick: 'GSE65682|GPL13667|whole_blood',
    corrBy: 'fam', fgroup: 'bulk_blood', hgroup: 'bulk_blood', hmeasure: 'loc', mix: 'mixing',
    task: 'sepsis_vs_healthy_control', scheme: 'leave_one_group_out', fold: 'G024',
    mtask: 'sepsis_vs_non_septic_control', etask: 'sepsis_vs_non_septic_control', sigScope: 'all_tissues',
  };

  // ---------- small helpers ----------
  const howto = t => `<p class="howto"><span>How to read it</span><span>${t}</span></p>`;
  const segBox = (id, options, value) => `<div class="seg" id="${id}">${segHTML(options, value)}</div>`;
  const widthOf = id => Math.max(300, Math.round(($(id) && $(id).clientWidth) || 860));
  const isSmall = w => w < 560;
  const cap = s => s ? s[0].toUpperCase() + s.slice(1) : s;
  const lower = s => /^[A-Z][a-z]/.test(s) ? s[0].toLowerCase() + s.slice(1) : s;
  const trunc = (s, n) => s.length > n ? s.slice(0, n - 1) + '…' : s;
  const CLASS_PLURAL = { sepsis: 'sepsis samples', healthy_control: 'healthy controls', non_septic_control: 'non-septic controls', other_condition: 'samples of other conditions', not_recorded: 'without a recorded class' };
  // text on a filled shape: dark or white, whichever reads better
  const onFill = c => { const [r, g, b] = hexRGB(c); return 0.299 * r + 0.587 * g + 0.114 * b > 150 ? '#1f1b17' : '#ffffff'; };
  // a color of the paper's figures, made readable on the dark page where it is too dark
  const lift = c => C.dark && sum(hexRGB(c)) < 330 ? mixHex(c, '#ffffff', 0.32) : c;
  const css = name => getComputedStyle($('resource')).getPropertyValue(name).trim();
  const plat = (study, platform) => platform ? `${study} ${platform}` : study;
  const tissueLab = t => lab('tissue', t);
  const halo = () => ({ stroke: C.bg, 'stroke-width': 3, 'paint-order': 'stroke', 'stroke-linejoin': 'round' });
  function hatch(s, id, color) {
    const defs = S('defs', {}, s), p = S('pattern', { id, width: 6, height: 6, patternUnits: 'userSpaceOnUse', patternTransform: 'rotate(45)' }, defs);
    S('rect', { width: 6, height: 6, fill: tint(color, 0.78) }, p);
    S('line', { x1: 0, y1: 0, x2: 0, y2: 6, stroke: color, 'stroke-width': 2.2 }, p);
    return `url(#${id})`;
  }
  function sec(id, nav, sub, title, body) {
    return `<section class="sec" id="${id}" data-nav="${nav}" data-sub="${sub}">
  <div class="sec-h"><span class="n"></span><h2>${title}</h2></div>
  ${body}
</section>`;
  }

  // ---------- generic charts ----------
  // Horizontal bars, one row per item, each a stack of parts. rows = [{label, parts: [{label, value, color, fill, tip, text, dash}], end}]
  // opts: share (each row as shares of its total), max, labelW, rowH, axis, tickF, axisLabel, aria
  function hbars(id, rows, opts = {}) {
    const w = widthOf(id), sm = isSmall(w), labelW = sm ? 0 : (opts.labelW == null ? 230 : opts.labelW);
    const rowH = opts.rowH || 22, gap = opts.gap == null ? 10 : opts.gap, lab_ = sm ? 17 : 0, per = rowH + gap + lab_;
    const endW = opts.endW == null ? (sm ? 46 : 70) : opts.endW, top = 4;
    const H = top + rows.length * per + (opts.axis ? 38 : 2);
    const s = svg(id, w, H, opts.aria || '');
    const max = opts.share ? 1 : (opts.max || Math.max(...rows.map(r => sum(r.parts.map(p => p.value)))));
    const x = lin(0, max, labelW, w - endW);
    if (opts.axis) {
      gridX(s, x, ticks(0, max, sm ? 3 : 5), top - 2, H - 34, opts.share ? (v => pct(v, 0)) : (opts.tickF || fmt));
      if (opts.axisLabel) axisTitle(s, (labelW + w - endW) / 2, H - 4, opts.axisLabel);
    }
    rows.forEach((r, i) => {
      const y = top + i * per + lab_, tot = sum(r.parts.map(p => p.value));
      if (sm) T(s, 0, y - 5, r.label, { 'font-size': 12, fill: C.ink, 'font-weight': r.bold ? 600 : 400 });
      else T(s, labelW - 8, y + rowH / 2 + 4, r.label, { 'text-anchor': 'end', 'font-size': 12.5, fill: C.ink, 'font-weight': r.bold ? 600 : 400 });
      let cx = labelW;
      r.parts.forEach(p => {
        const v = opts.share ? (tot ? p.value / tot : 0) : p.value, pw = x(v) - labelW;
        if (p.value <= 0) return;
        const rect = S('rect', { x: cx, y, width: Math.max(pw, 1), height: rowH, fill: p.fill || p.color, stroke: p.dash ? C.lineStrong : C.bg,
          'stroke-width': 1, 'stroke-dasharray': p.dash ? '3 2' : null }, s);
        if (p.text && pw > p.text.length * 6.6 + 8) T(s, cx + 6, y + rowH / 2 + 4, p.text, { 'font-size': 11, fill: p.fill ? C.ink : onFill(p.color), 'pointer-events': 'none' });
        if (p.tip) hover(rect, p.tip);
        cx += pw;
      });
      if (r.end) T(s, cx + 6, y + rowH / 2 + 4, r.end, { 'font-size': 11.5, fill: C.ink2 });
    });
    return s;
  }
  // A flow diagram (alluvial). cols = [[{key, label, color}]], links = [{c, s, t, v, color, tip}], c the column of the source.
  function sankey(id, cols, links, opts = {}) {
    const w = widthOf(id), sm = isSmall(w), H = opts.height || (sm ? 360 : 420), nodeW = sm ? 9 : 12, pad = opts.pad || (sm ? 7 : 10);
    const lw = sm ? (opts.leftS || 84) : (opts.left || 150), rw = sm ? (opts.rightS || 92) : (opts.right || 170);
    const s = svg(id, w, H, opts.aria || '');
    const xs = cols.map((_, i) => lw + i * ((w - lw - rw - nodeW) / Math.max(1, cols.length - 1)));
    const inOf = (i, k) => sum(links.filter(l => l.c === i - 1 && l.t === k).map(l => l.v));
    const outOf = (i, k) => sum(links.filter(l => l.c === i && l.s === k).map(l => l.v));
    const nodes = cols.map((col, i) => col.map(n => ({ ...n, i, v: Math.max(inOf(i, n.key), outOf(i, n.key)) })).filter(n => n.v > 0));
    const k = Math.min(...nodes.map(col => (H - 12 - pad * (col.length - 1)) / sum(col.map(n => n.v))));
    nodes.forEach(col => {
      const used = sum(col.map(n => n.v)) * k + pad * (col.length - 1); let y = (H - used) / 2;
      col.forEach(n => { n.y = y; n.h = Math.max(1, n.v * k); n.out = y; n.in = y; y += n.v * k + pad; });
    });
    const find = (i, key) => nodes[i].find(n => n.key === key);
    const order = (i, key) => nodes[i].findIndex(n => n.key === key);
    const sorted = links.filter(l => l.v > 0).slice().sort((a, b) => a.c - b.c || order(a.c, a.s) - order(b.c, b.s) || order(a.c + 1, a.t) - order(b.c + 1, b.t));
    const bands = [];
    sorted.forEach(l => {
      const a = find(l.c, l.s), b = find(l.c + 1, l.t); if (!a || !b) return;
      const h = l.v * k, x0 = xs[l.c] + nodeW, x1 = xs[l.c + 1], xm = (x0 + x1) / 2, y0 = a.out, y1 = b.in;
      a.out += h; b.in += h;
      const p = S('path', { d: `M${x0},${y0}C${xm},${y0} ${xm},${y1} ${x1},${y1}L${x1},${y1 + h}C${xm},${y1 + h} ${xm},${y0 + h} ${x0},${y0 + h}Z`,
        fill: l.color, 'fill-opacity': 0.42 }, s);
      p.dataset.s = `${l.c}:${l.s}`; p.dataset.t = `${l.c + 1}:${l.t}`;
      hover(p, l.tip || `<b>${a.label} → ${b.label}</b><br>${fmt(l.v)}`);
      p.addEventListener('mouseenter', () => p.setAttribute('fill-opacity', 0.75));
      p.addEventListener('mouseleave', () => p.setAttribute('fill-opacity', 0.42));
      bands.push(p);
    });
    nodes.forEach((col, i) => col.forEach(n => {
      const r = rrect(s, xs[i], n.y, nodeW, n.h, n.color, 2);
      const key = `${i}:${n.key}`, mine = p => p.dataset.s === key || p.dataset.t === key;
      hover(r, n.tip || `<b>${n.label}</b><br>${fmt(n.v)}${opts.unit ? ' ' + opts.unit : ''}`);
      r.addEventListener('mouseenter', () => bands.forEach(p => p.setAttribute('fill-opacity', mine(p) ? 0.75 : 0.1)));
      r.addEventListener('mouseleave', () => bands.forEach(p => p.setAttribute('fill-opacity', 0.42)));
      const last = i === cols.length - 1, first = i === 0, size = sm ? 10.5 : 12;
      const lx = first ? xs[i] - 6 : xs[i] + nodeW + 6, anchor = first ? 'end' : 'start';
      const label = sm && n.short ? n.short : n.label;
      if (n.h >= 7 || !opts.hideSmall) {
        T(s, lx, n.y + n.h / 2 + 4, sm && !first && !last ? label : `${label} ${fmt(n.v)}`, { 'text-anchor': anchor, 'font-size': size, fill: C.ink, ...((!first && !last) ? halo() : {}) });
      }
    }));
    return s;
  }

  // ======================================================================= the chapters
  function chapterHTML() {
    return [
      sec('r-studies', 'The studies', 'Searches and screening', 'The studies: searches and screening', `
  <p class="lead" id="rsLead"></p>
  <h3>From the searches to the included studies</h3>
  ${howto('Each bar is one step, as long as everything the step started with. The green part goes on to the next step; the gray parts were set aside, each for the first reason that applied; the hatched part met the criteria but had no data that could be used. Hover over a part to see its count.')}
  <div class="chart" id="rsFunnel"></div>
  <h3>Why candidates were excluded</h3>
  <div class="chart" id="rsReasons"></div>
  <p class="note" id="rsReasonsNote"></p>
  <h3>Studies by year of release and technology</h3>
  <div class="controls"><span>Show</span>${segBox('rsYearBy', [['platforms', 'Study platforms'], ['samples', 'Distinct samples'], ['added', 'Distinct samples, added up']], V.year)}</div>
  <div class="keys" id="rsYearKeys"></div>
  <div class="chart" id="rsYears"></div>
  ${howto('One bar for every year of release, stacked by technology. A study platform is the part of a study measured on one platform, so a study measured on two platforms counts twice; the samples count once each, in the year of their study platform. The last year runs only to the date of the searches. Hover over a bar to see its numbers.')}
  <p class="note" id="rsNote"></p>`),
      sec('r-samples', 'The samples', 'Classes, sample types, records', 'The samples and what their records say', `
  <p class="lead" id="raLead"></p>
  <h3>From technology to sample type to class</h3>
  ${howto('A flow diagram: the bars on the left are the technologies, those in the middle the material sampled and those on the right the curated classes. Every band is as wide as the distinct samples it carries, so the bands leaving a bar add up to the bar. Hover over a band or a bar to see its count.')}
  <div class="chart" id="raFlow"></div>
  <h3>What the records say</h3>
  <div class="controls"><span>Show</span>${segBox('raRecBy', [['age', 'Age'], ['clinical', 'Shock and outcome'], ['fields', 'Fields recorded'], ['basis', 'What a class rests on']], V.rec)}</div>
  <div class="keys" id="raRecKeys"></div>
  <div class="chart" id="raRec"></div>
  <p class="howto" id="raRecHow"><span>How to read it</span><span></span></p>
  <div class="exblock"><h4 class="exhead">Worked example: how the classes of GSE65682 were read <span class="m" id="raExMeta"></span></h4>
  <p class="readout" id="raExText"></p>
  <div class="tablewrap"><table class="t" id="raExTable"></table></div></div>
  <p class="note" id="raNote"></p>`),
      sec('r-repeats', 'Repeats and patients', 'Shared arrays, one per patient', 'Repeat deposits and patients', `
  <p class="lead" id="rrLead"></p>
  <h3>Study platforms that hold the same arrays or libraries</h3>
  ${howto('Every dot on the line is a study platform, colored by technology and grouped by cohort group. An arc joins two study platforms that hold records of the same physical arrays or libraries; it is as thick as the samples they share. Hover over an arc or a dot to see the numbers; click a dot to open its study.')}
  <div class="chart" id="rrArcs"></div>
  <h3>From records to samples to patients</h3>
  ${howto('Every bar is one count, drawn to one scale. The solid part of a bar goes on to the next bar; the dashed part is set aside in the count but kept in the resource and marked. Hover over a part to see its count.')}
  <div class="chart" id="rrSteps"></div>
  <div class="exblock"><h4 class="exhead">Worked example: the cohort group of GSE26440 <span class="m" id="rrExMeta"></span></h4>
  <p class="readout" id="rrExText"></p>
  <div class="chart" id="rrExample"></div></div>
  <p class="note" id="rrNote"></p>`),
      sec('r-genes', 'From probes to genes', 'Mapping, genes, scales', 'From probes to genes', `
  <p class="lead" id="rgLead"></p>
  <h3>How the rows of the tables were mapped to genes</h3>
  <div class="controls"><span>Show</span>${segBox('rgMapBy', [['share', 'Share of the rows'], ['rows', 'Number of rows']], V.map)}</div>
  <div class="keys" id="rgMapKeys"></div>
  <div class="chart" id="rgMap"></div>
  ${howto('One bar for every technology, split by what happened to its rows (probes, probe sets or gene identifiers) when they were mapped to HGNC genes. Only rows that map to one gene enter the gene-level tables. Hover over a part to see its count.')}
  <h3>The gene-level tables</h3>
  <div class="controls"><span>Show</span>${segBox('rgBy', [['tables', 'Genes per table'], ['shared', 'Genes shared by the tables'], ['scales', 'Value scales']], V.genes)}</div>
  <div class="keys" id="rgKeys"></div>
  <div class="chart" id="rgChart"></div>
  <p class="howto" id="rgHow"><span>How to read it</span><span></span></p>
  <p class="note" id="rgNote"></p>`),
      sec('r-checks', 'The checks', 'Sex, repeats, sepsis score', 'The checks inside single studies', `
  <p class="lead" id="rcLead"></p>
  <h3>Sex read from the genes against the recorded sex</h3>
  ${howto('In the scatter plot, every sample of a table that measures XIST is placed by its own table. Across, XIST from the mean of the samples read as male (0) to that of the samples read as female (1); up, the Y-chromosome genes from the female mean (0) to the male mean (1). Women should sit at the lower right and men at the upper left. The color is the recorded sex; a ring marks a recorded sex that disagrees with the genes. The bars give the share of recorded sexes that agree in every table with a recorded sex; orange bars are tables whose recorded sex is marked unreliable. Hover for details.')}
  <div class="keys" id="rcSexKeys"></div>
  <div class="grid2 r-grid"><div class="chart" id="rcSexPlane"></div><div class="chart" id="rcSexTables"></div></div>
  <p class="note" id="rcSexNote"></p>
  <h3>Two deposits of one sample</h3>
  <div class="controls"><span>Correlation</span>${segBox('rcRepBy', [['cc', 'After centering every gene'], ['sp', 'Spearman, as deposited']], V.rep)}</div>
  ${howto('Every dot is one pair of records of the same array or library, in the row of the evidence that links them, placed at the correlation of their two gene profiles. The bar shows the median and the middle half of the pairs. Centering every gene on its median in its own table removes the scale of each table, so that only the shape of the profile is compared. A cross marks a record that is not the closest match of its twin in the other table. Hover over a dot to see the pair.')}
  <div class="chart" id="rcRepeats"></div>
  <p class="note" id="rcRepNote"></p>
  <h3>The published sepsis genes and the sepsis score</h3>
  <div class="controls"><span>Show</span>${segBox('rcScoreBy', [['roc', 'ROC curves and AUROC'], ['genes', 'Sepsis genes, study by study']], V.score)}</div>
  <p class="howto" id="rcScoreHow"><span>How to read it</span><span></span></p>
  <div class="grid2 r-grid" id="rcScoreBox"><div class="chart" id="rcRoc"></div><div class="chart" id="rcAuroc"></div></div>
  <div class="chart" id="rcGenes"></div>
  <p class="readout" id="rcScorePick"></p>
  <p class="note" id="rcScoreNote"></p>
  <div class="result" id="rcResult"></div>`),
      sec('r-correction', 'The correction', 'Studies on one scale', 'The correction of the differences between studies', `
  <p class="lead" id="rxLead"></p>
  <h3>Before and after the correction (bulk blood)</h3>
  <div class="controls"><span>Color by</span>${segBox('rxBy', [['fam', 'Technology'], ['cls', 'Class'], ['tier', 'How the study was corrected']], V.corrBy)}</div>
  <div class="keys" id="rxKeys"></div>
  <div class="grid2 r-grid"><div class="chart"><canvas id="rxBefore" role="img" aria-label=""></canvas></div><div class="chart"><canvas id="rxAfter" role="img" aria-label=""></canvas></div></div>
  ${howto('Every dot is one sample of the corrected bulk-blood table, placed by the first two principal components of its most variable genes, before the correction and after it; the components are computed separately for each. Dots close together have similar profiles. Hover over a dot to see its sample.')}
  <h3>What the first two components follow</h3>
  <div class="controls"><span>Table</span>${segBox('rxFactorBy', [['bulk_blood', 'Bulk blood'], ['pbmc', 'PBMC']], V.fgroup)}</div>
  <div class="chart" id="rxFactors"></div>
  ${howto('For each factor, the share of the variation of a component that the factor explains: the open circle before the correction, the filled one after it, joined by an arrow, blue where the share falls and red where it rises. The last row counts the class over the studies corrected on their healthy controls alone, whose sepsis samples took no part in their correction.')}
  <h3>The held-out test</h3>
  <div class="controls"><span>Table</span>${segBox('rxHeldBy', [['bulk_blood', 'Bulk blood'], ['pbmc', 'PBMC']], V.hgroup)}<span>Measure</span>${segBox('rxHeldMeasure', [['loc', 'Location error'], ['scale', 'Scale error']], V.hmeasure)}</div>
  <p class="howto" id="rxHeldHow"><span>How to read it</span><span></span></p>
  <div class="chart" id="rxHeld"></div>
  <p class="readout" id="rxHeldText"></p>
  <h3>The checks of every corrected study, and where every record went</h3>
  <div class="controls"><span>Show</span>${segBox('rxMixBy', [['mixing', 'Study mixing'], ['ledger', 'Where every record went']], V.mix)}</div>
  <div class="keys" id="rxMixKeys"></div>
  <div class="chart" id="rxMix"></div>
  <p class="howto" id="rxMixHow"><span>How to read it</span><span></span></p>
  <p class="note" id="rxMixNote"></p>
  <div class="result" id="rxCautions"></div>`),
      sec('r-benchmark', 'Tasks and splits', 'Tasks, eligibility, folds', 'The benchmark tasks and splits', `
  <p class="lead" id="rbLead"></p>
  <div class="controls"><span>Task</span><div class="seg" id="rbTask"></div></div>
  <div class="cards" id="rbCards"></div>
  <h3>Which samples are eligible</h3>
  ${howto('Every record of the resource starts at the top. Each bar is what remains before a rule, and its dashed part is what the rule excludes, in the order in which the rules apply; the last bar splits the samples of the task into its positive and negative class. Hover over a part to see its count.')}
  <div class="chart" id="rbSteps"></div>
  <h3>The folds: which groups train and which are tested</h3>
  <div class="controls"><span>Scheme</span><div class="seg" id="rbScheme"></div><label id="rbFoldBox">Fold <select id="rbFold"></select></label></div>
  <div class="keys" id="rbKeys"></div>
  <div class="chart" id="rbFolds"></div>
  ${howto('Each block is one patient-safe group, as wide as its samples in the task: the darker part its positive class, the lighter part its negative class. A fold puts every group on one side only: it trains, it is tested, or it is held out but not scored, because it has fewer than the minimum samples of each class. Hover over a block to see the group and its studies.')}
  <p class="readout" id="rbFoldText"></p>
  <div class="result" id="rbWhy"></div>`),
      sec('r-evaluation', 'Benchmark', 'Twins, classifiers, signatures', 'What the benchmark measures', `
  <p class="lead" id="reLead"></p>
  <div class="cards" id="reCards"></div>
  <h3>Test records with a twin in training</h3>
  ${howto('One bar for every split of the records of the resource and of each task: the share of its test records that have a twin in the training data, by the kind of twin. The tick on a random split marks the share that chance predicts from the number of twins of each record. Hover over a bar to see its counts.')}
  <div class="keys" id="reLeakKeys"></div>
  <div class="chart" id="reLeak"></div>
  <p class="note" id="reLeakNote"></p>
  <h3>Baseline classifiers under each split</h3>
  <div class="controls"><span>Task</span><div class="seg" id="reModelTask"></div></div>
  <div class="keys" id="reModelKeys"></div>
  <div class="chart" id="reModels"></div>
  <p class="howto" id="reModelHow"><span>How to read it</span><span></span></p>
  <p class="readout" id="reModelText"></p>
  <h3>Where the higher scores come from</h3>
  <div class="controls"><span>Task</span><div class="seg" id="reEffectTask"></div></div>
  <div class="keys" id="reEffectKeys"></div>
  <div class="chart" id="reEffects"></div>
  <p class="howto" id="reEffectHow"><span>How to read it</span><span></span></p>
  <p class="readout" id="reEffectText"></p>
  <h3>What the study alone gives away</h3>
  <div class="grid2 r-grid"><div class="chart" id="reStudyOnly"></div><div class="chart" id="reIdent"></div></div>
  ${howto('On the left, a predictor that knows nothing but the study of a record, scored over the pooled test records of every random split and within every held-out study. On the right, how often the study platform of a sample is named correctly from its expression, from the values the classifiers read and from the corrected table, against guessing. Each dot is one split or one fold; hover for its numbers.')}
  <p class="note" id="reShortNote"></p>
  <h3>Published signatures, tested again within single studies</h3>
  <div class="controls"><span>Sample types</span>${segBox('reSigScope', [['all_tissues', 'All'], ['bulk_blood', 'Bulk blood']], V.sigScope)}</div>
  <div class="keys" id="reSigKeys"></div>
  <div class="chart" id="reSigs"></div>
  ${howto('One row for every signature and the comparison it was published for. Each dot is the AUROC of the signature in one study platform and sample type: filled where the stratum is independent of the studies that derived or validated the signature, hollow where it is not. The two ticks are the medians over all strata and over the independent ones. Hover over a dot for its study.')}
  <p class="readout" id="reSigText"></p>`),
      sec('r-data', 'Getting the data', 'Files, license, deposit', 'Getting the data', `
  <p class="lead" id="rdLead"></p>
  <div class="cards" id="rdCards"></div>
  <h3>The files of the data package</h3>
  ${howto('One bar for every file of the package, on a logarithmic scale: each grid line is a thousand times the one before (kB, MB, GB). An archive holds many tables, each a compressed CSV file. Hover over a bar for its contents.')}
  <div class="chart" id="rdFiles"></div>
  <h3>Which table answers which question</h3>
  <div class="tablewrap"><table class="t r-qa" id="rdQuestions"></table></div>
  <p class="note" id="rdNote"></p>`),
    ].join('\n');
  }

  // ======================================================================= 1. the studies
  function drawStudies() {
    const St = R.studies, src = N.studies_by_source;
    $('rsLead').innerHTML = `The searches of three public repositories, GEO, ArrayExpress and the European Nucleotide Archive (ENA), were run on ${fdate(St.searched)}, with terms for sepsis, septic shock, bacteremia, systemic inflammation and critical illness. They returned ${fmt(St.records_found)} records, which became ${fmt(St.candidates)} candidate studies once a study found in two repositories was counted once. Every candidate was screened against the eligibility rules: at least one transcriptome-wide profile of whole blood, leukocytes or peripheral blood mononuclear cells (PBMC), measured by microarray or RNA-seq, from a patient whom the study diagnosed with sepsis. ${fmt(St.included)} studies were included, ${fmt(src.geo)} from GEO and ${fmt(src.arrayexpress)} from ArrayExpress.`;
    drawFunnel(); drawReasons(); drawYears();
    const unav = St.excluded.filter(r => St.unavailable.includes(r[0])), sbY = R.studies.studies_by_year;
    const last5 = sbY.filter(r => r.year > N.last_year - 5), rs5 = sum(last5.map(r => r.rna_seq)), all5 = sum(last5.map(r => r.array + r.rna_seq));
    $('rsReasonsNote').textContent = `Each excluded candidate is counted under the first rule it failed. ${Spell(sum(unav.map(r => r[1])))} candidates met the criteria but no expression values could be obtained: ${andList(unav.map(([k, n]) => `${lower(St.reasons[k])} (${fmt(n)})`))}. They are listed, with every other excluded candidate and its reason, in the Studies tab.`;
    $('rsNote').textContent = `The ${fmt(St.included)} included studies were released from ${N.first_year} to ${N.last_year} and split into ${fmt(St.study_platforms)} study platforms, ${fmt(St.study_platforms_with_values)} of them with expression values: ${fmt(N.studies_by_technology.array)} studies measured RNA on microarrays and ${fmt(N.studies_by_technology.rna_seq)} by RNA-seq. Of the ${fmt(all5)} studies released from ${N.last_year - 4} on, ${fmt(rs5)} used RNA-seq.`;
  }
  function drawFunnel() {
    const St = R.studies, id = 'rsFunnel', w = widthOf(id), sm = isSmall(w);
    const labelW = sm ? 0 : 230, x0 = labelW, x1 = w - 6, barH = 30, zone = sm ? 54 : 44, labH = sm ? 34 : 0;
    const rows = [
      { big: St.records_found, cap: 'records returned by the searches', keep: St.candidates, keepCap: 'candidate studies',
        parts: St.not_candidate.map(([k, n]) => ({ k, n })) },
      { big: St.candidates, cap: 'candidate studies, screened', keep: St.included, keepCap: 'included studies',
        parts: [...St.excluded.filter(r => St.unavailable.includes(r[0])), ...St.excluded.filter(r => !St.unavailable.includes(r[0]))].map(([k, n]) => ({ k, n, unav: St.unavailable.includes(k) })) },
      { big: St.included, cap: 'included studies', parts: [] },
    ];
    const H = rows.length * (labH + barH) + (rows.length - 1) * zone + 8;
    const s = svg(id, w, H, `Funnel from ${fmt(St.records_found)} records to ${fmt(St.candidates)} candidate studies to ${fmt(St.included)} included studies, with the reasons for setting the others aside.`);
    const fill = hatch(s, 'rsHatch', C.muted), greys = [tint(C.muted, 0.5), tint(C.muted, 0.68)];
    let y = 4;
    rows.forEach((r, i) => {
      if (sm) { T(s, 0, y + 14, `${fmt(r.big)} ${r.cap}`, { 'font-size': 13, fill: C.ink, 'font-weight': 600 }); y += labH; }
      else {
        T(s, labelW - 12, y + 14, fmt(r.big), { 'text-anchor': 'end', 'font-size': 17, fill: C.ink, 'font-weight': 600, style: 'font-family: var(--serif)' });
        T(s, labelW - 12, y + 28, r.cap, { 'text-anchor': 'end', 'font-size': 12, fill: C.muted });
      }
      const x = lin(0, r.big, x0, x1);
      if (r.keep == null) {
        // the included studies: by repository, and by technology below
        const geo = N.studies_by_source.geo, ae = N.studies_by_source.arrayexpress;
        [[`GEO ${fmt(geo)}`, geo, C.accentInk, 0], [`ArrayExpress ${fmt(ae)}`, ae, C.accent, geo]].forEach(([t, n, col, at]) => {
          const rr = rrect(s, x(at), y, x(n) - x0, barH, col, 0, { stroke: C.bg });
          hover(rr, `<b>${fmt(n)} included studies</b> from ${t.split(' ')[0]}`);
          if (x(n) - x0 > t.length * 7) T(s, x(at) + 8, y + barH / 2 + 4, t, { 'font-size': 12, fill: onFill(col), 'pointer-events': 'none' });
        });
        y += barH; return;
      }
      const kw = x(r.keep) - x0;
      const kept = rrect(s, x0, y, kw, barH, C.accent, 0, { stroke: C.bg });
      hover(kept, `<b>${fmt(r.keep)} ${r.keepCap}</b><br>go on to the next step`);
      if (kw > 40) T(s, x0 + 6, y + barH / 2 + 4, fmt(r.keep), { 'font-size': 12, fill: onFill(C.accent), 'font-weight': 600, 'pointer-events': 'none' });
      let cx = x0 + kw;
      r.parts.forEach((p, j) => {
        const pw = Math.max(1, x(p.n) - x0), col = p.unav ? fill : greys[j % 2];
        const rr = S('rect', { x: cx, y, width: pw, height: barH, fill: col, stroke: C.bg, 'stroke-width': 1 }, s);
        hover(rr, `<b>${St.reasons[p.k]}</b><br>${fmt(p.n)} ${i === 0 ? 'records' : 'candidates'}${p.unav ? '<br><span class="m">met the criteria; no data that could be used</span>' : ''}`);
        if (pw > 120 && !p.unav) T(s, cx + 6, y + barH / 2 + 4, trunc(`${fmt(p.n)} ${lower(St.reasons[p.k])}`, Math.floor(pw / 6.4)), { 'font-size': 11, fill: C.ink, 'pointer-events': 'none' });
        cx += pw;
      });
      // the kept part opens out to the next bar
      const yb = y + barH, yt = yb + zone + (sm ? labH : 0);
      S('path', { d: `M${x0},${yb}L${x0 + kw},${yb}L${x1},${yt}L${x0},${yt}Z`, fill: C.accent, 'fill-opacity': 0.16 }, s);
      y = yb + zone;
    });
  }
  function drawReasons() {
    const St = R.studies;
    const rows = St.excluded.map(([k, n]) => ({ label: St.reasons[k], parts: [{ value: n, color: St.unavailable.includes(k) ? C.muted : tint(C.muted, 0.45), k,
      tip: `<b>${St.reasons[k]}</b><br>${fmt(n)} candidate ${n === 1 ? 'study' : 'studies'}${St.unavailable.includes(k) ? '<br><span class="m">met the criteria; no data that could be used</span>' : ''}` }], end: fmt(n) }));
    const s = hbars('rsReasons', rows, { labelW: 290, rowH: 16, gap: 6, axis: true, axisLabel: 'Candidate studies excluded', aria: 'Bar chart of the candidate studies excluded at screening, by the first reason that applied.' });
    const pat = hatch(s, 'rsHatch2', C.muted);
    s.querySelectorAll('rect').forEach(r => { if (r.getAttribute('fill') === C.muted) r.setAttribute('fill', pat); });
  }
  function drawYears() {
    const St = R.studies, id = 'rsYears', w = widthOf(id), sm = isSmall(w), fams = D.order.family;
    const years = [...new Set(St.by_year.map(r => r.year))].sort((a, b) => a - b);
    const val = {};
    years.forEach(y => { val[y] = {}; fams.forEach(f => { val[y][f] = 0; }); });
    const names = {};
    if (V.year === 'platforms') St.platforms.forEach(p => { val[p.year][p.family]++; (names[p.year] = names[p.year] || []).push(p.label); });
    else St.by_year.forEach(r => { val[r.year][r.family] = V.year === 'samples' ? r.samples : r.added_up; });
    const totals = years.map(y => sum(fams.map(f => val[y][f]))), vmax = niceMax(Math.max(...totals));
    const m = { l: sm ? 40 : 52, r: 8, t: 10, b: 36 }, H = sm ? 260 : 300;
    const s = svg(id, w, H, `Stacked bars of the ${V.year === 'platforms' ? 'study platforms' : 'distinct samples'} by year of release and technology.`);
    const x = lin(years[0] - 0.5, years[years.length - 1] + 0.5, m.l, w - m.r), y = lin(0, vmax, H - m.b, m.t), bw = (x(1) - x(0)) * 0.72;
    gridY(s, y, ticks(0, vmax, 4), m.l, w - m.r);
    const every = sm ? 5 : w < 800 ? 2 : 1;
    years.forEach((yr, i) => {
      if ((yr - years[0]) % every === 0 || i === years.length - 1) T(s, x(yr), H - m.b + 15, String(yr), { 'text-anchor': 'middle', 'font-size': 11, fill: C.muted });
      let base = 0;
      fams.forEach(f => {
        const v = val[yr][f]; if (!v) return;
        rrect(s, x(yr) - bw / 2, y(base + v), bw, y(base) - y(base + v), lift(C.fam[f]), 1.5, { stroke: C.bg, 'stroke-width': 0.6 });
        base += v;
      });
      const lines = fams.filter(f => val[yr][f] > 0).map(f => `${sw(lift(C.fam[f]))}${lab('family', f)}: ${fmt(val[yr][f])}`);
      const what = V.year === 'platforms' ? 'study platforms' : V.year === 'samples' ? 'distinct samples' : 'distinct samples, added up';
      const list = names[yr] ? `<br><span class="m">${esc(trunc(names[yr].join(', '), 220))}</span>` : '';
      hitRow(s, x(yr) - (x(1) - x(0)) / 2, m.t, x(1) - x(0), H - m.b - m.t,
        `<b>${yr}${yr === years[years.length - 1] ? `, to ${fdate(St.searched)}` : ''}</b>: ${fmt(totals[i])} ${what}${lines.length ? '<br>' + lines.join('<br>') : ''}${list}`);
    });
    axisTitle(s, sm ? 10 : 14, (m.t + H - m.b) / 2, V.year === 'platforms' ? 'Study platforms' : 'Distinct samples', true);
    $('rsYearKeys').innerHTML = keys(fams.map(f => [lab('family', f), lift(C.fam[f])]));
  }

  // ======================================================================= 2. the samples
  function drawSamples() {
    const A = R.samples, cls = D.order.class.filter(k => N.classes[k] > 0);
    $('raLead').innerHTML = `The included studies hold ${fmt(N.included_records)} sample records. The ${fmt(N.records)} with expression values make up the resource; the other ${fmt(N.records_without_values)} are listed with their reason in the data package. Some arrays and libraries were deposited in several studies, so the records with values are ${fmt(N.distinct)} distinct samples: ${andList(cls.map(k => `${fmt(N.classes[k])} ${CLASS_PLURAL[k]}`))}. Every sample was described from its records: its class, the material sampled, its age, sex and time point, and for patients the outcome and septic shock. A label follows the recorded wording. A sample is sepsis only where its record, or a quoted sentence of its study record, says sepsis, severe sepsis or septic shock; a control is a healthy control only where it is called healthy or normal, and other controls are non-septic controls. Diagnoses were not re-adjudicated, and a value that the records do not give is marked as not recorded.`;
    drawSampleFlow(); drawRecords(); drawLabelExample();
    const outcome = sum(A.clinical.filter(r => r[1] !== 'not_recorded').map(r => r[2])), shock = sum(A.clinical.filter(r => r[0] !== 'not_recorded').map(r => r[2]));
    const noAge = sum(Object.values(A.age_groups.not_recorded || {}));
    $('raNote').textContent = `The records leave gaps that limit some analyses: the outcome is recorded for ${fmt(outcome)} of the ${fmt(N.classes.sepsis)} distinct sepsis samples and septic shock for ${fmt(shock)}, the age of ${fmt(noAge)} samples is not recorded, and only ${spell(A.sepsis_3_studies)} studies state that they used the Sepsis-3 definition of sepsis.`;
  }
  function drawSampleFlow() {
    const A = R.samples, fams = D.order.family, cls = D.order.class;
    const SHORT = { affymetrix: 'Affymetrix', illumina_beadchip: 'Illumina', agilent: 'Agilent', other_array: 'Other array', rna_seq: 'RNA-seq',
      healthy_control: 'Healthy', non_septic_control: 'Non-septic', other_condition: 'Other', not_recorded: 'Not recorded', sepsis: 'Sepsis' };
    const cols = [fams.map(f => ({ key: f, label: lab('family', f), short: SHORT[f], color: lift(C.fam[f]) })),
      A.materials.map(m => ({ key: m, label: m, short: m.replace('Blood, unspecified', 'Blood, unsp.').replace('Nasopharyngeal swab', 'Swab'), color: mixHex(C.ink2, C.bg, 0.25) })),
      cls.map(k => ({ key: k, label: lab('class', k), short: SHORT[k], color: C.cls[k] }))];
    const links = [...A.tech_material.map(([f, m, v]) => ({ c: 0, s: f, t: m, v, color: lift(C.fam[f]), tip: `<b>${lab('family', f)} → ${m}</b><br>${fmt(v)} distinct samples` })),
      ...A.material_class.map(([m, k, v]) => ({ c: 1, s: m, t: k, v, color: C.cls[k], tip: `<b>${m} → ${lab('class', k)}</b><br>${fmt(v)} distinct samples` }))];
    sankey('raFlow', cols, links, { height: isSmall(widthOf('raFlow')) ? 400 : 460, unit: 'distinct samples', left: 140, right: 190, leftS: 74, rightS: 96,
      aria: `Flow diagram of the ${fmt(N.distinct)} distinct samples from their technology to the material sampled to their curated class.` });
  }
  function drawRecords() {
    const A = R.samples, how = $('raRecHow').lastElementChild, keysEl = $('raRecKeys'), cls = D.order.class;
    keysEl.innerHTML = '';
    if (V.rec === 'age') {
      const groups = ['adult', 'child', 'neonate', 'not_recorded'];
      const rows = groups.map(g => {
        const m = A.age_groups[g] || {}, tot = sum(Object.values(m));
        return { label: `${lab('age', g)}`, end: fmt(tot), parts: cls.filter(k => m[k]).map(k => ({ value: m[k], color: C.cls[k],
          tip: `<b>${lab('age', g)}, ${labLower('class', k)}</b><br>${fmt(m[k])} distinct samples (${pct(m[k] / tot)} of the age group)` })) };
      });
      const box = $('raRec'); box.innerHTML = '<div id="raRecA"></div><div id="raRecB"></div>';
      hbars('raRecA', rows, { labelW: 130, rowH: 22, axis: true, axisLabel: 'Distinct samples', aria: 'Stacked bars of the distinct samples by age group and class.' });
      drawAgeBoxes('raRecB');
      keysEl.innerHTML = keys(cls.filter(k => N.classes[k] > 0).map(k => [lab('class', k), C.cls[k]]));
      how.textContent = `Above, the distinct samples of every age group, split by class, as the records give them. Below, the recorded age in years of the samples of each class with an age: the box spans the middle half of the samples, the tick inside it is the median and the line runs from the youngest to the oldest. Hover over a bar or a box to see its numbers.`;
    } else if (V.rec === 'clinical') {
      const shock = [['yes', 'In septic shock'], ['no', 'Not in shock'], ['not_recorded', 'Shock not recorded']];
      const outc = [['survivor', 'Survived', C.ink2], ['non_survivor', 'Died', C.cls.sepsis], ['not_recorded', 'Outcome not recorded', C.cls.not_recorded]];
      const cols = [shock.map(([k, l]) => ({ key: k, label: l, color: k === 'not_recorded' ? C.cls.not_recorded : C.ink })),
        outc.map(([k, l, c]) => ({ key: k, label: l, color: c }))];
      const oc = Object.fromEntries(outc.map(([k, l, c]) => [k, [l, c]])), sl = Object.fromEntries(shock);
      const links = A.clinical.map(([sk, o, v]) => ({ c: 0, s: sk, t: o, v, color: oc[o][1], tip: `<b>${sl[sk]} → ${lower(oc[o][0])}</b><br>${fmt(v)} distinct sepsis samples` }));
      sankey('raRec', cols, links, { height: 300, unit: 'distinct sepsis samples', left: 160, right: 200, leftS: 110, rightS: 120,
        aria: 'Flow diagram of the distinct sepsis samples from their recorded septic shock to their recorded outcome.' });
      keysEl.innerHTML = keys(outc.map(([, l, c]) => [l, c]));
      how.textContent = `A flow diagram of the ${fmt(N.classes.sepsis)} distinct sepsis samples: on the left their recorded septic shock, on the right their recorded outcome, each band as wide as its samples. The gray bands are the point of the chart: many records state neither shock nor outcome. Outcomes follow each study's own definition, such as death within 28 days.`;
    } else if (V.rec === 'fields') {
      drawCoverage('raRec');
      how.textContent = `One column for every study platform, grouped by technology (the strip on top), and one row for every curated field, from the most recorded down. The shade is the share of the platform's sample records that record the field, from pale to dark; an empty cell means that no record of the platform gives it. On the right, the share of all ${fmt(N.records)} sample records that record the field. Hover over a cell to see its study platform.`;
      keysEl.innerHTML = keys(D.order.family.map(f => [lab('family', f), lift(C.fam[f])]));
    } else {
      const B = A.basis, order = ['sample_record', 'sample_record_and_study_description', 'study_description', 'not_recorded'], alpha = [1, 0.55, 0.28, 0.1];
      const rows = cls.filter(k => B[k]).map(k => {
        const tot = sum(order.map(o => B[k][o] || 0));
        return { label: lab('class', k), end: fmt(tot), parts: order.filter(o => B[k][o]).map((o, j) => ({ value: B[k][o], color: k === 'not_recorded' ? C.cls.not_recorded : mixHex(C.cls[k], C.bg, 1 - alpha[order.indexOf(o)]),
          tip: `<b>${lab('class', k)}</b>: ${fmt(B[k][o])} sample records (${pct(B[k][o] / tot)})<br>rest on ${lower(A.basis_labels[o])}` })) };
      });
      hbars('raRec', rows, { share: true, labelW: 150, rowH: 22, axis: true, axisLabel: 'Share of the sample records of the class', aria: 'Stacked bars of what the class of every sample record rests on.' });
      keysEl.innerHTML = keys(order.slice(0, 3).map((o, j) => [A.basis_labels[o], mixHex(C.ink2, C.bg, 1 - alpha[j])]));
      how.textContent = `One bar for every class, over all ${fmt(N.records)} sample records (here a repeat deposit counts once in every study that holds it). The darkest part rests on the sample's own record, the middle part on its record together with a quoted sentence of the study record, and the palest on the study record alone, for example where a record says only that the sample came from a patient and the study says which disease its patients had. Hover over a part to see its count.`;
    }
  }
  function drawAgeBoxes(id) {
    const A = R.samples, w = widthOf(id), sm = isSmall(w), labelW = sm ? 0 : 130, rowH = 18, gap = sm ? 22 : 10, rows = A.ages;
    const H = 14 + rows.length * (rowH + gap) + 34, s = svg(id, w, H, 'Box plots of the recorded age in years of the distinct samples of each class.');
    const x = lin(0, 100, labelW, w - (sm ? 46 : 70));
    gridX(s, x, [0, 18, 40, 60, 80, 100], 10, H - 32);
    S('rect', { x: x(0), y: 10, width: x(18) - x(0), height: H - 42, fill: C.band }, s);
    axisTitle(s, (labelW + w - 70) / 2, H - 4, 'Recorded age in years (the shaded band: under 18)');
    rows.forEach((r, i) => {
      const y = 14 + i * (rowH + gap) + (sm ? 14 : 0), col = C.cls[r.cls];
      if (sm) T(s, 0, y - 4, `${lab('class', r.cls)} (${fmt(r.n)} with an age)`, { 'font-size': 12, fill: C.ink });
      else T(s, labelW - 8, y + rowH / 2 + 4, lab('class', r.cls), { 'text-anchor': 'end', 'font-size': 12.5, fill: C.ink });
      S('line', { x1: x(r.min), x2: x(r.max), y1: y + rowH / 2, y2: y + rowH / 2, stroke: col, 'stroke-width': 1.5 }, s);
      rrect(s, x(r.q1), y, x(r.q3) - x(r.q1), rowH, tint(col, 0.45), 2, { stroke: col });
      S('line', { x1: x(r.median), x2: x(r.median), y1: y - 2, y2: y + rowH + 2, stroke: C.ink, 'stroke-width': 2 }, s);
      if (!sm) T(s, x(r.max) + 8, y + rowH / 2 + 4, `median ${fmt(r.median)}`, { 'font-size': 11.5, fill: C.ink2 });
      hitRow(s, 0, y - 4, w, rowH + 8, `<b>${lab('class', r.cls)}</b>: ${fmt(r.n)} distinct samples with a recorded age<br>median ${fx(r.median, 0)} years, middle half ${fx(r.q1, 0)}–${fx(r.q3, 0)}<br>youngest ${r.min < 1 ? 'under one year' : fx(r.min, 0) + ' years'}, oldest ${fx(r.max, 0)}`);
    });
  }
  function drawCoverage(id) {
    const A = R.samples, w = widthOf(id), sm = isSmall(w), fields = Object.keys(A.coverage_fields), fams = D.order.family;
    const share = f => sum(A.coverage.map(p => p[f] * p.records)) / sum(A.coverage.map(p => p.records));
    const rowsF = fields.slice().sort((a, b) => share(b) - share(a));
    const mean = p => mean_(rowsF.map(f => p[f]));
    const cols = A.coverage.slice().sort((a, b) => fams.indexOf(a.family) - fams.indexOf(b.family) || mean(b) - mean(a));
    const labelW = sm ? 96 : 140, rightW = sm ? 40 : 56, top = 12, rowH = sm ? 16 : 20, cw = (w - labelW - rightW) / cols.length;
    const H = top + rowsF.length * rowH + 8, s = svg(id, w, H, 'Heat map of the share of the sample records of every study platform that record each curated field.');
    cols.forEach((p, j) => S('rect', { x: labelW + j * cw, y: 2, width: Math.max(cw - 0.4, 0.6), height: 6, fill: lift(C.fam[p.family]) }, s));
    const shade = v => v <= 0 ? null : mixHex(C.bg, C.accentInk, 0.18 + 0.82 * v);
    rowsF.forEach((f, i) => {
      const y = top + i * rowH;
      T(s, labelW - 8, y + rowH / 2 + 4, A.coverage_fields[f], { 'text-anchor': 'end', 'font-size': sm ? 11 : 12.5, fill: C.ink });
      S('rect', { x: labelW, y, width: w - labelW - rightW, height: rowH - 2, fill: 'none', stroke: C.line }, s);
      cols.forEach((p, j) => { const c = shade(p[f]); if (c) S('rect', { x: labelW + j * cw, y, width: Math.max(cw - 0.4, 0.6), height: rowH - 2, fill: c }, s); });
      T(s, w - rightW + 6, y + rowH / 2 + 4, pct(share(f), 0), { 'font-size': 11.5, fill: C.ink2 });
    });
    const hit = S('rect', { x: labelW, y: 0, width: w - labelW - rightW, height: H, fill: 'transparent' }, s);
    hit.addEventListener('mousemove', e => {
      const r = s.getBoundingClientRect(), k = w / r.width, px = (e.clientX - r.left) * k, py = (e.clientY - r.top) * k;
      const j = clamp(Math.floor((px - labelW) / cw), 0, cols.length - 1), i = Math.floor((py - top) / rowH), p = cols[j];
      if (i < 0 || i >= rowsF.length) { tip(e, `<b>${esc(p.label)}</b><br>${lab('family', p.family)}, ${fmt(p.records)} sample records`); return; }
      const f = rowsF[i];
      tip(e, `<b>${esc(p.label)}</b> <span class="m">· ${lab('family', p.family)}</span><br>${A.coverage_fields[f]}: ${p[f] > 0 ? `recorded for ${pct(p[f], 0)} of its ${fmt(p.records)} sample records` : `not recorded in any of its ${fmt(p.records)} sample records`}`);
    });
    hit.addEventListener('mouseleave', untip);
  }
  const mean_ = a => a.length ? sum(a) / a.length : 0;
  function drawLabelExample() {
    const E = R.samples.example, rules = E.rules, tot = sum(rules.map(r => r.samples));
    const val = v => v === 'infection_sepsis_not_stated' ? 'Other condition (an infection; sepsis not stated)' : lab('class', v);
    const fields = [...new Set(rules.map(r => r.field))];
    const quoted = rules.filter(r => r.statement);
    $('raExMeta').textContent = `· ${E.study}, ${E.platform}`;
    $('raExText').innerHTML = `GSE65682, an intensive-care cohort, holds ${fmt(E.records)} sample records on one array. Its curation reads ${spell(fields.length)} fields of the records, and each row below is one recorded wording and the class it gives. The wordings ${andList(quoted.map(r => `<code>${esc(r.wording)}</code>`))} name only a group of the study's own analysis, so the study's sentence "${esc(quoted.length ? quoted[0].statement : '')}" is quoted for them. The rows add up to ${fmt(tot)} records, every record of the study.`;
    $('raExTable').innerHTML = `<thead><tr><th>Recorded field</th><th>Recorded wording</th><th>Curated class</th><th class="r">Records</th><th>Rests on</th></tr></thead><tbody>${rules.map(r => `<tr><td>${esc(r.field)}</td><td><code>${esc(r.wording.replace(/^\^/, ''))}${/^\^/.test(r.wording) ? ' …' : ''}</code></td><td>${sw(C.cls[r.value === 'infection_sepsis_not_stated' ? 'other_condition' : r.value] || C.muted)}${val(r.value)}</td><td class="r">${fmt(r.samples)}</td><td class="ui">${r.statement ? `its record, and the study's sentence "${esc(r.statement)}"` : 'its own record'}</td></tr>`).join('')}<tr class="total"><td colspan="3">All records of ${E.study}</td><td class="r">${fmt(tot)}</td><td></td></tr></tbody>`;
  }

  // ======================================================================= 3. repeat deposits and patients
  function drawRepeats() {
    const P = R.repeats, pat = P.patients;
    $('rrLead').innerHTML = `The same physical array or library is often deposited in more than one study, and one patient can give several samples. Nothing was removed. Records of one array or library were linked by what their raw files store: the scan time of an Affymetrix array, the barcode and position of an Illumina BeadChip, the slide of an Agilent array or identical RNA-seq counts, and by a near-identical profile only where studies are known to share samples. The first deposit by release date is the primary record, which identifies the distinct sample; ${fmt(P.repeated_samples)} samples have repeat deposits, and each keeps all its records. Patient identifiers were linked across deposits and studies, giving ${fmt(pat.identified)} identified patients, and one sample is marked for every patient.`;
    drawArcs(); drawSteps(); drawGroupExample();
    $('rrNote').textContent = `A sample without a patient identifier counts as its own patient, so the ${fmt(pat.one_per_patient.total)} samples marked one per patient are an upper bound on the number of patients; the records identify ${fmt(pat.identified)}, of whom ${fmt(pat.in_several_studies)} have samples in more than one study. Records of one sample can differ in class, because studies word their controls differently: this occurs for ${fmt(P.class_differs)} samples, and each record keeps its own label. Studies that share a sample or an identified patient form ${fmt(N.patient_safe_groups)} patient-safe groups, which the benchmark splits keep together.`;
  }
  function drawArcs() {
    const P = R.repeats, id = 'rrArcs', w = widthOf(id), sm = isSmall(w);
    const info = {}; P.platforms.forEach(p => { info[plat(p.study, p.platform)] = p; });
    const nodes = [...new Set(P.shared.flatMap(r => [r.a, r.b]))];
    nodes.forEach(n => { if (!info[n]) throw new Error(`repeats: no platform row for ${n}`); });
    nodes.sort((a, b) => info[a].group.localeCompare(info[b].group) || info[a].year - info[b].year || a.localeCompare(b));
    const studies = {}; nodes.forEach(n => { const s0 = info[n].study; studies[s0] = (studies[s0] || 0) + 1; });
    const label = n => studies[info[n].study] > 1 ? `${info[n].study} ${info[n].platform}` : info[n].study;
    const groups = [...new Set(nodes.map(n => info[n].group))], gapG = sm ? 8 : 18;
    const unit = (w - 24 - gapG * (groups.length - 1)) / Math.max(1, nodes.length - 1);
    const xs = {}; let cx = 12, last = null;
    nodes.forEach(n => { if (last && info[n].group !== info[last].group) cx += gapG; xs[n] = cx; cx += unit; last = n; });
    const maxSpan = Math.max(...P.shared.map(r => Math.abs(xs[r.a] - xs[r.b]))), lift_ = span => 18 + span * 0.62;
    const base = Math.round(lift_(maxSpan) * 0.75) + 26, H = base + (sm ? 92 : 112);
    const s = svg(id, w, H, `Arc diagram of the ${fmt(nodes.length)} study platforms that share arrays or libraries with another study platform.`);
    const maxS = Math.max(...P.shared.map(r => r.samples)), sw_ = v => 1.2 + Math.sqrt(v / maxS) * (sm ? 5 : 8);
    S('line', { x1: 6, x2: w - 6, y1: base, y2: base, stroke: C.lineStrong }, s);
    const arcs = P.shared.slice().sort((a, b) => b.samples - a.samples).map(r => {
      const xa = Math.min(xs[r.a], xs[r.b]), xb = Math.max(xs[r.a], xs[r.b]), h = lift_(xb - xa), col = lift(C.fam[info[r.a].family]);
      const p = S('path', { d: `M${xa},${base}C${xa},${base - h} ${xb},${base - h} ${xb},${base}`, fill: 'none', stroke: col, 'stroke-width': sw_(r.samples), 'stroke-opacity': 0.55, 'stroke-linecap': 'butt' }, s);
      p.dataset.a = r.a; p.dataset.b = r.b;
      hover(p, `<b>${esc(r.a)}</b> and <b>${esc(r.b)}</b><br>share ${fmt(r.samples)} arrays or libraries<br><span class="m">cohort group ${esc(r.group)}</span>`);
      p.addEventListener('mouseenter', () => p.setAttribute('stroke-opacity', 0.95)); p.addEventListener('mouseleave', () => p.setAttribute('stroke-opacity', 0.55));
      return p;
    });
    nodes.forEach(n => {
      const p = info[n], x = xs[n], partners = P.shared.filter(r => r.a === n || r.b === n);
      const c = S('circle', { cx: x, cy: base, r: sm ? 3.5 : 4.5, fill: lift(C.fam[p.family]), stroke: C.bg, 'stroke-width': 1.2 }, s);
      c.style.cursor = 'pointer';
      hover(c, `<b>${esc(n)}</b> <span class="m">· ${lab('family', p.family)}, released ${p.year}</span><br>${fmt(p.records)} sample records: ${fmt(p.primary)} primary, ${fmt(p.repeated)} repeat deposits<br>shares samples with ${partners.length === 1 ? 'one other study platform' : `${spell(partners.length)} other study platforms`}<br><span class="m">click to open the study</span>`);
      c.addEventListener('mouseenter', () => arcs.forEach(a => a.setAttribute('stroke-opacity', a.dataset.a === n || a.dataset.b === n ? 0.95 : 0.12)));
      c.addEventListener('mouseleave', () => arcs.forEach(a => a.setAttribute('stroke-opacity', 0.55)));
      c.addEventListener('click', () => { untip(); goTo(`studies/${p.study}`); });
      T(s, x + 3, base + 12, label(n), { 'font-size': sm ? 9.5 : 11, fill: C.ink2, transform: `rotate(55 ${x + 3} ${base + 12})` });
    });
  }
  function drawSteps() {
    const st = R.repeats.steps, max = st.included_sample_records, solid = C.accentInk, light = tint(C.muted, 0.72);
    const part = (v, color, t, dash, text) => ({ value: v, color, tip: t, dash, text, fill: dash ? tint(C.muted, 0.82) : null });
    const noValues = st.records_without_values, repeats = st.repeat_deposits, further = st.distinct_samples - st.one_sample_per_patient;
    const rows = [
      { label: 'Sample records', end: fmt(st.included_sample_records), parts: [part(st.sample_records_with_values, solid, `<b>${fmt(st.sample_records_with_values)} records with expression values</b>`), part(noValues, light, `<b>${fmt(noValues)} records without expression values</b><br>listed with the reason in a ledger`, true)] },
      { label: 'With expression values', end: fmt(st.sample_records_with_values), parts: [part(st.distinct_samples, solid, `<b>${fmt(st.distinct_samples)} distinct samples</b>: their primary records`), part(repeats, light, `<b>${fmt(repeats)} repeat deposits</b> of an array or library already counted<br>kept and linked to the primary record`, true)] },
      { label: 'Distinct samples', end: fmt(st.distinct_samples), parts: [part(st.one_sample_per_patient, solid, `<b>${fmt(st.one_sample_per_patient)} samples marked one per patient</b>`), part(further, light, `<b>${fmt(further)} further samples</b> of a patient already counted<br>kept: time series, repeat visits`, true)] },
      { label: 'One sample per patient', end: fmt(st.one_sample_per_patient), parts: [part(st.identified_patients, solid, `<b>${fmt(st.identified_patients)} identified patients</b>`, false, 'identified patients'), part(st.samples_without_identifier, tint(C.accentInk, 0.45), `<b>${fmt(st.samples_without_identifier)} samples without a patient identifier</b><br>each counted as its own patient`, false, 'no identifier')] },
    ];
    hbars('rrSteps', rows, { max, labelW: 190, rowH: 24, axis: true, axisLabel: 'Sample records of the included studies, and what they hold', aria: 'Bars from the sample records to the distinct samples to one sample per patient.' });
  }
  function drawGroupExample() {
    const P = R.repeats, E = P.example, rows = P.platforms.filter(p => p.group === E.group).sort((a, b) => a.year - b.year || a.study.localeCompare(b.study));
    const fam = rows[0] ? rows[0].family : 'affymetrix', col = lift(C.fam[fam]);
    $('rrExMeta').textContent = `· cohort group ${E.group}`;
    $('rrExText').innerHTML = `GSE26440 belongs to a group of ${spell(E.platforms)} studies of children with sepsis or septic shock on one Affymetrix array, released from ${rows[0].year} to ${rows[rows.length - 1].year}. Their arrays were deposited again and again: the group holds ${fmt(E.records)} sample records but only ${fmt(E.distinct)} distinct samples, ${fx(E.records / E.distinct, 1)} records for every array. Each bar is one study platform; the solid part is its primary records, the arrays it deposited first, and the fainter part its repeat deposits of arrays that an earlier study deposited.`;
    hbars('rrExample', rows.map(p => ({ label: `${p.study} (${p.year})`, end: fmt(p.records), parts: [
      { value: p.primary, color: col, tip: `<b>${p.study}</b>: ${fmt(p.primary)} primary records` },
      { value: p.repeated, color: tint(col, 0.62), tip: `<b>${p.study}</b>: ${fmt(p.repeated)} repeat deposits of arrays deposited earlier` }] })),
    { labelW: 150, rowH: 18, gap: 7, axis: true, axisLabel: 'Sample records', aria: `Bars of the primary records and repeat deposits of the study platforms of cohort group ${E.group}.` });
  }

  // ======================================================================= 4. from probes to genes
  function mapColor(k) { const c = (R.genes.mapping_colors || {})[k]; return k === 'no_gene_in_annotation' ? C.cls.not_recorded : lift(c || C.muted); }
  function drawGenes() {
    const G = R.genes;
    $('rgLead').innerHTML = `Arrays measure RNA with probes, short sequences that each match one transcript; RNA-seq counts the reads of every gene identifier of an annotation. So that the studies can be used together, every row of every table was mapped to one reference list of genes, that of the HUGO Gene Nomenclature Committee (HGNC), which gives every human gene one approved symbol and a fixed identifier. A row whose annotation names several genes was not assigned to any. Each study platform then has a gene-level table with one row per gene: on arrays the row of a gene with the highest mean expression is kept, and for RNA-seq the counts of the rows of a gene are summed. The row-level tables are released as well.`;
    drawMapping(); drawGeneChart();
  }
  function drawMapping() {
    const G = R.genes, fams = D.order.family.filter(f => G.mapping[f]), keys_ = Object.keys(G.mapping_labels);
    const total = {}; keys_.forEach(k => { total[k] = sum(fams.map(f => G.mapping[f][k])); });
    const row = (label, m, bold) => {
      const tot = sum(keys_.map(k => m[k]));
      return { label, bold, end: V.map === 'share' ? '' : `${fx(tot / 1e6, 2)} M`, parts: keys_.map(k => ({ value: m[k], color: mapColor(k),
        text: V.map === 'share' ? pct(m[k] / tot, 0) : '', tip: `<b>${label}</b><br>${G.mapping_labels[k]}: ${fmt(m[k])} rows (${pct(m[k] / tot, 1)})` })) };
    };
    const rows = [...fams.map(f => row(lab('family', f), G.mapping[f])), row('All rows', total, true)];
    hbars('rgMap', rows, { share: V.map === 'share', labelW: 160, rowH: 24, axis: true, endW: V.map === 'share' ? 10 : 64,
      tickF: v => v === 0 ? '0' : `${fmt(v / 1e6, 1)} M`, axisLabel: V.map === 'share' ? 'Share of the rows of the technology' : 'Rows of the tables (millions)',
      aria: 'Stacked bars of the rows of every technology by the outcome of their mapping to HGNC genes.' });
    $('rgMapKeys').innerHTML = keys(keys_.map(k => [G.mapping_labels[k], mapColor(k)]));
  }
  function drawGeneChart() {
    const G = R.genes, how = $('rgHow').lastElementChild, note = $('rgNote');
    const any = G.shared[0][1], all = G.shared[G.shared.length - 1][1], nT = G.tables.length;
    note.textContent = `Only ${fmt(all)} genes are in all ${fmt(nT)} gene-level tables, while ${fmt(any)} are in at least one; a table holds from ${fmt(G.genes.fewest)} to ${fmt(G.genes.most)} genes (median ${fmt(G.genes.median)}). This is why the corrected tables keep only the genes that nearly all their studies measure. The per-study tables are released as processed, each on its own scale, which the sample table names.`;
    $('rgKeys').innerHTML = '';
    if (V.genes === 'tables') {
      drawGenesPerTable('rgChart');
      $('rgKeys').innerHTML = keys(D.order.family.map(f => [lab('family', f), lift(C.fam[f])]));
      how.textContent = 'One dot for every gene-level table, in the row of its technology, placed at the number of genes it holds; the ring is the number of protein-coding genes among them. RNA-seq tables reach far to the right, because sequencing reads every annotated gene, while an array measures only the genes it has probes for. Hover over a dot to see its table; click it to open its study.';
    } else if (V.genes === 'shared') {
      drawShared('rgChart');
      how.textContent = `The line gives, for every number of tables, how many genes are held by at least that many of the ${fmt(nT)} gene-level tables. At the left are all the genes found in any table; at the right the genes that every table holds. Hover over the line to read it.`;
    } else {
      drawScales('rgChart');
      const sc = Object.keys(G.scale_labels).filter(k => G.scales.some(t => t.scale === k));
      $('rgKeys').innerHTML = keys(sc.map(k => [G.scale_labels[k], scaleColor(k)]));
      how.textContent = 'One column for every gene-level table, sorted by its median value within its scale: the thin line runs from the 5th to the 95th percentile of all its values, the box covers the middle half and the tick is the median. Tables on a log₂ scale share the left axis, in which one step is a doubling; the few tables that the authors deposited on another scale are on the right, on their own axis. Hover over a column to see its table.';
    }
  }
  function drawGenesPerTable(id) {
    const G = R.genes, w = widthOf(id), sm = isSmall(w), fams = D.order.family.filter(f => G.tables.some(t => t.family === f));
    const labelW = sm ? 0 : 150, rowH = sm ? 46 : 40, top = 6, H = top + fams.length * rowH + 36;
    const s = svg(id, w, H, 'Strip chart of the number of genes in every gene-level table, by technology.');
    const x = lin(0, niceMax(G.genes.most), labelW, w - 14);
    gridX(s, x, ticks(0, niceMax(G.genes.most), sm ? 3 : 5), top, H - 32);
    axisTitle(s, (labelW + w) / 2, H - 4, 'Genes in the table');
    fams.forEach((f, i) => {
      const y0 = top + i * rowH, cy = y0 + rowH / 2 + (sm ? 6 : 0), col = lift(C.fam[f]), list = G.tables.filter(t => t.family === f);
      if (sm) T(s, 0, y0 + 12, `${lab('family', f)} (${fmt(list.length)})`, { 'font-size': 12, fill: C.ink });
      else T(s, labelW - 10, cy + 4, `${lab('family', f)} (${fmt(list.length)})`, { 'text-anchor': 'end', 'font-size': 12.5, fill: C.ink });
      const r = rng(i + 7);
      list.forEach(t => {
        const jy = cy + (r() - 0.5) * (rowH - 18);
        S('circle', { cx: x(t.coding), cy: jy, r: 3, fill: 'none', stroke: col, 'stroke-width': 1 }, s);
        const c = S('circle', { cx: x(t.genes), cy: jy, r: 4, fill: col, 'fill-opacity': 0.8, stroke: C.bg, 'stroke-width': 0.8 }, s);
        c.style.cursor = 'pointer';
        hover(c, `<b>${esc(plat(t.study, t.platform))}</b> <span class="m">· ${lab('family', f)}</span><br>${fmt(t.genes)} genes, ${fmt(t.coding)} of them protein-coding<br><span class="m">click to open the study</span>`);
        c.addEventListener('click', () => { untip(); goTo(`studies/${t.study}`); });
      });
      const med = median(list.map(t => t.genes));
      S('line', { x1: x(med), x2: x(med), y1: y0 + (sm ? 16 : 6), y2: y0 + rowH - 4, stroke: C.ink, 'stroke-width': 1.5 }, s);
    });
  }
  function drawShared(id) {
    const G = R.genes, w = widthOf(id), sm = isSmall(w), m = { l: sm ? 50 : 62, r: 16, t: 12, b: 40 }, H = sm ? 260 : 300;
    const pts = G.shared, nT = pts[pts.length - 1][0], vmax = niceMax(pts[0][1]);
    const s = svg(id, w, H, 'Line chart of the number of genes held by at least a given number of gene-level tables.');
    const x = lin(1, nT, m.l, w - m.r), y = lin(0, vmax, H - m.b, m.t);
    gridY(s, y, ticks(0, vmax, 4), m.l, w - m.r); gridX(s, x, [1, ...ticks(0, nT, sm ? 4 : 7).filter(v => v >= 10), nT].filter((v, i, a) => a.indexOf(v) === i), H - m.b, H - m.b);
    axisTitle(s, (m.l + w - m.r) / 2, H - 4, 'Held by at least this many tables');
    axisTitle(s, 12, (m.t + H - m.b) / 2, 'Genes', true);
    S('path', { d: pts.map(([k, g], i) => `${i ? 'L' : 'M'}${x(k)},${y(g)}`).join('') + `L${x(nT)},${y(0)}L${x(1)},${y(0)}Z`, fill: C.accent, 'fill-opacity': 0.18 }, s);
    S('path', { d: pts.map(([k, g], i) => `${i ? 'L' : 'M'}${x(k)},${y(g)}`).join(''), fill: 'none', stroke: C.accentInk, 'stroke-width': 2 }, s);
    [[pts[0], 'start'], [pts[pts.length - 1], 'end']].forEach(([[k, g], a]) => {
      S('circle', { cx: x(k), cy: y(g), r: 4, fill: C.accentInk }, s);
      T(s, x(k) + (a === 'start' ? 8 : -8), y(g) - 8, `${fmt(g)} genes${k === nT ? ` in all ${fmt(nT)} tables` : ' in at least one table'}`, { 'text-anchor': a, 'font-size': 12, fill: C.ink, ...halo() });
    });
    const dot = S('circle', { r: 4, fill: C.ink, opacity: 0 }, s);
    const hit = S('rect', { x: m.l, y: m.t, width: w - m.l - m.r, height: H - m.t - m.b, fill: 'transparent' }, s);
    hit.addEventListener('mousemove', e => {
      const r = s.getBoundingClientRect(), px = (e.clientX - r.left) * w / r.width, k = clamp(Math.round(x.inv(px)), 1, nT), p = pts.find(q => q[0] === k) || pts[0];
      dot.setAttribute('cx', x(p[0])); dot.setAttribute('cy', y(p[1])); dot.setAttribute('opacity', 1);
      tip(e, `<b>${fmt(p[1])} genes</b><br>held by at least ${plural(p[0], 'table')} of ${fmt(nT)}`);
    });
    hit.addEventListener('mouseleave', () => { untip(); dot.setAttribute('opacity', 0); });
  }
  const scaleColor = k => ({ log2: lift(C.fam.affymetrix), log2_like: C.fam.agilent, 'log2 counts per million': lift(C.fam.rna_seq), linear: C.cls.sepsis, other: C.cls.other_condition })[k] || C.muted;
  function drawScales(id) {
    const G = R.genes, w = widthOf(id), sm = isSmall(w), logs = ['log2', 'log2_like', 'log2 counts per million'];
    const left = G.scales.filter(t => logs.includes(t.scale)).sort((a, b) => logs.indexOf(a.scale) - logs.indexOf(b.scale) || a.q[2] - b.q[2]);
    const right = G.scales.filter(t => !logs.includes(t.scale)).sort((a, b) => a.q[2] - b.q[2]);
    const m = { l: sm ? 34 : 44, t: 12, b: 42 }, H = sm ? 300 : 340, gapW = sm ? 46 : 64, rw = Math.max(60, right.length * (sm ? 14 : 18));
    const lw = w - m.l - gapW - rw - 8, s = svg(id, w, H, 'Range chart of the values of every gene-level table, on their own scales.');
    const lo = Math.floor(Math.min(...left.map(t => t.q[0]))), hi = Math.ceil(Math.max(...left.map(t => t.q[4])));
    const y = lin(lo, hi, H - m.b, m.t), cw = lw / left.length;
    gridY(s, y, ticks(lo, hi, 5), m.l, m.l + lw);
    axisTitle(s, 11, (m.t + H - m.b) / 2, 'Value (log₂ scale)', true);
    T(s, m.l, H - m.b + 16, `${fmt(left.length)} tables on a log₂ scale`, { 'font-size': 11.5, fill: C.muted });
    const col = (t, i, x0, cw_, yy) => {
      const cx = x0 + (i + 0.5) * cw_, c = scaleColor(t.scale), bw = Math.max(2, cw_ * 0.62);
      S('line', { x1: cx, x2: cx, y1: yy(t.q[0]), y2: yy(t.q[4]), stroke: c, 'stroke-width': 1 }, s);
      S('rect', { x: cx - bw / 2, y: yy(t.q[3]), width: bw, height: Math.max(1, yy(t.q[1]) - yy(t.q[3])), fill: tint(c, 0.35), stroke: c, 'stroke-width': 0.6 }, s);
      S('line', { x1: cx - bw / 2 - 1, x2: cx + bw / 2 + 1, y1: yy(t.q[2]), y2: yy(t.q[2]), stroke: C.ink, 'stroke-width': 1.4 }, s);
      hitRow(s, x0 + i * cw_, m.t, cw_, H - m.b - m.t, `<b>${esc(plat(t.study, t.platform))}</b> <span class="m">· ${G.scale_labels[t.scale]}</span><br>${fmt(t.samples)} samples<br>median ${fx(t.q[2], 2)}; middle half ${fx(t.q[1], 2)} to ${fx(t.q[3], 2)}<br>5th to 95th percentile ${fx(t.q[0], 2)} to ${fx(t.q[4], 2)}`);
    };
    left.forEach((t, i) => col(t, i, m.l, cw, y));
    // the tables on other scales: a signed logarithm, sign(v) log10(1 + |v|)
    const sl = v => Math.sign(v) * Math.log10(1 + Math.abs(v)), rmax = Math.max(...right.map(t => t.q[4])), rmin = Math.min(0, ...right.map(t => t.q[0]));
    const yr0 = lin(sl(rmin), sl(rmax) * 1.05, H - m.b, m.t), yr = v => yr0(sl(v)), x0r = m.l + lw + gapW;
    [0, 1, 10, 100, 1000, 10000, 100000].filter(v => v <= rmax * 1.05).forEach(v => {
      S('line', { x1: x0r, x2: x0r + rw, y1: yr(v), y2: yr(v), stroke: C.grid }, s);
      T(s, x0r - 5, yr(v) + 4, fmt(v), { 'text-anchor': 'end', 'font-size': 10.5, fill: C.muted });
    });
    right.forEach((t, i) => col(t, i, x0r, rw / right.length, yr));
    T(s, x0r, H - m.b + 16, `${spell(right.length)} on other scales`, { 'font-size': 11.5, fill: C.muted });
  }

  // ======================================================================= 5. the checks
  const SEX = () => ({ female: css('--r-female'), male: css('--r-male'), not_recorded: C.cls.not_recorded });
  function drawChecks() {
    const K = R.checks, sx = K.sex_receipt, rp = K.repeat_receipt, sc = K.score_receipt;
    $('rcLead').innerHTML = `The checks were made inside single studies, where no correction across studies can affect them. Each asks whether the data behave as known biology says they should: the sex read from the genes against the recorded sex, two deposits of one sample against each other, and published sepsis genes against the recorded classes. Such checks can expose a mislabeled table, a wrong assignment of files to records or a failed processing step. They cannot show that a label is clinically right, or that a score is useful at the bedside.`;
    drawSexPlane(); drawSexTables(); drawRepeatPairs(); drawScore();
    $('rcSexNote').textContent = `Over every tested table, including those without XIST, the sex read from the genes agreed with the recorded sex for ${fmt(sx.concordant)} of the ${fmt(sx.compared)} samples with both (${pct(sx.concordance)}). In ${spell(K.sex_matrices.unreliable.length)} tables (${andList(K.sex_matrices.unreliable)}) the agreement was low and their recorded sex is marked as unreliable; in the others it was ${pct(sx.concordance_outside_unreliable_matrices)}. The genes also give a sex to samples whose records state none.`;
    $('rcRepNote').textContent = `Across ${fmt(rp.pairs)} pairs, the Spearman correlation of the two deposits of a sample was at least ${fx(rp.spearman_lowest, 3)} (median ${fx(rp.spearman_median, 3)}), and ${fmt(rp.most_similar)} pairs were each other's closest match in the other table. Pairs linked because a record names a re-analysis correlate less, probably because the second study processed the data again in its own way.`;
    $('rcScoreNote').textContent = `In the ${fmt(sc.strata)} study platforms and tissues with at least ${spell(K.minimum_group)} sepsis samples and ${spell(K.minimum_group)} healthy controls, the published genes changed in their published direction in ${fmt(sc.as_published)} of ${fmt(sc.comparisons)} comparisons (${pct(sc.as_published / sc.comparisons)}), and the score separated sepsis from healthy controls with an AUROC from ${fx(sc.lowest, 2)} to ${fx(Math.max(...K.strata.map(r => r.auroc)), 2)} (median ${fx(sc.median, 3)}). The score was derived in earlier work partly from studies that are also in the resource; in the ${fmt(sc.independent_strata)} study platforms and tissues independent of its derivation, the median AUROC is ${fx(sc.independent_median, 3)}.`;
    $('rcResult').innerHTML = `<b>What the checks show, and what they do not.</b> They show that the tables carry the biology they should: sexes that separate, repeat deposits that agree and sepsis genes that move in their published direction, in study after study. They do not test the clinical labels themselves, which follow what the investigators recorded, and comparisons of sepsis against healthy controls test the plausibility of the data, not the clinical value of the score.`;
  }
  function drawSexPlane() {
    const K = R.checks, X = K.sex, id = 'rcSexPlane', box = $(id), n = X.x.length, col = SEX();
    if (!box.querySelector('canvas')) box.innerHTML = '<canvas role="img"></canvas>';
    const cv = box.querySelector('canvas'), w = box.clientWidth || 420, sm = w < 380;
    const fit = canvasFit(cv, { height: Math.round(Math.min(w, 520) * 0.92) }), { ctx } = fit, m = { l: 40, r: 10, t: 10, b: 38 };
    const xs = X.x, ys = X.y, xlo = Math.max(-0.6, Math.min(...xs)), xhi = Math.min(1.6, Math.max(...xs)), ylo = Math.max(-0.6, Math.min(...ys)), yhi = Math.min(1.6, Math.max(...ys));
    const sxf = lin(xlo, xhi, m.l, fit.w - m.r), syf = lin(ylo, yhi, fit.h - m.b, m.t);
    ctx.clearRect(0, 0, fit.w, fit.h);
    ctx.font = `11px ${getComputedStyle(document.body).fontFamily}`; ctx.lineWidth = 1;
    ticks(xlo, xhi, 4).forEach(v => { ctx.strokeStyle = C.grid; ctx.beginPath(); ctx.moveTo(sxf(v), m.t); ctx.lineTo(sxf(v), fit.h - m.b); ctx.stroke(); ctx.fillStyle = C.muted; ctx.textAlign = 'center'; ctx.fillText(fmt(v, 1), sxf(v), fit.h - m.b + 14); });
    ticks(ylo, yhi, 4).forEach(v => { ctx.strokeStyle = C.grid; ctx.beginPath(); ctx.moveTo(m.l, syf(v)); ctx.lineTo(fit.w - m.r, syf(v)); ctx.stroke(); ctx.fillStyle = C.muted; ctx.textAlign = 'right'; ctx.fillText(fmt(v, 1), m.l - 5, syf(v) + 4); });
    ctx.fillStyle = C.ink2; ctx.textAlign = 'center'; ctx.font = `12px ${getComputedStyle(document.body).fontFamily}`;
    ctx.fillText('XIST, from the male (0) to the female mean (1)', (m.l + fit.w - m.r) / 2, fit.h - 6);
    ctx.save(); ctx.translate(11, (m.t + fit.h - m.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(sm ? 'Y genes, female (0) to male (1)' : 'Y-chromosome genes, female (0) to male (1)', 0, 0); ctx.restore();
    const levels = X.levels.rec, rad = 2.1, sprites = levels.map(l => dotSprite(col[l], rad, l === 'not_recorded' ? 0.45 : 0.7));
    const px = new Float32Array(n), py = new Float32Array(n);
    // not recorded first, then the recorded sexes, then the disagreeing samples with a dark ring on top
    for (const pass of [2, 0, 1]) for (let i = 0; i < n; i++) {
      if (+X.rec[i] !== pass) continue;
      px[i] = sxf(clamp(xs[i], xlo, xhi)); py[i] = syf(clamp(ys[i], ylo, yhi)); stamp(ctx, sprites[pass], px[i], py[i]);
    }
    ctx.lineWidth = 1.3; ctx.strokeStyle = C.ink;
    for (let i = 0; i < n; i++) if (X.agr[i] === '1') { ctx.beginPath(); ctx.arc(px[i], py[i], rad + 1.6, 0, 2 * Math.PI); ctx.stroke(); }
    cv.setAttribute('aria-label', `Scatter plot of the ${fmt(n)} samples of tables that measure XIST, by their XIST and Y-chromosome genes, colored by recorded sex.`);
    cv.onmousemove = e => {
      const r = cv.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top; let best = -1, bd = 49;
      for (let i = 0; i < n; i++) { const d = (px[i] - x) ** 2 + (py[i] - y) ** 2; if (d < bd) { bd = d; best = i; } }
      if (best < 0) { untip(); return; }
      const rec = levels[+X.rec[best]], call = X.levels.call[+X.call[best]], agr = X.levels.agr[+X.agr[best]];
      tip(e, `<b>${esc(X.sample[best])}</b> <span class="m">· ${esc(X.tables[X.tab[best]])}</span><br>recorded: ${rec === 'not_recorded' ? 'not recorded' : rec}<br>read from the genes: ${call}${agr === 'discordant' ? '<br><b>disagrees with the record</b>' : ''}`);
    };
    cv.onmouseleave = untip;
    $('rcSexKeys').innerHTML = keys([['Recorded female', col.female, 'dot'], ['Recorded male', col.male, 'dot'], ['Sex not recorded', col.not_recorded, 'dot'], ['Disagrees with the genes', C.ink, 'ring']]);
  }
  function drawSexTables() {
    const K = R.checks, rows = K.sex_tables.slice().sort((a, b) => a.concordant / a.compared - b.concordant / b.compared);
    const id = 'rcSexTables', w = widthOf(id), labelW = 132, rowH = 11, top = 16, H = top + rows.length * rowH + 34;
    const s = svg(id, w, H, 'Bar chart of the share of recorded sexes that agree with the genes, in every table with a recorded sex.');
    const x = lin(0, 1, labelW, w - 44);
    gridX(s, x, [0, 0.25, 0.5, 0.75, 1], top - 4, H - 30, v => pct(v, 0));
    axisTitle(s, (labelW + w - 44) / 2, H - 4, 'Recorded sexes that agree');
    T(s, labelW, 10, `${fmt(rows.length)} tables with a recorded sex`, { 'font-size': 11, fill: C.muted });
    rows.forEach((r, i) => {
      const y = top + i * rowH, share = r.concordant / r.compared, c = r.marked ? C.tier.unanchored : C.accentInk;
      S('rect', { x: labelW, y: y + 1, width: Math.max(1, x(share) - labelW), height: rowH - 2, fill: r.marked ? c : tint(c, 0.25) }, s);
      T(s, labelW - 6, y + rowH - 2, r.label, { 'text-anchor': 'end', 'font-size': 9, fill: r.marked ? C.ink : C.muted });
      hitRow(s, 0, y, w, rowH, `<b>${esc(r.label)}</b>${r.marked ? ' <span class="m">· recorded sex marked unreliable</span>' : ''}<br>${fmt(r.concordant)} of ${fmt(r.compared)} recorded sexes agree (${pct(share)})`);
    });
  }
  function drawRepeatPairs() {
    const K = R.checks, Pr = K.repeats, id = 'rcRepeats', w = widthOf(id), sm = isSmall(w), vals = V.rep === 'cc' ? Pr.cc : Pr.sp;
    const order = Pr.links.slice().sort((a, b) => Pr.link.filter(v => Pr.links[v] === b).length - Pr.link.filter(v => Pr.links[v] === a).length);
    const labelW = sm ? 0 : 290, rowH = sm ? 50 : 40, top = 8, H = top + order.length * rowH + 38;
    const s = svg(id, w, H, 'Strip chart of the correlation of the two deposits of every repeated sample, by the evidence that links them.');
    const lo = V.rep === 'cc' ? Math.floor(Math.min(...vals) * 10) / 10 : Math.floor(Math.min(...vals) * 100) / 100;
    const x = lin(lo, 1, labelW + 8, w - 14);
    gridX(s, x, ticks(lo, 1, sm ? 3 : 6), top, H - 32, v => fx(v, V.rep === 'cc' ? 1 : 2));
    axisTitle(s, (labelW + w) / 2, H - 4, V.rep === 'cc' ? 'Correlation of the two profiles, every gene centered in its own table' : 'Spearman correlation of the two profiles');
    const r = rng(11);
    order.forEach((lk, i) => {
      const li = Pr.links.indexOf(lk), idx = Pr.link.map((v, j) => v === li ? j : -1).filter(j => j >= 0), y0 = top + i * rowH, cy = y0 + rowH / 2 + (sm ? 6 : 0);
      const label = `${Pr.link_labels[lk]} (${fmt(idx.length)})`;
      if (sm) T(s, 0, y0 + 12, label, { 'font-size': 11.5, fill: C.ink }); else T(s, labelW, cy + 4, label, { 'text-anchor': 'end', 'font-size': 12, fill: C.ink });
      const v = idx.map(j => vals[j]).sort((a, b) => a - b), q1 = quant(v, 0.25), q3 = quant(v, 0.75), md = quant(v, 0.5);
      idx.forEach(j => {
        const jy = cy + (r() - 0.5) * (rowH - 20), cx = x(vals[j]), not = Pr.best[j] === '0';
        const node = not ? S('path', { d: `M${cx - 4},${jy - 4}L${cx + 4},${jy + 4}M${cx - 4},${jy + 4}L${cx + 4},${jy - 4}`, stroke: C.cls.sepsis, 'stroke-width': 2 }, s)
          : S('circle', { cx, cy: jy, r: 2.3, fill: C.accentInk, 'fill-opacity': 0.4 }, s);
        hover(node, `<b>${esc(Pr.sample[j])}</b><br>${esc(Pr.a[j])} and ${esc(Pr.b[j])}<br>${V.rep === 'cc' ? 'centered correlation' : 'Spearman correlation'} ${fx(vals[j], 3)}${not ? '<br><b>not the closest match of its twin</b>' : ''}`);
      });
      rrect(s, x(q1), cy - 4, Math.max(2, x(q3) - x(q1)), 8, 'none', 2, { stroke: C.ink, 'stroke-width': 1.2 });
      S('line', { x1: x(md), x2: x(md), y1: cy - 7, y2: cy + 7, stroke: C.ink, 'stroke-width': 2.2 }, s);
    });
  }
  // the area under a ROC curve by the trapezoid rule
  const area = (fpr, tpr) => { let a = 0; for (let i = 1; i < fpr.length; i++) a += (fpr[i] - fpr[i - 1]) * (tpr[i] + tpr[i - 1]) / 2; return a; };
  function drawScore() {
    const K = R.checks, how = $('rcScoreHow').lastElementChild;
    if (!K.strata.some(r => r.key === V.pick)) V.pick = K.strata[0].key;
    const roc = V.score === 'roc';
    $('rcScoreBox').hidden = !roc; $('rcGenes').hidden = roc;
    if (roc) {
      how.textContent = `One ROC curve of the Sepsis MetaScore for every study platform and tissue with enough sepsis samples and healthy controls; the bold line is the median of the curves and the shade their middle half. A curve that climbs straight up the left edge separates sepsis from healthy controls perfectly; the diagonal is a coin toss. The dots give the AUROC of each, with its 95% interval. Click a curve or a dot to pick out one study.`;
      drawRoc(); drawAurocList();
    } else {
      how.textContent = `One row for each of the published genes (the eleven of the score and HLA-DRA), one column for every study platform and tissue, in the order of their AUROC. Each cell is the standardized difference between sepsis and healthy controls: red where the gene is higher in sepsis, blue where it is lower, the stronger the color the larger the difference. The strip on the left gives the direction in which each gene was published. Hover over a cell to see it.`;
      drawGeneMap();
    }
    const p = K.strata.find(r => r.key === V.pick);
    $('rcScorePick').innerHTML = `<b>${esc(plat(p.study, p.platform))}</b>, ${lower(tissueLab(p.tissue))}: ${fmt(p.sepsis)} sepsis samples and ${fmt(p.healthy)} healthy controls, AUROC ${fx(p.auroc, 3)}${p.low != null ? ` (95% interval ${fx(p.low, 3)} to ${fx(p.high, 3)})` : ''}; with one sample per patient, ${fx(p.one_per_patient, 3)}.${({ discovery: ' This study was one of those used to derive the score.', validation: ' This study was used to validate the score.', shared: ' This study shares samples with the studies used to derive the score.' })[p.derivation] || ''}`;
  }
  function drawRoc() {
    const K = R.checks, id = 'rcRoc', w = Math.min(widthOf(id), 520), m = { l: 44, r: 10, t: 10, b: 40 }, H = Math.round(w * 0.92);
    const s = svg(id, w, H, `ROC curves of the Sepsis MetaScore in ${fmt(K.strata.length)} study platforms and tissues.`);
    const x = lin(0, 1, m.l, w - m.r), y = lin(0, 1, H - m.b, m.t);
    gridX(s, x, [0, 0.25, 0.5, 0.75, 1], m.t, H - m.b, v => fx(v, 2)); gridY(s, y, [0, 0.25, 0.5, 0.75, 1], m.l, w - m.r, v => fx(v, 2));
    axisTitle(s, (m.l + w - m.r) / 2, H - 4, 'False-positive rate (healthy controls above the threshold)');
    axisTitle(s, 11, (m.t + H - m.b) / 2, 'True-positive rate (sepsis above it)', true);
    S('line', { x1: x(0), y1: y(0), x2: x(1), y2: y(1), stroke: C.lineStrong, 'stroke-dasharray': '4 4' }, s);
    // the median curve and the middle half: the highest true-positive rate at or below each false-positive rate
    const grid = Array.from({ length: 101 }, (_, i) => i / 100);
    const at = (r, f) => { let best = 0; for (let i = 0; i < r.fpr.length; i++) if (r.fpr[i] <= f + 1e-12) best = Math.max(best, r.tpr[i]); return best; };
    const q = grid.map(f => { const v = K.strata.map(r => at(r, f)).sort((a, b) => a - b); return [quant(v, 0.25), quant(v, 0.5), quant(v, 0.75)]; });
    S('path', { d: grid.map((f, i) => `${i ? 'L' : 'M'}${x(f)},${y(q[i][2])}`).join('') + grid.slice().reverse().map((f, i) => `L${x(f)},${y(q[100 - i][0])}`).join('') + 'Z', fill: C.ink, 'fill-opacity': 0.08 }, s);
    const curves = K.strata.map(r => {
      const p = S('path', { d: r.fpr.map((f, i) => `${i ? 'L' : 'M'}${x(f)},${y(r.tpr[i])}`).join(''), fill: 'none', stroke: C.lineStrong, 'stroke-width': 1, 'stroke-opacity': 0.8 }, s);
      p.style.cursor = 'pointer';
      hover(p, `<b>${esc(plat(r.study, r.platform))}</b>, ${lower(tissueLab(r.tissue))}<br>AUROC ${fx(r.auroc, 3)}<br>${fmt(r.sepsis)} sepsis, ${fmt(r.healthy)} healthy controls`);
      p.addEventListener('click', () => { V.pick = r.key; drawScore(); });
      return [r, p];
    });
    S('path', { d: grid.map((f, i) => `${i ? 'L' : 'M'}${x(f)},${y(q[i][1])}`).join(''), fill: 'none', stroke: C.ink, 'stroke-width': 2.4 }, s);
    const sel = curves.find(([r]) => r.key === V.pick);
    if (sel) { sel[1].setAttribute('stroke', C.cls.sepsis); sel[1].setAttribute('stroke-width', 2.6); sel[1].setAttribute('stroke-opacity', 1); s.appendChild(sel[1]); }
    T(s, x(0.98), y(0.06), 'bold: median curve', { 'text-anchor': 'end', 'font-size': 11, fill: C.ink2 });
  }
  function drawAurocList() {
    const K = R.checks, rows = K.strata.slice().sort((a, b) => b.auroc - a.auroc), id = 'rcAuroc', w = widthOf(id), sm = w < 380;
    const labelW = sm ? 118 : 170, rowH = 14, top = 16, H = top + rows.length * rowH + 36;
    const s = svg(id, w, H, 'Dot plot of the AUROC of the Sepsis MetaScore in every study platform and tissue, with 95% intervals.');
    const lo = Math.min(0.5, Math.floor(Math.min(...rows.map(r => r.low == null ? r.auroc : r.low)) * 10) / 10), x = lin(lo, 1, labelW, w - 12);
    gridX(s, x, ticks(lo, 1, 4), top - 4, H - 32, v => fx(v, 1));
    axisTitle(s, (labelW + w) / 2, H - 4, 'AUROC, sepsis against healthy controls');
    rows.forEach((r, i) => {
      const y = top + i * rowH + rowH / 2, on = r.key === V.pick;
      if (on) S('rect', { x: 0, y: y - rowH / 2, width: w, height: rowH, fill: C.mark }, s);
      T(s, labelW - 6, y + 3.5, `${r.study}${r.tissue !== 'whole_blood' ? ` · ${tissueLab(r.tissue).replace('Blood, unspecified', 'blood')}` : ''}`, { 'text-anchor': 'end', 'font-size': 10, fill: on ? C.ink : C.ink2, 'font-weight': on ? 600 : 400 });
      if (r.low != null) S('line', { x1: x(r.low), x2: x(r.high), y1: y, y2: y, stroke: C.lineStrong, 'stroke-width': 1.5 }, s);
      S('circle', { cx: x(r.auroc), cy: y, r: on ? 4.5 : 3.5, fill: on ? C.cls.sepsis : C.ink2 }, s);
      const h = hitRow(s, 0, y - rowH / 2, w, rowH, `<b>${esc(plat(r.study, r.platform))}</b>, ${lower(tissueLab(r.tissue))}<br>AUROC ${fx(r.auroc, 3)}${r.low != null ? `, 95% interval ${fx(r.low, 3)}–${fx(r.high, 3)}` : ', too few patients for an interval'}<br>${fmt(r.sepsis)} sepsis, ${fmt(r.healthy)} healthy controls<br><span class="m">click to pick it out</span>`);
      h.style.cursor = 'pointer'; h.addEventListener('click', () => { V.pick = r.key; drawScore(); });
    });
  }
  function drawGeneMap() {
    const K = R.checks, id = 'rcGenes', w = widthOf(id), sm = isSmall(w), cols = K.strata.slice().sort((a, b) => b.auroc - a.auroc);
    const order = K.genes.map((g, i) => i).sort((a, b) => (K.published[a] === K.published[b] ? 0 : K.published[a] === 'higher_in_sepsis' ? -1 : 1));
    const labelW = sm ? 70 : 92, top = 8, rowH = sm ? 17 : 21, cw = (w - labelW - 8) / cols.length, H = top + order.length * rowH + 74;
    const s = svg(id, w, H, 'Heat map of the standardized difference between sepsis and healthy controls of every published gene in every study platform and tissue.');
    const lim = 3, colr = v => v == null ? C.band : mixHex(C.bg, v > 0 ? C.cls.sepsis : C.cls.healthy_control, Math.min(1, Math.abs(v) / lim) * 0.92);
    order.forEach((gi, i) => {
      const y = top + i * rowH, up = K.published[gi] === 'higher_in_sepsis';
      T(s, labelW - 12, y + rowH / 2 + 4, K.genes[gi], { 'text-anchor': 'end', 'font-size': sm ? 10 : 11.5, fill: C.ink, 'font-style': 'italic' });
      S('rect', { x: labelW - 8, y, width: 5, height: rowH - 1, fill: up ? C.cls.sepsis : C.cls.healthy_control }, s);
      cols.forEach((c, j) => {
        const v = c.diff[gi], r = S('rect', { x: labelW + j * cw, y, width: Math.max(1, cw - 0.8), height: rowH - 1, fill: colr(v) }, s);
        hover(r, `<b><i>${K.genes[gi]}</i></b> in ${esc(plat(c.study, c.platform))}, ${lower(tissueLab(c.tissue))}<br>standardized difference ${v == null ? 'not measured' : fx(v, 2)}<br><span class="m">published ${up ? 'higher' : 'lower'} in sepsis</span>`);
      });
    });
    const y0 = top + order.length * rowH + 6;
    cols.forEach((c, j) => { if (cw >= 9) T(s, labelW + (j + 0.5) * cw + 3, y0 + 4, c.study, { 'font-size': 8.5, fill: C.muted, transform: `rotate(60 ${labelW + (j + 0.5) * cw + 3} ${y0 + 4})` }); });
    if (cw < 9) T(s, labelW, y0 + 14, `${fmt(cols.length)} study platforms and tissues, from the highest AUROC (left) to the lowest`, { 'font-size': 11, fill: C.muted });
  }

  // ======================================================================= 6. the correction
  function drawCorrection() {
    const X = R.correction, bulk = N.corrected.bulk_blood, pbmc = N.corrected.pbmc;
    $('rxLead').innerHTML = `The per-study tables are not on one scale: two studies can give the same gene different values for technical reasons alone, and such batch effects are large enough to dominate any pooled analysis. For analyses across studies, two corrected tables were built, one for bulk blood (${fmt(bulk.samples)} samples of ${fmt(bulk.studies)} studies) and one for PBMC (${fmt(pbmc.samples)} samples of ${fmt(pbmc.studies)} studies). First every sample was given the same distribution of values (quantile normalization). Then studies with at least ${spell(X.minimum_healthy)} healthy controls were corrected on those controls: every gene was shifted and scaled so that the study's healthy controls match the healthy controls of the other such studies, with estimates that borrow strength across genes (empirical Bayes). Studies without enough healthy controls were corrected by matched groups: the mean and spread of each of their classes were set to those of the same class in the studies corrected on healthy controls.`;
    drawBeforeAfter(); drawFactors(); drawHeldOut(); drawMix();
    const md = X.held_out_medians.bulk_blood, pooled = X.pooled.find(p => p.group === 'bulk_blood' && p.scope === 'units_with_healthy_controls' && p.check === 'sepsis_score_auroc');
    $('rxCautions').innerHTML = `<b>Cautions.</b> Because matched groups set the class means of a study to references, part of the separation of classes across those studies follows from the correction itself; the sample table records how each study was corrected, so that an analysis can keep the studies corrected on healthy controls alone. In those studies the sepsis samples took no part in the correction, and the AUROC of sepsis against healthy controls across them rose from ${fx(pooled.before, 3)} to ${fx(pooled.after, 3)}. A residual difference between studies remains: the location error of matched groups in the held-out test, ${fx(md.matched_groups.location_error, 2)}, is several times the sampling floor of ${fx(md.matched_groups.sampling_floor, 2)}, and a pooled analysis should allow for it, for example with a random effect of study. The corrected tables are not inputs for classifiers of sepsis: their correction used the healthy controls of every study, so a model tested on a held-out study would have seen its controls through the correction. The per-study tables remain the reference for any result that depends on the correction.`;
  }
  const BA = { geo: null };
  function drawBeforeAfter() {
    const H0 = D.hero, n = H0.n, by = V.corrBy;
    const levels = by === 'fam' ? D.order.family : by === 'cls' ? D.order.class : D.order.tier, codes = by === 'fam' ? H0.fam : by === 'cls' ? H0.cls : H0.tier;
    const color = k => by === 'fam' ? lift(C.fam[k]) : by === 'cls' ? C.cls[k] : lift(C.tier[k]), label = k => lab(by === 'fam' ? 'family' : by === 'cls' ? 'class' : 'tier', k);
    const counts = levels.map((_, i) => 0); for (let i = 0; i < n; i++) counts[+codes[i]]++;
    $('rxKeys').innerHTML = keys(levels.filter((_, i) => counts[i] > 0).map(k => [label(k), color(k), 'dot']));
    const sprites = levels.map(k => dotSprite(color(k), 1.7, 0.66));
    const draw = (id, xa, ya, title, share) => {
      const cv = $(id), w = cv.parentNode.clientWidth || 420, fit = canvasFit(cv, { height: Math.round(Math.min(w, 560) * 0.82) }), { ctx } = fit;
      const X = xa.map(v => v / H0.scale), Y = ya.map(v => v / H0.scale), m = { l: 10, r: 10, t: 30, b: 26 };
      const ext = a => [Math.min(...a), Math.max(...a)], [x0, x1] = ext(X), [y0, y1] = ext(Y), k = Math.min((fit.w - m.l - m.r) / (x1 - x0), (fit.h - m.t - m.b) / (y1 - y0));
      const cx = (m.l + fit.w - m.r) / 2, cy = (m.t + fit.h - m.b) / 2, xm = (x0 + x1) / 2, ym = (y0 + y1) / 2, sx = v => cx + (v - xm) * k, sy = v => cy - (v - ym) * k;
      ctx.clearRect(0, 0, fit.w, fit.h);
      ctx.font = `600 13px ${getComputedStyle(document.body).fontFamily}`; ctx.fillStyle = C.ink; ctx.textAlign = 'left'; ctx.fillText(title, 6, 18);
      ctx.font = `11.5px ${getComputedStyle(document.body).fontFamily}`; ctx.fillStyle = C.muted; ctx.textAlign = 'center';
      ctx.fillText(`Component 1 (${pct(share[0], 0)} of the variation) · component 2 (${pct(share[1], 0)})`, fit.w / 2, fit.h - 6);
      const px = new Float32Array(n), py = new Float32Array(n);
      for (let i = 0; i < n; i++) { px[i] = sx(X[i]); py[i] = sy(Y[i]); stamp(ctx, sprites[+codes[i]], px[i], py[i]); }
      cv.setAttribute('aria-label', `Scatter plot of the ${fmt(n)} samples of the corrected bulk-blood table ${title.toLowerCase()}, colored by ${by === 'fam' ? 'technology' : by === 'cls' ? 'class' : 'how the study was corrected'}.`);
      cv.onmousemove = e => {
        const r = cv.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top; let best = -1, bd = 49;
        for (let i = 0; i < n; i++) { const d = (px[i] - x) ** 2 + (py[i] - y) ** 2; if (d < bd) { bd = d; best = i; } }
        if (best < 0) { untip(); return; }
        const c = D.order.class[+H0.cls[best]], f = D.order.family[+H0.fam[best]], t = D.order.tier[+H0.tier[best]];
        tip(e, `<b>${esc(H0.sample[best])}</b> <span class="m">· ${esc(H0.studies[H0.study[best]])}</span><br>${sw(C.cls[c])}${lab('class', c)} &nbsp;${sw(lift(C.fam[f]))}${lab('family', f)}<br><span class="m">${esc(labLower('tier', t))}</span>`);
      };
      cv.onmouseleave = untip;
    };
    draw('rxBefore', H0.x0, H0.y0, 'Before the correction', H0.variance.before);
    draw('rxAfter', H0.x1, H0.y1, 'After the correction', H0.variance.after);
  }
  function drawFactors() {
    const X = R.correction, rows = X.factors.filter(r => r.group === V.fgroup), id = 'rxFactors', w = widthOf(id), sm = isSmall(w);
    const order = ['family', 'study_accession', 'sample_class', 'sample_class_anchored'];
    const name = { family: 'Technology', study_accession: 'Study', sample_class: 'Class', sample_class_anchored: 'Class, in studies with healthy controls' };
    const labelW = sm ? 0 : 250, panelGap = 30, pw = sm ? w - 12 : (w - labelW - panelGap - 10) / 2, rowH = 28, top = 26;
    const H = sm ? 2 * (top + order.length * (rowH + 14) + 30) : top + order.length * rowH + 40;
    const s = svg(id, w, H, 'Arrows of the share of the variation of the first two components explained by technology, study and class, before and after the correction.');
    [['c1', 'v1', 'Component 1'], ['c2', 'v2', 'Component 2']].forEach(([c, v, title], pi) => {
      const ox = sm ? 6 : labelW + pi * (pw + panelGap), oy = sm ? pi * (H / 2) : 0, x = lin(0, 1, ox, ox + pw), b = rows.find(r => r.stage === 'before'), a = rows.find(r => r.stage === 'after');
      T(s, ox, oy + 14, `${title}: ${pct(b[v], 0)} of the variation before, ${pct(a[v], 0)} after`, { 'font-size': 12, fill: C.ink, 'font-weight': 600 });
      gridX(s, x, [0, 0.25, 0.5, 0.75, 1], oy + top - 4, oy + top + order.length * (sm ? rowH + 14 : rowH), v_ => pct(v_, 0));
      order.forEach((f, i) => {
        const y = oy + top + i * (sm ? rowH + 14 : rowH) + rowH / 2 + (sm ? 10 : 0);
        const before = (rows.find(r => r.stage === 'before' && r.factor === f) || {})[c], after = (rows.find(r => r.stage === 'after' && r.factor === f) || {})[c];
        if (before == null) return;
        if (sm) T(s, ox, y - 10, name[f], { 'font-size': 11, fill: C.ink2 }); else if (pi === 0) T(s, labelW - 10, y + 4, name[f], { 'text-anchor': 'end', 'font-size': 12.5, fill: C.ink });
        const col = after < before ? C.cls.healthy_control : C.cls.sepsis, x0 = x(before), x1 = x(after), dir = x1 > x0 ? 1 : -1;
        if (Math.abs(x1 - x0) > 8) {
          S('line', { x1: x0, x2: x1 - dir * 6, y1: y, y2: y, stroke: col, 'stroke-width': 2 }, s);
          S('path', { d: `M${x1},${y}L${x1 - dir * 8},${y - 4.5}L${x1 - dir * 8},${y + 4.5}Z`, fill: col }, s);
        }
        S('circle', { cx: x0, cy: y, r: 5, fill: C.bg, stroke: col, 'stroke-width': 1.8 }, s);
        S('circle', { cx: x1, cy: y, r: 5, fill: col }, s);
        hitRow(s, ox - 4, y - rowH / 2, pw + 8, rowH, `<b>${name[f]}</b>, ${title.toLowerCase()}<br>explains ${pct(before, 0)} of its variation before the correction and ${pct(after, 0)} after`);
      });
    });
  }
  function drawHeldOut() {
    const X = R.correction, g = V.hgroup, meas = V.hmeasure, rows = X.held_out.filter(r => r.group === g), how = $('rxHeldHow').lastElementChild;
    const methods = X.methods.filter(mm => rows.some(r => r.method === mm)), units = [...new Set(rows.map(r => `${r.study}|${r.platform}`))];
    const key = meas === 'loc' ? 'loc' : 'scale', id = 'rxHeld', w = widthOf(id), sm = isSmall(w), m = { l: sm ? 42 : 56, r: sm ? 14 : 150, t: 14, b: 50 }, H = sm ? 320 : 360;
    const vmax = niceMax(Math.max(...rows.map(r => r[key]), ...(meas === 'loc' ? rows.map(r => r.disp) : [])) * 1.04);
    const s = svg(id, w, H, `Slope chart of the ${meas === 'loc' ? 'location' : 'scale'} error of the hidden healthy controls under each method, one line per held-out study.`);
    const x = lin(0, methods.length - 1, m.l + 30, w - m.r - 30), y = lin(0, vmax, H - m.b, m.t);
    gridY(s, y, ticks(0, vmax, 5), m.l, w - m.r, v => fx(v, 1));
    axisTitle(s, 12, (m.t + H - m.b) / 2, meas === 'loc' ? 'Location error (standard deviations)' : 'Scale error', true);
    methods.forEach((mm, i) => { const lines = X.method_labels[mm].split(' '); lines.forEach((t, j) => T(s, x(i), H - m.b + 16 + j * 13, t, { 'text-anchor': 'middle', 'font-size': 11.5, fill: C.ink })); });
    const md = X.held_out_medians[g];
    if (meas === 'loc') methods.forEach((mm, i) => {
      const d = md[mm].pool_dispersion, b = rrect(s, x(i) - 9, y(d), 18, Math.max(1, y(0) - y(d)), tint(C.muted, 0.55), 2);
      hover(b, `<b>${X.method_labels[mm]}</b><br>pool dispersion ${fx(d, 2)}: how far apart it leaves the healthy controls of the studies it learns from (median)`);
    });
    if (meas === 'loc') {
      const fl = md[methods[0]].sampling_floor;
      S('line', { x1: m.l, x2: w - m.r, y1: y(fl), y2: y(fl), stroke: C.ink2, 'stroke-dasharray': '2 3', 'stroke-width': 1.3 }, s);
      if (sm) T(s, w - m.r, y(fl) - 5, `sampling floor ${fx(fl, 2)}`, { 'text-anchor': 'end', 'font-size': 10.5, fill: C.ink2, ...halo() });
      else T(s, w - m.r + 4, y(fl) + 4, `sampling floor ${fx(fl, 2)}`, { 'font-size': 11, fill: C.ink2 });
    }
    const lines = units.map(u => {
      const pts = methods.map((mm, i) => { const r = rows.find(q => `${q.study}|${q.platform}` === u && q.method === mm); return r ? [i, r[key], r] : null; }).filter(Boolean);
      const p = S('path', { d: pts.map(([i, v], k) => `${k ? 'L' : 'M'}${x(i)},${y(v)}`).join(''), fill: 'none', stroke: C.lineStrong, 'stroke-width': 1.1, 'stroke-opacity': 0.85 }, s);
      const r0 = pts[0][2];
      const tipHtml = `<b>${esc(plat(r0.study, r0.platform))}</b> held out, ${fmt(r0.hidden)} healthy controls hidden<br>${pts.map(([i, v]) => `${X.method_labels[methods[i]]}: ${fx(v, 2)}`).join('<br>')}`;
      hover(p, tipHtml);
      p.addEventListener('mouseenter', () => { p.setAttribute('stroke', C.accentInk); p.setAttribute('stroke-width', 2.2); });
      p.addEventListener('mouseleave', () => { p.setAttribute('stroke', C.lineStrong); p.setAttribute('stroke-width', 1.1); });
      return p;
    });
    void lines;
    const mk = meas === 'loc' ? 'location_error' : 'scale_error';
    S('path', { d: methods.map((mm, i) => `${i ? 'L' : 'M'}${x(i)},${y(md[mm][mk])}`).join(''), fill: 'none', stroke: C.ink, 'stroke-width': 2.6 }, s);
    methods.forEach((mm, i) => {
      const v = md[mm][mk], c = S('circle', { cx: x(i), cy: y(v), r: 5, fill: C.ink }, s);
      hover(c, `<b>${X.method_labels[mm]}</b><br>median ${meas === 'loc' ? 'location' : 'scale'} error ${fx(v, 2)} over ${plural(md[mm].units, 'held-out study platform')}`);
      if (!sm && i === methods.length - 1) T(s, x(i) + 10, y(v) + 4, `median ${fx(v, 2)}`, { 'font-size': 11.5, fill: C.ink, 'font-weight': 600 });
    });
    how.textContent = meas === 'loc'
      ? `Each study with healthy controls was held out in turn, its controls hidden, and corrected from its other samples by each method. A thin line is one held-out study, joining the distance of its hidden controls from the healthy controls of the other studies under each method, in within-study standard deviations; lower is better. The bold line is the median. The gray bar at each method is the dispersion it leaves among the studies it learns from, and the dotted line is the sampling floor, the error that chance alone gives. Hover over a line to see its study.`
      : `The same test, for the spread of the hidden controls: over the genes, the root mean square of the log ratio of their standard deviation to that of the other studies' controls; lower is better. A thin line is one held-out study and the bold line the median. Hover over a line to see its study.`;
    const cm = X.combat_vs_matched, nu = md.none.units;
    $('rxHeldText').innerHTML = g === 'bulk_blood'
      ? `Over the ${fmt(nu)} bulk-blood study platforms held out, the median location error falls from ${fx(md.none.location_error, 2)} without correction to ${fx(md.plain_combat.location_error, 2)} for plain ComBat and ${fx(md.matched_groups.location_error, 2)} for matched groups, against a sampling floor of ${fx(md.matched_groups.sampling_floor, 2)}. Matched groups applied to ${fmt(md.matched_groups.units)} of them and brought the hidden controls closer than no correction in every one. Plain ComBat has the lowest median, though not significantly lower than matched groups when the studies are compared one by one (median difference ${sgn(cm.difference, 2)}, 95% interval ${sgn(cm.low, 2)} to ${sgn(cm.high, 2)}; p = ${fx(cm.p, 3)}), and it leaves the studies it learns from ${fx(md.plain_combat.pool_dispersion, 2)} standard deviations apart, against ${fx(md.matched_groups.pool_dispersion, 2)} for matched groups: making every study alike on average, whatever its mix of patients and controls, removes part of the difference between sepsis and controls as if it were a study effect. Flat genes, a method that learns the study effects from the genes that sepsis changes least, was kept only as a comparison.`
      : `Only ${spell(nu)} PBMC study platforms have enough healthy controls to be held out. Their median location error falls from ${fx(md.none.location_error, 2)} without correction to ${fx(md.plain_combat.location_error, 2)} for plain ComBat and ${fx(md.matched_groups.location_error, 2)} for matched groups, against a sampling floor of ${fx(md.matched_groups.sampling_floor, 2)}; with so few studies these numbers are uncertain.`;
  }
  function drawMix() {
    const X = R.correction, how = $('rxMixHow').lastElementChild, keysEl = $('rxMixKeys');
    if (V.mix === 'mixing') {
      drawMixing();
      keysEl.innerHTML = keys([[lab('tier', 'anchored'), lift(C.tier.anchored), 'dot'], [lab('tier', 'unanchored'), lift(C.tier.unanchored), 'dot'], ['Failed a check that does not decide for its tier (kept)', C.ink, 'ring']]);
      how.textContent = `Every corrected study and class is one bubble, as large as its samples. Across, how much its samples clustered with each other before the correction; up, after it: the share of the ten nearest neighbors of its samples that come from the same study, less the share expected by chance. Lower is better mixed, and points below the dashed diagonal mixed better after the correction. A square marks a study whose mixing worsened by more than ${fx(X.mixing_limit, 1)}; for studies corrected on healthy controls this check does not decide, because they are the reference, so they stay in. Hover over a bubble to see its study.`;
      const mx = X.mixing, unan = mx.filter(r => r.tier === 'unanchored'), worse = mx.filter(r => !r.passed);
      $('rxMixNote').textContent = `All ${fmt(unan.length)} mixing checks of the studies corrected by matched groups passed, one for each study and class, and every applicable check of the sex genes and of the sepsis score passed; mixing worsened beyond the limit in ${spell(worse.length)} classes of studies corrected on healthy controls, which were kept. These checks catch gross errors only: lining up class means mixes studies by itself, so the held-out test, not these checks, is the test of the methods.`;
    } else {
      drawLedger();
      keysEl.innerHTML = '';
      const rec = X.records;
      how.textContent = `A flow diagram of all ${fmt(rec.total)} sample records: from the group they belong to by their sample type (bulk blood, PBMC or none) to their outcome in the correction. Each band is as wide as its records. Hover over a band or a bar to see its count.`;
      $('rxMixNote').textContent = `${fmt(rec.by_outcome.in_corrected_table)} records are in a corrected table and ${fmt(rec.by_outcome.corrected_for_checks_only)} are repeat deposits, corrected only to compare the two deposits of a sample. The other ${fmt(rec.by_outcome.left_out)} were left out, each for one reason; every record left out remains in the per-study tables.`;
    }
  }
  function drawMixing() {
    const X = R.correction, rows = X.mixing, id = 'rxMix', w = widthOf(id), sm = isSmall(w), m = { l: sm ? 44 : 56, r: 14, t: 12, b: 42 }, H = Math.min(520, Math.round(w * (sm ? 0.95 : 0.62)));
    const lo = Math.min(0, ...rows.map(r => Math.min(r.before, r.after))), hi = Math.max(...rows.map(r => Math.max(r.before, r.after))) * 1.04;
    const s = svg(id, w, H, 'Bubble chart of the study mixing of every corrected study and class, before against after the correction.');
    const x = lin(lo, hi, m.l, w - m.r), y = lin(lo, hi, H - m.b, m.t);
    gridX(s, x, ticks(lo, hi, 5), m.t, H - m.b, v => fx(v, 1)); gridY(s, y, ticks(lo, hi, 5), m.l, w - m.r, v => fx(v, 1));
    axisTitle(s, (m.l + w - m.r) / 2, H - 6, 'Mixing excess before the correction');
    axisTitle(s, 12, (m.t + H - m.b) / 2, 'After the correction', true);
    S('line', { x1: x(lo), y1: y(lo), x2: x(hi), y2: y(hi), stroke: C.lineStrong, 'stroke-dasharray': '5 4' }, s);
    const maxS = Math.max(...rows.map(r => r.samples)), rad = v => 2.5 + Math.sqrt(v / maxS) * (sm ? 11 : 16);
    rows.slice().sort((a, b) => b.samples - a.samples).forEach(r => {
      const col = lift(C.tier[r.tier]), cx = x(r.before), cy = y(r.after), rr = rad(r.samples);
      const node = r.passed ? S('circle', { cx, cy, r: rr, fill: col, 'fill-opacity': 0.4, stroke: col, 'stroke-width': 1 }, s)
        : S('rect', { x: cx - rr, y: cy - rr, width: 2 * rr, height: 2 * rr, fill: col, 'fill-opacity': 0.4, stroke: C.ink, 'stroke-width': 1.6 }, s);
      hover(node, `<b>${esc(plat(r.study, r.platform))}</b>, ${labLower('class', r.cls)} <span class="m">· ${r.group === 'pbmc' ? 'PBMC' : 'bulk blood'}</span><br>${fmt(r.samples)} samples; ${esc(labLower('tier', r.tier))}<br>mixing excess ${fx(r.before, 2)} before, ${fx(r.after, 2)} after${r.passed ? '' : '<br><b>worsened beyond the limit</b>; kept, as a study corrected on healthy controls'}${r.role === 'check_only' ? '<br><span class="m">its records are repeat deposits, corrected for the checks only</span>' : ''}`);
    });
  }
  function drawLedger() {
    const X = R.correction, gl = { bulk_blood: 'Bulk blood', pbmc: 'PBMC', not_grouped: 'Not grouped' }, gc = k => lift((X.group_colors || {})[k] || C.muted);
    const outs = [...new Set(X.flows.filter(f => f[0] === 'group_to_outcome').map(f => f[2]))];
    const ocol = k => k === 'in_corrected_table' ? lift(C.tier.anchored) : k === 'repeat_deposit' ? C.cls.other_condition : tint(C.muted, 0.35);
    const rec = X.records.total;
    const cols = [[{ key: 'records', label: 'Sample records', short: 'Records', color: C.ink2 }],
      ['bulk_blood', 'pbmc', 'not_grouped'].map(k => ({ key: k, label: gl[k], color: gc(k) })),
      outs.sort((a, b) => sum(X.flows.filter(f => f[2] === b).map(f => f[3])) - sum(X.flows.filter(f => f[2] === a).map(f => f[3]))).map(k => ({ key: k, label: X.outcome_labels[k], short: trunc(X.outcome_labels[k], 22), color: ocol(k) }))];
    const links = X.flows.map(([step, a, b, v]) => step === 'records_to_group'
      ? { c: 0, s: 'records', t: b, v, color: gc(b), tip: `<b>${gl[b]}</b><br>${fmt(v)} of the ${fmt(rec)} sample records` }
      : { c: 1, s: a, t: b, v, color: ocol(b), tip: `<b>${gl[a]} → ${X.outcome_labels[b]}</b><br>${fmt(v)} sample records` });
    sankey('rxMix', cols, links, { height: isSmall(widthOf('rxMix')) ? 380 : 420, left: 130, right: 300, leftS: 70, rightS: 150, unit: 'sample records',
      aria: `Flow diagram of the ${fmt(rec)} sample records from their group to their outcome in the correction.` });
  }

  // ======================================================================= 7. the benchmark
  function benchTask() { return R.benchmark.tasks[V.task]; }
  function drawBenchmark() {
    const B = R.benchmark, T_ = benchTask(), tasks = Object.keys(B.tasks);
    $('rbLead').innerHTML = `A classifier, a rule learned by a computer from examples, takes the profile of a sample and returns a label or a score. It is learned on training samples and judged on test samples that it has not seen. The benchmark fixes ${spell(tasks.length)} tasks on the per-study tables and the splits of their samples into training and test, so that methods are judged by the same test and the test cannot be passed for the wrong reasons. The main danger is leakage: the same patient, or two deposits of one sample, on both sides of a split, or a classifier that learns to recognize the study instead of the illness. Every split therefore holds out whole patient-safe groups, studies joined wherever they share a sample or a patient.`;
    const tbox = $('rbTask'); if (!tbox.dataset.on) { tbox.innerHTML = segHTML(tasks.map(t => [t, B.task_labels[t]]), V.task); tbox.dataset.on = '1'; seg(tbox, v => { V.task = v; V.fold = null; drawBenchmark(); }); }
    setSeg(tbox, V.task);
    const folds = sum(Object.values(T_.folds_by_scheme));
    $('rbCards').innerHTML = `
      <div class="card"><div class="big">${fmt(T_.eligible)}</div><div class="lab">samples of the task, of ${fmt(T_.records)} records</div></div>
      <div class="card"><div class="big">${fmt(T_.positives)} : ${fmt(T_.negatives)}</div><div class="lab">${posLab(T_)} against ${negLab(T_)}</div></div>
      <div class="card"><div class="big">${fmt(T_.groups)}</div><div class="lab">patient-safe groups, ${fmt(T_.scored_groups)} with at least ${spell(B.minimum_per_class)} samples of each class</div></div>
      <div class="card"><div class="big">${fmt(folds)}</div><div class="lab">folds: ${andList(Object.entries(T_.folds_by_scheme).map(([k, n]) => `${fmt(n)} ${lower(B.scheme_labels[k])}`))}</div></div>`;
    drawEligibility(); drawFoldControls(); drawFolds();
    $('rbWhy').innerHTML = `<b>Why the per-study tables, not the corrected ones.</b> The correction used the labels: a study with healthy controls was shifted so that its controls match the controls of the other studies, and in a study corrected by matched groups the mean of every class, sepsis included, was set to that of the same class elsewhere. The groups that a fold holds out took part in both, so in a corrected table the values of a test sample already depend on the labels of its own study, and a classifier tested on them is given part of the answer. The per-study tables were processed within each study without the curated labels. Any scaling across studies has to be fitted on the training groups of a fold, or applied to each sample on its own. The correction used no survival outcome, so this reason does not concern the mortality task. Score a model within each held-out group, not over the held-out groups pooled: they differ in their share of each class, and a pooled score would reward recognizing the study.`;
  }
  const posLab = t => t.definition.positive === 'non_survivor' ? 'deaths' : 'sepsis samples';
  const negLab = t => t.definition.negative === 'survivor' ? 'survivors' : CLASS_PLURAL[t.definition.negative] || t.definition.negative;
  function drawEligibility() {
    const B = R.benchmark, T_ = benchTask(); let left = T_.records;
    const rows = [];
    T_.excluded.filter(([, n]) => n > 0).forEach(([k, n]) => {
      rows.push({ label: B.reason_labels[k], end: fmt(left - n), parts: [
        { value: left - n, color: C.accentInk, tip: `<b>${fmt(left - n)} records</b> remain after this rule` },
        { value: n, color: tint(C.muted, 0.72), fill: tint(C.muted, 0.82), dash: true, tip: `<b>${B.reason_labels[k]}</b><br>${fmt(n)} records excluded` }] });
      left -= n;
    });
    if (left !== T_.eligible) throw new Error(`benchmark ${V.task}: the exclusions leave ${left}, not ${T_.eligible}`);
    rows.push({ label: 'Samples of the task', bold: true, end: fmt(T_.eligible), parts: [
      { value: T_.positives, color: C.cls.sepsis, text: cap(posLab(T_)), tip: `<b>${fmt(T_.positives)} ${posLab(T_)}</b>, the positive class` },
      { value: T_.negatives, color: T_.definition.negative === 'survivor' ? C.ink2 : C.cls[T_.definition.negative], text: '', tip: `<b>${fmt(T_.negatives)} ${negLab(T_)}</b>, the negative class` }] });
    hbars('rbSteps', rows, { max: T_.records, labelW: 330, rowH: 18, gap: 8, axis: true, axisLabel: 'Sample records', aria: `Bars of the records excluded by each eligibility rule of the task ${B.task_labels[V.task]}.` });
  }
  function drawFoldControls() {
    const B = R.benchmark, T_ = benchTask(), schemes = Object.keys(T_.folds);
    if (!schemes.includes(V.scheme)) V.scheme = schemes[0];
    const sbox = $('rbScheme');
    if (!sbox.dataset.on) { sbox.dataset.on = '1'; seg(sbox, v => { V.scheme = v; V.fold = null; drawFoldControls(); drawFolds(); }); }
    sbox.innerHTML = segHTML(schemes.map(k => [k, B.scheme_labels[k]]), V.scheme);
    const folds = T_.folds[V.scheme], groups = Object.fromEntries(T_.group_rows.map(g => [g.g, g]));
    if (!folds.some(f => f.fold === V.fold)) V.fold = (folds.find(f => f.fold === 'G024') || folds[0]).fold;
    const sel = $('rbFold');
    sel.innerHTML = folds.map(f => `<option value="${esc(f.fold)}"${f.fold === V.fold ? ' selected' : ''}>${esc(foldName(f, groups))}</option>`).join('');
    if (!sel.dataset.on) { sel.dataset.on = '1'; sel.addEventListener('change', () => { V.fold = sel.value; drawFolds(); }); }
    $('rbFoldBox').hidden = folds.length < 2;
  }
  function foldName(f, groups) {
    const what = { array_to_rna_seq: 'arrays train, RNA-seq is tested', rna_seq_to_array: 'RNA-seq trains, arrays are tested', earlier_to_later: `groups released before ${fdate(R.benchmark.split_date)} train` };
    if (V.scheme !== 'leave_one_group_out') return `${R.benchmark.scheme_labels[V.scheme]}: ${what[f.fold] || f.fold.replace(/_/g, ' ')}`;
    const g = groups[f.test[0]]; return `${f.fold} · ${g ? trunc(g.studies.split(' ').join(', '), 54) : ''}`;
  }
  function foldRoles() {
    const T_ = benchTask(), f = T_.folds[V.scheme].find(q => q.fold === V.fold);
    const role = g => f.test.includes(g) ? 'test' : f.not_scored.includes(g) ? 'not_scored' : 'train';
    const byRole = { train: [], test: [], not_scored: [] };
    T_.group_rows.forEach(g => byRole[role(g.g)].push(g));
    return { f, byRole };
  }
  function drawFolds() {
    const B = R.benchmark, T_ = benchTask(), { f, byRole } = foldRoles(), id = 'rbFolds', w = widthOf(id), sm = isSmall(w);
    const roles = [['train', 'Training', C.accentInk], ['test', 'Test', C.cls.sepsis], ['not_scored', 'Held out, not scored', C.muted]];
    const labelW = sm ? 0 : 170, rowH = 34, gap = sm ? 30 : 16, top = 6, H = top + roles.length * (rowH + gap) + 4;
    const s = svg(id, w, H, `Blocks of the patient-safe groups of the task, by their role in the fold ${V.fold}.`);
    const x = lin(0, T_.eligible, labelW, w - 70);
    roles.forEach(([k, name, col], i) => {
      const y = top + i * (rowH + gap) + (sm ? 16 : 0), list = byRole[k].slice().sort((a, b) => a.first.localeCompare(b.first)), pos = sum(list.map(g => g.pos)), neg = sum(list.map(g => g.neg));
      const label = `${name}: ${plural(list.length, 'group')}`;
      if (sm) T(s, 0, y - 5, label, { 'font-size': 12, fill: C.ink }); else { T(s, labelW - 10, y + 14, name, { 'text-anchor': 'end', 'font-size': 13, fill: C.ink, 'font-weight': 600 }); T(s, labelW - 10, y + 29, plural(list.length, 'group'), { 'text-anchor': 'end', 'font-size': 11.5, fill: C.muted }); }
      S('rect', { x: labelW, y, width: x(T_.eligible) - labelW, height: rowH, fill: 'none', stroke: C.line, 'stroke-dasharray': '2 3' }, s);
      let cx = labelW;
      list.forEach(g => {
        const gw = x(g.pos + g.neg) - labelW; if (gw <= 0) return;
        const pw = x(g.pos) - labelW;
        S('rect', { x: cx, y, width: pw, height: rowH, fill: col }, s);
        S('rect', { x: cx + pw, y, width: gw - pw, height: rowH, fill: tint(col, 0.55) }, s);
        S('rect', { x: cx, y, width: gw, height: rowH, fill: 'none', stroke: C.bg, 'stroke-width': 1 }, s);
        hitRow(s, cx, y, gw, rowH, `<b>${esc(g.g)}</b> <span class="m">· ${g.tech === 'rna_seq' ? 'RNA-seq' : 'arrays'}, first released ${fdate(g.first)}</span><br>${fmt(g.pos)} ${posLab(T_)}, ${fmt(g.neg)} ${negLab(T_)}<br><span class="m">${esc(trunc(g.studies.split(' ').join(', '), 160))}</span>`);
        if (gw > 44 && list.length < 4) T(s, cx + 5, y + rowH / 2 + 4, g.g, { 'font-size': 11, fill: onFill(col), 'pointer-events': 'none' });
        cx += gw;
      });
      T(s, cx + 6, y + rowH / 2 + 4, `${fmt(pos)} : ${fmt(neg)}`, { 'font-size': 11.5, fill: C.ink2 });
    });
    $('rbKeys').innerHTML = keys([['Training', C.accentInk], ['Test', C.cls.sepsis], ['Held out, not scored', C.muted], [`Lighter part: ${negLab(T_)}`, tint(C.ink2, 0.55)]]);
    const tr = byRole.train, te = byRole.test, ns = byRole.not_scored, sp = a => sum(a.map(g => g.pos)), sn = a => sum(a.map(g => g.neg));
    let text = '';
    if (V.scheme === 'leave_one_group_out' && te.length === 1) {
      const g = te[0], st = g.studies.split(' ');
      text = `This fold holds out ${g.g}, the patient-safe group of ${andList(st)}. It gives the task ${fmt(g.pos + g.neg)} samples, ${fmt(g.pos)} ${posLab(T_)} and ${fmt(g.neg)} ${negLab(T_)}; the ${fmt(tr.length)} other groups train, with ${fmt(sp(tr))} and ${fmt(sn(tr))}. As a check, ${fmt(sp(tr))} + ${fmt(g.pos)} = ${fmt(sp(tr) + g.pos)} and ${fmt(sn(tr))} + ${fmt(g.neg)} = ${fmt(sn(tr) + g.neg)}, the totals of the task.${g.g === 'G024' ? ` Had the split been made by study rather than by group, the repeat deposits of a later study of this group could have trained while the same arrays were tested through their primary records.` : ''}`;
    } else {
      text = `In this fold ${plural(tr.length, 'group')} train (${fmt(sp(tr))} ${posLab(T_)}, ${fmt(sn(tr))} ${negLab(T_)}) and ${plural(te.length, 'group')} are tested (${fmt(sp(te))} and ${fmt(sn(te))}), each scored on its own; ${plural(ns.length, 'further group')} ${ns.length === 1 ? 'is' : 'are'} held out but not scored, with fewer than ${spell(B.minimum_per_class)} samples of a class.${V.scheme === 'temporal' ? ` A group is dated by the earliest release among its study platforms.` : V.scheme === 'technology_transfer' ? ' No patient-safe group holds both technologies, so the two sides never share a sample.' : ''}`;
    }
    $('rbFoldText').innerHTML = text;
    void f;
  }

  // ======================================================================= 8. what the benchmark measures
  const EV = () => R.evaluation;
  const TWIN = () => ({ array: css('--twin-array'), patient: css('--twin-patient'), both: css('--twin-both') });
  const MODEL = () => ({ logistic_regression: css('--r-model-lr'), gradient_boosting: css('--r-model-gb'), sepsis_metascore: C.muted });
  const taskOrder = () => Object.keys(R.benchmark.task_labels);
  const taskLower = t => lower(R.benchmark.task_labels[t]);
  const evMean = (task, model, scheme) => (EV().summary.find(r => r.task === task && r.model === model && r.scheme === scheme) || {}).mean;
  const evModels = task => Object.keys(EV().labels.models).filter(m => EV().summary.some(r => r.task === task && r.model === m));
  const jit = (k, h) => ((k * 0.6180339887) % 1 - 0.5) * h;
  const EV_SCHEMES = ['naive_record_split', 'naive_leave_one_study_out', 'leave_one_group_out', 'technology_transfer', 'temporal'];
  function taskSeg(id, key, redraw) {
    const box = $(id);
    if (!box.dataset.on) { box.innerHTML = segHTML(taskOrder().map(t => [t, R.benchmark.task_labels[t]]), V[key]); box.dataset.on = '1'; seg(box, v => { V[key] = v; redraw(); }); }
    setSeg(box, V[key]);
  }
  function drawEvaluation() {
    const E = EV(), rec = E.records, all = E.leakage.find(r => r.pop === 'all_records' && r.scheme === 'naive_record_split');
    const sigs = Object.keys(E.signatures.about).length, ident = E.identification;
    const totals = E.effects.filter(r => r.effect === 'total_inflation'), top = totals.slice().sort((a, b) => b.estimate - a.estimate)[0];
    $('reLead').innerHTML = `The tasks and splits of the previous chapter were used to measure what pooling public studies without the links of the resource gets wrong. The records of each task were split as a naive analysis would split them, at random or one study at a time, and by the patient-safe folds. Two baseline classifiers with fixed settings, logistic regression and gradient boosting on the ${fmt(E.genes)} genes that best separate the classes of each training set, were trained and tested under every split. A predictor that knows only the study was scored beside them, and ${spell(sigs)} published signatures were scored again by their published formulas, within single study platforms.`;
    $('reCards').innerHTML = `
      <div class="card"><div class="big">${pct(rec.share_with_a_linked_record)}</div><div class="lab">of the ${fmt(rec.total)} records have a twin: another record of the same array or another sample of the same patient</div></div>
      <div class="card"><div class="big">${pct(all.share)} → 0</div><div class="lab">of test records with a twin in training: a random split of all records, and the ${fmt(E.patient_safe_folds)} patient-safe folds</div></div>
      <div class="card"><div class="big">${sgn(top.estimate, 2)}</div><div class="lab">mean AUROC that a random split adds over leave one group out, at most (${taskLower(top.task)}, ${lower(E.labels.models[top.model])})</div></div>
      <div class="card"><div class="big">${fx(ident.balanced_accuracy.aligned, 3)}</div><div class="lab">balanced accuracy of naming the study platform of a sample from its expression, among ${fmt(ident.units)}; guessing gives ${fx(ident.chance, 3)}</div></div>`;
    drawTwinShare(); drawModels(); drawEffects(); drawShortcuts(); drawSignatures();
  }
  // A. the share of test records with a twin in training, by split and kind of twin
  function drawTwinShare() {
    const E = EV(), id = 'reLeak', w = widthOf(id), sm = isSmall(w), col = TWIN(), L = E.labels.schemes;
    const pops = [['all_records', '', 'All records of the resource']].concat(taskOrder().map(t => ['naive', t, R.benchmark.task_labels[t]]));
    const rows = [];
    pops.forEach(([pop, task, name]) => {
      rows.push({ head: name });
      ['naive_record_split', 'naive_leave_one_study_out'].forEach(sc => { const r = E.leakage.find(x => x.pop === pop && x.task === task && x.scheme === sc); if (r) rows.push({ r, name, label: L[sc] }); });
      const g = E.leakage.filter(x => x.pop === 'gts' && x.task === task);
      if (task && g.length) rows.push({ gts: g, name, label: L.leave_one_group_out });
    });
    const labelW = sm ? 0 : 190, endW = sm ? 96 : 190, rowH = 15, gap = 8, headH = 24, top = 2, lab_ = sm ? 15 : 0;
    const H = top + sum(rows.map(r => r.head ? headH : rowH + gap + lab_)) + 34;
    const s = svg(id, w, H, 'Bars of the share of test records with a twin in training under every split of the records of the resource and of each benchmark task.');
    const max = niceMax(Math.max(...E.leakage.map(r => Math.max(r.share || 0, r.expected || 0))));
    const x = lin(0, max, labelW, w - endW);
    gridX(s, x, ticks(0, max, sm ? 3 : 5), top, H - 34, v => pct(v, 0));
    axisTitle(s, (labelW + w - endW) / 2, H - 4, 'Share of the test records with a twin in training');
    let y = top;
    rows.forEach(row => {
      if (row.head) { T(s, 0, y + 17, row.head, { 'font-size': 12.5, fill: C.ink, 'font-weight': 600 }); y += headH; return; }
      if (sm) { T(s, 0, y + 10, row.label, { 'font-size': 11.5, fill: C.ink2, ...halo() }); y += lab_; }
      else T(s, labelW - 8, y + rowH / 2 + 4, row.label, { 'text-anchor': 'end', 'font-size': 12, fill: C.ink2 });
      if (row.gts) {
        const g = row.gts.find(q => q.scheme === 'leave_one_group_out') || row.gts[0];
        S('line', { x1: x(0), x2: x(0), y1: y, y2: y + rowH, stroke: C.ink, 'stroke-width': 2 }, s);
        T(s, x(0) + 6, y + rowH / 2 + 4, `${pct(g.share)} · ${fmt(g.exposed)} of ${fmt(g.test)} test samples, patient-safe`, { 'font-size': 11.5, fill: C.ink2 });
        hitRow(s, 0, y - lab_, w, rowH + lab_, `<b>${esc(row.name)} · the patient-safe folds</b><br>${row.gts.map(q => `${esc(L[q.scheme])}: ${fmt(q.exposed)} of ${fmt(q.test)} test samples, ${plural(q.folds, 'fold')}`).join('<br>')}`);
      } else {
        const r = row.r, parts = [['array', r.exposed - r.patient, 'the same array only'], ['both', r.array + r.patient - r.exposed, 'both'], ['patient', r.exposed - r.array, 'another sample of the patient only']];
        let cx = x(0);
        parts.forEach(([k, v]) => { const pw = x(v / r.test) - x(0); if (pw <= 0) return; S('rect', { x: cx, y, width: pw, height: rowH, fill: col[k], stroke: C.bg, 'stroke-width': 1 }, s); cx += pw; });
        if (r.expected != null) { const ex = x(r.expected); S('line', { x1: ex, x2: ex, y1: y - 3, y2: y + rowH + 3, stroke: C.ink, 'stroke-width': 1.6 }, s); }
        T(s, Math.max(cx, r.expected != null ? x(r.expected) : 0) + 6, y + rowH / 2 + 4, `${pct(r.share)} · ${fmt(r.exposed)} of ${fmt(r.test)}`, { 'font-size': 11.5, fill: C.ink2 });
        hitRow(s, 0, y - lab_, w, rowH + lab_, `<b>${esc(row.name)} · ${esc(row.label)}</b><br>${fmt(r.exposed)} of ${fmt(r.test)} test records (${pct(r.share)}) over ${plural(r.folds, 'fold')}<br>` +
          parts.map(([k, v, l]) => `${sw(col[k])}${l}: ${fmt(v)}`).join('<br>') + (r.expected != null ? `<br><span class="m">Expected by chance from the twins of each record: ${pct(r.expected)}</span>` : ''));
      }
      y += rowH + gap;
    });
    $('reLeakKeys').innerHTML = keys([['The same array', col.array], ['Both', col.both], ['Another sample of the patient', col.patient], ['Expected by chance', C.ink, 'line']]);
    const rec = E.records, loso = E.leakage.find(r => r.pop === 'all_records' && r.scheme === 'naive_leave_one_study_out'), all = E.leakage.find(r => r.pop === 'all_records' && r.scheme === 'naive_record_split');
    $('reLeakNote').innerHTML = `${fmt(rec.with_a_linked_record)} of the ${fmt(rec.total)} records of the resource have a twin, and ${fmt(rec.arrays_in_more_than_one_study)} arrays were deposited in more than one study, in ${fmt(rec.studies_sharing_an_array)} studies. A random split of all records leaves ${pct(all.share)} of its test records with a twin in training, as many as chance predicts. Leaving one study out lowers the share to ${pct(loso.share)} but does not remove it, because the twins of a record often sit in another study of its patient-safe group. The patient-safe folds hold out whole groups and leave none; the animation <a href="#story/twins/1">Twins and splits</a> in How it works shows the three splits on the records of one task.`;
  }
  // B. the baseline classifiers: the AUROC of every scored unit under every scheme
  const unitName = u => u.scheme === 'naive_record_split' ? `random split ${+u.fold.replace(/\D/g, '')}` : u.scheme === 'naive_leave_one_study_out' ? `held-out study ${u.unit}`
    : u.scheme === 'technology_transfer' ? `held-out group ${u.unit} (${u.fold.replace(/_/g, ' ').replace('rna seq', 'RNA-seq')})` : `held-out group ${u.unit}`;
  function drawModels() {
    taskSeg('reModelTask', 'mtask', drawModels);
    const E = EV(), task = V.mtask, id = 'reModels', w = widthOf(id), sm = isSmall(w), col = MODEL(), models = evModels(task), L = E.labels;
    const units = E.units.filter(u => u.task === task);
    const lo = Math.min(0.5, Math.floor(Math.min(...units.map(u => u.auroc)) * 10) / 10);
    const labelW = sm ? 0 : 170, endW = sm ? 44 : 64, laneH = 13, rowGap = 12, groupH = 22, top = 2, lab_ = sm ? 15 : 0;
    const groups = [['Naive splits', EV_SCHEMES.slice(0, 2)], ['Patient-safe splits', EV_SCHEMES.slice(2)]];
    const H = top + groups.length * groupH + EV_SCHEMES.length * (models.length * laneH + rowGap + lab_) + 34;
    const s = svg(id, w, H, `Dot plot of the AUROC of the baseline classifiers in every scored test unit, by split, for ${taskLower(task)}.`);
    const x = lin(lo, 1, labelW, w - endW);
    gridX(s, x, ticks(lo, 1, sm ? 4 : 6), top, H - 34, v => fx(v, 1));
    S('line', { x1: x(0.5), x2: x(0.5), y1: top, y2: H - 34, stroke: C.lineStrong, 'stroke-dasharray': '3 3' }, s);
    axisTitle(s, (labelW + w - endW) / 2, H - 4, `AUROC in each test unit, ${taskLower(task)}`);
    let y = top;
    groups.forEach(([gname, schemes]) => {
      T(s, 0, y + 16, gname, { 'font-size': 12.5, fill: C.ink, 'font-weight': 600 }); y += groupH;
      schemes.forEach(sc => {
        const sum_ = E.summary.filter(r => r.task === task && r.scheme === sc), first = sum_[0] || {};
        const label = `${L.schemes[sc]}`, sub = `${plural(first.scored || 0, sc === 'naive_record_split' ? 'split' : 'unit')} scored`;
        if (sm) { T(s, 0, y + 10, `${label} · ${sub}`, { 'font-size': 11.5, fill: C.ink2, ...halo() }); y += lab_; }
        else { T(s, labelW - 8, y + 10, label, { 'text-anchor': 'end', 'font-size': 12, fill: C.ink }); T(s, labelW - 8, y + 23, sub, { 'text-anchor': 'end', 'font-size': 10.5, fill: C.muted }); }
        models.forEach((m, k) => {
          const ly = y + k * laneH + laneH / 2, mine = units.filter(u => u.model === m && u.scheme === sc), mean = evMean(task, m, sc);
          mine.forEach((u, j) => {
            const c = S('circle', { cx: x(u.auroc), cy: ly + jit(j + 1, laneH * 0.55), r: sm ? 2.7 : 3.2, fill: col[m], 'fill-opacity': C.dark ? 0.75 : 0.6, stroke: C.bg, 'stroke-width': 0.8 }, s);
            hover(c, `<b>${esc(L.models[m])}</b> · ${esc(L.schemes[sc])}<br>${esc(unitName(u))}: AUROC ${fx(u.auroc, 3)}${u.low != null ? ` (95% interval ${fx(u.low, 3)}–${fx(u.high, 3)})` : ''}<br>${fmt(u.pos)} ${posLab(R.benchmark.tasks[task])}, ${fmt(u.neg)} ${negLab(R.benchmark.tasks[task])}`);
          });
          if (mean != null) {
            S('line', { x1: x(mean), x2: x(mean), y1: ly - laneH / 2 + 1, y2: ly + laneH / 2 - 1, stroke: m === 'sepsis_metascore' ? C.ink2 : col[m], 'stroke-width': 3, 'stroke-linecap': 'round' }, s);
            T(s, w - endW + 8, ly + 4, fx(mean, 3), { 'font-size': 11, fill: C.ink2 });
          }
        });
        hitRow(s, 0, y - lab_, labelW, models.length * laneH + lab_, `<b>${esc(L.schemes[sc])}</b><br>${models.map(m => `${sw(col[m])}${esc(L.models[m])}: mean AUROC ${fx(evMean(task, m, sc), 3)} over ${plural((sum_.find(r => r.model === m) || {}).scored || 0, 'unit')}`).join('<br>')}`);
        y += models.length * laneH + rowGap;
      });
    });
    $('reModelKeys').innerHTML = keys(models.map(m => [L.models[m] + (m === 'sepsis_metascore' ? ' (no training)' : ''), col[m], 'dot'])) + keys([['Mean of the units', C.ink2, 'line']]);
    $('reModelHow').lastElementChild.textContent = `One row for every split; each dot is one test unit: one of the ${spell(E.repeats)} random splits, scored over all its test records together; one held-out study; or one held-out patient-safe group. The bar is the mean of the units, also given on the right. Hover over a dot for its AUROC and interval.`;
    const m = (mm, sc) => fx(evMean(task, mm, sc), 3);
    $('reModelText').innerHTML = models.filter(mm => mm !== 'sepsis_metascore').map(mm => `${E.labels.models[mm]} averaged <b>${m(mm, 'naive_record_split')}</b> under the random split, ${m(mm, 'naive_leave_one_study_out')} leaving one study out and <b>${m(mm, 'leave_one_group_out')}</b> leaving one patient-safe group out; ${m(mm, 'technology_transfer')} across technologies and ${m(mm, 'temporal')} on the later groups.`).join(' ') +
      (models.includes('sepsis_metascore') ? ` The Sepsis MetaScore, which needs no training, gives ${m('sepsis_metascore', 'naive_record_split')} over the pooled test records of the random splits and ${m('sepsis_metascore', 'leave_one_group_out')} within the held-out groups.` : '');
  }
  // C. the total difference beside the four effects, each measured on its own comparison
  function drawEffects() {
    taskSeg('reEffectTask', 'etask', drawEffects);
    const E = EV(), task = V.etask, id = 'reEffects', w = widthOf(id), sm = isSmall(w), col = MODEL(), L = E.labels;
    const effects = Object.keys(L.effects), rows = E.effects.filter(r => r.task === task), models = evModels(task).filter(m => rows.some(r => r.model === m));
    const vals = rows.flatMap(r => [r.estimate, r.low, r.high]).filter(v => v != null).concat([0]);
    const span = Math.max(...vals) - Math.min(...vals) || 0.01, lo = Math.min(...vals) - span * 0.08, hi = Math.max(...vals) + span * 0.08;
    const labelW = sm ? 0 : 190, endW = sm ? 64 : 84, laneH = 15, rowGap = 14, top = 6, lab_ = sm ? 15 : 0;
    const H = top + effects.length * (models.length * laneH + rowGap + lab_) + 10 + 34;
    const s = svg(id, w, H, `The total difference in mean AUROC between the random split and leave one group out, and four separately measured effects with their 95% intervals, for ${taskLower(task)}.`);
    const x = lin(lo, hi, labelW, w - endW);
    gridX(s, x, ticks(lo, hi, sm ? 4 : 6), top, H - 34, v => sgn(v, Math.abs(hi - lo) < 0.02 ? 3 : 2));
    S('line', { x1: x(0), x2: x(0), y1: top, y2: H - 34, stroke: C.ink2, 'stroke-width': 1.2 }, s);
    axisTitle(s, (labelW + w - endW) / 2, H - 4, `Difference in AUROC, ${taskLower(task)}`);
    let y = top;
    effects.forEach((ef, i) => {
      const [name, about] = L.effects[ef];
      if (i === 1) { S('line', { x1: 0, x2: w, y1: y - rowGap / 2 + 2, y2: y - rowGap / 2 + 2, stroke: C.lineStrong }, s); y += 10; }
      if (sm) { T(s, 0, y + 10, name, { 'font-size': 11.5, fill: C.ink, 'font-weight': i ? 400 : 600, ...halo() }); y += lab_; }
      else T(s, labelW - 8, y + models.length * laneH / 2 + 4, name, { 'text-anchor': 'end', 'font-size': 12, fill: C.ink, 'font-weight': i ? 400 : 600 });
      models.forEach((m, k) => {
        const r = rows.find(q => q.effect === ef && q.model === m); if (!r) return;
        const ly = y + k * laneH + laneH / 2;
        if (r.low != null) S('line', { x1: x(r.low), x2: x(r.high), y1: ly, y2: ly, stroke: col[m], 'stroke-width': 2, 'stroke-linecap': 'round' }, s);
        if (ef === 'total_inflation') S('path', { d: `M${x(r.estimate)},${ly - 5.5}l5.5,5.5l-5.5,5.5l-5.5,-5.5z`, fill: col[m], stroke: C.bg, 'stroke-width': 1.5 }, s);
        else S('circle', { cx: x(r.estimate), cy: ly, r: 4.4, fill: col[m], stroke: C.bg, 'stroke-width': 1.5 }, s);
        const n = rows.filter(q => q.effect === ef).map(q => q.studies);
        if (k === 0 || n.some(v => v !== n[0])) T(s, w - endW + 8, n.every(v => v === n[0]) ? y + models.length * laneH / 2 + 4 : ly + 4, plural(r.studies, 'study', 'studies'), { 'font-size': 11, fill: C.muted });
        hitRow(s, labelW, ly - laneH / 2, w - labelW, laneH, `<b>${esc(name)}</b> · ${esc(L.models[m])}<br><span class="m">${esc(about)}</span><br>${fx(r.first, 3)} − ${fx(r.second, 3)} = <b>${sgn(r.estimate, 3)}</b>` +
          (r.low != null ? `<br>${pct(E.level, 0)} interval ${sgn(r.low, 3)} to ${sgn(r.high, 3)}` : '') + `<br>${plural(r.studies, 'study', 'studies')}${r.matched ? ', the same studies on both sides' : ''}`);
      });
      y += models.length * laneH + rowGap;
    });
    $('reEffectKeys').innerHTML = keys(models.map(m => [L.models[m], col[m], 'dot']));
    $('reEffectHow').lastElementChild.textContent = `The first row is the total difference: the mean AUROC under the random split minus that under leave one group out. The rows below it are ${spell(effects.length - 1)} effects, each measured on its own comparison and its own studies, with its ${pct(E.level, 0)} interval: they are not parts of the total and do not add up to it. Hover over a mark for the comparison and the number of studies.`;
    const r = (ef, m) => rows.find(q => q.effect === ef && q.model === m), both = ef => andList(models.map(m => `${sgn(r(ef, m).estimate, 3)}`));
    $('reEffectText').innerHTML = `Under the random split the mean AUROC is <b>${both('total_inflation')}</b> higher than leaving one group out (${andList(models.map(m => lower(L.models[m])))}). Measured one at a time, each on the studies where it can be measured: ${effects.slice(1).map(ef => `${lower(L.effects[ef][0])} ${both(ef)}`).join('; ')}.`;
  }
  // D. the study-only predictor and the study identification
  function drawShortcuts() {
    const E = EV(), so = E.study_only, tasks = taskOrder().filter(t => so[t]), id = 'reStudyOnly', w = widthOf(id), sm = w < 420;
    const labelW = sm ? 0 : 128, laneH = 14, rowGap = 12, top = 24, lab_ = sm ? 15 : 0, endW = 44;
    const H = top + tasks.length * (2 * laneH + rowGap + 16) + 34;
    const s = svg(id, w, H, 'The AUROC of a predictor that knows only the study, over the pooled test records of the random splits and within every held-out study.');
    T(s, 0, 13, 'A predictor that knows only the study', { 'font-size': 12.5, fill: C.ink, 'font-weight': 600 });
    const lo = Math.min(0.4, Math.floor(Math.min(...E.pooled.map(p => p.auroc)) * 10) / 10), x = lin(lo, 1, labelW, w - endW);
    gridX(s, x, ticks(lo, 1, sm ? 3 : 4), top, H - 34, v => fx(v, 1));
    axisTitle(s, (labelW + w - endW) / 2, H - 4, 'AUROC');
    let y = top;
    tasks.forEach(t => {
      T(s, 0, y + 11, R.benchmark.task_labels[t], { 'font-size': 11.5, fill: C.ink, 'font-weight': 600 }); y += 16;
      [['pooled', 'Pooled', E.pooled.filter(p => p.task === t), so[t].random_split.pooled_auroc_mean],
       ['within', 'Within each study', E.within.filter(p => p.task === t), so[t].held_out_study.within_study_auroc_mean]].forEach(([k, name, pts, mean], j) => {
        const ly = y + j * laneH + laneH / 2;
        if (!sm) T(s, labelW - 8, ly + 4, name, { 'text-anchor': 'end', 'font-size': 11.5, fill: C.ink2 });
        else T(s, x(mean) + 8, ly + 4, name, { 'font-size': 10.5, fill: C.muted, ...halo() });
        pts.forEach((p, q) => {
          const c = S('circle', { cx: x(p.auroc), cy: ly + jit(q + 1, laneH * 0.6), r: 3, fill: k === 'pooled' ? C.ink2 : C.muted, 'fill-opacity': 0.55, stroke: C.bg, 'stroke-width': 0.8 }, s);
          hover(c, k === 'pooled' ? `<b>Random split ${+p.fold.replace(/\D/g, '')}</b>, ${esc(taskLower(t))}<br>AUROC ${fx(p.auroc, 3)} over its ${fmt(p.n)} test records` : `<b>Held-out study ${esc(p.study)}</b>, ${esc(taskLower(t))}<br>AUROC ${fx(p.auroc, 3)}: all its ${fmt(p.n)} records get one score`);
        });
        S('line', { x1: x(mean), x2: x(mean), y1: ly - laneH / 2 + 1, y2: ly + laneH / 2 - 1, stroke: C.ink, 'stroke-width': 3, 'stroke-linecap': 'round' }, s);
        T(s, w - endW + 6, ly + 4, fx(mean, 3), { 'font-size': 11, fill: C.ink2 });
      });
      y += 2 * laneH + rowGap;
    });
    // the study identification
    const I = E.identification, id2 = 'reIdent', w2 = widthOf(id2), sm2 = w2 < 420, tables = [['aligned', 'Values the classifiers read'], ['corrected', 'Corrected table']];
    const labelW2 = sm2 ? 0 : 150, rowH = 18, top2 = 24, lab2 = sm2 ? 15 : 0, H2 = top2 + tables.length * (rowH + 18 + lab2) + 44;
    const s2 = svg(id2, w2, H2, 'Bars of the balanced accuracy of naming the study platform of a sample from its aligned and from its corrected expression values, against guessing.');
    T(s2, 0, 13, 'Naming the study platform from the expression', { 'font-size': 12.5, fill: C.ink, 'font-weight': 600 });
    const x2 = lin(0, 1, labelW2, w2 - 48);
    gridX(s2, x2, [0, 0.25, 0.5, 0.75, 1], top2, H2 - 34, v => fx(v, 2));
    axisTitle(s2, (labelW2 + w2 - 48) / 2, H2 - 4, `Balanced accuracy, ${fmt(I.units)} study platforms`);
    let y2 = top2 + 6;
    tables.forEach(([k, name]) => {
      if (sm2) { T(s2, x2(I.chance) + 6, y2 + 10, name, { 'font-size': 11.5, fill: C.ink2 }); y2 += lab2; }
      else T(s2, labelW2 - 8, y2 + rowH / 2 + 4, name, { 'text-anchor': 'end', 'font-size': 11.5, fill: C.ink2 });
      const v = I.balanced_accuracy[k], folds = I.folds.filter(f => f.table === k && f.fold !== 'all');
      rrect(s2, x2(0), y2, x2(v) - x2(0), rowH, k === 'aligned' ? C.ink2 : tint(C.ink2, 0.45), 3);
      folds.forEach((f, q) => { const c = S('circle', { cx: x2(f.balanced), cy: y2 + rowH / 2 + jit(q + 1, rowH * 0.5), r: 3, fill: C.bg, stroke: C.ink, 'stroke-width': 1.2 }, s2); hover(c, `<b>${esc(name)}</b>, fold ${+f.fold.replace(/\D/g, '')}<br>balanced accuracy ${fx(f.balanced, 3)} over ${fmt(f.test)} test samples`); });
      T(s2, x2(Math.max(v, ...folds.map(f => f.balanced))) + 8, y2 + rowH / 2 + 4, fx(v, 3), { 'font-size': 11.5, fill: C.ink2, ...halo() });
      hitRow(s2, 0, y2 - lab2, labelW2 || w2, rowH + lab2, `<b>${esc(name)}</b><br>balanced accuracy ${fx(v, 3)} over ${fmt(I.samples)} samples in ${spell(folds.length)} folds; guessing ${fx(I.chance, 3)}`);
      y2 += rowH + 18;
    });
    S('line', { x1: x2(I.chance), x2: x2(I.chance), y1: top2, y2: H2 - 34, stroke: C.ink2, 'stroke-width': 1.5, 'stroke-dasharray': '4 3' }, s2);
    T(s2, x2(I.chance) + 5, H2 - 36, `guessing ${fx(I.chance, 3)}`, { 'font-size': 11, fill: C.ink2 });
    const t0 = tasks.map(t => `${fx(so[t].random_split.pooled_auroc_mean, 3)} (${taskLower(t)})`);
    $('reShortNote').textContent = `The predictor gives every test record the share of sepsis, or of deaths, among the training records of its study, and knows nothing else. Over the pooled test records of the random splits it reaches ${andList(t0)}, because the studies differ in their share of each class; within each held-out study every record gets one score, and the AUROC is ${fx(so[tasks[0]].held_out_study.within_study_auroc_mean, 1)}. A classifier can also tell the studies apart: from the aligned values that the baseline classifiers read, the study platform of a sample is named with a balanced accuracy of ${fx(I.balanced_accuracy.aligned, 3)}, and after the correction with ${fx(I.balanced_accuracy.corrected, 3)}, against ${fx(I.chance, 3)} by guessing. Scoring within each held-out study or group, as the benchmark does, removes this shortcut from the score.`;
  }
  // E. the published signatures, all strata against the independent ones
  function drawSignatures() {
    const E = EV(), G = E.signatures, scope = V.sigScope, id = 'reSigs', w = widthOf(id), sm = isSmall(w), L = E.labels;
    const comps = [...new Set(G.order.map(([, c]) => c))];
    const strata = G.strata.filter(r => scope === 'all_tissues' || r.scope);
    const med = (sig, c, set) => G.summary.find(r => r.sig === sig && r.comparison === c && r.scope === scope && r.set === set);
    const lo = Math.min(0.3, Math.floor(Math.min(...G.strata.map(r => r.auroc)) * 10) / 10);
    const labelW = sm ? 0 : 190, endW = sm ? 0 : 170, rowH = 26, headH = 24, top = 2, lab_ = sm ? 30 : 0;
    const H = top + comps.length * headH + G.order.length * (rowH + lab_ + 6) + 34;
    const s = svg(id, w, H, 'Dot plot of the within-study AUROC of every tested stratum of each published signature, with the medians over all strata and over the independent strata.');
    const x = lin(lo, 1, labelW, w - endW - 8);
    gridX(s, x, ticks(lo, 1, sm ? 4 : 7), top, H - 34, v => fx(v, 1));
    S('line', { x1: x(0.5), x2: x(0.5), y1: top, y2: H - 34, stroke: C.lineStrong, 'stroke-dasharray': '3 3' }, s);
    axisTitle(s, (labelW + w - endW) / 2, H - 4, 'AUROC within each study platform and sample type');
    let y = top;
    comps.forEach(c => {
      T(s, 0, y + 17, R.benchmark.task_labels[c], { 'font-size': 12.5, fill: C.ink, 'font-weight': 600 }); y += headH;
      const tk = R.benchmark.tasks[c];
      G.order.filter(([, cc]) => cc === c).forEach(([sig]) => {
        const about = G.about[sig], all = med(sig, c, 'all'), ind = med(sig, c, 'independent'), mine = strata.filter(r => r.sig === sig && r.comparison === c);
        const text = `${fx(all.median, 3)} (${fmt(all.tested)}) → ${fx(ind.median, 3)} (${fmt(ind.tested)})`;
        if (sm) { T(s, 0, y + 10, `${L.signatures[sig]} · ${about.cite}`, { 'font-size': 11.5, fill: C.ink, ...halo() }); T(s, 0, y + 24, `median ${text}`, { 'font-size': 10.5, fill: C.muted, ...halo() }); y += lab_; }
        else {
          T(s, labelW - 8, y + 11, L.signatures[sig], { 'text-anchor': 'end', 'font-size': 12, fill: C.ink });
          T(s, labelW - 8, y + 24, about.cite, { 'text-anchor': 'end', 'font-size': 10.5, fill: C.muted });
          T(s, w - endW + 4, y + rowH / 2 + 4, text, { 'font-size': 11.5, fill: C.ink2 });
        }
        const cy = y + rowH / 2;
        mine.forEach((r, q) => {
          const not = r.role === 'discovery' ? 'one of the studies that derived it' : r.role === 'validation' ? 'one of the studies that validated it' : r.shared ? `shares ${plural(r.shared, 'sample')} with those studies` : r.repeats ? `holds ${plural(r.repeats, 'repeat deposit')}` : 'not independent';
          const c_ = S('circle', { cx: x(r.auroc), cy: cy + jit(q + 1, rowH * 0.55), r: sm ? 2.8 : 3.3, fill: r.independent ? C.ink2 : C.bg, 'fill-opacity': r.independent ? 0.6 : 1, stroke: C.ink2, 'stroke-width': r.independent ? 0.6 : 1.2 }, s);
          hover(c_, `<b>${esc(plat(r.study, r.platform))}</b>, ${esc(lower(tissueLab(r.tissue)))}<br>${esc(L.signatures[sig])}: AUROC ${fx(r.auroc, 3)}${r.low != null ? ` (95% interval ${fx(r.low, 3)}–${fx(r.high, 3)})` : ''}<br>${fmt(r.pos)} ${posLab(tk)}, ${fmt(r.neg)} ${negLab(tk)}<br>${r.independent ? 'Independent stratum' : `<span class="m">Not independent: ${esc(not)}</span>`}`);
        });
        S('line', { x1: x(all.median), x2: x(all.median), y1: cy - 11, y2: cy + 11, stroke: C.ink, 'stroke-width': 2, 'stroke-linecap': 'round' }, s);
        S('line', { x1: x(ind.median), x2: x(ind.median), y1: cy - 12, y2: cy + 12, stroke: C.accentInk, 'stroke-width': 3.5, 'stroke-linecap': 'round' }, s);
        const skip = G.skipped.filter(([s_, c2]) => s_ === sig && c2 === c);
        hitRow(s, 0, y - lab_, labelW || w, rowH + lab_, `<b>${esc(about.label)}</b> <span class="m">· ${esc(about.cite)}, ${plural(about.genes, 'gene')}</span><br>Median AUROC ${fx(all.median, 3)} over ${fmt(all.tested)} strata, ${fx(ind.median, 3)} over the ${fmt(ind.tested)} independent ones` +
          (skip.length ? `<br><span class="m">Skipped over all sample types: ${skip.map(([, , k, n]) => `${lower(L.reasons[k])} (${fmt(n)})`).join('; ')}</span>` : ''));
        y += rowH + 6;
      });
    });
    $('reSigKeys').innerHTML = keys([['Independent stratum', C.ink2, 'dot'], ['Not independent', C.ink2, 'ring'], ['Median, all strata', C.ink, 'line'], ['Median, independent strata', C.accentInk, 'line']]);
    const rows = G.order.map(([sig, c]) => ({ sig, c, all: med(sig, c, 'all'), ind: med(sig, c, 'independent') }));
    const drop = rows.slice().sort((a, b) => (b.all.median - b.ind.median) - (a.all.median - a.ind.median))[0];
    $('reSigText').innerHTML = `Over all its strata the ${L.signatures[drop.sig]} has a median AUROC of <b>${fx(drop.all.median, 3)}</b> in ${fmt(drop.all.tested)} strata (${taskLower(drop.c)}); over the ${fmt(drop.ind.tested)} independent strata, <b>${fx(drop.ind.median, 3)}</b>. ` +
      `The other signatures, all strata against independent ones: ${rows.filter(r => r !== drop).map(r => `${L.signatures[r.sig]}, ${taskLower(r.c)}, ${fx(r.all.median, 3)} and ${fx(r.ind.median, 3)}`).join('; ')}.` +
      (scope === 'all_tissues' ? '' : ` Bulk blood only: ${andList(E.scope_tissues.map(t => lower(tissueLab(t))))}.`);
  }

  // ======================================================================= 9. getting the data
  const bytes = b => b >= 1e9 ? `${fx(b / 1e9, 1)} GB` : b >= 1e6 ? `${fx(b / 1e6, b >= 1e8 ? 0 : 1)} MB` : b >= 1e3 ? `${fx(b / 1e3, 0)} kB` : `${fmt(b)} bytes`;
  function drawData() {
    const P = R.package, total = sum(P.files.map(f => f.bytes));
    $('rdLead').innerHTML = `The data package holds every table of the resource as compressed CSV files inside archives, with a data dictionary that defines every column of every table. It is organized so that samples are chosen first and analyzed second: the sample table says which samples to keep, and their study and sample accessions lead to their columns in the expression tables of their study. Every distinct sample is identified by its primary record, so keeping the primary records keeps one copy of every array or library, and keeping the samples marked one per patient keeps one sample of every patient. The per-study tables are released as processed, each on its own scale; the corrected tables are separate, with their own sample tables.`;
    $('rdCards').innerHTML = `
      <div class="card"><div class="lab">Repository</div><div class="big">${esc(P.repository)}</div><div class="lab">DOI: ${P.doi ? `<a href="https://doi.org/${esc(P.doi)}">${esc(P.doi)}</a>` : '<span id="rdDoi">to be assigned</span>'}</div></div>
      <div class="card"><div class="lab">License of the data</div><div class="big"><a href="${esc(P.license_url)}">${esc(P.license.replace('CC0-1.0', 'CC0 1.0'))}</a></div><div class="lab">public domain: use it for any purpose, without asking</div></div>
      <div class="card"><div class="lab">Version ${esc(P.version)}</div><div class="big">${bytes(total)}</div><div class="lab">in ${plural(P.files.length, 'file')}; ${fmt(P.counts.sample_records_with_values)} sample records of ${fmt(P.counts.included_studies)} studies</div></div>
      <div class="card"><div class="lab">Code (MIT License)</div><div class="big r-code"><a href="${esc(P.code)}">GitHub</a></div><div class="lab">${esc(P.code.replace('https://', ''))}</div></div>`;
    drawFiles();
    $('rdQuestions').innerHTML = `<thead><tr><th>Question</th><th>Where to look</th><th>What it holds</th></tr></thead><tbody>${P.questions.map(q => `<tr><td>${esc(q.q)}</td><td><code>${esc(q.archive)}</code>${q.member ? `<br><code class="m">${esc(q.member)}</code>` : ''}${q.rows ? `<br><span class="ui">${fmt(q.rows)} rows</span>` : ''}</td><td>${esc(q.about)}</td></tr>`).join('')}</tbody>`;
    $('rdNote').textContent = `For most analyses, start from the sample table: choose the classes, sample types and age groups of interest, keep the primary records and, where patients matter, the one sample marked for each patient, with or without the samples that carry a quality flag. Analyses within studies use the per-study tables, and their results can be combined across studies by meta-analysis, which leaves each study on its own scale. Pooled analyses use the corrected tables, preferably the studies corrected on healthy controls, with study as a random effect. The expression values were deposited by the investigators of the source studies: please cite the studies whose samples you use, together with the resource.`;
  }
  function drawFiles() {
    const P = R.package, rows = P.files.slice().sort((a, b) => b.bytes - a.bytes), id = 'rdFiles', w = widthOf(id), sm = isSmall(w);
    const labelW = sm ? 0 : 190, rowH = 18, gap = sm ? 22 : 7, top = 4, H = top + rows.length * (rowH + gap) + 34;
    const s = svg(id, w, H, 'Bar chart of the sizes of the files of the data package, on a logarithmic scale.');
    const x = lg(100, 1e10, labelW, w - (sm ? 60 : 74));
    [1e3, 1e6, 1e9].forEach(v => { S('line', { x1: x(v), x2: x(v), y1: top, y2: H - 30, stroke: C.grid }, s); T(s, x(v), H - 14, v === 1e3 ? '1 kB' : v === 1e6 ? '1 MB' : '1 GB', { 'text-anchor': 'middle', 'font-size': 11, fill: C.muted }); });
    rows.forEach((f, i) => {
      const y = top + i * (rowH + gap) + (sm ? 16 : 0), archive = /\.tar(\.gz)?$/.test(f.name);
      if (sm) T(s, 0, y - 4, f.name, { 'font-size': 11.5, fill: C.ink, 'font-family': 'ui-monospace, Consolas, monospace' });
      else T(s, labelW - 8, y + rowH / 2 + 4, f.name, { 'text-anchor': 'end', 'font-size': 12, fill: C.ink, 'font-family': 'ui-monospace, Consolas, monospace' });
      rrect(s, x(100), y, x(f.bytes) - x(100), rowH, archive ? C.accentInk : tint(C.accentInk, 0.5), 2);
      T(s, x(f.bytes) + 6, y + rowH / 2 + 4, `${bytes(f.bytes)}${f.members ? ` · ${fmt(f.members)} tables or files` : ''}`, { 'font-size': 11, fill: C.ink2 });
      hitRow(s, 0, y - (sm ? 14 : 0), w, rowH + (sm ? 14 : 0), `<b>${esc(f.name)}</b> <span class="m">· ${bytes(f.bytes)}${f.members ? `, ${fmt(f.members)} members` : ''}</span><br>${esc(f.about)}`);
    });
  }

  // ======================================================================= the view
  const DRAW = { 'r-studies': drawStudies, 'r-samples': drawSamples, 'r-repeats': drawRepeats, 'r-genes': drawGenes, 'r-checks': drawChecks,
    'r-correction': drawCorrection, 'r-benchmark': drawBenchmark, 'r-evaluation': drawEvaluation, 'r-data': drawData };
  function init() {
    R = part('resource');
    const pages = $('resource').querySelector('.pages');
    pages.querySelectorAll('section.sec').forEach(s => s.remove());
    pages.querySelector('.pager').insertAdjacentHTML('beforebegin', chapterHTML());
    const wire = (id, key, redraw) => seg(id, v => { V[key] = v; redraw(); });
    wire('rsYearBy', 'year', drawYears);
    wire('raRecBy', 'rec', drawRecords);
    wire('rgMapBy', 'map', drawMapping);
    wire('rgBy', 'genes', drawGeneChart);
    wire('rcRepBy', 'rep', drawRepeatPairs);
    wire('rcScoreBy', 'score', drawScore);
    wire('rxBy', 'corrBy', drawBeforeAfter);
    wire('rxFactorBy', 'fgroup', drawFactors);
    wire('rxHeldBy', 'hgroup', drawHeldOut);
    wire('rxHeldMeasure', 'hmeasure', drawHeldOut);
    wire('rxMixBy', 'mix', drawMix);
    wire('reSigScope', 'sigScope', drawSignatures);
  }
  registerView('resource', {
    part: 'resource',
    init,
    show(route) { const f = DRAW[route.chapter]; if (f) f(); },
  });

  // ---------- checks (?selftest=1) ----------
  const ok = (cond, detail) => ({ ok: !!cond, detail });
  selftest('resource', 'the funnel of the studies adds up, and ends at the included studies of the core', () => {
    const St = part('resource').studies, nc = sum(St.not_candidate.map(r => r[1])), ex = sum(St.excluded.map(r => r[1]));
    return ok(St.records_found - nc === St.candidates && St.candidates - ex === St.included && St.included === N.studies && St.candidates === N.candidates && St.records_found === N.census_records,
      `${fmt(St.records_found)} − ${fmt(nc)} = ${fmt(St.records_found - nc)} candidates; − ${fmt(ex)} = ${fmt(St.candidates - ex)} included, core ${fmt(N.studies)}`);
  });
  selftest('resource', 'the study platforms by year hold every study platform and sample record, and the studies by year every study', () => {
    const St = part('resource').studies, recs = sum(St.platforms.map(p => p.records)), st = sum(St.studies_by_year.map(r => r.array + r.rna_seq));
    const lastYear = Math.max(...St.by_year.map(r => r.year)), added = sum(St.by_year.filter(r => r.year === lastYear).map(r => r.added_up)), samp = sum(St.by_year.map(r => r.samples));
    return ok(St.platforms.length === N.study_platforms && recs === N.records && st === N.studies && added === N.distinct && samp === N.distinct,
      `${fmt(St.platforms.length)} platforms, ${fmt(recs)} records, ${fmt(st)} studies, ${fmt(samp)} samples by year`);
  });
  selftest('resource', 'the flow from technology to sample type to class holds every distinct sample, class by class', () => {
    const A = part('resource').samples, a = sum(A.tech_material.map(r => r[2])), b = sum(A.material_class.map(r => r[2]));
    const bad = D.order.class.filter(k => sum(A.material_class.filter(r => r[1] === k).map(r => r[2])) !== N.classes[k]);
    return ok(a === N.distinct && b === N.distinct && !bad.length, `${fmt(a)} and ${fmt(b)} against ${fmt(N.distinct)}${bad.length ? `; classes differ: ${bad.join(', ')}` : ''}`);
  });
  selftest('resource', 'the age groups and the sexes add up to the distinct samples; what a class rests on adds up to the records', () => {
    const A = part('resource').samples, ages = sum(Object.values(A.age_groups).flatMap(m => Object.values(m))), sexes = sum(Object.values(A.sex));
    const basis = sum(Object.values(A.basis).flatMap(m => Object.values(m)));
    return ok(ages === N.distinct && sexes === N.distinct && basis === N.records, `ages ${fmt(ages)}, sexes ${fmt(sexes)}, bases ${fmt(basis)}`);
  });
  selftest('resource', 'shock and outcome cover every distinct sepsis sample, and the rules of GSE65682 every one of its records', () => {
    const A = part('resource').samples, c = sum(A.clinical.map(r => r[2])), e = sum(A.example.rules.map(r => r.samples));
    return ok(c === N.classes.sepsis && e === A.example.records, `${fmt(c)} against ${fmt(N.classes.sepsis)}; ${fmt(e)} against ${fmt(A.example.records)}`);
  });
  selftest('resource', 'records, distinct samples and patients add up, and match the core', () => {
    const P = part('resource').repeats, st = P.steps;
    return ok(st.included_sample_records === N.included_records && st.records_without_values === N.records_without_values && st.sample_records_with_values === N.records
      && st.included_sample_records - st.records_without_values === st.sample_records_with_values && st.sample_records_with_values - st.repeat_deposits === st.distinct_samples
      && st.distinct_samples === N.distinct && st.one_sample_per_patient === st.identified_patients + st.samples_without_identifier
      && P.patients.identified === N.patients && P.patients.one_per_patient.total === N.one_per_patient && st.one_sample_per_patient === N.one_per_patient && P.repeated_samples === N.repeated_samples,
      `${fmt(st.included_sample_records)} − ${fmt(st.records_without_values)} − ${fmt(st.repeat_deposits)} = ${fmt(st.distinct_samples)}; ${fmt(st.one_sample_per_patient)} = ${fmt(st.identified_patients)} + ${fmt(st.samples_without_identifier)}`);
  });
  selftest('resource', 'the cohort group of the example holds the records and distinct samples that the text gives', () => {
    const P = part('resource').repeats, rows = P.platforms.filter(p => p.group === P.example.group), recs = sum(rows.map(p => p.records)), prim = sum(rows.map(p => p.primary));
    return ok(rows.length === P.example.platforms && recs === P.example.records && prim === P.example.distinct, `${rows.length} platforms, ${fmt(recs)} records, ${fmt(prim)} primary records against ${fmt(P.example.distinct)} distinct samples`);
  });
  selftest('resource', 'every row of every table is mapped once, and there is one gene-level table per table of the core', () => {
    const G = part('resource').genes, rows = sum(Object.values(G.mapping).flatMap(m => Object.values(m)));
    return ok(rows === G.rows && G.tables.length === N.tables && G.scales.length === N.tables && G.shared.length === N.tables && Math.max(...G.tables.map(t => t.genes)) === G.genes.most,
      `${fmt(rows)} rows against ${fmt(G.rows)}; ${G.tables.length} tables`);
  });
  selftest('resource', 'the sex tables and the sex plane agree with the validation stage', () => {
    const K = part('resource').checks, c = sum(K.sex_tables.map(r => r.concordant)), n = sum(K.sex_tables.map(r => r.compared)), marked = K.sex_tables.filter(r => r.marked).length;
    const plane = K.sex.x.length === K.sex.y.length && K.sex.rec.length === K.sex.x.length && K.sex.agr.length === K.sex.x.length;
    return ok(c === K.sex_receipt.concordant && n === K.sex_receipt.compared && n === N.sex_compared && marked === N.sex_unreliable_tables && plane && Math.abs(c / n - N.sex_concordance) < 1e-4,
      `${fmt(c)} of ${fmt(n)} agree; ${marked} tables marked`);
  });
  selftest('resource', 'the repeat pairs are those of the validation stage, with their closest matches', () => {
    const K = part('resource').checks, Pr = K.repeats, best = [...Pr.best].filter(v => v === '1').length, sp = Pr.sp.slice().sort((a, b) => a - b);
    return ok(Pr.sp.length === K.repeat_receipt.pairs && best === K.repeat_receipt.most_similar && Math.abs(sp[0] - K.repeat_receipt.spearman_lowest) < 6e-4 && Math.abs(median(sp) - K.repeat_receipt.spearman_median) < 6e-4,
      `${fmt(Pr.sp.length)} pairs, ${fmt(best)} closest matches, lowest ${fx(sp[0], 4)}, median ${fx(median(sp), 4)}`);
  });
  selftest('resource', 'the area under every ROC curve equals the AUROC of the score table', () => {
    const K = part('resource').checks, bad = K.strata.filter(r => Math.abs(area(r.fpr, r.tpr) - r.table_auroc) > 2e-3 || Math.abs(r.auroc - r.table_auroc) > 1e-4);
    return ok(K.strata.length === K.score_receipt.strata && !bad.length, bad.length ? bad.slice(0, 3).map(r => `${r.key}: ${fx(area(r.fpr, r.tpr), 4)} against ${fx(r.table_auroc, 4)}`).join('; ') : `${K.strata.length} curves`);
  });
  selftest('resource', 'the median AUROC and the published directions match the validation stage', () => {
    const K = part('resource').checks, md = median(K.strata.map(r => r.auroc)), flags = K.strata.map(r => r.pub).join('');
    const cells = flags.replace(/-/g, '').length, asPub = [...flags].filter(v => v === '1').length;
    // the sign of every difference agrees with its flag, where it does not round to zero
    const signs = K.strata.every(r => r.diff.every((v, i) => v == null ? r.pub[i] === '-' : Math.abs(v) < 0.01 || ((v > 0) === (K.published[i] === 'higher_in_sepsis')) === (r.pub[i] === '1')));
    return ok(Math.abs(md - K.score_receipt.median) < 6e-4 && cells === K.score_receipt.comparisons && asPub === K.score_receipt.as_published && signs,
      `median ${fx(md, 4)}; ${asPub} of ${cells} in the published direction`);
  });
  selftest('resource', 'the held-out medians recomputed from the figure points equal those of the held-out table', () => {
    const X = part('resource').correction, bad = [];
    Object.entries(X.held_out_medians).forEach(([g, by]) => Object.entries(by).forEach(([mm, v]) => {
      const rows = X.held_out.filter(r => r.group === g && r.method === mm);
      [['loc', 'location_error'], ['scale', 'scale_error'], ['disp', 'pool_dispersion'], ['floor', 'sampling_floor']].forEach(([a, b]) => {
        const m_ = median(rows.map(r => r[a])); if (Math.abs(m_ - v[b]) > 2e-4 || rows.length !== v.units) bad.push(`${g} ${mm} ${b}: ${fx(m_, 4)} against ${fx(v[b], 4)}`);
      });
    }));
    return ok(!bad.length, bad.slice(0, 3).join('; ') || `${Object.keys(X.held_out_medians).length} groups`);
  });
  selftest('resource', 'every record has one outcome in the correction, and the corrected tables hold the samples of the core', () => {
    const X = part('resource').correction, out = sum(X.outcomes.map(r => r[1])), flows = sum(X.flows.filter(f => f[0] === 'records_to_group').map(f => f[3]));
    const inT = sum(X.outcomes.filter(r => /^in_table/.test(r[0])).map(r => r[1]));
    return ok(out === N.records && flows === N.records && inT === N.corrected.bulk_blood.samples + N.corrected.pbmc.samples && inT === X.records.by_outcome.in_corrected_table,
      `${fmt(out)} outcomes, ${fmt(flows)} in the flow, ${fmt(inT)} in the tables`);
  });
  selftest('resource', 'every fold puts every group on one side, and its sides add up to the totals of the task', () => {
    const B = part('resource').benchmark, bad = []; let folds = 0;
    Object.entries(B.tasks).forEach(([t, T_]) => {
      const pos = sum(T_.group_rows.map(g => g.pos)), neg = sum(T_.group_rows.map(g => g.neg));
      if (pos !== T_.positives || neg !== T_.negatives || T_.group_rows.length !== T_.groups || pos + neg !== T_.eligible) bad.push(`${t}: groups give ${pos} and ${neg}`);
      if (T_.records - sum(T_.excluded.map(r => r[1])) !== T_.eligible) bad.push(`${t}: exclusions do not add up`);
      Object.entries(T_.folds).forEach(([sc, list]) => list.forEach(f => {
        folds++;
        const names = new Set(T_.group_rows.map(g => g.g)), side = [...f.test, ...f.not_scored];
        if (new Set(side).size !== side.length || side.some(g => !names.has(g)) || !f.test.length) bad.push(`${t} ${sc} ${f.fold}`);
      }));
      Object.entries(T_.folds_by_scheme).forEach(([sc, n]) => { if ((T_.folds[sc] || []).length !== n) bad.push(`${t} ${sc}: ${n} folds expected`); });
    });
    return ok(!bad.length && folds === N.bench_folds, bad.slice(0, 3).join('; ') || `${folds} folds`);
  });
  selftest('resource', 'the worked fold of G024 gives the counts of the guide\'s example: training and test add up to the task', () => {
    const B = part('resource').benchmark, T_ = B.tasks.sepsis_vs_healthy_control, f = T_.folds.leave_one_group_out.find(q => q.fold === 'G024');
    if (!f) return ok(false, 'no fold G024');
    const test = T_.group_rows.filter(g => f.test.includes(g.g)), train = T_.group_rows.filter(g => !f.test.includes(g.g) && !f.not_scored.includes(g.g));
    return ok(sum(test.map(g => g.pos)) + sum(train.map(g => g.pos)) === T_.positives && sum(test.map(g => g.neg)) + sum(train.map(g => g.neg)) === T_.negatives,
      `test ${sum(test.map(g => g.pos))} : ${sum(test.map(g => g.neg))}, ${train.length} training groups ${sum(train.map(g => g.pos))} : ${sum(train.map(g => g.neg))}`);
  });
  // the benchmark of evaluation designs
  selftest('resource', 'benchmark: the kinds of twin add up to the exposed test records of every split, and no patient-safe fold has one', () => {
    const E = part('resource').evaluation, bad = [];
    E.leakage.forEach(r => {
      const parts = [r.exposed - r.patient, r.array + r.patient - r.exposed, r.exposed - r.array];
      if (parts.some(v => v < 0) || sum(parts) !== r.exposed || Math.abs(r.exposed / r.test - r.share) > 1e-6) bad.push(`${r.pop} ${r.task} ${r.scheme}`);
      if (r.pop === 'gts' && r.exposed) bad.push(`${r.task} ${r.scheme} exposes ${r.exposed}`);
    });
    const folds = sum(E.leakage.filter(r => r.pop === 'gts').map(r => r.folds));
    return ok(!bad.length && folds === E.patient_safe_folds, bad.slice(0, 3).join('; ') || `${E.leakage.length} splits; ${folds} patient-safe folds with no twin in training`);
  });
  selftest('resource', 'benchmark: the mean AUROC of every classifier and split, recomputed from its scored units, equals the summary', () => {
    const E = part('resource').evaluation, bad = E.summary.filter(r => {
      const u = E.units.filter(x => x.task === r.task && x.model === r.model && x.scheme === r.scheme);
      return u.length !== r.scored || Math.abs(mean(u.map(x => x.auroc)) - r.mean) > 2e-6;
    });
    return ok(!bad.length && E.units.length === sum(E.summary.map(r => r.scored)), bad.slice(0, 3).map(r => `${r.task} ${r.model} ${r.scheme}`).join('; ') || `${E.summary.length} means from ${fmt(E.units.length)} units`);
  });
  selftest('resource', 'benchmark: the total difference is the random split minus leave one group out, and every effect is drawn on its own row', () => {
    const E = part('resource').evaluation, bad = E.effects.filter(r => r.effect === 'total_inflation' &&
      Math.abs(r.estimate - (evMean(r.task, r.model, 'naive_record_split') - evMean(r.task, r.model, 'leave_one_group_out'))) > 2e-6);
    const intervals = E.effects.filter(r => r.effect !== 'total_inflation' && r.status === 'estimable' && !(r.low <= r.high && r.studies >= E.minimum_studies));
    const drawn = $('reEffects').querySelectorAll('circle, path').length, want = E.effects.filter(r => r.task === V.etask).length;
    return ok(!bad.length && !intervals.length && drawn === want, `${bad.length} totals and ${intervals.length} intervals wrong; ${drawn} marks drawn for ${want} effects of ${V.etask}`);
  });
  selftest('resource', 'benchmark: the study-only predictor scores 0.5 in every held-out study, and its pooled means and the study identification match the stage', () => {
    const E = part('resource').evaluation, bad = [];
    Object.entries(E.study_only).forEach(([t, so]) => {
      const p = E.pooled.filter(x => x.task === t), w = E.within.filter(x => x.task === t);
      if (p.length !== so.random_split.repeats || Math.abs(mean(p.map(x => x.auroc)) - so.random_split.pooled_auroc_mean) > 2e-6) bad.push(`${t} pooled`);
      if (w.length !== so.held_out_study.studies_scored || w.some(x => x.auroc !== so.held_out_study.within_study_auroc_mean)) bad.push(`${t} within`);
    });
    const I = E.identification, tested = ['aligned', 'corrected'].map(k => sum(I.folds.filter(f => f.table === k && f.fold !== 'all').map(f => f.test)));
    if (tested.some(v => v !== I.samples) || Math.abs(I.chance - 1 / I.units) > 1e-6) bad.push('identification');
    return ok(!bad.length, bad.join('; ') || `${Object.keys(E.study_only).length} tasks; ${fmt(I.samples)} samples in every table`);
  });
  selftest('resource', 'benchmark: the medians of every signature, recomputed from its strata, equal the summary for all and independent strata and both scopes', () => {
    const G = part('resource').evaluation.signatures, bad = G.summary.filter(r => {
      const m = G.strata.filter(x => x.sig === r.sig && x.comparison === r.comparison && (r.scope === 'all_tissues' || x.scope) && (r.set === 'all' || x.independent));
      return m.length !== r.tested || Math.abs(median(m.map(x => x.auroc)) - r.median) > 2e-6;
    });
    return ok(!bad.length && G.summary.length === 4 * G.order.length, bad.slice(0, 3).map(r => `${r.sig} ${r.comparison} ${r.scope} ${r.set}`).join('; ') || `${G.summary.length} medians from ${G.strata.length} strata`);
  });
  selftest('resource', 'benchmark: the colors of the twins and of the two classifiers are set in both themes and stand out from the page', () => {
    const before = themePref, out = [];
    try {
      ['light', 'dark'].forEach(theme => {
        applyTheme(theme, { persist: false, redraw: false });
        const cols = [...Object.values(TWIN()), css('--r-model-lr'), css('--r-model-gb')];
        if (cols.some(c => !/^#[0-9a-f]{6}$/i.test(c)) || new Set(cols).size !== cols.length) out.push(`${theme}: ${cols.join(' ')}`);
        cols.forEach(c => { if (/^#/.test(c) && Math.hypot(...hexRGB(c).map((v, k) => v - hexRGB(C.bg)[k])) < 90) out.push(`${theme}: ${c} is close to the page`); });
      });
    } finally { applyTheme(before, { persist: false, redraw: false }); }
    return ok(!out.length, out.join('; ') || 'five colors in each theme');
  });
  selftest('resource', 'the data package lists every file with its size, and shows the DOI as to be assigned until it is', () => {
    const P = part('resource').package, files = P.files.length, sized = P.files.every(f => f.bytes > 0 && f.about);
    const shown = P.doi ? true : !!$('rdDoi') && $('rdDoi').textContent === 'to be assigned';
    return ok(files >= 10 && sized && shown && P.counts.sample_records_with_values === N.records && P.counts.included_sample_records === N.included_records
      && P.counts.distinct_samples === N.distinct && P.questions.length > 0,
      `${files} files${P.doi ? `, DOI ${P.doi}` : ', DOI to be assigned'}`);
  });
})();
