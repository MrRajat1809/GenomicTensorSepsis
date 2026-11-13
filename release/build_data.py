"""Write the data of the results page (release/index.html) from the outputs of the pipeline.

The page shows only numbers written here. data.js is the core: the main numbers, the labels and colors of classes,
technologies and sample types, the glossary and the start page. Each view has its own part, which the page loads when
the view is first opened: data/story.js, data/resource.js, data/methods.js, and data/studies/index.js with one file
data/studies/<study>.js per study. Every file assigns window.DATA_PARTS["<part>"].

Sources: the tables in data/processed/, the review tables, receipts and QC reports in outputs/, the points of the
paper's figures in outputs/visualization/figure_data.csv and the stage configurations in src/. A missing file, key or
column stops the build with a message naming it.

Run from the repository root, on the host or in the pipeline's image (pandas and pyarrow are the only packages used):
    python release/build_data.py                every part, then stamp the page's scripts with content versions
    python release/build_data.py core methods   the named parts only
"""
import argparse
import hashlib
import json
import math
import os
import random
import re
import tempfile
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROCESSED = ROOT / 'data' / 'processed'
OUTPUTS = ROOT / 'outputs'
SRC = ROOT / 'src'

# Short stage names (as in data/processed/) with their folder in src/, their folder in outputs/ and their configuration.
STAGES = {
    'acquisition': ('01_data_acquisition', 'data_acquisition', 'acquisition_config.json'),
    'expression': ('02_expression_processing', 'expression_processing', 'expression_config.json'),
    'annotation': ('03_gene_annotation', 'gene_annotation', 'annotation_config.json'),
    'cohort': ('04_cohort_curation', 'cohort_curation', 'cohort_config.json'),
    'validation': ('05_technical_validation', 'technical_validation', 'validation_config.json'),
    'correction': ('06_batch_correction', 'batch_correction', 'correction_config.json'),
    'benchmark': ('07_benchmark', 'benchmark', 'benchmark_config.json'),
    'visualization': ('08_visualization', 'visualization', 'visualization_config.json'),
}

# Size limits of the page's files, in bytes: the core is loaded with the page, a part when its view is opened.
CORE_LIMIT, PART_LIMIT = 1_000_000, 1_500_000


class BuildError(Exception):
    """A source of the page data is missing or does not hold what the page needs."""


# ---------- readers ----------
def _stage(stage):
    if stage not in STAGES:
        raise BuildError(f'unknown stage {stage!r}; expected one of {", ".join(STAGES)}')
    return STAGES[stage]


def _existing(path):
    if not path.is_file():
        raise BuildError(f'{path.relative_to(ROOT).as_posix()} is missing; run the pipeline stage that writes it')
    return path


def _columns_present(frame, columns, where):
    missing = [c for c in columns or [] if c not in frame.columns]
    if missing:
        raise BuildError(f'{where} has no column {", ".join(missing)}')


def table(name, columns=None):
    """A table of data/processed/, e.g. table('cohort/cohort_samples', ['sample_id', 'sample_class'])."""
    path = _existing(PROCESSED / f'{name}.parquet')
    present = pq.read_schema(path).names
    missing = [c for c in columns or [] if c not in present]
    if missing:
        raise BuildError(f'{path.relative_to(ROOT).as_posix()} has no column {", ".join(missing)}')
    return pd.read_parquet(path, columns=columns)


def review(stage, name, columns=None):
    """A review table of outputs/, e.g. review('cohort', 'benchmark') reads outputs/cohort_curation/benchmark_review.csv."""
    path = _existing(OUTPUTS / _stage(stage)[1] / f'{name}_review.csv')
    frame = pd.read_csv(path, low_memory=False)
    _columns_present(frame, columns, path.relative_to(ROOT).as_posix())
    return frame[columns] if columns else frame


def _json(path):
    return json.loads(_existing(path).read_text(encoding='utf-8'))


def receipt(stage, number):
    """The receipt of one script, e.g. receipt('cohort', 5) reads outputs/cohort_curation/05.json; it must be complete."""
    path = OUTPUTS / _stage(stage)[1] / f'{number:02d}.json'
    data = _json(path)
    if data.get('status') != 'COMPLETE':
        raise BuildError(f'{path.relative_to(ROOT).as_posix()} has status {data.get("status")!r}, not COMPLETE')
    return data


def quality(stage):
    """The QC report of a stage (outputs/<stage>/qc.json); it must have passed."""
    path = OUTPUTS / _stage(stage)[1] / 'qc.json'
    data = _json(path)
    if data.get('status') != 'PASS':
        raise BuildError(f'{path.relative_to(ROOT).as_posix()} has status {data.get("status")!r}, not PASS')
    return data


def config(stage):
    """The configuration of a stage, e.g. config('cohort') reads src/04_cohort_curation/cohort_config.json."""
    folder, _, name = _stage(stage)
    return _json(SRC / folder / name)


def register(stage):
    """The register of a stage, which its configuration names, e.g. register('cohort') reads src/04_cohort_curation/curation_register.json:
    the entries for single studies and study platforms."""
    folder, _, name = _stage(stage)
    return _json(SRC / folder / pick(config(stage), 'register', name))


_FIGURE_DATA = None


def figure_data(figure, panel, columns=None):
    """The points of one panel of the paper's figures (outputs/visualization/figure_data.csv), with the columns that the
    panel uses. Named columns must exist and hold a value in at least one row of the panel."""
    global _FIGURE_DATA
    if _FIGURE_DATA is None:
        _FIGURE_DATA = pd.read_csv(_existing(OUTPUTS / 'visualization' / 'figure_data.csv'), low_memory=False)
    rows = _FIGURE_DATA[(_FIGURE_DATA['figure'] == figure) & (_FIGURE_DATA['panel'] == panel)]
    if rows.empty:
        raise BuildError(f'figure_data.csv has no rows for {figure} panel {panel}')
    rows = rows.dropna(axis=1, how='all')
    _columns_present(rows, columns, f'figure_data.csv {figure} panel {panel}')
    return rows[columns].reset_index(drop=True) if columns else rows.reset_index(drop=True)


def pick(data, path, source='data'):
    """A nested value by dotted path, e.g. pick(receipt, 'groups.bulk_blood.samples'), or a clear error."""
    value = data
    for key in path.split('.'):
        if not isinstance(value, dict) or key not in value:
            raise BuildError(f'{source} has no {path}')
        value = value[key]
    return value


# ---------- numbers ----------
def rnd(value, digits=0):
    """Round to the digits the page shows; None stays None (also for NaN)."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    value = round(float(value), digits)
    return int(value) if digits <= 0 else value


def sig(value, digits=3):
    """Round to significant digits."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    v = float(value)
    return 0 if v == 0 else round(v, digits - 1 - int(math.floor(math.log10(abs(v)))))


def clean(x):
    """Plain Python values for JSON: numpy and pandas scalars become numbers, NaN becomes None."""
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if hasattr(x, 'item') and not isinstance(x, (str, bytes)):
        x = x.item()
    if isinstance(x, float) and math.isnan(x):
        return None
    return x


def pack(records):
    """Records that share their fields, written once as {"__f": fields, "__r": [values]}; common.js unpacks them into
    the same list of objects when the part is loaded."""
    records = list(records)
    if not records:
        return []
    fields = list(records[0])
    if any(list(r) != fields for r in records):
        raise BuildError('pack() needs records with the same fields in the same order')
    return {'__f': fields, '__r': [[r[f] for f in fields] for r in records]}


# ---------- writing ----------
PART_NAME = re.compile(r'^[a-z]+(?:/[A-Za-z0-9._-]+)?$')


def part_file(key):
    """The page file of a part: core -> data.js, studies -> data/studies/index.js, studies/<id> -> data/studies/<id>.js,
    any other -> data/<part>.js (the same rule as partFile in common.js)."""
    if not PART_NAME.match(key):
        raise BuildError(f'invalid part name {key!r}')
    if key == 'core':
        return HERE / 'data.js'
    if key == 'studies':
        return HERE / 'data' / 'studies' / 'index.js'
    return HERE / 'data' / f'{key}.js'


def serialize(key, payload):
    text = json.dumps(clean(payload), separators=(',', ':'), ensure_ascii=False, allow_nan=False)
    text = text.replace('\u2028', '\\u2028').replace('\u2029', '\\u2029')
    return f'(window.DATA_PARTS = window.DATA_PARTS || {{}})[{json.dumps(key)}] = {text};\n'.encode('utf-8')


def write_atomic(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(dir=path.parent, prefix=f'.{path.name}.', suffix='.tmp')
    try:
        with os.fdopen(handle, 'wb') as stream:
            stream.write(content)
        os.replace(temporary, path)
    except BaseException:
        if os.path.exists(temporary):
            os.remove(temporary)
        raise


def write_part(name, files):
    """Write the files of one part after checking their names and sizes; files of the part that it no longer writes
    (a study that left the resource) are removed."""
    if not isinstance(files, dict) or name not in files:
        raise BuildError(f'{name}() must return a dict with the key {name!r}')
    content = {}
    for key, payload in files.items():
        if key != name and not key.startswith(f'{name}/'):
            raise BuildError(f'{name}() returned the part {key!r}, which belongs to another part')
        content[key] = serialize(key, payload)
        limit = CORE_LIMIT if key == 'core' else PART_LIMIT
        if len(content[key]) > limit:
            raise BuildError(f'{part_file(key).relative_to(HERE).as_posix()} would be {len(content[key]) / 1e3:.0f} kB, '
                             f'over its limit of {limit / 1e3:.0f} kB')
    for key, data in content.items():
        write_atomic(part_file(key), data)
    written = {part_file(k) for k in content}
    folder = HERE / 'data' / name
    if name != 'core' and folder.is_dir():
        for old in folder.glob('*.js'):
            if old not in written:
                old.unlink()
    sizes = sorted(len(v) for v in content.values())
    detail = f'{sizes[0] / 1e3:.0f} kB' if len(sizes) == 1 else \
        f'{len(sizes)} files, {sum(sizes) / 1e3:.0f} kB, largest {sizes[-1] / 1e3:.0f} kB'
    print(f'Wrote {name}: {part_file(name).relative_to(HERE).as_posix()} ({detail})')


PAGE_SCRIPTS = ['data.js', 'common.js', 'story.js', 'resource.js', 'studies.js', 'methods.js', 'main.js']


def stamp():
    """Give every script of index.html and every part file a version from its content (src="main.js?v=..." and
    window.PART_VERSIONS), so that a browser never runs a cached old file together with new ones."""
    page = HERE / 'index.html'
    html = _existing(page).read_text(encoding='utf-8')
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()[:10]

    def script(m):
        return f'src="{m.group(1)}?v={digest(_existing(HERE / m.group(1)))}"'
    pattern = r'src="(' + '|'.join(re.escape(s) for s in PAGE_SCRIPTS) + r')(?:\?v=[0-9a-f]+)?"'
    html, n = re.subn(pattern, script, html)
    if n != len(PAGE_SCRIPTS):
        raise BuildError(f'index.html: expected {len(PAGE_SCRIPTS)} page scripts to stamp, found {n}')
    versions = {p.relative_to(HERE).as_posix(): digest(p) for p in sorted((HERE / 'data').rglob('*.js'))}
    line = f'<script id="part-versions">window.PART_VERSIONS = {json.dumps(versions, separators=(",", ":"))};</script>'
    html, n = re.subn(r'<script id="part-versions">.*?</script>', lambda m: line, html, flags=re.S)
    if n != 1:
        raise BuildError('index.html has no <script id="part-versions"> to fill')
    write_atomic(page, html.encode('utf-8'))
    print(f'Stamped the {len(PAGE_SCRIPTS)} scripts of index.html and {len(versions)} part files with content versions')


# ---------- core: data.js ----------
GLOSSARY = [
    ('Sepsis', 'Life-threatening organ dysfunction caused by a dysregulated response of the body to an infection (the Sepsis-3 '
               'definition). Earlier definitions drew the line differently, so a study\'s "sepsis" means what that study meant by it.'),
    ('Septic shock', 'Sepsis in which the circulation and the metabolism of the cells fail, so that blood pressure stays low '
                     'despite fluids.'),
    ('Gene, RNA and expression', 'A gene is a stretch of DNA with the instructions for one product. A cell that uses a gene '
                                 'copies it into messenger RNA; the more copies it holds, the more highly the gene is expressed.'),
    ('Transcriptome', 'All the RNA of a sample, measured as one number for each of thousands of genes. A profile is the '
                      'transcriptome of one sample; a table holds one row per gene and one column per sample.'),
    ('Microarray', 'A chip with millions of short DNA probes at known spots. RNA of the sample sticks to the probes that match '
                   'it, and the light of each spot measures how much there was.'),
    ('RNA-seq', 'RNA sequencing: the RNA of a sample is broken into fragments (its library), the fragments are read, and the '
                'reads of every gene are counted.'),
    ('Logarithmic scale', 'Expression is shown as log\u2082 of the amount: each step of one is a doubling, so a difference of '
                          'two is a fourfold change.'),
    ('GEO, ArrayExpress and ENA', 'Public repositories of gene expression and sequencing data. Accessions such as GSE65682 or '
                                  'E-MTAB-4421 name a study; GSM accessions name a sample record.'),
    ('HGNC', 'The HUGO Gene Nomenclature Committee, which gives every human gene one approved symbol, such as XIST, and a fixed '
             'identifier. GTS uses these identifiers as the key of every gene.'),
    ('Study platform', 'The samples of one study measured on one platform (one array design or one sequencer): the unit in '
                       'which values were processed.'),
    ('Sample record', 'One deposited profile of a sample in one study.'),
    ('Distinct sample', 'One physical array or library, however many times it was deposited. Its first deposit by release '
                        'date is its primary record.'),
    ('Repeat deposit', 'A record of an array or library that is already deposited in another study. It is kept and linked to '
                       'the primary record, so that the sample is counted once.'),
    ('Sample class', 'Sepsis, healthy control, non-septic control or other condition, as the records of the sample describe '
                     'it. A healthy control is called healthy or normal by its study; a non-septic control is another control, '
                     'often a patient without infection.'),
    ('Whole blood, leukocytes and PBMC', 'Whole blood keeps every cell; leukocytes are the white cells; peripheral blood '
                                         'mononuclear cells (PBMC) keep only lymphocytes and monocytes. Bulk blood means whole '
                                         'blood, blood not further specified and leukocytes.'),
    ('Batch effect', 'A difference between groups of samples caused by how they were measured or processed, not by biology.'),
    ('Correction', 'Removing the differences between studies so that samples of many studies share one scale. Studies with '
                   'enough healthy controls were corrected on them; the others by matched groups, which line up the averages '
                   'of their classes with those of the first.'),
    ('Principal components', 'Axes that summarize thousands of genes at once: the first is the direction in which the samples '
                             'differ most, the second the next, at right angles to it. Samples close together have similar '
                             'profiles.'),
    ('Patient-safe group', 'Studies that share a sample or an identified patient, joined into one group. A split by these '
                           'groups keeps every deposit of a sample and every sample of a patient on one side.'),
    ('Patient-safe split', 'A split into training and test data that moves whole patient-safe groups, so that no test record '
                           'has a twin in the training data. All the folds of the benchmark are patient-safe.'),
    ('Twin', 'Another record of the same physical array, deposited in another study, or another sample of the same patient, '
             'often from another day. A twin carries much of the same information as the record itself.'),
    ('Leakage', 'Information about the test data reaching a model during training, such as a twin of a test record in the '
                'training data. It makes the model look better than it is.'),
]

# Labels of the sample types; the visualization configuration names the groups, and sorted cells are named one by one.
SORTED_CELLS = {'neutrophils': 'Neutrophils', 'monocytes': 'Monocytes', 'lymphocytes': 'Lymphocytes', 'b_cells': 'B cells',
                'cd4_t_cells': 'CD4 T cells', 'cd8_t_cells': 'CD8 T cells', 'nk_cells': 'NK cells'}


def core():
    """Main numbers, labels and colors, the glossary and the start page's view of the corrected bulk-blood table."""
    vis, cohort_config = config('visualization'), config('cohort')
    acquisition, acquisition_qc = receipt('acquisition', 3), quality('acquisition')
    expression_qc = quality('expression')
    links, flags, cohort, bench = receipt('cohort', 2), receipt('cohort', 4), receipt('cohort', 5), receipt('cohort', 6)
    cohort_qc = quality('cohort')
    sex = receipt('validation', 2)
    corrected = receipt('correction', 7)
    for stage in ('validation', 'correction', 'visualization'):
        quality(stage)

    classes = vis['class_order']
    families = vis['family_order']
    tissues = cohort_config['curation']['vocabularies']['tissue']
    tiers = list(vis['correction']['tier_labels'])
    group_of = {t: label for label, members in vis['tissue_groups'].items() for t in members}
    missing = [t for t in tissues if t not in group_of]
    if missing:
        raise BuildError(f'visualization_config.json tissue_groups does not name {", ".join(missing)}')
    tissue_labels = {t: SORTED_CELLS.get(t, group_of[t]) for t in tissues}

    # main numbers
    summary = table('cohort/cohort_summary', ['study_accession', 'release_date'])
    years = pd.to_datetime(summary['release_date']).dt.year
    distinct = pick(cohort, 'distinct_samples', '05.json')
    # The included studies hold the records with expression values, which the resource holds, and those without, which the
    # cohort stage lists in its ledger with their reason.
    included_records = pick(acquisition, 'included.sample_records', '03.json')
    held = pick(receipt('cohort', 1), 'records', 'cohort 01.json')
    records, without = held['with_values'], pick(held, 'without_values.total', 'cohort 01.json')
    if held['included'] != included_records or pick(expression_qc, 'summary.samples', 'expression qc.json') != included_records:
        raise BuildError(f'cohort 01.json and expression qc.json do not count the {included_records} sample records of the included studies')
    if records + without != included_records:
        raise BuildError('cohort 01.json: the records with and without expression values do not add up to the included records')
    for other, source in [(pick(cohort, 'records.total', '05.json'), 'cohort 05.json records.total'),
                          (pick(cohort_qc, 'summary.samples', 'cohort qc.json'), 'cohort qc.json summary.samples'),
                          (pick(expression_qc, 'summary.in_a_matrix', 'expression qc.json'), 'expression qc.json summary.in_a_matrix')]:
        if other != records:
            raise BuildError(f'{source} is {other}, but cohort 01.json counts {records} sample records with expression values')
    if sum(distinct[c] for c in classes) != distinct['total']:
        raise BuildError('the distinct samples by class in cohort 05.json do not add up to their total')
    tasks = pick(bench, 'tasks', '06.json')
    groups = pick(corrected, 'groups', 'correction 07.json')
    numbers = {
        'census_records': pick(acquisition_qc, 'summary.census_records', 'acquisition qc.json'),
        'candidates': pick(acquisition_qc, 'summary.candidates', 'acquisition qc.json'),
        'studies': pick(acquisition, 'included.studies', '03.json'),
        'studies_by_source': pick(acquisition, 'included.by_source', '03.json'),
        'studies_by_technology': {'array': pick(acquisition, 'included.by_assay.array', '03.json'),
                                  'rna_seq': pick(acquisition, 'included.by_assay.rna-seq', '03.json')},
        'first_year': int(years.min()), 'last_year': int(years.max()),
        'records': records, 'included_records': included_records, 'records_without_values': without,
        'distinct': distinct['total'],
        'classes': {c: distinct[c] for c in classes},
        'study_platforms': pick(cohort, 'cohorts.study_platforms', '05.json'),
        'tables': pick(expression_qc, 'summary.matrices', 'expression qc.json'),
        'repeated_samples': pick(links, 'samples.deposited_more_than_once', 'cohort 02.json'),
        'flagged_records': pick(flags, 'flagged.samples', 'cohort 04.json'),
        'patients': pick(cohort, 'patients.identified', '05.json'),
        'one_per_patient': pick(cohort, 'patients.one_per_patient.total', '05.json'),
        'patient_safe_groups': pick(cohort, 'cohorts.patient_safe_groups', '05.json'),
        'sex_compared': pick(sex, 'samples.compared', 'validation 02.json'),
        'sex_concordance': rnd(pick(sex, 'samples.concordance', 'validation 02.json'), 4),
        'sex_concordance_reliable': rnd(pick(sex, 'samples.concordance_outside_unreliable_matrices', 'validation 02.json'), 4),
        'sex_unreliable_tables': len(pick(sex, 'unreliable_matrices', 'validation 02.json')),
        'bench_tasks': len(pick(cohort_config, 'benchmark.tasks', 'cohort_config.json')),
        'bench_folds': sum(sum(pick(t, 'folds', '06.json').values()) for t in tasks.values()),
        'corrected': {g: {'samples': v['samples'], 'studies': v['studies'], 'genes': v['genes'],
                          'by_method': v['by_method'], 'by_class': v['by_class'], 'by_tissue': v['by_tissue']}
                      for g, v in groups.items()},
    }
    if len(tasks) != numbers['bench_tasks']:
        raise BuildError('the benchmark receipt (cohort 06.json) and cohort_config.json name different numbers of tasks')

    # the start page: every sample of the corrected bulk-blood table, before and after the correction (Fig6 A-C)
    view = figure_data('Fig6', 'views_bulk_blood', ['sample_id', 'study_accession', 'platform', 'sample_class', 'family',
                                                    'tissue', 'unit_tier', 'pc1_before', 'pc2_before', 'pc1_after', 'pc2_after'])
    samples = table('correction/corrected/bulk_blood_samples', ['sample_id'])
    if len(view) != len(samples) or set(view['sample_id']) != set(samples['sample_id']):
        raise BuildError(f'Fig6 views_bulk_blood holds {len(view)} samples, the corrected bulk-blood table {len(samples)}; '
                         'rerun stage 08')
    if len(view) != groups['bulk_blood']['samples']:
        raise BuildError('the corrected bulk-blood table and correction 07.json differ in their number of samples')
    for column, allowed in [('sample_class', classes), ('family', families), ('tissue', tissues), ('unit_tier', tiers)]:
        unknown = sorted(set(view[column]) - set(allowed))
        if unknown:
            raise BuildError(f'Fig6 views_bulk_blood: {column} values {unknown} are not in the configuration')
    if view[['pc1_before', 'pc2_before', 'pc1_after', 'pc2_after']].isna().any().any():
        raise BuildError('Fig6 views_bulk_blood has samples without a position')
    order = list(range(len(view)))
    random.Random(vis['seed']).shuffle(order)          # a fixed order, so that no class is always drawn on top
    view = view.iloc[order].reset_index(drop=True)
    studies = sorted(view['study_accession'].unique())
    platforms = sorted(view['platform'].fillna('').unique())
    accession = view['sample_id'].str.split(':', n=1).str[1]
    if (view['sample_id'].str.split(':', n=1).str[0] != view['study_accession']).any():
        raise BuildError('Fig6 views_bulk_blood: a sample_id does not start with its study')
    code = lambda column, levels: ''.join(str(levels.index(v)) for v in view[column])
    factors = figure_data('Fig6', 'factors', ['group', 'stage', 'factor', 'component_1', 'component_2', 'variance_1',
                                              'variance_2', 'genes_used'])
    factors = factors[factors['group'] == 'bulk_blood']
    if set(factors['stage']) != {'before', 'after'}:
        raise BuildError('Fig6 factors has no before and after rows for bulk blood')
    hero = {
        'figure': 'Fig6', 'n': len(view), 'scale': 10,
        'x0': [rnd(10 * v) for v in view['pc1_before']], 'y0': [rnd(10 * v) for v in view['pc2_before']],
        'x1': [rnd(10 * v) for v in view['pc1_after']], 'y1': [rnd(10 * v) for v in view['pc2_after']],
        'cls': code('sample_class', classes), 'fam': code('family', families), 'tis': code('tissue', tissues),
        'tier': code('unit_tier', tiers),
        'studies': studies, 'study': [studies.index(s) for s in view['study_accession']],
        'platforms': platforms, 'platform': [platforms.index(p) for p in view['platform'].fillna('')],
        'sample': accession.tolist(),
        'genes': int(factors['genes_used'].iloc[0]),
        'variance': {s: [rnd(r['variance_1'], 4), rnd(r['variance_2'], 4)]
                     for s, r in factors.drop_duplicates('stage').set_index('stage').iterrows()},
        'explained': {s: {r['factor']: [rnd(r['component_1'], 2), rnd(r['component_2'], 2)] for _, r in rows.iterrows()}
                      for s, rows in factors.groupby('stage')},
    }

    return {'core': {
        'numbers': numbers,
        'order': {'class': classes, 'family': families, 'tissue': tissues, 'tier': tiers},
        'labels': {'class': vis['class_labels'], 'family': vis['family_labels'], 'tissue': tissue_labels,
                   'tissue_group': group_of, 'tier': vis['correction']['tier_labels'], 'age': vis['age_labels']},
        'colors': {'class': vis['class_colors'], 'family': vis['family_colors'], 'tier': vis['correction']['tier_colors']},
        'glossary': [list(entry) for entry in GLOSSARY],
        'hero': hero,
    }}


# ---------- the parts of the other views ----------
# ---------- story: data/story.js ----------
# The animation follows the cohort group of one study, the worked example of the guide; the sex step uses a study whose
# records give the sex of every sample, since the records of that group rarely do; the last step holds that group out.
STORY_STUDY = 'GSE26440'
STORY_SEX_STUDY = 'GSE65682'
STORY_TASK, STORY_SCHEME = 'sepsis_vs_healthy_control', 'leave_one_group_out'
# Short names of the screening reasons, in three groups for the page; the full rule of each is the acquisition config's.
STORY_REASON_GROUPS = [
    ('scope', 'Not a study of blood RNA in sepsis', [
        ('not_a_dataset', 'Not a dataset'), ('not_human', 'No human samples'),
        ('not_transcriptome_wide_mrna', 'Not genome-wide mRNA'), ('single_cell', 'Single cells only'),
        ('not_blood', 'No blood sample'), ('sorted_cells', 'Sorted cells or plasma only'),
        ('ex_vivo_stimulated', 'Blood treated outside the body'), ('no_sepsis_patients', 'No patient with sepsis')]),
    ('overlap', 'Samples already in another study', [
        ('third_party_reanalysis', 'Re-analysis of another study'),
        ('superseries_of_included_subseries', 'Collection of included studies'),
        ('overlaps_included_study', 'Overlaps an included study')]),
    ('access', 'Met the criteria, but no expression values', [
        ('controlled_access', 'Controlled access'), ('no_public_reads', 'Reads not public'),
        ('raw_reads_only', 'Raw reads only'), ('no_expression_files', 'No expression file'),
        ('no_usable_expression_values', 'Files not readable as values')]),
]
STORY_SOURCES = {'geo': 'GEO', 'arrayexpress': 'ArrayExpress', 'ena': 'ENA', 'ega': 'EGA', 'dbgap': 'dbGaP', 'node': 'NODE'}
STORY_GLOSSARY = [
    ('Candidate study', 'A study that the searches found, counted once however many repositories hold it. It is not yet '
                        'judged eligible.'),
    ('Screening', 'Checking every candidate study against the eligibility rules, in a fixed order; the first rule that a '
                  'study fails is its reason for exclusion.'),
    ('Cohort group', 'Study platforms that share at least one sample, joined so that every deposit of an array or library '
                     'lies in one group.'),
    ('Primary record', 'The first deposit of an array or library, from the study released first. Keeping only primary '
                       'records keeps one copy of every array or library.'),
    ('Recorded wording', 'The words of a sample record, or of a sentence of its study, from which a curated value such as '
                         'the class was read. The resource keeps them beside the value.'),
    ('Patient identifier', 'The label of a patient in the records of a study, such as 01_0024. Samples with one identifier '
                           'are samples of one patient.'),
    ('One sample per patient', 'The sample marked for each patient: preferably without a quality flag, and the '
                               'earliest. A sample without a patient identifier is marked as its own patient.'),
    ('XIST and the Y-chromosome genes', 'XIST is active in the blood cells of women and almost silent in those of men; '
                                        'genes of the Y chromosome, such as RPS4Y1, are active only in men. Together they '
                                        'show the sex of a sample.'),
    ('Training and test data', 'A classifier is fitted on the training data and judged on the test data, which it has not '
                               'seen. A fold is one such division of the samples.'),
    ('Random split', 'A split that puts each record in the test data by chance, here a fifth of the records, whatever study '
                     'or patient it belongs to.'),
    ('Leave one study out', 'A set of splits that tests each study in turn, with all the other studies as training data.'),
]


def _story_nat(text):
    """A sort key that orders accessions by their numbers: GSM9 before GSM10."""
    return tuple((0, int(p), '') if p.isdigit() else (1, 0, p) for p in re.split(r'(\d+)', str(text)) if p)


def _story_screen(acquisition_qc):
    reasons = pick(config('acquisition'), 'screening.reasons', 'acquisition_config.json')
    decisions = table('acquisition/study_decisions', ['source', 'accession', 'title', 'assay', 'release_date', 'decision', 'reason'])
    candidates = pick(acquisition_qc, 'summary.candidates', 'acquisition qc.json')
    if len(decisions) != candidates or decisions['accession'].duplicated().any():
        raise BuildError(f'study_decisions holds {len(decisions)} rows for {candidates} candidate studies')
    named = [code for _, _, rows in STORY_REASON_GROUPS for code, _ in rows]
    unknown = sorted(set(decisions['reason']) - set(named) - {'eligible'})
    if unknown or set(named) | {'eligible'} != set(reasons):
        raise BuildError(f'the screening reasons of the page and of acquisition_config.json differ: {unknown or sorted(set(reasons) ^ set(named))}')
    included = decisions['decision'] == 'include'
    if (included != (decisions['reason'] == 'eligible')).any():
        raise BuildError('study_decisions: a study is included without the reason eligible, or excluded with it')
    exclusions = pick(acquisition_qc, 'summary.exclusions', 'acquisition qc.json')
    counts = decisions['reason'].value_counts()
    rows = [{'code': 'eligible', 'label': 'Included', 'group': 'included', 'rule': reasons['eligible'], 'n': int(counts.get('eligible', 0))}]
    for group, _, members in STORY_REASON_GROUPS:
        for code, label in members:
            n = int(counts.get(code, 0))
            if n != exclusions.get(code, 0):
                raise BuildError(f'{code}: {n} studies in study_decisions, {exclusions.get(code, 0)} in acquisition qc.json')
            if n:
                rows.append({'code': code, 'label': label, 'group': group, 'rule': reasons[code], 'n': n})
    index = {r['code']: i for i, r in enumerate(rows)}
    sources = [s for s in STORY_SOURCES if s in set(decisions['source'])]
    if set(decisions['source']) - set(sources):
        raise BuildError(f'study_decisions: unknown sources {sorted(set(decisions["source"]) - set(sources))}')
    decisions = decisions.assign(order=decisions['reason'].map(index), key=decisions['accession'].map(_story_nat))
    decisions = decisions.sort_values(['order', 'key'], key=lambda c: c if c.name == 'order' else c.map(tuple))
    title = lambda t: (t[:117].rstrip() + '…') if len(str(t)) > 120 else str(t)
    year = lambda d: int(str(d)[:4]) if re.match(r'^\d{4}', str(d or '')) else None
    return {
        'census_records': pick(acquisition_qc, 'summary.census_records', 'acquisition qc.json'), 'candidates': candidates,
        'reasons': pack(rows), 'groups': [[g, label] for g, label, _ in STORY_REASON_GROUPS],
        'sources': [STORY_SOURCES[s] for s in sources],
        'studies': pack({'a': r.accession, 't': title(r.title), 'y': year(r.release_date), 's': sources.index(r.source),
                         'k': r.assay or '', 'r': int(r.order)} for r in decisions.itertuples()),
    }


def _story_group():
    groups = table('cohort/cohort_groups', ['cohort_group', 'study_accession', 'platform', 'release_date', 'records',
                                            'primary_records', 'repeated_records'])
    found = groups.loc[groups['study_accession'] == STORY_STUDY, 'cohort_group'].unique()
    if len(found) != 1:
        raise BuildError(f'{STORY_STUDY} must lie in one cohort group, not {len(found)}')
    gid = found[0]
    members = groups[groups['cohort_group'] == gid].sort_values(['release_date', 'study_accession', 'platform']).reset_index(drop=True)
    cols = ['study_accession', 'sample_accession', 'platform', 'cohort_group', 'patient_safe_group', 'sample_id',
            'is_primary_record', 'records_of_sample', 'linked_by', 'sample_class', 'class_differs_between_records', 'sex',
            'patient_id', 'patient_key', 'patient_identified', 'one_per_patient', 'time_point', 'time_order', 'time_unit',
            'any_flag', 'quality_flags']
    records = table('cohort/cohort_samples', cols)
    records = records[records['cohort_group'] == gid].copy()
    summary = table('cohort/cohort_group_summary')
    summary = summary[summary['cohort_group'] == gid]
    if len(summary) != 1:
        raise BuildError(f'cohort_group_summary has no single row for {gid}')
    summary = summary.iloc[0]
    if len(records) != members['records'].sum():
        raise BuildError(f'{gid}: {len(records)} records in cohort_samples, {members["records"].sum()} in cohort_groups')
    classes = config('visualization')['class_order']
    unit = lambda s, p: f'{s}|{p or ""}'
    units = [unit(s, p) for s, p in zip(members['study_accession'], members['platform'])]
    records['unit'] = [unit(s, p) for s, p in zip(records['study_accession'], records['platform'])]
    primary = records[records['is_primary_record']].copy()
    if len(primary) != records['sample_id'].nunique() or len(primary) != summary['distinct_samples']:
        raise BuildError(f'{gid}: the primary records, distinct sample identifiers and cohort_group_summary disagree')
    # distinct samples in the order of their primary study platform, then of their accession
    primary['col'] = primary['unit'].map(units.index)
    primary = primary.assign(key=primary['sample_accession'].map(_story_nat).map(tuple)).sort_values(['col', 'key'])
    sample_index = {sid: i for i, sid in enumerate(primary['sample_id'])}
    records['sample'] = records['sample_id'].map(sample_index)
    records['col'] = records['unit'].map(units.index)
    if records['sample'].isna().any():
        raise BuildError(f'{gid}: a record belongs to a sample without a primary record in the group')
    # in each column: its primary records first, in the order of the samples, then its repeat deposits by their sample
    records = records.assign(rep=~records['is_primary_record']).sort_values(['col', 'rep', 'sample'])
    # the recorded wording of each primary record's class, and the sentence of its study where one was used
    labels = table('cohort/label_sources', ['study_accession', 'sample_accession', 'platform', 'variable', 'fields', 'recorded',
                                            'basis', 'study_statement'])
    labels = labels[(labels['variable'] == 'sample_class') & labels['study_accession'].isin(members['study_accession'])]
    labels = labels.assign(unit=[unit(s, p) for s, p in zip(labels['study_accession'], labels['platform'])])
    labels = labels.drop_duplicates(['unit', 'sample_accession']).set_index(['unit', 'sample_accession'])
    statements = []

    def wording(row):
        key = (row.unit, row.sample_accession)
        if key not in labels.index:
            raise BuildError(f'label_sources has no class source for {row.study_accession} {row.sample_accession}')
        lab = labels.loc[key]
        sentence = str(lab['study_statement'] or '')
        if sentence and sentence not in statements:
            statements.append(sentence)
        return [str(lab['fields'] or ''), str(lab['recorded'] or ''), statements.index(sentence) if sentence else -1]
    # patients: identified ones by their key, in the order of their first sample
    identified = primary[primary['patient_identified']]
    keys = list(dict.fromkeys(identified['patient_key']))
    patient_index = {k: i for i, k in enumerate(keys)}
    shown_id = identified.groupby('patient_key')['patient_id'].first()
    samples = []
    for row in primary.itertuples():
        f, w, s = wording(row)
        if row.sample_class not in classes:
            raise BuildError(f'{row.sample_id}: class {row.sample_class!r} is not in the configuration')
        samples.append({'col': int(row.col), 'acc': row.sample_accession, 'cls': classes.index(row.sample_class), 'field': f,
                        'said': w, 'stmt': s, 'differs': bool(row.class_differs_between_records), 'nrec': int(row.records_of_sample),
                        'pat': patient_index[row.patient_key] if row.patient_identified else -1, 'marked': bool(row.one_per_patient),
                        'time': str(row.time_point or ''), 'order': rnd(row.time_order, 2), 'unit': str(row.time_unit or ''),
                        'flag': str(row.quality_flags or ''), 'sex': str(row.sex or 'not_recorded')})
    # worked examples: a sepsis sample of the study itself, and the sample of the group with the most records whose records
    # differ in class (the array of the guide)
    own = primary[(primary['study_accession'] == STORY_STUDY) & (primary['sample_class'] == 'sepsis')]
    differ = primary[primary['class_differs_between_records']].sort_values(['records_of_sample', 'col', 'key'], ascending=[False, True, True])
    if own.empty or differ.empty:
        raise BuildError(f'{gid}: no sepsis sample of {STORY_STUDY}, or no sample whose records differ in class')
    twin = differ.iloc[0]
    twin_records = records[records['sample_id'] == twin['sample_id']]
    twin_rows = []
    for row in twin_records.itertuples():
        lab = labels.loc[(row.unit, row.sample_accession)]
        twin_rows.append({'col': int(row.col), 'acc': row.sample_accession, 'cls': classes.index(row.sample_class),
                          'field': str(lab['fields'] or ''), 'said': str(lab['recorded'] or ''), 'primary': bool(row.is_primary_record)})
    # a patient with two samples of different class in one study, else the first patient with several samples
    several = identified.groupby('patient_key').filter(lambda p: len(p) >= 2)
    changing = several.groupby('patient_key').filter(lambda p: p['sample_class'].nunique() > 1 and p['unit'].nunique() == 1)
    example_patient = (changing if len(changing) else several)['patient_key']
    if example_patient.empty:
        raise BuildError(f'{gid}: no patient with several samples')
    # how the repeat deposits of the study itself were first deposited
    focus = records[records['study_accession'] == STORY_STUDY]
    first_in = focus[~focus['is_primary_record']].merge(primary[['sample_id', 'col']].rename(columns={'col': 'first'}), on='sample_id')
    per_sample = records.groupby('sample_id').size()
    return gid, primary, {
        'id': gid, 'study': STORY_STUDY, 'patient_safe_group': str(primary['patient_safe_group'].iloc[0]),
        'columns': pack({'study': r.study_accession, 'platform': r.platform or '', 'date': r.release_date, 'records': int(r.records),
                         'primary': int(r.primary_records)} for r in members.itertuples()),
        'rec_col': records['col'].astype(int).tolist(), 'rec_acc': records['sample_accession'].tolist(),
        'rec_sample': records['sample'].astype(int).tolist(),
        'samples': pack(samples), 'statements': statements,
        'patients': [str(shown_id[k]) for k in keys],
        'stated': {'records': int(members['records'].sum()), 'distinct': int(summary['distinct_samples']),
                   'classes': {c: int(summary[c]) for c in classes}, 'patients': int(summary['patients_identified']),
                   'without_patient': int(summary['samples_without_patient']), 'one_per_patient': int(summary['one_per_patient'])},
        'linked_by': str(records.loc[records['linked_by'] != '', 'linked_by'].mode().iloc[0]) if (records['linked_by'] != '').any() else '',
        'most_records': int(per_sample.max()), 'samples_with_most': int((per_sample == per_sample.max()).sum()),
        'focus_first_in': {int(c): int(n) for c, n in first_in['first'].value_counts().sort_index().items()},
        'example_sepsis': sample_index[own.iloc[0]['sample_id']], 'example_twin': sample_index[twin['sample_id']],
        'twin_records': pack(twin_rows), 'example_patient': patient_index[example_patient.iloc[0]],
    }


def _story_sex(primary):
    plane = figure_data('Fig5', 'sex_plane', ['study_accession', 'platform', 'sample_accession', 'recorded', 'expression_sex',
                                              'agreement', 'x', 'y'])
    plane = plane[plane['study_accession'] == STORY_SEX_STUDY]
    matrices = table('validation/sex_check_matrices')
    matrix = matrices[matrices['study_accession'] == STORY_SEX_STUDY]
    if len(matrix) != 1 or plane.empty:
        raise BuildError(f'{STORY_SEX_STUDY} must have one tested matrix in sex_check_matrices and points in Fig5 sex_plane')
    m = matrix.iloc[0]
    called = plane['expression_sex'].value_counts()
    agree = plane['agreement'].value_counts()
    if len(plane) != m['samples'] or called.get('female', 0) != m['called_female'] or called.get('male', 0) != m['called_male'] \
            or agree.get('concordant', 0) != m['concordant'] or agree.get('discordant', 0) != m['discordant']:
        raise BuildError(f'Fig5 sex_plane of {STORY_SEX_STUDY} does not hold the calls of its matrix; rerun stage 08')
    plane = plane.assign(key=plane['sample_accession'].map(_story_nat).map(tuple)).sort_values('key')
    values = table('validation/sex_check_samples', ['study_accession', 'sample_accession', 'recorded_sex', 'female_gene_value',
                                                    'male_genes_value', 'difference', 'expression_sex'])
    values = values[values['study_accession'] == STORY_SEX_STUDY]
    values = values.assign(key=values['sample_accession'].map(_story_nat).map(tuple)).sort_values('key')
    examples = []
    for sex in ('male', 'female'):
        hit = values[(values['recorded_sex'] == sex) & (values['expression_sex'] == sex)]
        if hit.empty:
            raise BuildError(f'{STORY_SEX_STUDY}: no sample recorded and read as {sex}')
        r = hit.iloc[0]
        examples.append({'acc': r['sample_accession'], 'sex': sex, 'y': rnd(r['male_genes_value'], 2),
                         'xist': rnd(r['female_gene_value'], 2), 'd': rnd(r['difference'], 1)})
    rec, call = ['female', 'male', 'not_recorded'], ['female', 'male', 'ambiguous']
    bad = sorted(set(plane['recorded']) - set(rec)) + sorted(set(plane['expression_sex']) - set(call))
    if bad:
        raise BuildError(f'Fig5 sex_plane: unexpected values {bad}')
    sex = config('validation')['sex']
    return {
        'study': STORY_SEX_STUDY, 'platform': str(m['platform'] or ''), 'acc': plane['sample_accession'].tolist(),
        'x': [rnd(v, 3) for v in plane['x']], 'y': [rnd(v, 3) for v in plane['y']],
        'rec': ''.join(str(rec.index(v)) for v in plane['recorded']), 'call': ''.join(str(call.index(v)) for v in plane['expression_sex']),
        'matrix': {k: clean(m[k]) for k in ('samples', 'called_female', 'called_male', 'ambiguous', 'with_recorded_sex',
                                            'concordant', 'discordant')} | {'separation': rnd(m['separation'], 1)},
        'examples': examples, 'y_genes': list(sex['male_genes']), 'female_gene': sex['female_gene'],
        'group_sex': {k: int((primary['sex'] == k).sum()) for k in rec},
    }


def _story_correction(gid):
    corrected = table('correction/corrected/bulk_blood_samples', ['sample_id', 'study_accession', 'cohort_group', 'unit_method'])
    mine = corrected[corrected['cohort_group'] == gid]
    methods = config('correction')
    by_study = mine.groupby(['study_accession', 'unit_method']).size().reset_index(name='n')
    return {'samples': sorted(mine['sample_id'], key=_story_nat), 'minimum_healthy': pick(methods, 'inputs.minimum_healthy_controls', 'correction_config.json'),
            'by_study': pack({'study': r.study_accession, 'method': r.unit_method, 'n': int(r.n)} for r in by_study.itertuples())}


def _story_bench(gid, primary):
    tasks_config = {t['task']: t for t in pick(config('cohort'), 'benchmark.tasks', 'cohort_config.json')}
    if STORY_TASK not in tasks_config:
        raise BuildError(f'cohort_config.json has no benchmark task {STORY_TASK}')
    group = str(primary['patient_safe_group'].iloc[0])
    splits = table('cohort/benchmark_splits')
    of_task = splits[(splits['task'] == STORY_TASK) & (splits['scheme'] == STORY_SCHEME)]
    fold = of_task[of_task['fold'] == group]
    if fold.empty or list(fold.loc[fold['role'] != 'train', 'patient_safe_group']) != [group]:
        raise BuildError(f'{STORY_TASK} {STORY_SCHEME} has no fold that holds out {group} alone')
    tasks = table('cohort/benchmark_tasks', ['task', 'patient_safe_group', 'eligible', 'reason', 'label'])
    tasks = tasks[tasks['task'] == STORY_TASK]
    eligible = tasks[tasks['eligible']]
    totals = {'pos': int((eligible['label'] == 1).sum()), 'neg': int((eligible['label'] == 0).sum())}
    if fold['positives'].sum() != totals['pos'] or fold['negatives'].sum() != totals['neg']:
        raise BuildError(f'the fold that holds out {group} does not hold every eligible sample of {STORY_TASK}')
    samples = table('cohort/cohort_samples', ['patient_safe_group', 'study_accession'])
    studies = samples.groupby('patient_safe_group')['study_accession'].agg(lambda s: sorted(set(s), key=_story_nat))
    fold = fold.assign(size=fold['positives'] + fold['negatives']).sort_values(['size', 'patient_safe_group'], ascending=[False, True])
    held = tasks[tasks['patient_safe_group'] == group]
    reasons = held.loc[~held['eligible'], 'reason'].value_counts()
    t = tasks_config[STORY_TASK]
    return {
        'task': STORY_TASK, 'scheme': STORY_SCHEME, 'fold': group, 'positive': t['positive'], 'negative': t['negative'],
        'groups': pack({'g': r.patient_safe_group, 'role': r.role, 'tech': r.technology, 'pos': int(r.positives), 'neg': int(r.negatives),
                        'studies': studies.get(r.patient_safe_group, [])} for r in fold.itertuples()),
        'totals': totals, 'folds': int(of_task['fold'].nunique()),
        'held': {'records': int(len(held)), 'pos': int(((held['eligible']) & (held['label'] == 1)).sum()),
                 'neg': int(((held['eligible']) & (held['label'] == 0)).sum()), 'reasons': {k: int(v) for k, v in reasons.items()}},
        'minimum_per_class': pick(config('cohort'), 'benchmark.minimum_per_class', 'cohort_config.json'),
    }


# The whole resource, from the searches to the classes of the samples with expression values (the flow of Figure 1).
STORY_NOT_CANDIDATE = {'no_human_rnaseq_runs': 'ENA records without a human RNA-seq run',
                       'cross_source_duplicate': 'Copies of a study found in another repository',
                       'study_type_not_expression': 'ArrayExpress records that are not expression studies'}
STORY_NO_VALUES = {'absent_from_table': 'Not named in the table of its study', 'no_own_data_file': 'No data file of its own',
                   'not_processed': 'Study platform not processed', 'not_in_ncbi_counts': 'Missing from the NCBI counts',
                   'incomplete_array_file': 'Raw array file incomplete'}


def _story_whole(classes):
    """One small integer per sample record of the included studies, in bands by technology and study platform: its class (bits
    0-2), whether it is a primary record (3), has expression values (4), is marked as one per patient (5), and the reason it has
    none (6-8). The records with values are those of the sample table; those without come from the ledger of the cohort stage,
    each its own record, with no class (0)."""
    flow = figure_data('Fig1', 'flow', ['step', 'count'])
    F = {s: int(c) for s, c in zip(flow['step'], flow['count'])}
    needed = ['records_found', 'candidate_studies', 'included_studies', 'included_study_platforms', 'study_platforms_with_values',
              'included_sample_records', 'records_without_values', 'sample_records_with_values', 'repeat_deposits', 'distinct_samples',
              'one_sample_per_patient', 'identified_patients', 'samples_without_identifier'] + [f'class_{c}' for c in classes]
    missing = [s for s in needed if s not in F]
    if missing:
        raise BuildError(f'figure_data.csv Fig1 flow has no step {", ".join(missing)}')
    group = lambda prefix: {s[len(prefix):]: n for s, n in F.items() if s.startswith(prefix)}
    not_candidate, excluded, no_values = group('not_candidate_'), group('excluded_'), group('without_values_')
    unknown = sorted(set(not_candidate) - set(STORY_NOT_CANDIDATE)) + sorted(set(no_values) - set(STORY_NO_VALUES))
    if unknown:
        raise BuildError(f'Fig1 flow: the page has no label for {", ".join(unknown)}')
    keys = ['study_accession', 'sample_accession', 'platform']
    records = table('cohort/cohort_samples', keys + ['sample_id', 'is_primary_record', 'sample_class', 'one_per_patient'])
    ledger = table('cohort/records_without_values', keys + ['reason'])
    index = table('annotation/sample_index', keys + ['in_matrix', 'reason'])
    for frame in (records, ledger, index):
        frame['platform'] = frame['platform'].fillna('')
    # every included record is in the sample table, with values, or in the ledger, without them and with the reason of the expression stage
    listed = pd.concat([records[keys].assign(values=True, why=''), ledger[keys].assign(values=False, why=ledger['reason'])], ignore_index=True)
    joined = index.merge(listed, on=keys, how='outer', validate='one_to_one', indicator=True)
    if ((joined['_merge'] != 'both').any() or (joined['in_matrix'] != joined['values']).any()
            or (joined['reason'].fillna('') != joined['why']).any()):
        raise BuildError('cohort_samples and the ledger records_without_values do not hold the records of annotation/sample_index with their reasons')
    reasons = list(STORY_NO_VALUES)
    if not set(ledger['reason']) <= set(reasons):
        raise BuildError(f'records_without_values: the page has no label for {sorted(set(ledger["reason"]) - set(reasons))}')
    gone = ledger.assign(sample_id='', is_primary_record=True, sample_class=classes[0], one_per_patient=False, values=False,
                         why=ledger['reason'].map(lambda r: reasons.index(r) + 1))
    records = pd.concat([records.assign(values=True, why=0), gone.drop(columns='reason')], ignore_index=True)
    units = table('annotation/matrix_index', ['accession', 'platform', 'family']).rename(columns={'accession': 'study_accession'})
    released = table('acquisition/included_studies', ['accession', 'release_date']).rename(columns={'accession': 'study_accession'})
    units['platform'] = units['platform'].fillna('')
    units = units.merge(released, on='study_accession', how='left', validate='many_to_one')
    order = config('visualization')['family_order']
    if units['family'].isna().any() or not set(units['family']) <= set(order) or units['release_date'].isna().any():
        raise BuildError('a study platform of annotation/matrix_index has no technology family or no included study')
    units = units.assign(f=units['family'].map(order.index), key=units['study_accession'].map(_story_nat))
    units = units.sort_values(['f', 'release_date', 'key', 'platform']).reset_index(drop=True)
    units['u'] = range(len(units))
    records = records.merge(units[['study_accession', 'platform', 'u']], on=['study_accession', 'platform'], how='left', validate='many_to_one')
    if records['u'].isna().any():
        raise BuildError('a sample record lies on a study platform that annotation/matrix_index does not hold')
    records = records.assign(key=records['sample_accession'].map(_story_nat)).sort_values(['u', 'key']).reset_index(drop=True)
    units['records'] = records.groupby('u').size().reindex(units['u'], fill_value=0).values
    units['with_values'] = records[records['values']].groupby('u').size().reindex(units['u'], fill_value=0).values
    cls = records['sample_class'].map(classes.index)
    if cls.isna().any():
        raise BuildError('a record has a class outside the configuration')
    code = (cls + 8 * records['is_primary_record'].astype(int) + 16 * records['values'].astype(int)
            + 32 * records['one_per_patient'].astype(int) + 64 * records['why'])
    first = {sid: i for i, sid in enumerate(records['sample_id']) if sid and records['is_primary_record'].iat[i]}
    repeats = [(i, first[sid]) for i, sid in enumerate(records['sample_id']) if not records['is_primary_record'].iat[i]]
    # the page draws the flow from these records; they must give the counts of Figure 1
    held = records[records['values']]
    shown = held[held['is_primary_record']]
    checks = {'included_sample_records': len(records), 'records_without_values': len(gone), 'sample_records_with_values': len(held),
              'repeat_deposits': len(repeats), 'distinct_samples': len(shown), 'one_sample_per_patient': int(held['one_per_patient'].sum()),
              'included_study_platforms': len(units), 'study_platforms_with_values': int((units['with_values'] > 0).sum())}
    checks.update({f'class_{c}': int((shown['sample_class'] == c).sum()) for c in classes})
    checks.update({f'without_values_{r}': int((gone['why'] == reasons.index(r) + 1).sum()) for r in reasons if r in no_values})
    wrong = {k: (v, F.get(k)) for k, v in checks.items() if v != F.get(k)}
    if wrong or sum(no_values.values()) != len(gone) or len(first) != len(shown):
        raise BuildError(f'the sample records do not give the counts of Fig1 flow: {wrong}; rerun stage 08')
    excluded_groups = []
    for g, label, members in STORY_REASON_GROUPS:
        excluded_groups.append([g, label, sum(excluded.get(code_, 0) for code_, _ in members)])
    if sum(n for _, _, n in excluded_groups) != sum(excluded.values()):
        raise BuildError('Fig1 flow names an exclusion reason that the page does not group')
    return {
        'flow': {k: F[k] for k in needed},
        'not_candidate': [[k, STORY_NOT_CANDIDATE[k], not_candidate[k]] for k in STORY_NOT_CANDIDATE if k in not_candidate],
        'excluded': excluded_groups,
        'no_values': [[r, STORY_NO_VALUES[r], no_values.get(r, 0)] for r in reasons],
        'units': pack({'study': r.study_accession, 'platform': r.platform, 'fam': int(r.f), 'n': int(r.records)} for r in units.itertuples()),
        'code': [int(v) for v in code], 'rep': [v for pair in repeats for v in pair],
    }


# Twins and splits: the records of one task as a naive analysis pools them, under one real random split, leave one study
# out and the patient-safe folds of the resource (the leakage audit of the benchmark stage, script 01).
LEAK_TASK, LEAK_REPEAT, LEAK_SCHEME = 'sepsis_vs_non_septic_control', 1, 'leave_one_group_out'


def _leak_marks(frame, test, train):
    """The test records of a split that have a twin in training: another record of the same array, another array of the
    same identified patient, and either (the rule of the leakage audit)."""
    held = frame[train]
    arrays = held.groupby('sample_id').size()
    named = held[held['patient_key'] != '']
    patients, pairs = named.groupby('patient_key').size(), named.groupby(['patient_key', 'sample_id']).size()
    key = frame['patient_key']
    same_array = test & frame['sample_id'].map(arrays).fillna(0).gt(0)
    other = key.map(patients).fillna(0) - pd.Series([pairs.get((p, s), 0) for p, s in zip(key, frame['sample_id'])], index=frame.index)
    same_patient = test & (key != '') & other.gt(0)
    return {'exposed': int((same_array | same_patient).sum()), 'array': int(same_array.sum()),
            'patient': int(same_patient.sum()), 'test': int(test.sum())}, same_array | same_patient


def _story_leak():
    """Every record of the naive population of LEAK_TASK in blocks by patient-safe group and study, with its twins, the test
    records of one random split, and the folds of leave one group out; every count is recomputed here and must equal the
    leakage audit."""
    bench, audit = config('benchmark'), receipt('benchmark', 1)
    policy, tasks = bench['splits'], {t['task']: t for t in pick(config('cohort'), 'benchmark.tasks', 'cohort_config.json')}
    if LEAK_TASK not in tasks:
        raise BuildError(f'cohort_config.json has no benchmark task {LEAK_TASK}')
    keys = ['source', 'study_accession', 'sample_accession', 'platform']
    people = table('benchmark/split_populations', ['population', 'task', *keys, 'sample_id', 'patient_key', 'patient_safe_group',
                                                   'linked_records', 'same_array_links', 'same_patient_links'])
    people['platform'] = people['platform'].fillna('')
    frame = people[(people['population'] == 'naive') & (people['task'] == LEAK_TASK)].reset_index(drop=True)
    safe = people[(people['population'] == 'gts') & (people['task'] == LEAK_TASK)]
    eligible = pd.MultiIndex.from_frame(frame[keys]).isin(pd.MultiIndex.from_frame(safe[keys]))
    if eligible.sum() != len(safe) or frame.duplicated(keys).any():
        raise BuildError(f'split_populations: the samples of {LEAK_TASK} are not all records of its naive population')
    leak = table('benchmark/split_leakage', ['population', 'task', 'scheme', 'fold', 'test_records', 'exposed', 'exposed_same_array',
                                             'exposed_same_patient', 'share_exposed', 'expected_share'])
    rows = leak[leak['task'] == LEAK_TASK]
    review_rows = review('benchmark', 'split_leakage')
    totals = {s: r for s, r in review_rows[review_rows['task'] == LEAK_TASK].set_index('scheme').iterrows()}

    def same(found, row, where):
        want = {'exposed': row['exposed'], 'array': row['exposed_same_array'], 'patient': row['exposed_same_patient'], 'test': row['test_records']}
        if any(found[k] != int(want[k]) for k in want):
            raise BuildError(f'{where}: the page counts {found}, the leakage audit {want}')
    # the twins of every record, by the same rule as the audit
    per_array = frame.groupby('sample_id')['sample_id'].transform('size')
    named = frame['patient_key'] != ''
    arrays_of_patient = frame[named].groupby('patient_key')['sample_id'].transform('nunique').reindex(frame.index, fill_value=1)
    twin_array, twin_patient = per_array.gt(1), named & arrays_of_patient.gt(1)
    if (twin_array != frame['same_array_links'].gt(0)).any() or (twin_patient != frame['same_patient_links'].gt(0)).any() \
            or ((twin_array | twin_patient) != frame['linked_records'].gt(0)).any():
        raise BuildError(f'split_populations: the twins of {LEAK_TASK} recomputed here differ from the audit')
    for column in ('sample_id', 'patient_key'):
        spread = frame[frame[column] != ''].groupby(column)['patient_safe_group'].nunique()
        if (spread > 1).any():
            raise BuildError(f'split_populations: records of one {column} lie in two patient-safe groups ({spread[spread > 1].index[0]})')
    # one real random split: the first test_fraction of the records by SHA-256 of seed, repeat and record
    digest = [hashlib.sha256(f'{policy["seed"]}:{LEAK_REPEAT}:{s}:{a}'.encode('utf-8')).hexdigest()
              for s, a in zip(frame['study_accession'], frame['sample_accession'])]
    test = pd.Series(False, index=frame.index)
    test.iloc[sorted(range(len(frame)), key=digest.__getitem__)[:math.floor(policy['test_fraction'] * len(frame) + 0.5)]] = True
    fold = f'repeat_{LEAK_REPEAT:02d}'
    random_row = rows[(rows['population'] == 'naive') & (rows['scheme'] == 'naive_record_split') & (rows['fold'] == fold)]
    if len(random_row) != 1:
        raise BuildError(f'split_leakage has no single row for {LEAK_TASK} {fold}')
    random_marks, _ = _leak_marks(frame, test, ~test)
    same(random_marks, random_row.iloc[0], f'{LEAK_TASK} random split {fold}')
    # leave one study out, every study once
    loso = {'exposed': 0, 'array': 0, 'patient': 0, 'test': 0}
    for study in sorted(set(frame['study_accession'])):
        held_out = frame['study_accession'] == study
        found, _ = _leak_marks(frame, held_out, ~held_out)
        loso = {k: loso[k] + found[k] for k in loso}
    same(loso, totals['naive_leave_one_study_out'], f'{LEAK_TASK} leave one study out')
    # the patient-safe folds: the samples of the task in the held-out groups against every record of the training groups
    splits = table('cohort/benchmark_splits', ['task', 'scheme', 'fold', 'patient_safe_group', 'role'])
    splits = splits[(splits['task'] == LEAK_TASK) & (splits['scheme'] == LEAK_SCHEME)]
    # blocks: patient-safe groups by their records, studies within a group by their records, then the records of a patient
    # and of an array side by side
    size = frame.groupby('patient_safe_group').size()
    group_order = sorted(size.index, key=lambda g: (-size[g], _story_nat(g)))
    unknown = sorted(set(splits['patient_safe_group']) - set(group_order))
    if unknown:
        raise BuildError(f'benchmark_splits names groups without records of {LEAK_TASK}: {unknown[:3]}')
    study_size = frame.groupby('study_accession').size()
    frame = frame.assign(eligible=eligible, test=test.values, g=frame['patient_safe_group'].map(group_order.index),
                         s=frame['study_accession'].map(lambda s: -study_size[s]),
                         p=[k or f'~{a}' for k, a in zip(frame['patient_key'], frame['sample_accession'])],
                         key=frame['sample_accession'].map(_story_nat))
    frame = frame.sort_values(['g', 's', 'study_accession', 'p', 'key']).reset_index(drop=True)
    # a band is the records of one study in one group; a study measured on two platforms can lie in two groups
    studies = list(dict.fromkeys(frame['study_accession']))
    bands = list(dict.fromkeys(zip(frame['g'], frame['study_accession'])))
    folds, gts = [], {'exposed': 0, 'array': 0, 'patient': 0, 'test': 0}
    gts_rows = rows[(rows['population'] == 'gts') & (rows['scheme'] == LEAK_SCHEME)].set_index('fold')
    for name, f in splits.groupby('fold', sort=False):
        out = sorted(f.loc[f['role'] != 'train', 'patient_safe_group'], key=group_order.index)
        train = sorted(f.loc[f['role'] == 'train', 'patient_safe_group'], key=group_order.index)
        found, _ = _leak_marks(frame, frame['eligible'] & frame['patient_safe_group'].isin(out), frame['patient_safe_group'].isin(train))
        if name not in gts_rows.index or found['exposed'] or found['test'] != gts_rows.loc[name, 'test_records']:
            raise BuildError(f'{LEAK_TASK} fold {name}: the page counts {found}, the leakage audit {dict(gts_rows.loc[name]) if name in gts_rows.index else "nothing"}')
        gts = {k: gts[k] + found[k] for k in gts}
        folds.append({'fold': name, 'out': [group_order.index(g) for g in out], 'train': [group_order.index(g) for g in train]})
    same(gts, totals[LEAK_SCHEME], f'{LEAK_TASK} {LEAK_SCHEME}')
    folds.sort(key=lambda f: f['out'], reverse=True)      # the largest group last, where the animation ends
    # sets of twins: arrays with two or more records here, identified patients with two or more records here
    def sets(column):
        counts = frame[frame[column] != ''].groupby(column).size()
        index = {k: i for i, k in enumerate(dict.fromkeys(k for k in frame[column] if k and counts.get(k, 0) > 1))}
        return [index.get(k, -1) for k in frame[column]], index
    array_set, arrays = sets('sample_id')
    patient_set, patients = sets('patient_key')
    # the examples: the array with the most records, and the identified patient with the most arrays
    deposits = frame[frame['sample_id'].isin(arrays)].groupby('sample_id')['study_accession'].nunique()
    top_array = max(arrays, key=lambda k: (deposits[k], -arrays[k]))
    spread = frame[frame['patient_key'].isin(patients)].groupby('patient_key')['sample_id'].nunique()
    top_patient = max(patients, key=lambda k: (spread[k], -patients[k]))
    cohort = table('cohort/cohort_samples', ['study_accession', 'sample_accession', 'patient_id', 'time_point', 'time_unit'])
    mine = frame[frame['patient_key'] == top_patient].merge(cohort, on=['study_accession', 'sample_accession'], how='left', validate='one_to_one')
    t, scores = tasks[LEAK_TASK], {}
    summary = review('benchmark', 'model_summary', ['task', 'model', 'scheme', 'auroc_mean'])
    for model, part_ in summary[summary['task'] == LEAK_TASK].groupby('model'):
        scores[model] = {s: rnd(v, 3) for s, v in zip(part_['scheme'], part_['auroc_mean']) if s in ('naive_record_split', LEAK_SCHEME)}
    row = lambda r: {'test': int(r['test_records']), 'exposed': int(r['exposed']), 'array': int(r['exposed_same_array']),
                     'patient': int(r['exposed_same_patient']), 'folds': int(r['folds'])}
    return {
        'task': LEAK_TASK, 'positive': t['positive'], 'negative': t['negative'], 'repeat': LEAK_REPEAT, 'repeats': policy['repeats'],
        'test_fraction': policy['test_fraction'], 'scheme': LEAK_SCHEME,
        'groups': group_order, 'studies': studies, 'bands': pack({'s': studies.index(s), 'g': int(g)} for g, s in bands),
        'rec_band': [bands.index(k) for k in zip(frame['g'], frame['study_accession'])], 'rec_acc': frame['sample_accession'].tolist(),
        'rec_array': array_set, 'rec_patient': patient_set,
        'rec_flag': ''.join(str(int(e) + 2 * int(x)) for e, x in zip(frame['eligible'], frame['test'])),
        'folds': pack(folds),
        'example_array': arrays[top_array], 'example_patient': patients[top_patient],
        'example_patient_id': str(mine['patient_id'].dropna().iloc[0]) if mine['patient_id'].notna().any() else '',
        'example_patient_times': [f'{u} {p}'.strip() if p else '' for p, u in zip(mine['time_point'].fillna(''), mine['time_unit'].fillna(''))],
        'audit': {
            'records': len(frame), 'linked': int((twin_array | twin_patient).sum()), 'linked_array': int(twin_array.sum()),
            'linked_patient': int(twin_patient.sum()), 'samples': int(len(safe)),
            'random': random_marks | {'share': rnd(random_row.iloc[0]['share_exposed'], 6), 'expected': rnd(random_row.iloc[0]['expected_share'], 6)},
            'random_all': row(totals['naive_record_split']) | {'expected': rnd(totals['naive_record_split']['expected_share'], 6)},
            'loso': row(totals['naive_leave_one_study_out']), 'gts': row(totals[LEAK_SCHEME]),
            'patient_safe_folds': pick(audit, 'patient_safe_folds.folds', 'benchmark 01.json'),
            'all_records': pick(audit, 'records', 'benchmark 01.json'),
        },
        'scores': scores, 'models': {m: RESOURCE_MODELS[m].lower() for m in RESOURCE_MODELS if m in scores},
    }


def story():
    """How it works: three animations of real samples. The whole resource, from the searches to the final counts of the
    samples with expression values; one cohort group close up, from its records to a patient-safe benchmark split; and the
    twins of the records of one task across a random split, leave one study out and the patient-safe folds."""
    acquisition_qc = quality('acquisition')
    quality('cohort'), quality('validation'), quality('correction'), quality('benchmark'), quality('visualization')
    gid, primary, group = _story_group()
    return {'story': {
        'whole': _story_whole(config('visualization')['class_order']),
        'screen': _story_screen(acquisition_qc), 'group': group, 'sex': _story_sex(primary),
        'correction': _story_correction(gid), 'bench': _story_bench(gid, primary), 'leak': _story_leak(),
        'glossary': [list(entry) for entry in STORY_GLOSSARY],
    }}


# ---------- the resource: data/resource.js ----------
# Plain-language names of the keys that the resource's chapters show; a key of the data without a name here stops the build.
RESOURCE_REASONS = {
    # Fig1 flow: records of the searches that did not become candidates, and candidates excluded at screening
    'no_human_rnaseq_runs': 'ENA record without human RNA-seq runs',
    'cross_source_duplicate': 'The same study, found in a second repository',
    'study_type_not_expression': 'ArrayExpress record that is not an expression study',
    'not_blood': 'No blood sample',
    'no_sepsis_patients': 'No patient with sepsis',
    'sorted_cells': 'Sorted cells only',
    'single_cell': 'Single-cell profiles only',
    'not_transcriptome_wide_mrna': 'Not a transcriptome-wide mRNA profile',
    'raw_reads_only': 'Deposited only as raw reads',
    'ex_vivo_stimulated': 'Samples stimulated outside the body',
    'overlaps_included_study': 'Holds the samples of an included study',
    'controlled_access': 'Data under controlled access',
    'no_public_reads': 'No public reads',
    'not_a_dataset': 'Not a data set',
    'third_party_reanalysis': 'Re-analysis of other studies by a third party',
    'not_human': 'Not human',
    'superseries_of_included_subseries': 'Umbrella record of an included study',
    'no_expression_files': 'No expression files deposited',
    'no_usable_expression_values': 'Deposited files not readable as comparable values',
}
# The curated fields of SuppFig2, the bases of a class (Fig3 A), the evidence of a link (Fig5 B), the mapping outcomes and scales
# (Fig4), the methods and outcomes of the correction (Fig6), and the tasks, schemes and exclusions of the benchmark
RESOURCE_COVERAGE = {'sample_class': 'Class', 'tissue': 'Sample type', 'sex': 'Sex', 'age_group': 'Age group',
                     'age_years': 'Age in years', 'patient': 'Patient identifier', 'time_point': 'Time point',
                     'outcome': 'Outcome', 'septic_shock': 'Septic shock', 'score': 'Severity score'}
RESOURCE_BASIS = {'sample_record': 'The sample\'s own record',
                  'sample_record_and_study_description': 'Its record and a quoted sentence of the study record',
                  'study_description': 'The study record alone', 'not_recorded': 'Nothing recorded'}
RESOURCE_LINKS = {'chip_type_scan_time_and_median_intensity': 'Scan time of an Affymetrix array',
                  'beadchip_barcode_and_position': 'Barcode and position of an Illumina BeadChip',
                  'slide_identifier_and_median_intensities': 'Slide of an Agilent array',
                  'identical_counts': 'Identical RNA-seq counts',
                  'expression_identity': 'Near-identical profile, studies known to share samples',
                  'record_names_a_re_analysis': 'A record that names a re-analysis',
                  'expression_identity; identical_counts': 'Several kinds of evidence together'}
RESOURCE_MAPPING = {'one_gene': 'One gene', 'several_genes': 'Several genes (not assigned)',
                    'no_gene_in_annotation': 'No gene in the annotation', 'gene_not_in_hgnc': 'A gene not approved by HGNC',
                    'accession_not_resolved': 'An accession that could not be resolved'}
RESOURCE_SCALES = {'log2': 'log₂', 'log2_like': 'log₂-like (authors\' values)',
                   'log2 counts per million': 'log₂ counts per million', 'linear': 'Linear (authors\' values)',
                   'other': 'Other (authors\' values)'}
RESOURCE_METHODS = {'none': 'No correction', 'plain_combat': 'Plain ComBat', 'flat_genes': 'Flat genes',
                    'matched_groups': 'Matched groups'}
RESOURCE_OUTCOMES = {'in_corrected_table': 'In a corrected table', 'in_table_bulk_blood': 'In the bulk-blood table',
                     'in_table_pbmc': 'In the PBMC table',
                     'repeat_deposit': 'Repeat deposit, corrected only to compare the deposits',
                     'sample_type_not_corrected': 'Swabs, plasma, sorted cells or type not recorded',
                     'newborn_study_not_corrected': 'Newborn study, not corrected',
                     'no_class_to_line_up': 'No class to line up', 'value_scale_not_comparable': 'Table not on a log scale',
                     'study_has_too_few_samples': 'Too few samples to estimate the study\'s effects',
                     'no_check_applicable': 'No check applied'}
RESOURCE_TASKS = {'sepsis_vs_healthy_control': 'Sepsis against healthy controls',
                  'sepsis_vs_non_septic_control': 'Sepsis against non-septic controls',
                  'mortality_in_sepsis': 'Death against survival in sepsis'}
RESOURCE_SCHEMES = {'leave_one_group_out': 'Leave one group out', 'technology_transfer': 'Technology transfer',
                    'temporal': 'Temporal'}
RESOURCE_BENCH_REASONS = {
    'not_primary_record': 'Not the primary record of its sample',
    'not_one_per_patient': 'Not the one sample marked for its patient', 'quality_flag': 'Carries a quality flag',
    'tissue_outside_task': 'Not whole blood, blood not further specified or leukocytes',
    'phase_outside_task': 'Drawn before the onset of sepsis or in recovery',
    'class_outside_task': 'Of a class outside the task', 'label_not_recorded': 'No recorded outcome'}
# Which file of the data package answers which question (the archives and members named here must exist in the manifest)
RESOURCE_QUESTIONS = [
    ('Which samples are there, and what do their records say?', 'cohort.tar', 'cohort_samples.csv.gz',
     'One row for every sample record with expression values: class, sample type, sex, age, time point, outcome, septic shock, '
     'the sample it belongs to, its patient, its patient-safe group and its quality flags. Start here to choose samples.'),
    ('Which records have no expression values, and why?', 'provenance.tar.gz', 'provenance_records_without_values.csv.gz',
     'The sample records of the included studies that have no expression values, each with its reason and its recorded '
     'fields. They take no part in the sample tables.'),
    ('What are the expression values of one study?', 'gene_expression.tar', None,
     'One gene-level table for every study platform, genes as rows and samples as columns, on the scale its study platform '
     'names. These per-study tables are the reference.'),
    ('What are the read counts of an RNA-seq study?', 'gene_counts.tar', None,
     'Gene-level counts of every RNA-seq study platform, for methods that start from counts.'),
    ('Which records are the same array or library, and which samples share a patient?', 'cohort.tar',
     'cohort_sample_links.csv.gz',
     'The primary record of every sample and how its deposits were linked; with sample_id and patient_key of the sample '
     'table, keep one copy of every sample or one sample of every patient.'),
    ('Where does a label come from?', 'cohort.tar', 'cohort_label_sources.csv.gz',
     'For every curated value, the record field or study sentence it was read from; the recorded metadata are kept verbatim '
     'in cohort_sample_metadata.'),
    ('Which gene does a probe measure?', 'gene_mapping.tar', 'gene_mapping_genes.csv.gz',
     'The HGNC gene set, and for every study platform the gene of every row and the rows behind every gene.'),
    ('Can I compare samples of many studies on one scale?', 'corrected_expression.tar',
     'corrected_expression_bulk_blood_samples.csv.gz',
     'The bulk-blood and PBMC tables with the differences between studies removed (corrected_expression_bulk_blood_expression '
     'and corrected_expression_pbmc_expression), their sample tables with how every study was corrected, and the outcome of '
     'every record in the correction.'),
    ('How do I test a classifier without leakage?', 'cohort.tar', 'cohort_benchmark_splits.csv.gz',
     'The samples of the three tasks with the reason for every exclusion (cohort_benchmark_tasks), and the folds with the '
     'role of every patient-safe group.'),
    ('What do naive splits change, and how do published signatures fare?', 'benchmark.tar', None,
     'The twins of every record and their share in every split, the scores of the baseline classifiers in every test unit '
     'with the measured effects, the study-only predictor and the study identification, and every signature stratum with '
     'its AUROC or the reason it was skipped.'),
    ('What does a column mean?', 'data_dictionary.csv', None,
     'The definition, type, unit and missing-value code of every column of every table.'),
    ('How many records reached each step?', 'cohort_flow.csv', None,
     'Counts from the searches to the samples, for the resource and for every study platform.'),
]
RESOURCE_GLOSSARY = [
    ('Candidate study', 'A study that the searches found, counted once even when it was found in several repositories. '
                        'Every candidate was screened against the eligibility rules and either included or excluded for '
                        'the first rule it failed.'),
    ('Primary record', 'The first deposit of a sample, by release date. Its accession identifies the distinct sample, so '
                       'keeping the primary records keeps one copy of every array or library.'),
    ('Cohort group', 'Study platforms that share at least one sample.'),
    ('One sample per patient', 'One sample marked for every patient. A sample without a patient identifier counts as its '
                               'own patient, so the number of these samples is an upper bound on the number of patients.'),
    ('Probe and probe set', 'A short sequence on an array that binds one transcript; a probe set is a group of probes read '
                            'as one value. The rows of an array table are probes or probe sets.'),
    ('Gene-level table', 'A table of one study platform with one row for every HGNC gene and one column for every sample. '
                         'On arrays it keeps, for each gene, the row with the highest mean; for RNA-seq it sums the counts.'),
    ('XIST and the Y-chromosome genes', 'XIST is expressed from the second X chromosome, so it is high in women; genes '
                                        'of the Y chromosome, such as RPS4Y1, are expressed only in men. Together they '
                                        'show the sex of a sample.'),
    ('Spearman correlation', 'A number from −1 to 1 that says how far two lists of values put their genes in the '
                             'same order; 1 is the same order.'),
    ('Sepsis MetaScore', 'A published score of eleven genes, six higher and five lower in sepsis: the mean of the higher '
                         'genes minus the mean of the lower ones.'),
    ('ROC curve and AUROC', 'A ROC curve shows how well a score separates two groups at every threshold. The area under '
                            'it (AUROC) is the chance that a sample of one group scores above a sample of the other: 1 is '
                            'perfect separation, 0.5 none.'),
    ('Standardized difference', 'The difference of two group means divided by their pooled standard deviation, so that '
                                'genes measured on different scales can be compared.'),
    ('Matched groups', 'The correction of a study without enough healthy controls: the mean and spread of each of its '
                       'classes are set to those of the same class in the studies corrected on healthy controls.'),
    ('ComBat', 'A widely used method that removes batch effects gene by gene. Plain ComBat, with the study as the batch, '
               'makes every study alike on average, whatever its mix of patients and controls.'),
    ('Held-out test', 'A test of the correction where the answer is known: the healthy controls of one study are hidden, '
                      'the study is corrected from its other samples, and the hidden controls should land among the '
                      'healthy controls of the other studies.'),
    ('Sampling floor', 'The error that chance alone gives in the held-out test, because a few hidden controls are only a '
                       'sample: no method can be expected to come much closer.'),
    ('Study mixing', 'A check of every corrected study: among the nearest neighbors of its samples, the share that come '
                     'from the same study, less the share expected by chance. Lower is better mixed.'),
    ('Classifier', 'A rule, usually learned by a computer from examples, that takes the profile of a sample and returns a '
                   'label or a score, such as the chance that the patient has sepsis.'),
    ('Fold', 'One division of the samples into training and test. A benchmark scheme is a set of folds.'),
    ('Benchmark', 'Fixed tasks, samples and splits into training and test, so that different methods are judged by the '
                  'same test.'),
    ('Baseline classifier', 'A simple classifier with fixed settings, here logistic regression and gradient boosting on the '
                            'genes that best separate the classes of the training data, against which other methods can be '
                            'compared.'),
    ('Study shortcut', 'Scoring the test records of many studies together, so that a classifier gains from telling the '
                       'studies apart, since studies differ in their share of each class.'),
    ('Seen study', 'A test study whose other samples were in the training data, so that a classifier has already learned '
                   'how that study measures its samples.'),
    ('Independent stratum', 'A study platform and sample type that took no part in deriving or validating a published '
                            'signature: not one of its studies, sharing no sample with them, and holding no repeat deposit.'),
    ('Balanced accuracy', 'The share of correct answers in each class, averaged over the classes, so that large classes do '
                          'not dominate. Guessing gives one divided by the number of classes.'),
    ('CC0 1.0', 'A public-domain dedication: anyone may copy, change and use the data for any purpose, without asking.'),
    ('Zenodo', 'A general-purpose open data repository, run by CERN, which gives every deposited data set a DOI.'),
    ('DOI', 'Digital object identifier: a permanent address of a published object, such as a data set.'),
]
# The benchmark of evaluation designs (stage 07): models, schemes, effects, signatures and the reasons a stratum was skipped
RESOURCE_MODELS = {'logistic_regression': 'Logistic regression', 'gradient_boosting': 'Gradient boosting',
                   'sepsis_metascore': 'Sepsis MetaScore'}
RESOURCE_EVAL_SCHEMES = {'naive_record_split': 'Random split', 'naive_leave_one_study_out': 'Leave one study out',
                         'leave_one_group_out': 'Leave one group out', 'technology_transfer': 'Technology transfer',
                         'temporal': 'Temporal'}
RESOURCE_EFFECTS = {
    'total_inflation': ('Total difference', 'The random split, scored over its pooled test records, minus leave one group '
                                            'out, scored within each held-out group'),
    'study_shortcut': ('Study shortcut', 'The random split scored over its pooled test records, minus the same predictions '
                                         'scored within each study'),
    'seen_study': ('Seen study', 'Within the studies of the random split, without the test records that have a twin in '
                                 'training, minus leave one group out'),
    'twins': ('Twins in the same study', 'Within the studies of the random split, with minus without the test records that '
                                         'have a twin in training'),
    'twins_across_studies': ('Twins in other studies', 'Leave one study out, with minus without the test records that have '
                                                       'a twin in another study'),
}
RESOURCE_SIGNATURES = {'sepsis_metascore': 'Sepsis MetaScore', 'hla_dra': 'HLA-DRA', 'faim3_plac8_ratio': 'FCMR:PLAC8',
                       'septicyte_lab': 'SeptiCyte LAB', 'sweeney2018_mortality_12gene': 'Stanford 12-gene score'}
RESOURCE_SIGNATURE_REASONS = {
    'definition_not_verified': 'Its definition could not be checked against the publication',
    'scale_not_comparable': 'The values of the platform are not on a log scale',
    'fewer_than_minimum_of_a_class': 'Too few samples of a class',
    'gene_coverage_below_threshold': 'Too few of its genes measured on the platform',
    'too_few_genes_on_a_side': 'Too few genes measured on one side of the score',
    'fewer_than_minimum_with_a_score': 'Too few samples with a score'}
RESOURCE_CODE_URL = 'https://github.com/MrRajat1809/GenomicTensorSepsis'
RESOURCE_RELEASE = ROOT / 'data_release'


def resource_named(keys, names, where):
    """Fail when a key of the data has no plain-language name in the page."""
    missing = sorted({str(k) for k in keys} - set(names))
    if missing:
        raise BuildError(f'{where}: no name for {", ".join(missing)}; add it to build_data.py')


def resource_bool(value):
    """A flag of figure_data.csv, which a mixed column holds as the text True or False."""
    return str(value) in ('True', '1', '1.0')


def resource_studies(acquisition):
    flow = figure_data('Fig1', 'flow', ['step', 'count'])
    steps = dict(zip(flow['step'], flow['count'].astype(int)))
    prefix = lambda p: [(s[len(p):], n) for s, n in steps.items() if s.startswith(p)]
    not_candidate, excluded = prefix('not_candidate_'), prefix('excluded_')
    resource_named([k for k, _ in not_candidate + excluded], RESOURCE_REASONS, 'Fig1 flow')
    unavailable = list(pick(acquisition, 'eligible_but_unavailable', 'acquisition 03.json'))
    for k in unavailable:
        if k not in dict(excluded):
            raise BuildError(f'acquisition 03.json eligible_but_unavailable names {k}, which Fig1 flow does not')
    if steps['records_found'] - sum(n for _, n in not_candidate) != steps['candidate_studies']:
        raise BuildError('Fig1 flow: the records that did not become candidates do not add up')
    if steps['candidate_studies'] - sum(n for _, n in excluded) != steps['included_studies']:
        raise BuildError('Fig1 flow: the excluded candidates do not add up')
    atlas = figure_data('Fig2', 'atlas', ['label', 'year', 'family', 'records'])
    by_year = figure_data('Fig2', 'samples_by_year', ['year', 'family', 'samples', 'added_up'])
    included = table('acquisition/included_studies', ['accession', 'assay', 'release_date', 'source'])
    included['year'] = pd.to_datetime(included['release_date']).dt.year
    searched = str(pick(receipt('acquisition', 1), 'searched_at_utc', 'acquisition 01.json'))[:10]
    return {
        'searched': searched,
        'records_found': steps['records_found'], 'candidates': steps['candidate_studies'],
        'included': steps['included_studies'], 'study_platforms': steps['included_study_platforms'],
        'study_platforms_with_values': steps['study_platforms_with_values'],
        'not_candidate': [[k, n] for k, n in sorted(not_candidate, key=lambda r: -r[1])],
        'excluded': [[k, n] for k, n in sorted(excluded, key=lambda r: -r[1])],
        'unavailable': unavailable,
        'reasons': {k: RESOURCE_REASONS[k] for k, _ in not_candidate + excluded},
        'platforms': pack({'label': r['label'], 'year': int(r['year']), 'family': r['family'], 'records': int(r['records'])}
                          for _, r in atlas.iterrows()),
        'by_year': pack({'year': int(r['year']), 'family': r['family'], 'samples': int(r['samples']),
                         'added_up': int(r['added_up'])} for _, r in by_year.iterrows()),
        'studies_by_year': pack({'year': int(y), 'array': int((g['assay'] == 'array').sum()),
                                 'rna_seq': int((g['assay'] == 'rna-seq').sum())}
                                for y, g in included.groupby('year')),
    }


def resource_samples(cohort):
    flows = figure_data('Fig2', 'flows', ['step', 'source', 'target', 'value'])
    tm = flows[flows['step'] == 'technology_to_material']
    mc = flows[flows['step'] == 'material_to_class']
    if tm.empty or mc.empty or tm['value'].sum() != mc['value'].sum():
        raise BuildError('Fig2 flows: the bands from technology to material and from material to class differ in total')
    basis = figure_data('Fig3', 'label_basis')
    resource_named([c for c in basis.columns if c not in ('sample_class', 'figure', 'panel')], RESOURCE_BASIS, 'Fig3 label_basis')
    ages = figure_data('Fig3', 'ages', ['sample_class', 'count', 'mean', 'min', '25%', '50%', '75%', 'max'])
    clinical = figure_data('Fig3', 'clinical_flow', ['samples', 'outcome', 'shock'])
    coverage = figure_data('SuppFig2', 'coverage', ['label', 'family', 'records', *RESOURCE_COVERAGE])
    samples = table('cohort/cohort_samples', ['study_accession', 'is_primary_record', 'sample_class', 'age_group', 'sex',
                                              'tissue', 'sepsis_definition'])
    primary = samples[samples['is_primary_record']]
    age = primary.groupby(['age_group', 'sample_class']).size()
    sex = primary['sex'].value_counts()
    rules = review('cohort', 'label_rules', ['study_accession', 'platform', 'variable', 'rule', 'fields', 'matched', 'value',
                                             'samples', 'basis', 'study_statement', 'examples'])
    example = rules[(rules['study_accession'] == 'GSE65682') & (rules['variable'] == 'sample_class')]
    if example.empty:
        raise BuildError('label_rules_review.csv has no class rules for GSE65682, the example of the samples chapter')
    order = ['Whole blood', 'Blood, unspecified', 'Leukocytes', 'PBMC', 'Sorted cells', 'Plasma', 'Nasopharyngeal swab',
             'Not recorded']
    resource_named(set(tm['target']) | set(mc['source']), order, 'Fig2 flows materials')
    return {
        'records': pick(cohort, 'records', 'cohort 05.json'),
        'materials': order,
        'tech_material': [[r['source'], r['target'], int(r['value'])] for _, r in tm.iterrows()],
        'material_class': [[r['source'], r['target'], int(r['value'])] for _, r in mc.iterrows()],
        'basis': {r['sample_class']: {k: int(r[k]) for k in RESOURCE_BASIS if k in basis.columns} for _, r in basis.iterrows()},
        'basis_labels': RESOURCE_BASIS,
        'ages': pack({'cls': r['sample_class'], 'n': int(r['count']), 'mean': rnd(r['mean'], 1), 'min': rnd(r['min'], 2),
                      'q1': rnd(r['25%'], 1), 'median': rnd(r['50%'], 1), 'q3': rnd(r['75%'], 1), 'max': rnd(r['max'], 1)}
                     for _, r in ages.iterrows()),
        'age_groups': {g: {c: int(n) for c, n in age[g].items()} for g in age.index.get_level_values(0).unique()},
        'sex': {k: int(v) for k, v in sex.items()},
        'clinical': [[r['shock'], r['outcome'], int(r['samples'])] for _, r in clinical.iterrows()],
        'coverage_fields': RESOURCE_COVERAGE,
        'coverage': pack({'label': re.sub(r'\s+', ' ', r['label']).strip(), 'family': r['family'], 'records': int(r['records']),
                          **{k: rnd(r[k], 3) for k in RESOURCE_COVERAGE}} for _, r in coverage.iterrows()),
        'sepsis_3_studies': int(samples.loc[samples['sepsis_definition'] == 'sepsis-3', 'study_accession'].nunique()),
        'example': {'study': 'GSE65682', 'platform': example['platform'].iloc[0],
                    'records': int((samples['study_accession'] == 'GSE65682').sum()),
                    'rules': pack({'field': r['fields'], 'wording': r['matched'], 'value': r['value'], 'samples': int(r['samples']),
                                   'basis': r['basis'], 'statement': r['study_statement'] if isinstance(r['study_statement'], str) else ''}
                                  for _, r in example.iterrows())},
    }


_RESOURCE_FAMILIES = None


def resource_family(study, platform):
    """The technology of a study platform, from the sectors of Fig2 C (labeled by study, and by platform where a study
    has several)."""
    global _RESOURCE_FAMILIES
    if _RESOURCE_FAMILIES is None:
        atlas = figure_data('Fig2', 'atlas', ['label', 'family'])
        _RESOURCE_FAMILIES = dict(zip(atlas['label'].str.replace(r'\s+', ' ', regex=True).str.strip(), atlas['family']))
    platform = platform if isinstance(platform, str) else ''
    for key in (f'{study} {platform}'.strip(), study):
        if key in _RESOURCE_FAMILIES:
            return _RESOURCE_FAMILIES[key]
    raise BuildError(f'Fig2 atlas has no sector for {study} {platform}')


def resource_repeats(cohort, links):
    flow = figure_data('Fig1', 'flow', ['step', 'count'])
    steps = dict(zip(flow['step'], flow['count'].astype(int)))
    shared = figure_data('Fig2', 'shared_samples', ['samples', 'first', 'second', 'cohort_group'])
    summary = table('cohort/cohort_summary', ['cohort_group', 'patient_safe_group', 'study_accession', 'platform', 'records',
                                              'primary_records', 'repeated_records', 'release_date'])
    groups = summary.groupby('cohort_group')['study_accession'].transform('size')
    multi = summary[groups > 1].sort_values(['cohort_group', 'release_date', 'study_accession'])
    samples = table('cohort/cohort_samples', ['cohort_group', 'sample_id', 'study_accession', 'platform'])
    g024 = samples[samples['cohort_group'] == 'G024']
    if g024.empty or 'GSE26440' not in set(g024['study_accession']):
        raise BuildError('cohort_samples: G024 is no longer the cohort group of GSE26440, the example of the repeats chapter')
    return {
        'shared': pack({'a': r['first'].replace('|', ' '), 'b': r['second'].replace('|', ' '), 'samples': int(r['samples']),
                        'group': r['cohort_group']} for _, r in shared.iterrows()),
        'steps': {k: steps[k] for k in ['included_sample_records', 'records_without_values', 'sample_records_with_values',
                                        'repeat_deposits', 'distinct_samples', 'one_sample_per_patient', 'identified_patients',
                                        'samples_without_identifier']},
        'repeated_samples': pick(links, 'samples.deposited_more_than_once', 'cohort 02.json'),
        'patients': pick(cohort, 'patients', 'cohort 05.json'),
        'class_differs': pick(cohort, 'distinct_samples.whose_records_differ_in_class', 'cohort 05.json'),
        'cohort_groups': pick(cohort, 'cohorts.cohort_groups', 'cohort 05.json'),
        'platforms': pack({'study': r['study_accession'], 'platform': r['platform'], 'group': r['cohort_group'],
                           'family': resource_family(r['study_accession'], r['platform']),
                           'year': int(str(r['release_date'])[:4]), 'records': int(r['records']),
                           'primary': int(r['primary_records']), 'repeated': int(r['repeated_records'])}
                          for _, r in multi.iterrows()),
        'example': {'group': 'G024', 'platforms': int(g024.groupby(['study_accession', 'platform']).ngroups),
                    'records': len(g024), 'distinct': int(g024['sample_id'].nunique())},
    }


def resource_genes(annotation):
    mapping = figure_data('Fig4', 'mapping')
    keys = list(RESOURCE_MAPPING)
    _columns_present(mapping, [f'{k}_rows' for k in keys], 'figure_data.csv Fig4 mapping')
    extra = [c[:-5] for c in mapping.columns if c.endswith('_rows') and c[:-5] not in keys and mapping[c].sum() > 0]
    if extra:
        raise BuildError(f'Fig4 mapping: no name for {", ".join(extra)}')
    tables = figure_data('Fig4', 'genes_per_table', ['family', 'study_accession', 'platform', 'genes', 'protein_coding_genes'])
    shared = figure_data('Fig4', 'shared_genes', ['genes', 'tables'])
    scales = figure_data('Fig4', 'scales', ['samples', 'study_accession', 'platform', 'value_scale', 'q05', 'q25', 'q50', 'q75', 'q95'])
    resource_named(scales['value_scale'], RESOURCE_SCALES, 'Fig4 scales')
    plat = lambda p: p if isinstance(p, str) else ''
    family = dict(zip(tables['study_accession'] + ' ' + tables['platform'].fillna(''), tables['family']))
    rows = pick(annotation, 'rows.total', 'annotation 04.json')
    if int(mapping[[f'{k}_rows' for k in keys]].to_numpy().sum()) != rows:
        raise BuildError('Fig4 mapping and annotation 04.json differ in the number of rows')
    return {
        'rows': rows, 'mapping_labels': RESOURCE_MAPPING,
        'mapping_colors': {k: v for k, v in pick(config('visualization'), 'mapping_colors', 'visualization_config.json').items()
                           if k in RESOURCE_MAPPING},
        'mapping': {r['family']: {k: int(r[f'{k}_rows']) for k in keys} for _, r in mapping.iterrows()},
        'tables': pack({'study': r['study_accession'], 'platform': plat(r['platform']), 'family': r['family'],
                        'genes': int(r['genes']), 'coding': int(r['protein_coding_genes'])} for _, r in tables.iterrows()),
        'shared': [[int(r['tables']), int(r['genes'])] for _, r in shared.sort_values('tables').iterrows()],
        'scale_labels': RESOURCE_SCALES,
        'scales': pack({'study': r['study_accession'], 'platform': plat(r['platform']),
                        'family': family.get(f'{r["study_accession"]} {plat(r["platform"])}', ''),
                        'scale': r['value_scale'], 'samples': int(r['samples']),
                        'q': [sig(r[c], 4) for c in ['q05', 'q25', 'q50', 'q75', 'q95']]} for _, r in scales.iterrows()),
        'genes': pick(annotation, 'genes', 'annotation 04.json'),
    }


def resource_checks(sex_receipt, repeat_receipt, score_receipt):
    plane = figure_data('Fig5', 'sex_plane', ['study_accession', 'platform', 'sample_accession', 'recorded', 'expression_sex',
                                              'agreement', 'x', 'y'])
    plane['table'] = plane['study_accession'] + ' ' + plane['platform'].fillna('')
    tables = sorted(plane['table'].unique())
    rec, sexes, agree = ['female', 'male', 'not_recorded'], ['female', 'male', 'ambiguous'], \
        ['concordant', 'discordant', 'sex_not_recorded', 'ambiguous']
    for column, allowed in [('recorded', rec), ('expression_sex', sexes), ('agreement', agree)]:
        if set(plane[column]) - set(allowed):
            raise BuildError(f'Fig5 sex_plane: unexpected {column} {sorted(set(plane[column]) - set(allowed))}')
    by_table = figure_data('Fig5', 'sex_by_table', ['label', 'concordant', 'discordant', 'compared', 'concordance', 'marked'])
    repeats = figure_data('Fig5', 'repeats', ['sample_id', 'study_a', 'study_b', 'linked_by', 'spearman',
                                              'centred_correlation', 'most_similar'])
    resource_named(repeats['linked_by'], RESOURCE_LINKS, 'Fig5 repeats')
    links = [k for k in RESOURCE_LINKS if k in set(repeats['linked_by'])]
    score = figure_data('Fig5', 'score', ['sepsis', 'healthy_control', 'tissue', 'study_accession', 'platform', 'auroc',
                                          'auroc_low', 'auroc_high', 'auroc_one_per_patient', 'role'])
    stored = table('validation/sepsis_score_by_study', ['study_accession', 'platform', 'tissue', 'auroc'])
    stored_auroc = {f'{r["study_accession"]}|{r["platform"]}|{r["tissue"]}': float(r['auroc']) for _, r in stored.iterrows()}
    roc = figure_data('Fig5', 'roc', ['stratum', 'fpr', 'tpr'])
    curves = {s: [[rnd(v, 4) for v in g['fpr']], [rnd(v, 4) for v in g['tpr']]] for s, g in roc.groupby('stratum', sort=False)}
    directions = figure_data('Fig5', 'directions', ['tissue', 'study_accession', 'platform', 'role', 'published_symbol',
                                                    'standardised_difference', 'as_published'])
    directions['stratum'] = directions['study_accession'] + '|' + directions['platform'].fillna('') + '|' + directions['tissue']
    genes = list(dict.fromkeys(directions['published_symbol']))
    role = dict(zip(directions['published_symbol'], directions['role']))
    plat = lambda p: p if isinstance(p, str) else ''
    strata = []
    for _, r in score.iterrows():
        key = f'{r["study_accession"]}|{plat(r["platform"])}|{r["tissue"]}'
        if key not in curves or key not in stored_auroc:
            raise BuildError(f'Fig5: the stratum {key} has no ROC curve or no stored AUROC')
        rows = directions[directions['stratum'] == key].set_index('published_symbol')
        d, flag = rows['standardised_difference'], rows['as_published']
        strata.append({'key': key, 'study': r['study_accession'], 'platform': plat(r['platform']),
                       'tissue': r['tissue'], 'sepsis': int(r['sepsis']), 'healthy': int(r['healthy_control']),
                       'auroc': rnd(r['auroc'], 4), 'low': rnd(r['auroc_low'], 3), 'high': rnd(r['auroc_high'], 3),
                       'one_per_patient': rnd(r['auroc_one_per_patient'], 4),
                       'derivation': r['role'] if isinstance(r['role'], str) else '',
                       'table_auroc': rnd(stored_auroc[key], 4),
                       'fpr': curves[key][0], 'tpr': curves[key][1],
                       'diff': [rnd(d.get(g), 2) for g in genes],
                       'pub': ''.join('-' if pd.isna(d.get(g)) else ('1' if resource_bool(flag.get(g)) else '0') for g in genes)})
    rename = {'C9orf95': 'NMRK1', 'KIAA1370': 'ATOSA'}        # the approved symbols, as in the caption of Fig5 D
    return {
        'sex': {'x': [rnd(v, 3) for v in plane['x']], 'y': [rnd(v, 3) for v in plane['y']],
                'rec': ''.join(str(rec.index(v)) for v in plane['recorded']),
                'call': ''.join(str(sexes.index(v)) for v in plane['expression_sex']),
                'agr': ''.join(str(agree.index(v)) for v in plane['agreement']),
                'tab': [tables.index(t) for t in plane['table']], 'tables': tables,
                'sample': plane['sample_accession'].tolist(), 'levels': {'rec': rec, 'call': sexes, 'agr': agree}},
        'sex_tables': pack({'label': re.sub(r'\s+', ' ', r['label']).strip(), 'concordant': int(r['concordant']),
                            'compared': int(r['compared']), 'marked': resource_bool(r['marked'])} for _, r in by_table.iterrows()),
        'sex_receipt': {k: pick(sex_receipt, f'samples.{k}', 'validation 02.json')
                        for k in ['total', 'called', 'ambiguous', 'compared', 'concordant', 'concordance',
                                  'compared_outside_unreliable_matrices', 'concordance_outside_unreliable_matrices']},
        'sex_matrices': {'total': pick(sex_receipt, 'matrices.total', 'validation 02.json'),
                         'with_recorded_sex': pick(sex_receipt, 'matrices.tested_with_recorded_sex', 'validation 02.json'),
                         'unreliable': [t.replace('|', ' ') for t in pick(sex_receipt, 'unreliable_matrices', 'validation 02.json')]},
        'repeats': {'link': [links.index(v) for v in repeats['linked_by']], 'links': links,
                    'link_labels': {k: RESOURCE_LINKS[k] for k in links},
                    'sp': [rnd(v, 4) for v in repeats['spearman']], 'cc': [rnd(v, 3) for v in repeats['centred_correlation']],
                    'best': ''.join('1' if resource_bool(v) else '0' for v in repeats['most_similar']),
                    'sample': repeats['sample_id'].tolist(), 'a': repeats['study_a'].tolist(), 'b': repeats['study_b'].tolist()},
        'repeat_receipt': {'pairs': pick(repeat_receipt, 'pairs.total', 'validation 03.json'),
                           'spearman_lowest': rnd(pick(repeat_receipt, 'spearman.lowest', 'validation 03.json'), 3),
                           'spearman_median': rnd(pick(repeat_receipt, 'spearman.median', 'validation 03.json'), 3),
                           'centred_median': rnd(pick(repeat_receipt, 'centred_correlation.median', 'validation 03.json'), 3),
                           'most_similar': pick(repeat_receipt, 'most_similar.pairs', 'validation 03.json')},
        'genes': [rename.get(g, g) for g in genes], 'published': [role[g] for g in genes],
        'strata': pack(strata),
        'score_receipt': {'strata': pick(score_receipt, 'strata.total', 'validation 04.json'),
                          'comparisons': pick(score_receipt, 'genes.comparisons', 'validation 04.json'),
                          'as_published': pick(score_receipt, 'genes.as_published', 'validation 04.json'),
                          'median': rnd(pick(score_receipt, 'score.auroc_median', 'validation 04.json'), 3),
                          'lowest': rnd(pick(score_receipt, 'score.auroc_lowest', 'validation 04.json'), 2),
                          'independent_strata': pick(score_receipt, 'score.independent.strata', 'validation 04.json'),
                          'independent_median': rnd(pick(score_receipt, 'score.independent.auroc_median', 'validation 04.json'), 3)},
        'minimum_group': pick(config('validation'), 'biology.minimum_group', 'validation_config.json'),
    }


def resource_correction(selected, corrected):
    factors = figure_data('Fig6', 'factors', ['group', 'stage', 'factor', 'component_1', 'component_2', 'variance_1', 'variance_2'])
    held = figure_data('Fig6', 'held_out', ['group', 'study_accession', 'platform', 'method', 'hidden_healthy',
                                            'location_error', 'scale_error', 'pool_dispersion', 'sampling_floor'])
    resource_named(held['method'], RESOURCE_METHODS, 'Fig6 held_out')
    test = table('correction/held_out_test', ['group', 'method', 'location_error', 'scale_error', 'pool_dispersion',
                                              'sampling_floor'])
    if len(test) != len(held):
        raise BuildError('Fig6 held_out and correction/held_out_test differ in their rows; rerun stage 08')
    measures = ['location_error', 'scale_error', 'pool_dispersion', 'sampling_floor']
    medians = {g: {m: {c: rnd(rows[c].median(), 4) for c in measures} | {'units': len(rows)} for m, rows in grp.groupby('method')}
               for g, grp in test.groupby('group')}
    comparison = table('correction/held_out_comparison', ['group', 'method', 'baseline', 'units', 'method_closer',
                                                          'median_difference', 'difference_low', 'difference_high', 'wilcoxon_p'])
    pair = comparison[(comparison['group'] == 'bulk_blood') & (comparison['method'] == 'plain_combat') &
                      (comparison['baseline'] == 'matched_groups')]
    if len(pair) != 1:
        raise BuildError('correction/held_out_comparison has no bulk-blood row of plain ComBat against matched groups')
    pair = pair.iloc[0]
    mixing = figure_data('Fig6', 'mixing', ['samples', 'group', 'sample_class', 'study_accession', 'platform', 'role', 'tier',
                                            'gate', 'before', 'after', 'limit', 'passed'])
    outcomes = figure_data('Fig6', 'outcomes', ['records', 'outcome'])
    flows = figure_data('Fig6', 'ledger_flows', ['step', 'source', 'target', 'value'])
    resource_named(set(outcomes['outcome']) | set(flows.loc[flows['step'] == 'group_to_outcome', 'target'].dropna()),
                   RESOURCE_OUTCOMES, 'Fig6 outcomes')
    pooled = table('correction/pooled_checks', ['group', 'scope', 'check', 'samples', 'before', 'after'])
    plat = lambda p: p if isinstance(p, str) else ''
    return {
        'factors': pack({'group': r['group'], 'stage': r['stage'], 'factor': r['factor'], 'c1': rnd(r['component_1'], 3),
                         'c2': rnd(r['component_2'], 3), 'v1': rnd(r['variance_1'], 3), 'v2': rnd(r['variance_2'], 3)}
                        for _, r in factors.iterrows()),
        'methods': list(RESOURCE_METHODS), 'method_labels': RESOURCE_METHODS,
        'held_out': pack({'group': r['group'], 'study': r['study_accession'], 'platform': plat(r['platform']),
                          'method': r['method'], 'hidden': int(r['hidden_healthy']), 'loc': rnd(r['location_error'], 4),
                          'scale': rnd(r['scale_error'], 4), 'disp': rnd(r['pool_dispersion'], 4),
                          'floor': rnd(r['sampling_floor'], 4)} for _, r in held.iterrows()),
        'held_out_medians': medians,
        'combat_vs_matched': {'units': int(pair['units']), 'closer': int(pair['method_closer']),
                              'difference': rnd(pair['median_difference'], 2), 'low': rnd(pair['difference_low'], 2),
                              'high': rnd(pair['difference_high'], 2), 'p': rnd(pair['wilcoxon_p'], 3)},
        'mixing': pack({'group': r['group'], 'cls': r['sample_class'], 'study': r['study_accession'], 'platform': plat(r['platform']),
                        'tier': r['tier'], 'role': r['role'], 'gate': resource_bool(r['gate']), 'samples': int(r['samples']),
                        'before': rnd(r['before'], 3), 'after': rnd(r['after'], 3), 'passed': resource_bool(r['passed'])}
                       for _, r in mixing.iterrows()),
        'mixing_limit': float(mixing['limit'].iloc[0]),
        'outcomes': [[r['outcome'], int(r['records'])] for _, r in outcomes.iterrows()],
        'flows': [[r['step'], r['source'] if isinstance(r['source'], str) else 'not_grouped',
                   r['target'] if isinstance(r['target'], str) else 'not_grouped', int(r['value'])] for _, r in flows.iterrows()],
        'outcome_labels': RESOURCE_OUTCOMES,
        'group_colors': {k or 'not_grouped': v for k, v in
                         pick(config('visualization'), 'correction.group_colors', 'visualization_config.json').items()},
        'records': pick(corrected, 'records', 'correction 07.json'),
        'pooled': pack({'group': r['group'], 'scope': r['scope'], 'check': r['check'], 'samples': int(r['samples']),
                        'before': rnd(r['before'], 3), 'after': rnd(r['after'], 3)} for _, r in pooled.iterrows()),
        'units': {g: {'units': pick(selected, f'groups.{g}.units', 'correction 06.json'),
                      'by_method': pick(selected, f'groups.{g}.by_method', 'correction 06.json'),
                      'tables_units': pick(corrected, f'groups.{g}.units', 'correction 07.json')}
                  for g in pick(corrected, 'groups', 'correction 07.json')},
        'minimum_healthy': pick(config('correction'), 'inputs.minimum_healthy_controls', 'correction_config.json'),
    }


def resource_benchmark(bench):
    tasks = pick(bench, 'tasks', 'cohort 06.json')
    resource_named(tasks, RESOURCE_TASKS, 'cohort 06.json tasks')
    reasons = pick(bench, 'reasons', 'cohort 06.json')
    resource_named(reasons, RESOURCE_BENCH_REASONS, 'cohort 06.json reasons')
    splits = table('cohort/benchmark_splits', ['task', 'scheme', 'fold', 'patient_safe_group', 'role', 'technology',
                                               'first_release', 'positives', 'negatives'])
    resource_named(splits['scheme'], RESOURCE_SCHEMES, 'benchmark_splits schemes')
    unknown = set(splits['role']) - {'train', 'test', 'held_out_not_scored'}
    if unknown:
        raise BuildError(f'benchmark_splits: unexpected roles {sorted(unknown)}')
    summary = table('cohort/cohort_summary', ['patient_safe_group', 'study_accession'])
    natural = lambda a: [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', a)]
    members = {g: sorted(set(rows['study_accession']), key=natural) for g, rows in summary.groupby('patient_safe_group')}
    studies_of = lambda g: sorted({s for part_ in g.split('+') for s in members.get(part_, [])}, key=natural)
    out = {}
    for task, rows in splits.groupby('task', sort=False):
        groups = rows.drop_duplicates('patient_safe_group').sort_values('patient_safe_group')
        folds = {}
        for (scheme, fold), f in rows.groupby(['scheme', 'fold'], sort=False):
            folds.setdefault(scheme, []).append({'fold': fold, 'test': sorted(f.loc[f['role'] == 'test', 'patient_safe_group']),
                                                 'not_scored': sorted(f.loc[f['role'] == 'held_out_not_scored', 'patient_safe_group'])})
        t = tasks[task]
        out[task] = {
            'eligible': t['eligible'], 'positives': t['positives'], 'negatives': t['negatives'],
            'groups': t['patient_safe_groups'], 'scored_groups': t['groups_with_the_minimum_of_each_class'],
            'excluded': [[r, t['excluded'][r]] for r in reasons if r in t['excluded']], 'records': t['records'],
            'folds_by_scheme': t['folds'], 'definition': t['definition'],
            'group_rows': pack({'g': r['patient_safe_group'], 'tech': r['technology'], 'first': str(r['first_release'])[:10],
                                'pos': int(r['positives']), 'neg': int(r['negatives']),
                                'studies': ' '.join(studies_of(r['patient_safe_group']))} for _, r in groups.iterrows()),
            'folds': folds,
        }
    return {'tasks': out, 'task_labels': RESOURCE_TASKS, 'scheme_labels': RESOURCE_SCHEMES,
            'reason_labels': RESOURCE_BENCH_REASONS,
            'split_date': pick(bench, 'policy.temporal_split_date', 'cohort 06.json'),
            'minimum_per_class': pick(bench, 'policy.minimum_per_class', 'cohort 06.json')}


def resource_evaluation():
    """The benchmark of evaluation designs (stage 07): the share of test records with a twin in training under every split,
    the baseline classifiers in every scored unit, the total difference beside the separately measured effects, the
    study-only predictor and the study identification, and every tested stratum of the published signatures. The summaries
    are checked against the units, folds and strata they come from."""
    quality('benchmark')
    audit, models, shortcuts, signatures = (receipt('benchmark', n) for n in (1, 2, 3, 4))
    bench = config('benchmark')
    close = lambda a, b, tol=1e-6: a is not None and b is not None and abs(float(a) - float(b)) <= tol
    # 1. twins across the splits
    leak = review('benchmark', 'split_leakage')
    leak['task'] = leak['task'].fillna('')
    resource_named(leak['scheme'], RESOURCE_EVAL_SCHEMES, 'split_leakage_review schemes')
    resource_named(set(leak['task']) - {''}, RESOURCE_TASKS, 'split_leakage_review tasks')
    for r in leak.itertuples():
        found = pick(audit, f'leakage.{r.task or r.population}.{f"{r.population}:{r.scheme}" if r.task else r.scheme}', 'benchmark 01.json')
        if found['exposed'] != r.exposed or found['test_records'] != r.test_records or not close(found['share_exposed'], r.share_exposed):
            raise BuildError(f'split_leakage_review and benchmark 01.json differ for {r.population} {r.task} {r.scheme}')
    # 2. the baseline classifiers: every scored unit, and the means of the summary recomputed from them
    summary = review('benchmark', 'model_summary')
    units = table('benchmark/model_units', ['task', 'model', 'scheme', 'fold', 'unit', 'records', 'positives', 'negatives', 'status',
                                            'auroc', 'auroc_low', 'auroc_high'])
    scored = units[units['status'] == 'scored']
    resource_named(summary['model'], RESOURCE_MODELS, 'model_summary models')
    resource_named(summary['scheme'], RESOURCE_EVAL_SCHEMES, 'model_summary schemes')
    for r in summary.itertuples():
        mine = scored[(scored['task'] == r.task) & (scored['model'] == r.model) & (scored['scheme'] == r.scheme)]
        if len(mine) != r.units_scored or not close(mine['auroc'].mean(), r.auroc_mean):
            raise BuildError(f'model_summary {r.task} {r.model} {r.scheme}: {len(mine)} scored units with mean {mine["auroc"].mean()}, '
                             f'the summary {r.units_scored} with {r.auroc_mean}')
    # 3. the total difference and the four effects, each with its own interval and studies; they are not parts of the total
    effects = review('benchmark', 'model_inflation_effects')
    resource_named(effects['effect'], RESOURCE_EFFECTS, 'model_inflation_effects effects')
    mean = lambda task, model, scheme: summary.loc[(summary['task'] == task) & (summary['model'] == model) & (summary['scheme'] == scheme), 'auroc_mean'].iloc[0]
    for r in effects[effects['effect'] == 'total_inflation'].itertuples():
        if not close(r.estimate, mean(r.task, r.model, 'naive_record_split') - mean(r.task, r.model, bench['metrics']['inflation_reference'])):
            raise BuildError(f'model_inflation_effects: the total difference of {r.task} {r.model} is not that of model_summary')
    for r in effects.itertuples():
        effect = pick(models, f'inflation_effects.{r.task}.{r.model}.{r.effect}', 'benchmark 02.json')
        if not close(effect['estimate'], r.estimate) or effect['studies'] != r.studies:
            raise BuildError(f'model_inflation_effects and benchmark 02.json differ for {r.task} {r.model} {r.effect}')
    # 4. the study-only predictor and the study identification
    folds = review('benchmark', 'shortcut_folds')
    pooled = folds[(folds['scheme'] == 'naive_record_split') & (folds['scoring'] == 'pooled')]
    within = folds[(folds['scheme'] == 'naive_leave_one_study_out') & (folds['scoring'] == 'within_study') & (folds['status'] == 'scored')]
    study_only = pick(shortcuts, 'study_only', 'benchmark 03.json')
    for task, s in study_only.items():
        p, w = pooled[pooled['task'] == task], within[within['task'] == task]
        if len(p) != s['random_split']['repeats'] or not close(p['auroc'].mean(), s['random_split']['pooled_auroc_mean']) \
                or len(w) != s['held_out_study']['studies_scored'] or not ((w['auroc'] - 0.5).abs() < 1e-12).all():
            raise BuildError(f'shortcut_folds_review and benchmark 03.json differ for the study-only predictor of {task}')
    ident = review('benchmark', 'study_identification')
    found = pick(shortcuts, 'study_identification', 'benchmark 03.json')
    for name in ('aligned', 'corrected'):
        everything = ident[(ident['table'] == name) & (ident['fold'] == 'all')]
        if len(everything) != 1 or not close(everything['balanced_accuracy'].iloc[0], found['balanced_accuracy'][name]):
            raise BuildError(f'study_identification_review and benchmark 03.json differ for the {name} table')
    if not close(found['chance'], 1 / found['units']):
        raise BuildError('benchmark 03.json: the chance of the study identification is not one over its study platforms')
    # 5. the published signatures: every tested stratum of the primary variant, and the medians recomputed from them
    sig_summary = review('benchmark', 'signature_summary')
    sig_summary = sig_summary[sig_summary['primary_variant']]
    tests = table('benchmark/signature_tests', ['signature', 'variant', 'primary_variant', 'comparison', 'study_accession', 'platform',
                                                'tissue', 'status', 'reason', 'positives', 'negatives', 'auroc', 'auroc_low', 'auroc_high',
                                                'derivation_role', 'samples_shared_with_derivation', 'repeat_deposits',
                                                'in_scope_group', 'independent'])
    tests = tests[tests['primary_variant']].copy()
    tests['platform'] = tests['platform'].fillna('')
    resource_named(tests['signature'], RESOURCE_SIGNATURES, 'signature_tests signatures')
    resource_named(set(tests['reason'].fillna('')) - {''}, RESOURCE_SIGNATURE_REASONS, 'signature_tests reasons')
    resource_named(tests['comparison'], RESOURCE_TASKS, 'signature_tests comparisons')
    tested = tests[tests['status'] == 'tested']
    for r in sig_summary.itertuples():
        mine = tested[(tested['signature'] == r.signature) & (tested['comparison'] == r.comparison)]
        if r.scope == bench['signatures']['scope_group']:
            mine = mine[mine['in_scope_group']]
        if r.strata_set == 'independent':
            mine = mine[mine['independent']]
        if len(mine) != r.tested or not close(mine['auroc'].median(), r.auroc_median, 1e-9):
            raise BuildError(f'signature_summary {r.signature} {r.comparison} {r.scope} {r.strata_set}: {len(mine)} strata with median '
                             f'{mine["auroc"].median()}, the summary {r.tested} with {r.auroc_median}')
    definitions = {d['name']: d for d in bench['signatures']['definitions']}
    cite = lambda d: (lambda m: f'{m.group(1)} et al., {m.group(2)}' if m else '')(re.match(r'([A-Z][A-Za-z-]+)\b.*?(\d{4})\.\s*$', d['citation']))
    order = [(s, c) for s in definitions for c in definitions[s]['intended_comparison'] if ((tests['signature'] == s) & (tests['comparison'] == c)).any()]
    order.sort(key=lambda sc: (list(RESOURCE_TASKS).index(sc[1]), list(RESOURCE_SIGNATURES).index(sc[0])))
    skipped = tests[tests['status'] != 'tested'].groupby(['signature', 'comparison', 'reason']).size()
    r6 = lambda v: rnd(v, 6)
    return {
        'records': pick(audit, 'records', 'benchmark 01.json'), 'patient_safe_folds': pick(audit, 'patient_safe_folds.folds', 'benchmark 01.json'),
        'test_fraction': bench['splits']['test_fraction'], 'repeats': bench['splits']['repeats'],
        'leakage': pack({'pop': r.population, 'task': r.task, 'scheme': r.scheme, 'folds': int(r.folds), 'test': int(r.test_records),
                         'exposed': int(r.exposed), 'array': int(r.exposed_same_array), 'patient': int(r.exposed_same_patient),
                         'share': r6(r.share_exposed), 'expected': r6(r.expected_share)} for r in leak.itertuples()),
        'summary': pack({'task': r.task, 'model': r.model, 'scheme': r.scheme, 'units': int(r.units), 'scored': int(r.units_scored),
                         'records': int(r.records_scored), 'mean': r6(r.auroc_mean), 'median': r6(r.auroc_median),
                         'inflation': r6(r.inflation_auroc)} for r in summary.itertuples()),
        'units': pack({'task': r.task, 'model': r.model, 'scheme': r.scheme, 'fold': r.fold, 'unit': r.unit, 'n': int(r.records), 'pos': int(r.positives),
                       'neg': int(r.negatives), 'auroc': r6(r.auroc), 'low': r6(r.auroc_low), 'high': r6(r.auroc_high)}
                      for r in scored.sort_values(['task', 'model', 'scheme', 'auroc', 'unit']).itertuples()),
        'genes': pick(models, 'features.genes_selected', 'benchmark 02.json'),
        'effects': pack({'task': r.task, 'model': r.model, 'effect': r.effect, 'matched': bool(r.matched), 'status': r.status,
                         'first': r6(r.first), 'second': r6(r.second), 'estimate': r6(r.estimate), 'low': r6(r.low), 'high': r6(r.high),
                         'studies': int(r.studies)} for r in effects.itertuples()),
        'minimum_studies': bench['metrics']['inflation_effects']['minimum_studies'],
        'level': bench['metrics']['inflation_effects']['bootstrap']['level'],
        'study_only': study_only,
        'pooled': pack({'task': r.task, 'fold': r.fold, 'auroc': r6(r.auroc), 'n': int(r.records)} for r in pooled.itertuples()),
        'within': pack({'task': r.task, 'study': r.unit, 'n': int(r.records), 'pos': int(r.positives), 'neg': int(r.negatives),
                        'auroc': r6(r.auroc)} for r in within.itertuples()),
        'identification': found | {'folds': pack({'table': r.table, 'fold': r.fold, 'test': int(r.test_samples),
                                                  'balanced': r6(r.balanced_accuracy)} for r in ident.itertuples())},
        'signatures': {
            'order': [list(sc) for sc in order], 'scope_group': bench['signatures']['scope_group'],
            'about': {s: {'label': d['label'], 'cite': cite(d), 'doi': d.get('doi', ''), 'genes': len(d['genes_up']) + len(d['genes_down'])}
                      for s, d in definitions.items() if s in RESOURCE_SIGNATURES},
            'summary': pack({'sig': r.signature, 'comparison': r.comparison, 'scope': r.scope, 'set': r.strata_set, 'strata': int(r.strata),
                             'tested': int(r.tested), 'median': r6(r.auroc_median), 'q1': r6(r.auroc_q1), 'q3': r6(r.auroc_q3)}
                            for r in sig_summary.itertuples()),
            'strata': pack({'sig': r.signature, 'comparison': r.comparison, 'study': r.study_accession, 'platform': r.platform,
                            'tissue': r.tissue, 'pos': int(r.positives), 'neg': int(r.negatives), 'auroc': r6(r.auroc),
                            'low': r6(r.auroc_low), 'high': r6(r.auroc_high), 'role': r.derivation_role or '',
                            'shared': int(r.samples_shared_with_derivation), 'repeats': int(r.repeat_deposits),
                            'scope': bool(r.in_scope_group), 'independent': bool(r.independent)}
                           for r in tested.sort_values(['signature', 'comparison', 'auroc']).itertuples()),
            'skipped': [[s, c, k, int(n)] for (s, c, k), n in skipped.items()],
        },
        'labels': {'models': RESOURCE_MODELS, 'schemes': RESOURCE_EVAL_SCHEMES, 'signatures': RESOURCE_SIGNATURES,
                   'effects': {k: list(v) for k, v in RESOURCE_EFFECTS.items()}, 'reasons': RESOURCE_SIGNATURE_REASONS},
        'scope_tissues': pick(signatures, 'scope.tissues', 'benchmark 04.json'),
        'minimum_group': bench['signatures']['minimum_group'],
    }


def resource_package():
    """The files of the data package (read-only, from data_release/): sizes and members from the manifest, descriptions
    from the readme."""
    manifest = _json(RESOURCE_RELEASE / 'release_manifest.json')
    readme = _existing(RESOURCE_RELEASE / 'readme.txt').read_text(encoding='utf-8')
    block = readme.split('Files (name, then description)', 1)
    if len(block) != 2:
        raise BuildError('data_release/readme.txt has no "Files (name, then description)" list')
    described = dict(line.split('\t', 1) for line in block[1].split('\n\n', 1)[0].strip().splitlines() if '\t' in line)
    files, archives = pick(manifest, 'files', 'release_manifest.json'), pick(manifest, 'archives', 'release_manifest.json')
    member_rows = pick(manifest, 'members', 'release_manifest.json')
    out = []
    for name, info in files.items():
        if name not in described:
            raise BuildError(f'data_release/readme.txt does not describe {name}')
        out.append({'name': name, 'bytes': int(info['size_bytes']), 'members': archives.get(name), 'about': described[name]})
    # the manifest and the checksum list do not list themselves; their sizes are those of the files
    for name in described:
        if name not in files:
            out.append({'name': name, 'bytes': _existing(RESOURCE_RELEASE / name).stat().st_size, 'members': None,
                        'about': described[name]})
    questions = []
    for q, archive, member, about in RESOURCE_QUESTIONS:
        if archive not in files:
            raise BuildError(f'the data package has no {archive}')
        rows = None
        if member:
            key = f'{archive}/{member}'
            if key not in member_rows:
                raise BuildError(f'release_manifest.json has no member {key}')
            rows = member_rows[key].get('rows')
        questions.append({'q': q, 'archive': archive, 'member': member or '', 'rows': rows, 'about': about})
    return {'version': pick(manifest, 'dataset_version', 'release_manifest.json'),
            'license': pick(manifest, 'license.data', 'release_manifest.json'),
            'license_url': pick(manifest, 'license.data_url', 'release_manifest.json'),
            'counts': pick(manifest, 'counts', 'release_manifest.json'),
            'files': pack(out), 'questions': pack(questions)}


def resource():
    """The resource: the data of its chapters, from the studies over time to the files of the data package. Besides the
    pipeline outputs, the last chapter reads the file list of the data package from data_release/ (the manifest and the
    readme), without changing it."""
    ZENODO_DOI = '10.5281/zenodo.24205193'   # the DOI of the data set in Zenodo; without one the page shows "to be assigned"
    for stage in ('cohort', 'validation', 'correction', 'benchmark', 'visualization'):
        quality(stage)
    cohort, links, bench = receipt('cohort', 5), receipt('cohort', 2), receipt('cohort', 6)
    payload = {
        'glossary': [list(entry) for entry in RESOURCE_GLOSSARY],
        'studies': resource_studies(receipt('acquisition', 3)),
        'samples': resource_samples(cohort),
        'repeats': resource_repeats(cohort, links),
        'genes': resource_genes(receipt('annotation', 4)),
        'checks': resource_checks(receipt('validation', 2), receipt('validation', 3), receipt('validation', 4)),
        'correction': resource_correction(receipt('correction', 6), receipt('correction', 7)),
        'benchmark': resource_benchmark(bench),
        'evaluation': resource_evaluation(),
        'package': resource_package() | {'doi': ZENODO_DOI, 'code': RESOURCE_CODE_URL, 'repository': 'Zenodo'},
    }
    return {'resource': payload}


# ---------- studies: data/studies/index.js and data/studies/<id>.js ----------
STUDIES_REASONS = {
    # short labels of the screening reasons (acquisition_config.json screening.reasons gives the full wording)
    'not_a_dataset': 'Not a study',
    'third_party_reanalysis': 'Re-analysis without own samples',
    'not_human': 'No human samples',
    'not_transcriptome_wide_mrna': 'No transcriptome-wide mRNA',
    'single_cell': 'Single-cell profiles only',
    'not_blood': 'No whole blood, leukocyte or PBMC sample',
    'sorted_cells': 'Only sorted cells, vesicles, plasma or serum',
    'ex_vivo_stimulated': 'Only blood treated outside the body',
    'no_sepsis_patients': 'No patient the study diagnoses with sepsis',
    'superseries_of_included_subseries': 'SuperSeries of included studies',
    'overlaps_included_study': 'Shares samples with an included study',
    'controlled_access': 'Controlled access',
    'no_public_reads': 'Registered without public reads',
    'raw_reads_only': 'Raw reads only',
    'no_expression_files': 'No expression file',
    'no_usable_expression_values': 'Files not readable as values',
}
STUDIES_AVAILABILITY = ['controlled_access', 'no_public_reads', 'raw_reads_only', 'no_expression_files', 'no_usable_expression_values']
STUDIES_CHAPTERS = ['overview', 'samples', 'platform', 'genes', 'checks', 'correction', 'benchmark']
STUDIES_GLOSSARY = [
    ('Processing route', 'The lowest data level from which a study platform was read: raw array files, unnormalized intensities, '
                         'read counts or values normalized by the authors. Each study platform has one route, never a mix.'),
    ('Value scale', 'What the numbers of a table mean: log₂ values (each step of one a doubling), log₂-like values '
                    'normalized by the authors, log₂ counts per million, or linear values such as FPKM.'),
    ('Series matrix', "A table of the authors' own processed values that GEO keeps with a study. It was used to check that every "
                      'column of the reprocessed table belongs to its record.'),
    ('Sex plane', 'A chart of the samples of one table by their XIST and their Y-chromosome genes, each scaled from the mean of '
                  'the samples read as one sex to that of the other, so that tables of every platform look alike.'),
    ('Spearman correlation', 'A correlation of ranks: it compares the order of the genes, from low to high, in two profiles. '
                             'It is 1 when the order is the same, whatever the scale of the values.'),
    ('Study mixing', 'How far the samples of one study keep to themselves among their nearest neighbors of the same class, beyond '
                     "what the study's share of the samples would give. It should fall with the correction."),
    ('Split scheme', 'A way of dividing the benchmark into training and test: leave one group out, technology transfer (arrays '
                     'against RNA-seq) or temporal (earlier studies against later ones).'),
    ('Quality flag', 'A mark on a sample whose library is small, whose values span a narrow range, that its investigators reported '
                     'as failed or whose identity is uncertain. Flagged samples are kept.'),
]


def _studies_year(value):
    m = re.match(r'^(\d{4})', str(value or ''))
    return int(m.group(1)) if m else None


def _studies_date(value):
    m = re.match(r'^(\d{4})-(\d{2})-(\d{2})', str(value or ''))
    return m.group(0) if m else None


def _studies_counts(series):
    """{value: count} of a column, in order of frequency and then of value, without empty values."""
    s = series.fillna('').astype(str)
    s = s[s != '']
    # Ties are ordered by value: value_counts orders them differently between pandas versions.
    return {k: int(v) for k, v in sorted(s.value_counts().items(), key=lambda kv: (-kv[1], kv[0]))}


def _studies_split(series, sep='; '):
    out = {}
    for v in series.fillna('').astype(str):
        for k in [x.strip() for x in v.split(sep) if x.strip()]:
            out[k] = out.get(k, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def _studies_quartiles(values):
    v = pd.Series(values).dropna().astype(float)
    if v.empty:
        return None
    return [int(len(v)), rnd(v.min(), 1), rnd(v.quantile(0.25), 1), rnd(v.median(), 1), rnd(v.quantile(0.75), 1), rnd(v.max(), 1)]


def _studies_key(acc, platform):
    return f'{acc}|{platform or ""}'


def studies():
    """Studies: the list of included studies and excluded candidates (data/studies/index.js, with the study accessions in
    'ids'), and one file per study with the data of its chapters (data/studies/<id>.js)."""
    vis, cohort_config, validation_config = config('visualization'), config('cohort'), config('validation')
    correction_config, acquisition_config = config('correction'), config('acquisition')
    acquisition, cohort05 = receipt('acquisition', 3), receipt('cohort', 5)
    classes, families = vis['class_order'], vis['family_order']

    # ---- the studies and their screening ----
    decisions = table('acquisition/study_decisions', ['source', 'accession', 'title', 'assay', 'release_date', 'url', 'decision',
                                                     'reason', 'basis', 'flag', 'overlaps_with', 'recorded_samples', 'also_in',
                                                     'superseries_of', 'subseries_of'])
    reasons_config = pick(acquisition_config, 'screening.reasons', 'acquisition_config.json')
    unknown = sorted(set(decisions['reason']) - set(reasons_config))
    if unknown:
        raise BuildError(f'study_decisions: reasons {unknown} are not in acquisition_config.json')
    missing_labels = sorted(set(reasons_config) - set(STUDIES_REASONS) - {'eligible'})
    if missing_labels:
        raise BuildError(f'studies(): no short label for the screening reasons {missing_labels}')
    included = decisions[decisions['decision'] == 'include'].copy()
    excluded = decisions[decisions['decision'] != 'include'].copy()
    if len(included) != pick(acquisition, 'included.studies', '03.json'):
        raise BuildError('study_decisions and acquisition 03.json differ in the number of included studies')
    meta = table('acquisition/studies', ['accession', 'summary', 'design', 'pubmed_ids', 'platforms'])
    meta = meta.drop_duplicates('accession').set_index('accession')

    # ---- samples ----
    cs = table('cohort/cohort_samples', ['source', 'study_accession', 'sample_accession', 'platform', 'cohort_group',
                                         'patient_safe_group', 'sample_id', 'is_primary_record', 'records_of_sample',
                                         'linked_by', 'sample_class', 'infection_sepsis_not_stated', 'class_basis',
                                         'class_differs_between_records', 'sepsis_phase', 'sepsis_definition', 'tissue', 'sex',
                                         'age_years', 'age_group', 'patient_identified', 'patient_key', 'samples_of_patient',
                                         'one_per_patient', 'time_point', 'time_order', 'time_unit', 'outcome',
                                         'outcome_definition', 'septic_shock', 'quality_flags', 'any_flag'])
    # the records without expression values are not in the sample table; the ledger of the cohort stage lists them with their reason
    without = table('cohort/records_without_values', ['study_accession', 'sample_accession', 'platform', 'reason'])
    without['platform'] = without['platform'].fillna('')
    ids = sorted(included['accession'])
    if sorted(set(cs['study_accession']) | set(without['study_accession'])) != ids:
        raise BuildError('cohort_samples, records_without_values and the included studies of study_decisions name different studies')
    for column, allowed in [('sample_class', classes)]:
        bad = sorted(set(cs[column]) - set(allowed))
        if bad:
            raise BuildError(f'cohort_samples: {column} values {bad} are not in the configuration')
    labels =table('cohort/label_sources', ['study_accession', 'sample_accession', 'platform', 'variable', 'value', 'basis',
                                            'fields', 'recorded', 'study_statement', 'note'])
    routes = table('expression/processing_routes', ['accession', 'platform', 'family', 'platform_title', 'samples', 'tier',
                                                    'data_source', 'covered', 'uncovered', 'value_kind', 'flag', 'basis'])
    mindex = table('annotation/matrix_index', ['accession', 'platform', 'method', 'value_scale', 'in_matrix', 'features',
                                               'identifier_shape', 'identifier_example', 'reason'])
    gms = table('annotation/gene_matrix_summary', ['accession', 'platform', 'value_scale', 'rule', 'samples', 'features',
                                                   'rows_used', 'rows_left_out', 'genes', 'protein_coding_genes',
                                                   'genes_with_several_rows'])
    fms = table('annotation/feature_map_summary', ['accession', 'platform', 'annotation_source', 'annotation_key', 'features',
                                                   'one_gene', 'several_genes', 'no_gene_in_annotation', 'gene_not_in_hgnc',
                                                   'accession_not_resolved', 'not_in_annotation', 'genes', 'protein_coding_genes'])
    affy = table('expression/affymetrix_studies', ['accession', 'platform', 'chip_type', 'package', 'pm_probes', 'features',
                                                   'first_scan', 'last_scan'])
    single = table('expression/single_channel_studies', ['accession', 'platform', 'layout', 'software', 'features', 'probes',
                                                         'first_scan', 'last_scan'])
    inten = table('expression/intensity_table_studies', ['accession', 'platform', 'correction', 'probes', 'tables'])
    counts = table('expression/count_table_studies', ['accession', 'platform', 'data_source', 'tables', 'genes', 'count_kind',
                                                      'library_size_minimum', 'library_size_median', 'library_size_maximum'])
    author = table('expression/author_matrix_studies', ['accession', 'platform', 'values', 'scale', 'features',
                                                        'value_minimum', 'value_median', 'value_maximum'])
    value_check = review('expression', 'value_check', ['accession', 'platform', 'status', 'samples_compared',
                                                       'own_column_most_similar', 'exception'])
    for frame in (routes, mindex, gms, fms, affy, single, inten, counts, author, value_check):
        frame['platform'] = frame['platform'].fillna('')
    platform_titles = table('expression/platforms', ['platform', 'family', 'title'])

    # ---- validation ----
    sex_tables = table('validation/sex_check_matrices', ['study_accession', 'platform', 'samples', 'status', 'female_gene_read',
                                                         'male_genes_read', 'called_female', 'called_male', 'ambiguous',
                                                         'with_recorded_sex', 'concordant', 'discordant', 'concordance',
                                                         'recorded_sex_unreliable'])
    sex_samples = table('validation/sex_check_samples', ['study_accession', 'sample_accession', 'platform', 'recorded_sex',
                                                         'difference', 'expression_sex', 'agreement'])
    splits = figure_data('SuppFig3', 'splits', ['study_accession', 'platform', 'samples', 'lower_group', 'higher_group', 'gap',
                                                'separation', 'midpoint'])
    split_rows = figure_data('SuppFig3', 'splits')
    splits['tissue'] = split_rows['tissue'] if 'tissue' in split_rows else None
    plane = figure_data('Fig5', 'sex_plane', ['study_accession', 'platform', 'sample_accession', 'x', 'y'])
    for frame in (plane, splits):
        frame['platform'] = frame['platform'].fillna('').astype(str)
    score_strata = table('validation/sepsis_score_by_study', ['study_accession', 'platform', 'tissue', 'value_scale', 'sepsis',
                                                              'healthy_control', 'sepsis_left_out_by_phase', 'genes_present',
                                                              'genes_as_published', 'score_computed', 'auroc', 'auroc_low',
                                                              'auroc_high', 'patients_sepsis', 'patients_healthy_control',
                                                              'auroc_one_per_patient', 'sepsis_one_per_patient',
                                                              'healthy_control_one_per_patient', 'score_derivation'])
    score_samples = table('validation/sepsis_score_samples', ['study_accession', 'sample_accession', 'platform', 'tissue',
                                                              'sample_class', 'one_per_patient', 'score'])
    genes_by_study = table('validation/sepsis_genes_by_study', ['study_accession', 'platform', 'tissue', 'test', 'role',
                                                                'published_symbol', 'present', 'standardised_difference',
                                                                'as_published'])
    vmat = table('validation/validation_matrices', ['study_accession', 'platform', 'value_scale', 'comparable_scale'])
    vsam = table('validation/validation_samples', ['study_accession', 'sample_accession', 'platform', 'sample_class',
                                                   'sepsis_phase', 'tissue'])
    repeats = table('validation/repeat_pairs', ['study_a', 'platform_a', 'study_b', 'platform_b', 'sample_a', 'sample_b',
                                                'linked_by', 'spearman', 'most_similar'])
    markers = table('validation/cell_markers_by_study', ['study_accession', 'platform', 'cell_type', 'markers_read', 'samples',
                                                         'other_samples', 'auroc'])
    minimum_group = pick(validation_config, 'biology.minimum_group', 'validation_config.json')
    score_phases = pick(validation_config, 'biology.sepsis_phases', 'validation_config.json')

    # ---- correction ----
    units = table('correction/correction_units', ['group', 'study_accession', 'platform', 'family', 'value_scale', 'role',
                                                  'records', 'primary_records', 'estimation_samples', 'healthy_anchors', 'sepsis',
                                                  'non_septic_control', 'other_condition', 'neonatal', 'tier', 'genes',
                                                  'genes_in_group_set', 'status', 'reason'])
    methods_table = table('correction/unit_methods', ['group', 'study_accession', 'platform', 'method_tried', 'method', 'basis',
                                                      'status', 'reason'])
    checks = table('correction/unit_checks', ['group', 'study_accession', 'platform', 'check', 'sample_class', 'gate',
                                              'applicable', 'samples', 'before', 'after', 'limit', 'passed'])
    held = table('correction/held_out_test', ['group', 'study_accession', 'platform', 'method', 'hidden_healthy',
                                              'estimation_samples', 'location_error', 'scale_error', 'pool_dispersion',
                                              'sampling_floor'])
    ledger = table('correction/correction_ledger', ['study_accession', 'sample_accession', 'platform', 'group', 'unit_method',
                                                    'outcome', 'reason'])
    for frame in (units, methods_table, checks, held, ledger):
        frame['platform'] = frame['platform'].fillna('')
        frame['group'] = frame['group'].fillna('')
    corrected = {g: table(f'correction/corrected/{g}_samples', ['study_accession', 'sample_accession', 'platform', 'unit_method'])
                 for g in correction_config['groups']}
    minimum_healthy = pick(correction_config, 'inputs.minimum_healthy_controls', 'correction_config.json')

    # ---- benchmark ----
    btasks = table('cohort/benchmark_tasks', ['study_accession', 'sample_accession', 'platform', 'task', 'patient_safe_group',
                                              'eligible', 'reason', 'label'])
    bsplits = table('cohort/benchmark_splits', ['task', 'scheme', 'fold', 'patient_safe_group', 'role', 'technology',
                                                'first_release', 'positives', 'negatives'])
    task_order = [t['task'] for t in pick(cohort_config, 'benchmark.tasks', 'cohort_config.json')]
    unknown = sorted(set(btasks['task']) - set(task_order))
    if unknown:
        raise BuildError(f'benchmark_tasks: tasks {unknown} are not in cohort_config.json')

    # ---- groups and links ----
    groups = table('cohort/cohort_groups', ['cohort_group', 'study_accession', 'platform', 'release_date', 'records',
                                            'primary_records', 'repeated_records'])
    links = table('cohort/sample_links', ['study_accession', 'sample_accession', 'platform', 'sample_id', 'is_primary_record',
                                          'linked_by', 'primary_study_accession'])

    # ======================================================== the list
    fam_of_platform = {}
    for _, r in routes.iterrows():
        fam_of_platform[_studies_key(r['accession'], r['platform'])] = r['family']
    unknown = sorted(set(routes['family']) - set(families))
    if unknown:
        raise BuildError(f'processing_routes: families {unknown} are not in the visualization configuration')
    by_study = {acc: part for acc, part in cs.groupby('study_accession')}
    without_by_study = {acc: part for acc, part in without.groupby('study_accession')}
    inc = included.set_index('accession')
    rows = []
    for acc in ids:
        part, gone = by_study.get(acc, cs.iloc[0:0]), without_by_study.get(acc, without.iloc[0:0])
        r = inc.loc[acc]
        plats = routes[routes['accession'] == acc]
        rec_by_family = {}
        for _, p in plats.iterrows():
            n = int((part['platform'].fillna('') == p['platform']).sum()) + int((gone['platform'] == p['platform']).sum())
            rec_by_family[p['family']] = rec_by_family.get(p['family'], 0) + n
        fams = sorted(rec_by_family, key=lambda f: (-rec_by_family[f], families.index(f)))
        if not fams:
            raise BuildError(f'{acc} has no study platform in processing_routes')
        # records: the sample records with expression values, which the sample table holds; without: those of the ledger
        rows.append({
            'id': acc, 'src': r['source'], 'title': r['title'], 'date': _studies_date(r['release_date']),
            'year': _studies_year(r['release_date']), 'assay': r['assay'],
            'fam': [families.index(f) for f in fams],
            'plats': [p for p in plats['platform']],
            'records': int(len(part)), 'distinct': int(part['is_primary_record'].sum()), 'without': int(len(gone)),
            'cls': [int((part['sample_class'] == c).sum()) for c in classes],
            'clsd': [int(((part['sample_class'] == c) & part['is_primary_record']).sum()) for c in classes],
            'group': sorted(part['patient_safe_group'].unique().tolist()),
            'flag': r['flag'] or '',
        })
    if sum(x['records'] + x['without'] for x in rows) != pick(acquisition, 'included.sample_records', '03.json'):
        raise BuildError('the records of the included studies do not add up to acquisition 03.json included.sample_records')
    if sum(x['records'] for x in rows) != pick(cohort05, 'records.total', 'cohort 05.json'):
        raise BuildError('the records with values of the included studies do not add up to the records of cohort 05.json')
    if sum(x['distinct'] for x in rows) != pick(cohort05, 'distinct_samples.total', 'cohort 05.json'):
        raise BuildError('the primary records of the included studies do not add up to the distinct samples of cohort 05.json')

    excl = []
    for _, r in excluded.sort_values(['reason', 'accession']).iterrows():
        excl.append({'id': r['accession'], 'src': r['source'], 'title': r['title'] or '', 'year': _studies_year(r['release_date']),
                     'assay': r['assay'] or '', 'reason': r['reason'], 'basis': r['basis'] or '', 'url': r['url'] or '',
                     'overlaps': r['overlaps_with'] or '', 'flag': r['flag'] or ''})
    order = [k for k in reasons_config if k != 'eligible']

    # ======================================================== comparisons with the whole resource
    scales = figure_data('Fig4', 'scales', ['study_accession', 'platform', 'value_scale', 'q05', 'q25', 'q50', 'q75', 'q95'])
    gpt = figure_data('Fig4', 'genes_per_table', ['family', 'study_accession', 'platform', 'genes', 'protein_coding_genes'])
    mapping = figure_data('Fig4', 'mapping')
    mapping_order = vis['mapping_order']
    _columns_present(mapping, ['family'] + [f'{k}_share' for k in mapping_order], 'figure_data.csv Fig4 mapping')
    ages = figure_data('Fig3', 'ages', ['sample_class', 'count', 'min', '25%', '50%', '75%', 'max'])
    score_all = figure_data('Fig5', 'score', ['study_accession', 'platform', 'tissue', 'sepsis', 'healthy_control', 'auroc',
                                              'auroc_low', 'auroc_high'])
    pbmc_view = figure_data('SuppFig6', 'views_pbmc', ['sample_id', 'study_accession', 'sample_class', 'pc1_before', 'pc2_before',
                                                       'pc1_after', 'pc2_after'])
    factors = figure_data('Fig6', 'factors', ['group', 'stage', 'variance_1', 'variance_2'])
    held_medians = {}
    for (g, m), part in held.groupby(['group', 'method']):
        held_medians.setdefault(g, {})[m] = {'location': rnd(part['location_error'].median(), 3),
                                              'dispersion': rnd(part['pool_dispersion'].median(), 3),
                                              'floor': rnd(part['sampling_floor'].median(), 3), 'units': int(len(part))}
    pv_studies = sorted(pbmc_view['study_accession'].unique())
    pv = {'n': len(pbmc_view), 'scale': 10,
          'x0': [rnd(10 * v) for v in pbmc_view['pc1_before']], 'y0': [rnd(10 * v) for v in pbmc_view['pc2_before']],
          'x1': [rnd(10 * v) for v in pbmc_view['pc1_after']], 'y1': [rnd(10 * v) for v in pbmc_view['pc2_after']],
          'cls': ''.join(str(classes.index(c)) for c in pbmc_view['sample_class']),
          'studies': pv_studies, 'study': [pv_studies.index(s) for s in pbmc_view['study_accession']],
          'sample': pbmc_view['sample_id'].str.split(':', n=1).str[1].tolist(),
          'variance': {r['stage']: [rnd(r['variance_1'], 4), rnd(r['variance_2'], 4)]
                       for _, r in factors[factors['group'] == 'pbmc'].drop_duplicates('stage').iterrows()}}
    if len(pbmc_view) != corrected['pbmc'].shape[0]:
        raise BuildError('SuppFig6 views_pbmc and the corrected PBMC table differ in their number of samples')
    coverage_cols = {'sex': lambda d: d['sex'] != 'not_recorded', 'age': lambda d: d['age_group'] != 'not_recorded',
                     'patient': lambda d: d['patient_identified'].astype(bool),
                     'time': lambda d: d['time_unit'].fillna('') != '',
                     'outcome': lambda d: d['outcome'] != 'not_recorded', 'shock': lambda d: d['septic_shock'] != 'not_recorded'}
    resource_cov = {k: rnd(f(cs).mean(), 4) for k, f in coverage_cols.items()}
    sex_conc = sex_tables[sex_tables['concordance'].notna()]
    index = {
        'ids': ids,
        'chapters': STUDIES_CHAPTERS,
        'list': pack(rows),
        'excluded': pack(excl),
        'reasons': {'order': order, 'label': {k: STUDIES_REASONS[k] for k in order}, 'text': {k: reasons_config[k] for k in order},
                    'availability': STUDIES_AVAILABILITY},
        'candidates': int(len(decisions)),
        'mapping': {'order': mapping_order, 'label': vis['mapping_labels'], 'color': vis['mapping_colors'],
                    'by_family': {r['family']: [rnd(r[f'{k}_share'], 4) for k in mapping_order] for _, r in mapping.iterrows()}},
        'links': {k: v.replace(chr(10), ' ') for k, v in vis['link_labels'].items()},
        'flags': {k: v.format(**cohort_config['quality']) for k, v in vis['flag_labels'].items()},
        'basis': vis['basis_labels'],
        'methods': pick(vis, 'correction.method_order', 'visualization_config.json'),
        'resource': {
            'coverage': resource_cov,
            'ages': {r['sample_class']: [int(r['count']), rnd(r['min'], 1), rnd(r['25%'], 1), rnd(r['50%'], 1), rnd(r['75%'], 1),
                                         rnd(r['max'], 1)] for _, r in ages.iterrows()},
            'scales': pack([{'id': r['study_accession'], 'p': r['platform'] if isinstance(r['platform'], str) else '',
                             'scale': r['value_scale'], 'q': [rnd(r[c], 2) for c in ('q05', 'q25', 'q50', 'q75', 'q95')]}
                            for _, r in scales.iterrows()]),
            'genes': pack([{'id': r['study_accession'], 'p': r['platform'] if isinstance(r['platform'], str) else '',
                            'fam': families.index(r['family']), 'genes': int(r['genes']), 'coding': int(r['protein_coding_genes'])}
                           for _, r in gpt.iterrows()]),
            'sex': pack([{'id': r['study_accession'], 'p': r['platform'] or '', 'compared': int(r['concordant'] + r['discordant']),
                          'conc': rnd(r['concordance'], 4), 'unreliable': bool(r['recorded_sex_unreliable'])}
                         for _, r in sex_conc.iterrows()]),
            'scores': pack([{'id': r['study_accession'], 'p': r['platform'] if isinstance(r['platform'], str) else '',
                             'tissue': r['tissue'], 'sepsis': int(r['sepsis']), 'healthy': int(r['healthy_control']),
                             'auroc': rnd(r['auroc'], 3), 'low': rnd(r['auroc_low'], 3), 'high': rnd(r['auroc_high'], 3)}
                            for _, r in score_all.iterrows()]),
            'held_out': held_medians,
            'score_minimum': minimum_group, 'minimum_healthy': minimum_healthy,
            'sex_rules': {k: pick(validation_config, f'sex.{k}', 'validation_config.json')
                          for k in ('minimum_samples', 'minimum_separation', 'minimum_group', 'ambiguous_zone', 'unreliable_below',
                                    'minimum_compared')},
            'minimum_unit': pick(correction_config, 'inputs.minimum_unit_samples', 'correction_config.json'),
            'bench_minimum': pick(cohort_config, 'benchmark.minimum_per_class', 'cohort_config.json'),
            'bench_split_date': pick(cohort_config, 'benchmark.temporal_split_date', 'cohort_config.json'),
            'folds': {t: {s: int(part['fold'].nunique()) for s, part in bsplits[bsplits['task'] == t].groupby('scheme')}
                      for t in task_order},
            'tasks': task_order,
        },
        'pbmc_view': pv,
        'glossary': [list(e) for e in STUDIES_GLOSSARY],
    }

    # ======================================================== one file per study
    files = {'studies': index}
    lab_vars = ['sample_class', 'tissue', 'sex', 'age', 'time', 'outcome', 'septic_shock', 'patient']
    keep = {'sample_class': 8, 'tissue': 3, 'sex': 3, 'age': 2, 'time': 4, 'outcome': 4, 'septic_shock': 4, 'patient': 2}
    labels_by_study = {acc: part for acc, part in labels.groupby('study_accession')}
    links_by_sample = links.groupby('sample_id')
    sex_by_study = {acc: part for acc, part in sex_samples.groupby('study_accession')}
    plane_by = {k: part for k, part in plane.groupby(['study_accession', 'platform'])}
    for acc in ids:
        part, gone = by_study.get(acc, cs.iloc[0:0]), without_by_study.get(acc, without.iloc[0:0])
        r = inc.loc[acc]
        m = meta.loc[acc] if acc in meta.index else None
        cohort_groups = sorted(part['cohort_group'].unique())
        psg = sorted(part['patient_safe_group'].unique())
        part = part.assign(platform=part['platform'].fillna(''))

        # -- overview: links, shares with other studies --
        shares = {}
        multi = part[part['records_of_sample'] > 1]
        mine = set(zip(part['study_accession'], part['sample_accession']))
        for sid in multi['sample_id'].unique():
            others = links_by_sample.get_group(sid)
            here = others[others['study_accession'] == acc]
            primary_here = bool(here['is_primary_record'].any())
            primary_study = others['primary_study_accession'].iloc[0]
            for _, o in others[others['study_accession'] != acc].iterrows():
                s = shares.setdefault(o['study_accession'], {'id': o['study_accession'], 'samples': 0, 'first_here': 0,
                                                             'first_there': 0, 'first_elsewhere': 0, 'by': {}})
                s['samples'] += 1
                if primary_here:
                    s['first_here'] += 1
                elif primary_study == o['study_accession']:
                    s['first_there'] += 1
                else:
                    s['first_elsewhere'] += 1
                for b in [x.strip() for x in str(o['linked_by'] or '').split(';') if x.strip()]:
                    s['by'][b] = s['by'].get(b, 0) + 1
        share_rows = sorted(shares.values(), key=lambda s: (-s['samples'], s['id']))
        lk = links[links['study_accession'] == acc]
        repeat_other = int(((~lk['is_primary_record'].astype(bool)) & (lk['primary_study_accession'] != acc)).sum())
        repeat_same = int(((~lk['is_primary_record'].astype(bool)) & (lk['primary_study_accession'] == acc)).sum())
        members = groups[groups['cohort_group'].isin(cohort_groups)]
        group_rows = [{'group': g['cohort_group'], 'id': g['study_accession'], 'p': g['platform'] or '', 'records': int(g['records']),
                       'primary': int(g['primary_records'])} for _, g in members.sort_values(['cohort_group', 'release_date']).iterrows()]

        # -- samples --
        sepsis = part[part['sample_class'] == 'sepsis']
        tp = part[part['time_unit'].fillna('') != '']
        times = []
        for (order_, unit), tpart in tp.groupby(['time_order', 'time_unit'], dropna=False):
            times.append({'order': rnd(order_, 2) if order_ == order_ else None, 'unit': unit, 'n': int(len(tpart)),
                          'words': list(_studies_counts(tpart['time_point']).keys())[:3]})
        times.sort(key=lambda t: (t['order'] is None, t['order'] if t['order'] is not None else 0, t['unit']))
        wording = {}
        lab = labels_by_study.get(acc)
        if lab is not None:
            lab = lab.fillna('')
            for var in lab_vars:
                v = lab[lab['variable'] == var]
                n_not = int((v['basis'] == 'not_recorded').sum())
                v = v[~v['basis'].isin(['not_recorded', 'from_sample_class'])]
                if var in ('age', 'patient'):
                    grp = v.groupby(['fields', 'basis'], sort=False)
                    items = []
                    for (fields, basis), gp in grp:
                        ex = gp.iloc[[0, len(gp) // 2]].drop_duplicates('sample_accession') if len(gp) > 1 else gp.iloc[[0]]
                        for _, e in ex.iterrows():
                            items.append({'field': fields, 'recorded': e['recorded'], 'value': e['value'], 'n': int(len(gp)),
                                          'basis': basis, 'statement': e['study_statement'], 'note': e['note'],
                                          'ex': e['sample_accession']})
                    items = items[:keep[var]]
                else:
                    g = v.groupby(['fields', 'recorded', 'value', 'basis', 'study_statement', 'note'], sort=False)
                    items = [{'field': k[0], 'recorded': k[1], 'value': k[2], 'basis': k[3], 'statement': k[4], 'note': k[5],
                              'n': int(len(gp)), 'ex': gp['sample_accession'].iloc[0]} for k, gp in g]
                    items.sort(key=lambda x: -x['n'])
                    total_groups = len(items)
                    items = items[:keep[var]]
                    if total_groups > len(items):
                        items.append({'more': total_groups - len(items),
                                      'n': int(len(v)) - sum(x['n'] for x in items)})
                wording[var] = {'items': items, 'not_recorded': n_not, 'records': int(len(lab[lab['variable'] == var]))}
        ages_by_class = {}
        for c in classes:
            q = _studies_quartiles(part.loc[part['sample_class'] == c, 'age_years'])
            if q:
                ages_by_class[c] = q
        samples = {
            'tissue': _studies_counts(part['tissue']), 'age_group': _studies_counts(part['age_group']),
            'sex': _studies_counts(part['sex']), 'ages': ages_by_class,
            'coverage': {k: rnd(f(part).mean(), 4) for k, f in coverage_cols.items()},
            'times': times, 'time_records': int(len(tp)), 'time_patients': int(tp.loc[tp['patient_identified'].astype(bool), 'patient_key'].nunique()),
            'sepsis': int(len(sepsis)),
            'outcome': _studies_counts(sepsis['outcome']), 'outcome_all': _studies_counts(part['outcome']),
            'outcome_definition': _studies_counts(part.loc[part['outcome'] != 'not_recorded', 'outcome_definition']),
            'shock': _studies_counts(sepsis['septic_shock']),
            'phase': _studies_counts(sepsis['sepsis_phase']),
            'sepsis3': int((part['sepsis_definition'] == 'sepsis-3').sum()),
            'patients': int(part.loc[part['patient_identified'].astype(bool), 'patient_key'].nunique()),
            'with_patient': int(part['patient_identified'].astype(bool).sum()),
            'one_per_patient': int(part['one_per_patient'].astype(bool).sum()),
            'one_per_patient_primary': int((part['one_per_patient'].astype(bool) & part['is_primary_record']).sum()),
            'flags': _studies_split(part['quality_flags']), 'flagged': int(part['any_flag'].astype(bool).sum()),
            'basis': _studies_counts(part['class_basis']),
            'infection_not_stated': int(part['infection_sepsis_not_stated'].astype(bool).sum()),
            'class_differs': int(part['class_differs_between_records'].astype(bool).sum()),
            'wording': wording,
        }

        # -- platforms, processing and genes --
        plats = []
        for _, p in routes[routes['accession'] == acc].iterrows():
            key = (acc, p['platform'])
            pp, pg = part[part['platform'] == p['platform']], gone[gone['platform'] == p['platform']]
            sel = lambda f: f[(f['accession'] == acc) & (f['platform'] == p['platform'])]
            mi, gm, fm = sel(mindex), sel(gms), sel(fms)
            title = p['platform_title'] or ''
            item = {
                'p': p['platform'], 'title': title, 'fam': families.index(p['family']),
                'records': int(len(pp)), 'distinct': int(pp['is_primary_record'].sum()), 'without': int(len(pg)),
                'cls': [int((pp['sample_class'] == c).sum()) for c in classes],
                'tier': p['tier'], 'source': p['data_source'], 'covered': int(p['covered']), 'uncovered': int(p['uncovered']),
                'value_kind': p['value_kind'] or '', 'flag': p['flag'] or '', 'basis': p['basis'] or '',
                'method': mi['method'].iloc[0] if len(mi) else '', 'scale': mi['value_scale'].iloc[0] if len(mi) else '',
                'features': int(mi['features'].iloc[0]) if len(mi) and mi['features'].iloc[0] == mi['features'].iloc[0] else None,
                'id_example': mi['identifier_example'].iloc[0] if len(mi) else '',
                'no_values': _studies_counts(pg['reason']),
                'tissue': _studies_counts(pp['tissue']),
            }
            a = sel(affy)
            if len(a):
                a = a.iloc[0]
                item['affy'] = {'chip': a['chip_type'], 'package': a['package'], 'probes': int(a['pm_probes']),
                                'sets': int(a['features']), 'first': _studies_date(a['first_scan']), 'last': _studies_date(a['last_scan'])}
            s1 = sel(single)
            if len(s1):
                s1 = s1.iloc[0]
                item['single'] = {'layout': s1['layout'], 'software': s1['software'] or '', 'features': int(s1['features']),
                                  'probes': int(s1['probes']), 'first': _studies_date(s1['first_scan']), 'last': _studies_date(s1['last_scan'])}
            it = sel(inten)
            if len(it):
                it = it.iloc[0]
                item['inten'] = {'correction': it['correction'], 'probes': int(it['probes'])}
            ct = sel(counts)
            if len(ct):
                ct = ct.iloc[0]
                item['counts'] = {'source': ct['data_source'], 'kind': ct['count_kind'], 'genes': int(ct['genes']),
                                  'lib': [rnd(ct['library_size_minimum']), rnd(ct['library_size_median']), rnd(ct['library_size_maximum'])]}
            au = sel(author)
            if len(au):
                au = au.iloc[0]
                item['author'] = {'values': au['values'] or '', 'scale': au['scale'] or ''}
            vc = sel(value_check)
            if len(vc):
                vc = vc.iloc[0]
                item['check'] = {'status': vc['status'], 'compared': int(vc['samples_compared']),
                                 'own': int(vc['own_column_most_similar']), 'exception': bool(vc['exception'])}
            if len(gm):
                g0, f0 = gm.iloc[0], fm.iloc[0] if len(fm) else None
                item['genes'] = {'rule': g0['rule'], 'features': int(g0['features']), 'rows_used': int(g0['rows_used']),
                                 'rows_left_out': int(g0['rows_left_out']), 'genes': int(g0['genes']),
                                 'coding': int(g0['protein_coding_genes']), 'several_rows': int(g0['genes_with_several_rows'])}
                if f0 is not None:
                    item['genes']['source'] = f0['annotation_source']
                    item['genes']['key'] = f0['annotation_key']
                    item['genes']['map'] = [int(f0[k]) for k in mapping_order]
            plats.append(item)

        # -- checks --
        sx = []
        ssam = sex_by_study.get(acc)
        for _, t in sex_tables[sex_tables['study_accession'] == acc].iterrows():
            pf = t['platform'] or ''
            pts = ssam[ssam['platform'].fillna('') == pf] if ssam is not None else ssam
            pl = plane_by.get((acc, pf))
            xy = {}
            if pl is not None:
                xy = {s: (rnd(x, 3), rnd(y, 3)) for s, x, y in zip(pl['sample_accession'], pl['x'], pl['y'])}
            sp = splits[(splits['study_accession'] == acc) & (splits['platform'] == pf)]
            pts_rows = []
            by_tissue = int((sp['gap'] == sp['gap']).sum()) > 1 or t['status'] == 'tested_by_tissue'
            tissue_of = dict(zip(part.loc[part['platform'] == pf, 'sample_accession'], part.loc[part['platform'] == pf, 'tissue']))
            if pts is not None:
                for _, s in pts.iterrows():
                    if s['difference'] != s['difference']:
                        continue
                    pos = xy.get(s['sample_accession'])
                    row = {'s': s['sample_accession'], 'r': s['recorded_sex'], 'e': s['expression_sex'],
                           'a': s['agreement'], 'd': rnd(s['difference'], 3),
                           'x': pos[0] if pos else None, 'y': pos[1] if pos else None}
                    if by_tissue:
                        row['t'] = tissue_of.get(s['sample_accession'], '')
                    pts_rows.append(row)
            sx.append({
                'p': pf, 'status': t['status'], 'samples': int(t['samples']), 'xist': bool(t['female_gene_read']),
                'y_genes': int(t['male_genes_read']), 'female': int(t['called_female']), 'male': int(t['called_male']),
                'ambiguous': int(t['ambiguous']), 'recorded': int(t['with_recorded_sex']), 'concordant': int(t['concordant']),
                'discordant': int(t['discordant']), 'conc': rnd(t['concordance'], 4), 'unreliable': bool(t['recorded_sex_unreliable']),
                'splits': [{'tissue': (s['tissue'] if isinstance(s['tissue'], str) else ''), 'samples': int(s['samples']),
                            'lower': int(s['lower_group']), 'higher': int(s['higher_group']), 'gap': rnd(s['gap'], 3),
                            'separation': rnd(s['separation'], 2), 'midpoint': rnd(s['midpoint'], 3)}
                           for _, s in sp.iterrows() if s['gap'] == s['gap']],
                'plane': pl is not None,
                'points': pack(pts_rows) if pts_rows else [],
            })
        strata = []
        for _, t in score_strata[score_strata['study_accession'] == acc].iterrows():
            ss = score_samples[(score_samples['study_accession'] == acc) & (score_samples['platform'].fillna('') == (t['platform'] or ''))
                               & (score_samples['tissue'] == t['tissue'])]
            gb = genes_by_study[(genes_by_study['study_accession'] == acc) & (genes_by_study['platform'].fillna('') == (t['platform'] or ''))
                                & (genes_by_study['tissue'] == t['tissue'])]
            strata.append({
                'p': t['platform'] or '', 'tissue': t['tissue'], 'sepsis': int(t['sepsis']), 'healthy': int(t['healthy_control']),
                'left_out': int(t['sepsis_left_out_by_phase']), 'genes': int(t['genes_present']), 'as_published': int(t['genes_as_published']),
                'auroc': rnd(t['auroc'], 4), 'low': rnd(t['auroc_low'], 3), 'high': rnd(t['auroc_high'], 3),
                'auroc_pp': rnd(t['auroc_one_per_patient'], 4), 'sepsis_pp': int(t['sepsis_one_per_patient']),
                'healthy_pp': int(t['healthy_control_one_per_patient']), 'role': t['score_derivation'] or '',
                'pos': [rnd(v, 4) for v in ss.loc[ss['sample_class'] == 'sepsis', 'score']],
                'neg': [rnd(v, 4) for v in ss.loc[ss['sample_class'] == 'healthy_control', 'score']],
                'pos_s': ss.loc[ss['sample_class'] == 'sepsis', 'sample_accession'].tolist(),
                'neg_s': ss.loc[ss['sample_class'] == 'healthy_control', 'sample_accession'].tolist(),
                'dirs': pack([{'g': g['published_symbol'], 'role': g['role'], 'test': g['test'], 'present': bool(g['present']),
                               'z': rnd(g['standardised_difference'], 2), 'ok': bool(g['as_published']) if g['present'] else None}
                              for _, g in gb.iterrows()]) if len(gb) else [],
            })
        # where no score was computed: the sepsis samples (phases of the score) and healthy controls in a table, by table and tissue
        cand = vsam[vsam['study_accession'] == acc].copy()
        cand['platform'] = cand['platform'].fillna('')
        comparable = {(v['platform'] or ''): (bool(v['comparable_scale']), v['value_scale'])
                      for _, v in vmat[vmat['study_accession'] == acc].iterrows()}
        not_scored = []
        for (pf, tissue), cp in cand.groupby(['platform', 'tissue']):
            if any(s['p'] == pf and s['tissue'] == tissue for s in strata):
                continue
            ns = int(((cp['sample_class'] == 'sepsis') & cp['sepsis_phase'].isin(score_phases)).sum())
            nh = int((cp['sample_class'] == 'healthy_control').sum())
            comp, scale = comparable.get(pf, (False, ''))
            not_scored.append({'p': pf, 'tissue': tissue, 'sepsis': ns, 'healthy': nh, 'comparable': comp, 'scale': scale or ''})
        rp = repeats[(repeats['study_a'] == acc) | (repeats['study_b'] == acc)]
        rep_rows = []
        if len(rp):
            other = rp['study_b'].where(rp['study_a'] == acc, rp['study_a'])
            for o, op in rp.groupby(other):
                rep_rows.append({'id': o, 'pairs': int(len(op)), 'lowest': rnd(op['spearman'].min(), 4),
                                 'median': rnd(op['spearman'].median(), 4), 'best': int(op['most_similar'].sum()),
                                 'by': _studies_split(op['linked_by'])})
            rep_rows.sort(key=lambda x: -x['pairs'])
        mk = [{'p': m_['platform'] or '', 'cell': m_['cell_type'], 'markers': int(m_['markers_read']), 'samples': int(m_['samples']),
               'others': int(m_['other_samples']), 'auroc': rnd(m_['auroc'], 3)} for _, m_ in markers[markers['study_accession'] == acc].iterrows()]

        # -- correction --
        cu = []
        for _, u in units[units['study_accession'] == acc].iterrows():
            sel = lambda f: f[(f['study_accession'] == acc) & (f['platform'] == u['platform']) & (f['group'] == u['group'])]
            mt, ch, ho = sel(methods_table), sel(checks), sel(held)
            mt = mt.iloc[0] if len(mt) else None
            cu.append({
                'group': u['group'], 'p': u['platform'], 'role': u['role'], 'tier': u['tier'], 'status': u['status'],
                'reason': u['reason'] or '', 'records': int(u['records']), 'primary': int(u['primary_records']),
                'estimation': int(u['estimation_samples']), 'anchors': int(u['healthy_anchors']), 'sepsis': int(u['sepsis']),
                'neonatal': bool(u['neonatal']), 'genes': int(u['genes_in_group_set']), 'scale': u['value_scale'],
                'method': (mt['method'] or '') if mt is not None else '', 'tried': (mt['method_tried'] or '') if mt is not None else '',
                'basis': (mt['basis'] or '') if mt is not None else '', 'mstatus': (mt['status'] or '') if mt is not None else '',
                'mreason': (mt['reason'] or '') if mt is not None else '',
                'checks': [{'check': c['check'], 'cls': c['sample_class'] or '', 'gate': bool(c['gate']), 'applicable': bool(c['applicable']),
                            'samples': int(c['samples']), 'before': rnd(c['before'], 3), 'after': rnd(c['after'], 3),
                            'limit': rnd(c['limit'], 3), 'passed': bool(c['passed'])} for _, c in ch.iterrows()],
                'held': [{'method': h['method'], 'hidden': int(h['hidden_healthy']), 'estimation': int(h['estimation_samples']),
                          'location': rnd(h['location_error'], 3), 'scale': rnd(h['scale_error'], 3),
                          'dispersion': rnd(h['pool_dispersion'], 3), 'floor': rnd(h['sampling_floor'], 3)} for _, h in ho.iterrows()],
            })
        lg = ledger[ledger['study_accession'] == acc].fillna('')
        ledger_rows = [{'group': k[0], 'method': k[1], 'outcome': k[2], 'reason': k[3], 'n': int(len(gp))}
                       for k, gp in lg.groupby(['group', 'unit_method', 'outcome', 'reason'])]
        ledger_rows.sort(key=lambda x: -x['n'])
        in_tables = {g: int((t['study_accession'] == acc).sum()) for g, t in corrected.items()}

        # -- benchmark --
        bt = btasks[btasks['study_accession'] == acc]
        bench = []
        for task in task_order:
            tp_ = bt[bt['task'] == task]
            el = tp_[tp_['eligible'].astype(bool)]
            gsel = bsplits[(bsplits['task'] == task) & bsplits['patient_safe_group'].isin(psg)]
            folds = [{'group': f['patient_safe_group'], 'scheme': f['scheme'], 'fold': f['fold'], 'role': f['role']}
                     for _, f in gsel.iterrows()]
            # roles summarized by group and scheme: train in n folds, test in m, held out but not scored in k
            roles = {}
            for f in folds:
                key = f'{f["group"]}|{f["scheme"]}'
                d = roles.setdefault(key, {'group': f['group'], 'scheme': f['scheme'], 'train': 0, 'test': 0, 'held_out_not_scored': 0,
                                           'test_folds': []})
                d[f['role']] = d.get(f['role'], 0) + 1
                if f['role'] != 'train':
                    d['test_folds'].append(f['fold'])
            gt = gsel.drop_duplicates('patient_safe_group')
            bench.append({
                'task': task, 'records': int(len(tp_)), 'eligible': int(len(el)),
                'pos': int((el['label'] == 1).sum()), 'neg': int((el['label'] == 0).sum()),
                'reasons': _studies_counts(tp_.loc[~tp_['eligible'].astype(bool), 'reason']),
                'groups': [{'group': g['patient_safe_group'], 'tech': g['technology'], 'first': g['first_release'],
                            'pos': int(g['positives']), 'neg': int(g['negatives'])} for _, g in gt.iterrows()],
                'roles': list(roles.values()),
            })
        group_studies = sorted(set(cs.loc[cs['patient_safe_group'].isin(psg), 'study_accession']) - {acc})

        files[f'studies/{acc}'] = {
            'id': acc, 'src': r['source'], 'title': r['title'], 'url': r['url'] or '', 'date': _studies_date(r['release_date']),
            'assay': r['assay'], 'also_in': [x.strip() for x in str(r['also_in'] or '').split(';') if x.strip()],
            'superseries': [x.strip() for x in str(r['superseries_of'] or '').split(';') if x.strip()],
            'subseries_of': [x.strip() for x in str(r['subseries_of'] or '').split(';') if x.strip()],
            'basis': r['basis'] or '', 'flag': r['flag'] or '',
            'summary': (str(m['summary'])[:3000] if m is not None and m['summary'] else ''),
            'design': (str(m['design'])[:1500] if m is not None and m['design'] else ''),
            'pubmed': [x for x in str(m['pubmed_ids'] if m is not None else '').split(';') if x.strip()][:12],
            'records': int(len(part)), 'distinct': int(part['is_primary_record'].sum()), 'without': int(len(gone)),
            'repeat_other': repeat_other, 'repeat_same': repeat_same,
            'cls': [int((part['sample_class'] == c).sum()) for c in classes],
            'clsd': [int(((part['sample_class'] == c) & part['is_primary_record']).sum()) for c in classes],
            'cohort_groups': cohort_groups, 'group_members': group_rows, 'psg': psg, 'psg_studies': group_studies,
            'shares': share_rows,
            'samples': samples, 'platforms': plats,
            'checks': {'sex': sx, 'scores': strata, 'not_scored': not_scored, 'repeats': rep_rows, 'markers': mk},
            'correction': {'units': cu, 'ledger': ledger_rows, 'in_tables': in_tables},
            'bench': bench,
        }
    return files


# ---------- methods: data/methods.js ----------
METHODS_GLOSSARY = [
    ('Candidate study', 'A study that matched at least one search term, counted once across the repositories. Screening then decides '
                        'whether it is included.'),
    ('Recall', 'Of the studies known to belong to the field, the share that the search found.'),
    ('Screening reason', 'The reason a candidate study was left out. Reasons are checked in a fixed order and the first that applies '
                         'is recorded, so every study is counted once.'),
    ('CEL file', 'The file an Affymetrix scanner writes: the light of every cell of the chip, with the chip type and the moment of '
                 'the scan.'),
    ('Probe and probe set', 'A probe is a short DNA sequence on an array that matches part of one messenger RNA; the probes for one '
                            'transcript form a probe set, whose values are combined into one number.'),
    ('Background correction', 'Removing the light that a spot of an array gives when nothing has bound to it, so that what remains '
                              'estimates the signal. The normexp model keeps the result small but positive near the background.'),
    ('Quantile normalization', 'Giving every sample of a group the same distribution of values: the value at each rank is replaced '
                               'by the average of the values at that rank across the samples, so each sample keeps the order of '
                               'its own values.'),
    ('Median polish', 'Combining the probes of a probe set by removing row and column medians in turn, so that one misbehaving probe '
                      'does not move the result.'),
    ('RMA', 'Robust multi-array average: background correction, quantile normalization and median polish of Affymetrix arrays, '
            'giving log₂ values.'),
    ('Library size', 'The total count of reads of one RNA-seq sample; a count means something only relative to it.'),
    ('TMM', 'Trimmed mean of M-values: a factor for each RNA-seq library that corrects its size when a few genes take a large share '
            'of the reads. The library size times the factor is the effective library size.'),
    ('Counts per million (CPM)', 'A count divided by the effective library size and multiplied by a million. Its log₂ is the '
                                 'value of the RNA-seq tables.'),
    ('Prior count', 'A small amount added to every RNA-seq count before taking the logarithm, so that zero counts get a value and '
                    'small counts do not swing widely.'),
    ('Expressed-gene flag', 'A mark on a gene of a count table that has enough reads in enough samples of its study (the rule of '
                            'filterByExpr); it removes nothing.'),
    ('Feature', 'A row of a matrix as the technology measures it: a probe, a probe set, a transcript cluster or a gene.'),
    ('Gene-level table', 'The table beside each matrix with one row per HGNC gene, built only from rows that map to exactly one gene.'),
    ('Curation rule', 'A written rule that reads one field of a sample record and gives a value from a fixed list. The rules of one '
                      'variable are tried in order and the first that gives a value decides.'),
    ('Primary record', 'Of the records of one array or library, the record of the study released first; its study and accession '
                       'are the sample identifier.'),
    ('Cohort group', 'Study platforms that share at least one array or library.'),
    ('One sample per patient', 'The one sample marked for each patient: without a quality flag, earliest in time, of whole blood, '
                               'and lowest in accession, in that order of preference.'),
    ('Two-means split', 'Dividing the samples of a table into two groups along one score: two centers start at the lowest and '
                        'highest value, each sample joins the nearer center, the centers move to their group means, and this '
                        'repeats until no sample changes group.'),
    ('Separation', 'The gap between the means of two groups divided by their pooled standard deviation, so that it does not depend '
                   'on the scale of the platform.'),
    ('Geometric mean', 'The n-th root of the product of n positive values. It lies between the smallest and largest value and is '
                       'pulled less than the ordinary mean by one large value.'),
    ('Sepsis MetaScore', 'The geometric mean of six genes reported higher in sepsis minus that of five genes reported lower, '
                         'published by Sweeney and colleagues in 2015.'),
    ('ROC curve and AUROC', 'The receiver operating characteristic (ROC) curve draws, for every threshold of a score, the share of '
                            'sepsis samples above it (sensitivity) against the share of controls above it (false-positive rate). '
                            'The area under it (AUROC) is the share of all pairs of a sepsis sample and a control in which the '
                            'sepsis sample has the higher score, a tie counting one half: 0.5 is chance and 1 is perfect order.'),
    ('Bootstrap interval', 'The range of a result over many repeats in which patients are drawn at random with replacement; the '
                           'middle 95% of the repeats is the interval.'),
    ('Study unit', 'The samples of one study platform in one corrected group (bulk blood or PBMC): the unit that is corrected.'),
    ('Anchor', 'The healthy controls of all units with enough of them. Each gene has an anchor mean and an anchor standard '
               'deviation, the spread of controls around their own unit\'s mean.'),
    ('Location and scale', 'How far a unit reads a gene above or below the anchor, and how much more or less it spreads it, both in '
                           'anchor standard deviations.'),
    ('Empirical Bayes shrinkage', 'Pulling each gene\'s estimate toward the middle of the estimates of all genes of the unit, more '
                                  'strongly the less information it carries; the distribution pulled toward is fitted from the '
                                  'data.'),
    ('Matched groups', 'The correction of a unit without enough healthy controls: the mean and spread of each of its classes are '
                       'lined up with those of the same class and age group in the units corrected on healthy controls.'),
    ('Held-out test', 'Hiding the healthy controls of a unit that has them, correcting the unit from its other samples, and '
                      'measuring how close the hidden controls land to the controls of the other units.'),
    ('Location error', 'In the held-out test, the root mean square over genes of the distance between the hidden controls and the '
                       'other units\' controls, in within-study standard deviations.'),
    ('Sampling floor', 'The location error that chance alone gives for a group of hidden controls of that size: half the root mean '
                       'square distance between the means of two random halves of them.'),
    ('Fold', 'One division of the samples into training and test.'),
]


def methods_check(condition, message):
    if not condition:
        raise BuildError(f'methods: {message}')


def methods_one(frame, message):
    methods_check(len(frame) == 1, message)
    return frame.iloc[0]


def methods_stem(study, platform):
    return f'{study}__{platform}'


def methods_platform_title(title):
    return title.split('] ', 1)[1] if str(title).startswith('[') else title


# Short names of the screening reasons, in plain words; build_data.py stops if the configuration names another reason.
METHODS_REASONS = {
    'not_a_dataset': 'Not a dataset', 'third_party_reanalysis': 'Re-analysis owning no samples', 'not_human': 'Not human',
    'not_transcriptome_wide_mrna': 'No transcriptome-wide mRNA', 'single_cell': 'Single-cell profiles only',
    'not_blood': 'No whole blood, leukocyte or PBMC sample', 'sorted_cells': 'Only sorted cells, vesicles, plasma or serum',
    'ex_vivo_stimulated': 'Every blood sample treated outside the body',
    'no_sepsis_patients': 'No blood sample from a patient the study diagnoses with sepsis',
    'superseries_of_included_subseries': 'SuperSeries of included studies', 'overlaps_included_study': 'Shares samples with an included study',
    'controlled_access': 'Controlled access', 'no_public_reads': 'Registered without public reads', 'raw_reads_only': 'Raw reads only',
    'no_expression_files': 'No expression file', 'no_usable_expression_values': 'Files not readable as values'}


def methods_search():
    """Chapter 1: the searches, duplicates across repositories, the screening reasons in their order, the funnel, and the
    worked examples GSE103119 (excluded) and GSE65682 (included)."""
    acq = config('acquisition')
    terms = pick(acq, 'search.terms', 'acquisition_config.json')
    sources = pick(acq, 'deduplication.source_priority', 'acquisition_config.json')
    reasons = [r for r in pick(acq, 'screening.reasons', 'acquisition_config.json') if r != 'eligible']
    unknown = [r for r in reasons if r not in METHODS_REASONS]
    methods_check(not unknown, f'screening reasons without a short name: {unknown}')
    census, screen, qc = receipt('acquisition', 1), receipt('acquisition', 3), quality('acquisition')['summary']
    availability = list(pick(screen, 'eligible_but_unavailable', 'acquisition 03.json'))
    methods_check(reasons[-len(availability):] == availability, 'the availability reasons are not the last screening reasons')
    hits = [[None] * len(terms) for _ in sources]
    for q in census['queries']:
        hits[sources.index(q['source'])][terms.index(q['term'])] = q['hits']
    methods_check(all(v is not None for row in hits for v in row), 'a search of one term in one repository is missing')
    attrition = {a['source']: a for a in census['attrition']}
    by_source = [{'source': s, 'records': attrition[s]['records'], 'candidates': attrition[s].get('candidate', 0),
                  'duplicates': attrition[s].get('cross_source_duplicate', 0),
                  'no_runs': attrition[s].get('no_human_rnaseq_runs', 0), 'not_expression': attrition[s].get('study_type_not_expression', 0)}
                 for s in sources]
    methods_check(sum(s['records'] for s in by_source) == qc['census_records'], 'the census records by repository do not add up')
    methods_check(sum(s['candidates'] for s in by_source) == qc['candidates'], 'the candidates by repository do not add up')
    # duplicates and the repository of the record that each duplicates
    census_table = table('acquisition/census', ['source', 'accession', 'census_status', 'duplicate_of'])
    owners = census_table[census_table['census_status'] == 'candidate'].set_index('accession')['source']
    dups = census_table[census_table['census_status'] == 'cross_source_duplicate']
    owner = owners.reindex(dups['duplicate_of']).to_numpy()
    methods_check(not pd.isna(owner).any(), 'a duplicate names no candidate')
    pairs = pd.DataFrame({'source': dups['source'].to_numpy(), 'owner': owner}).value_counts()
    duplicates = [{'source': s, 'of': o, 'records': int(n)} for (s, o), n in sorted(pairs.items(), key=lambda x: (sources.index(x[0][0]), sources.index(x[0][1])))]
    # every candidate, with the terms that reached it and its screening reason (the term simulation)
    candidates = table('acquisition/candidates', ['source', 'accession', 'matched_terms', 'release_date', 'platform', 'also_in'])
    decisions = table('acquisition/study_decisions', ['accession', 'decision', 'reason', 'samples', 'basis'])
    merged = candidates.merge(decisions, on='accession', validate='one_to_one').sort_values('accession').reset_index(drop=True)
    methods_check(len(merged) == qc['candidates'], 'every candidate needs one decision')
    order = ['eligible'] + reasons
    methods_check(set(merged['reason']) <= set(order), 'a decision names an unknown reason')
    rows = []
    for r in merged.itertuples(index=False):
        matched = [t for t in str(r.matched_terms).split(';') if t]
        # an ENA record can reach the census without a recorded term (ENA matches words inside words); it counts as reached
        methods_check(set(matched) <= set(terms) and (matched or r.source == 'ena'), f'{r.accession} matched terms that are not search terms')
        rows.append([r.accession, sources.index(r.source), sum(1 << terms.index(t) for t in matched), order.index(r.reason), int(r.samples)])
    exclusions = qc['exclusions']
    methods_check(sum(exclusions.values()) + qc['included_studies'] == qc['candidates'], 'the screening decisions do not add up')
    funnel = [{'step': a['step'], 'removed': a['removed'], 'remaining': a['remaining'], 'records': a.get('sample_records')}
              for a in screen['attrition']]
    methods_check([a['step'] for a in funnel[1:]] == reasons, 'the steps of the screening receipt are not the screening reasons in order')
    recall = pick(census, 'recall', 'acquisition 01.json')
    found = {s['accession'] for s in recall['studies'] if s['found']}

    # worked examples
    titles = table('annotation/matrix_index', ['platform', 'platform_title']).drop_duplicates('platform').set_index('platform')['platform_title']
    summary = table('cohort/cohort_summary', ['study_accession', 'platform', 'records', 'sepsis', 'healthy_control',
                                              'non_septic_control', 'other_condition', 'outcome_recorded'])
    attributes = table('acquisition/sample_attributes', ['study_accession', 'key', 'value'])

    def example(accession):
        row = methods_one(merged[merged['accession'] == accession], f'{accession} is not a candidate')
        matched = [t for t in terms if t in str(row['matched_terms']).split(';')]
        return {'accession': accession, 'source': row['source'], 'year': str(row['release_date'])[:4], 'records': int(row['samples']),
                'platform': row['platform'], 'platform_title': methods_platform_title(titles.get(row['platform'], '')),
                'terms': matched, 'reason': row['reason'], 'basis': row['basis'], 'copy': row['also_in'] or '',
                'recall': accession in found}
    excluded, included = example('GSE103119'), example('GSE65682')
    methods_check(excluded['reason'] == 'no_sepsis_patients' and included['reason'] == 'eligible', 'the worked examples changed decision')
    condition = attributes[(attributes['study_accession'] == 'GSE103119') & (attributes['key'] == 'condition')]['value'].value_counts()
    methods_check(condition.sum() == excluded['records'], 'every record of GSE103119 needs a condition')
    excluded['conditions'] = [[k, int(v)] for k, v in condition.items()]
    own = summary[summary['study_accession'] == 'GSE65682']
    methods_check(len(own) == 1 and int(own['records'].iloc[0]) == included['records'], 'GSE65682 is no longer one study platform')
    included['study_platform'] = f"GSE65682|{own['platform'].iloc[0]}"
    included['classes'] = {c: int(own[c].iloc[0]) for c in ('sepsis', 'healthy_control', 'non_septic_control', 'other_condition')}
    included['outcome_recorded'] = int(own['outcome_recorded'].iloc[0])
    methods_check(included['outcome_recorded'] <= included['records'], 'GSE65682 records more outcomes than records')
    return {
        'terms': terms, 'sources': sources, 'hits': hits, 'by_source': by_source, 'duplicates': duplicates,
        'census_records': qc['census_records'], 'candidates_total': qc['candidates'], 'included': qc['included_studies'],
        'included_records': qc['included_sample_records'], 'harvested_records': qc['sample_records'],
        'reasons': [{'key': r, 'label': METHODS_REASONS[r], 'description': acq['screening']['reasons'][r],
                     'group': 'availability' if r in availability else 'scientific', 'studies': exclusions[r]} for r in reasons],
        'funnel': funnel, 'searched': str(census['searched_at_utc'])[:10],
        'recall': {'checked': recall['checked'], 'found': recall['found'],
                   'outside_scope': len(pick(acq, 'recall.outside_scope.studies', 'acquisition_config.json'))},
        'candidates': {'__f': ['accession', 'source', 'terms', 'reason', 'records'], '__r': rows},
        'excluded': excluded, 'included_example': included,
    }


def methods_day(text):
    """'2004-03-23T...' -> '2004-03-23' (the page writes dates in words)."""
    return str(text)[:10]


def methods_arrays():
    """Chapter 2: a real Affymetrix array of GSE26440 through RMA, and the background of a GSE54514 BeadChip (normexp)."""
    study_id, platform, array_id = 'GSE26440', 'GPL570', 'GSM648590'
    studies = table('expression/affymetrix_studies')
    study = methods_one(studies[(studies['accession'] == study_id) & (studies['platform'] == platform)], f'{study_id} is no longer one Affymetrix study platform')
    arrays = table('expression/affymetrix_arrays', ['study_accession', 'platform', 'sample_accession', 'chip_type', 'member', 'scan_time', 'raw_pm_median'])
    arrays = arrays[(arrays['study_accession'] == study_id) & (arrays['platform'] == platform)]
    array = methods_one(arrays[arrays['sample_accession'] == array_id], f'{array_id} is no longer an array of {study_id}')
    methods_check(array['scan_time'] == arrays['scan_time'].min(), f'{array_id} is no longer the first array of {study_id} to be scanned')
    matrix = table(f'expression/study_matrices/{methods_stem(study_id, platform)}').set_index('feature_id')
    methods_check(len(arrays) == study['samples'] == matrix.shape[1] and len(matrix) == study['features'], f'the {study_id} matrix differs from its record')
    medians = matrix.median()
    arrays = arrays.sort_values(['scan_time', 'sample_accession']).reset_index(drop=True)
    rma = {
        'study': study_id, 'platform': platform, 'array': array_id, 'arrays': int(study['samples']), 'chip': study['chip_type'],
        'package': study['package'], 'first_year': str(study['first_scan'])[:4], 'last_year': str(study['last_scan'])[:4],
        'pm_probes': int(study['pm_probes']), 'probe_sets': int(study['features']), 'cel': array['member'], 'scan_day': methods_day(array['scan_time']),
        'raw_median': rnd(array['raw_pm_median'], 2), 'array_median': rnd(medians[array_id], 4),
        # every array of the study: its accession, the median of its raw probe intensities and the median of its RMA values
        'each': {'__f': ['accession', 'raw', 'rma'], '__r': [[a.sample_accession, rnd(a.raw_pm_median, 2), rnd(medians[a.sample_accession], 4)]
                                                          for a in arrays.itertuples(index=False)]},
    }
    # the background of the first array of GSE54514, from its detection p-values (limma's normexp)
    expr = config('expression')
    offset, detected = pick(expr, 'intensity_tables.offset', 'expression_config.json'), pick(expr, 'intensity_tables.detected_below', 'expression_config.json')
    bead = table('expression/intensity_table_arrays')
    row = methods_one(bead[(bead['study_accession'] == 'GSE54514') & (bead['sample_accession'] == 'GSM1317896')], 'GSM1317896 is no longer an array of GSE54514')
    methods_check(row['correction'] == 'detection_p_values', 'GSE54514 is no longer corrected from its detection p-values')
    mu, sigma, alpha = (float(row[n]) for n in ('background_mean', 'background_sd', 'signal_mean'))

    def normexp(x):
        shifted = x - mu - sigma ** 2 / alpha
        z = shifted / sigma
        return shifted, shifted + sigma * (math.exp(-z * z / 2) / math.sqrt(2 * math.pi)) / (0.5 * (1 + math.erf(z / math.sqrt(2))))
    probes = []
    for name, x in (('median', float(row['signal_median'])), ('bright', float(row['signal_percentile_99'])), ('dim', 100.0)):
        m, s = normexp(x)
        probes.append({'probe': name, 'x': rnd(x, 2), 'subtracted': rnd(x - mu, 2), 'm': rnd(m, 2), 'signal': rnd(s, 3), 'log': rnd(math.log2(s + offset), 3)})
    background = {'study': 'GSE54514', 'array': 'GSM1317896', 'label': row['label'], 'mu': rnd(mu, 4), 'sigma': rnd(sigma, 4), 'alpha': rnd(alpha, 4),
                  'detected': rnd(row['detected_share'], 4), 'detected_below': detected, 'offset': offset, 'probes': probes,
                  'invented': ['dim']}
    routes = {t['tier']: t['study_platforms'] for t in pick(quality('expression'), 'summary.by_tier', 'expression qc.json')}
    return {'rma': rma, 'background': background, 'routes': routes,
            'affymetrix': {'study_platforms': len(studies), 'arrays': int(studies['samples'].sum())}}


def methods_counts():
    """Chapter 3: log2 counts per million of HLA-DRA in GSM4683511 of GSE154918, the value of a zero count in every
    normalization group, and the expressed-gene flag of GSE154918."""
    study_id, platform, sample_id, gene = 'GSE154918', 'GPL20301', 'GSM4683511', 'HLA-DRA'
    prior = pick(config('expression'), 'count_tables.prior_count', 'expression_config.json')
    samples = table('expression/count_table_samples', ['study_accession', 'platform', 'sample_accession', 'normalization_group', 'library_size', 'norm_factor', 'log2_cpm_of_zero'])
    study = samples[(samples['study_accession'] == study_id) & (samples['platform'] == platform)]
    methods_check(set(study['normalization_group']) == {''}, f'{study_id} is now normalized in groups')
    sample = methods_one(study[study['sample_accession'] == sample_id], f'{sample_id} is no longer a library of {study_id}')
    effective = study['library_size'] * study['norm_factor']
    mean = float(effective.mean())
    methods_check(abs(sum(math.log(f) for f in study['norm_factor'])) < 1e-6 * len(study), f'the TMM factors of {study_id} no longer multiply to one')
    size, factor = float(sample['library_size']), float(sample['norm_factor'])
    own = size * factor
    scaled = prior * own / mean
    stem = methods_stem(study_id, platform)
    genes = table(f'annotation/feature_maps/{stem}', ['feature_id', 'symbol', 'status'])
    feature = methods_one(genes[(genes['symbol'] == gene) & (genes['status'] == 'one_gene')], f'{gene} no longer has one row in {study_id}')['feature_id']
    counts = table(f'expression/study_counts/{stem}', ['feature_id', sample_id]).set_index('feature_id')[sample_id]
    values = table(f'expression/study_matrices/{stem}', ['feature_id', sample_id]).set_index('feature_id')[sample_id]
    count = float(counts[feature])
    result = math.log2((count + scaled) / (own + 2 * scaled) * 1e6)
    zero = math.log2(scaled / (own + 2 * scaled) * 1e6)
    methods_check(abs(result - float(values[feature])) < 1e-3 and abs(zero - float(sample['log2_cpm_of_zero'])) < 1e-3
                  and abs(float(counts.sum()) - size) < 0.5, f'the log2 counts per million of {study_id} no longer reproduce the stored values')
    cpm = {'study': study_id, 'platform': platform, 'sample': sample_id, 'gene': gene, 'feature': str(feature), 'libraries': len(study),
           'count': rnd(count, 2), 'library': int(round(size)), 'factor': rnd(factor, 6), 'mean_effective': rnd(mean, 1), 'prior': prior,
           'value': rnd(float(values[feature]), 4), 'zero': rnd(float(sample['log2_cpm_of_zero']), 4)}
    # the value of a zero count, c / (mean + 2c) x 10^6, in every normalization group
    groups = samples.assign(effective=samples['library_size'] * samples['norm_factor']).groupby(
        ['study_accession', 'platform', 'normalization_group'], sort=True).agg(zero=('log2_cpm_of_zero', 'first'), lo=('log2_cpm_of_zero', 'min'),
                                                                               hi=('log2_cpm_of_zero', 'max'), mean=('effective', 'mean'), n=('effective', 'size'))
    methods_check(((groups['hi'] - groups['lo']).abs() < 1e-6).all(), 'the log2 value of a zero count differs within a normalization group')
    zeros = [[s, p, g, rnd(r['mean'] / 1e6, 4), rnd(r['zero'], 4), int(r['n'])] for (s, p, g), r in groups.iterrows()]
    rc = receipt('expression', 7)
    sources = {k: v['study_platforms'] for k, v in pick(rc, 'by_data_source', 'expression 07.json').items()}
    # the expressed-gene flag (edgeR filterByExpr) of the same study
    policy = pick(config('annotation'), 'gene_matrices.expressed', 'annotation_config.json')
    n = len(study)
    methods_check(n > policy['large_n'], f'{study_id} no longer has more libraries than the large-n limit')
    needed = policy['large_n'] + (n - policy['large_n']) * policy['minimum_proportion']
    median = float(effective.median())
    cutoff = policy['minimum_count'] / median * 1e6
    rows = table(f'annotation/gene_rows/{stem}', ['symbol', 'expressed']).set_index('symbol')['expressed']
    gene_counts = table(f'annotation/gene_counts/{stem}').drop(columns='hgnc_id').set_index('symbol')
    eff = pd.Series(effective.to_numpy(), index=study['sample_accession'].to_numpy())

    def above(symbol):
        return int((gene_counts.loc[symbol, eff.index].astype(float) / eff * 1e6 >= cutoff).sum())
    cohort = table('cohort/cohort_samples', ['study_accession', 'platform', 'sex'])
    sexes = cohort[(cohort['study_accession'] == study_id) & (cohort['platform'] == platform)]['sex'].value_counts()
    checked = []
    for symbol in (gene, 'XIST', 'IL10'):
        checked.append({'gene': symbol, 'samples': above(symbol), 'total': rnd(float(gene_counts.loc[symbol].sum()), 1), 'flagged': bool(rows[symbol])})
    methods_check(checked[0]['flagged'] and checked[0]['samples'] == n and not checked[1]['flagged'] and not checked[2]['flagged'],
                  'the expressed flags of HLA-DRA, XIST and IL10 in GSE154918 have changed')
    expressed = {'large_n': policy['large_n'], 'proportion': policy['minimum_proportion'], 'minimum_count': policy['minimum_count'],
                 'minimum_total': policy['minimum_total_count'], 'median_effective': rnd(median, 1), 'genes': len(rows), 'flagged': int(rows.sum()),
                 'women': int(sexes.get('female', 0)), 'checked': checked}
    return {'cpm': cpm, 'zeros': {'__f': ['study', 'platform', 'group', 'mean_millions', 'zero', 'libraries'], '__r': zeros},
            'study_platforms': rc['study_platforms']['processed'], 'samples': rc['samples']['in_matrix'], 'sources': sources, 'expressed': expressed}


def methods_genes():
    """Chapter 4: the rows of GSE26440 by mapping status, two rows with several candidate genes, the two probe sets of
    HLA-DRA and the one kept, matching through previous symbols, and a count gene whose identical copy counts once."""
    study_id, platform, array_id, gene = 'GSE26440', 'GPL570', 'GSM648590', 'HLA-DRA'
    stem = methods_stem(study_id, platform)
    maps = table(f'annotation/feature_maps/{stem}')
    status = maps['status'].value_counts()
    names = ['one_gene', 'several_genes', 'no_gene_in_annotation', 'gene_not_in_hgnc']
    methods_check(set(status.index) <= set(names), f'the map of {study_id} holds another status')
    reference = table('annotation/genes', ['hgnc_id', 'symbol'])
    symbols = dict(zip(reference['hgnc_id'], reference['symbol']))
    several = []
    for feature in ('209312_x_at', '205033_s_at'):
        row = methods_one(maps[maps['feature_id'] == feature], f'{feature} is no longer a row of {study_id}')
        methods_check(row['status'] == 'several_genes', f'{feature} no longer maps to several genes')
        several.append({'feature': feature, 'genes': sorted(symbols[g] for g in row['candidate_genes'].split('|'))})
    rows = table(f'annotation/gene_rows/{stem}')
    own = methods_one(rows[rows['symbol'] == gene], f'{gene} is no longer one gene of {study_id}')
    features = own['feature_ids'].split('|')
    methods_check(own['rule'] == 'highest_mean' and len(features) == own['rows'] == 2, f'{gene} is no longer kept from two rows by the highest mean')
    matrix = table(f'expression/study_matrices/{stem}').set_index('feature_id')
    arrays = table('expression/affymetrix_arrays', ['study_accession', 'platform', 'sample_accession', 'scan_time'])
    order = arrays[(arrays['study_accession'] == study_id) & (arrays['platform'] == platform)].sort_values(['scan_time', 'sample_accession'])['sample_accession'].tolist()
    means = matrix.loc[features].mean(axis=1).sort_values(ascending=False)
    methods_check(means.index[0] == own['kept_feature_id'], f'the kept row of {gene} is no longer the one with the higher mean')
    level = table(f'annotation/gene_matrices/{stem}', ['symbol', array_id]).set_index('symbol')[array_id]
    methods_check(float(level[gene]) == float(matrix.loc[own['kept_feature_id'], array_id]), f'the gene-level value of {gene} is not its kept row')
    maprow = maps.set_index('feature_id').loc[features]
    probe_sets = [{'feature': f, 'matched_by': maprow.loc[f, 'matched_by'], 'values': [rnd(v, 4) for v in matrix.loc[f, order]],
                   'mean': rnd(means[f], 4)} for f in means.index]
    summary = table('annotation/gene_matrix_summary')
    srow = methods_one(summary[(summary['accession'] == study_id) & (summary['platform'] == platform)], f'{study_id} is no longer one gene-level table')
    hla = {'study': study_id, 'platform': platform, 'array': array_id, 'gene': gene, 'arrays': order, 'probe_sets': probe_sets,
           'kept': own['kept_feature_id'], 'value': rnd(float(level[gene]), 4)}
    # rows reached through a symbol that HGNC no longer approves, over every map of the resource
    matched = {'previous_symbol': 0, 'alias_symbol': 0}
    for path in table('annotation/feature_map_summary', ['map'])['map']:
        name = str(path).replace('data/processed/', '').removesuffix('.parquet')
        found = table(name, ['matched_by'])['matched_by'].value_counts()
        for k in matched:
            matched[k] += int(found.get(k, 0))
    renamed = table('annotation/feature_maps/GSE310929__GPL18573', ['feature_id', 'symbol', 'matched_by'])
    renamed = renamed[renamed['feature_id'].isin(['IL1RL2', 'IL1F10'])]
    methods_check(set(renamed['matched_by']) == {'previous_symbol'} and len(renamed) == 2, 'the renamed interleukin genes of GSE310929 have changed')
    # a count gene with an identical copy on the Y chromosome
    copy_study, copy_platform, copy_sample, copy_gene = 'GSE171741', 'GPL24676', 'GSM5233186', 'CSF2RA'
    cstem = methods_stem(copy_study, copy_platform)
    crows = table(f'annotation/gene_rows/{cstem}')
    crow = methods_one(crows[crows['symbol'] == copy_gene], f'{copy_gene} is no longer one gene of {copy_study}')
    cfeatures = crow['feature_ids'].split('|')
    methods_check(len(cfeatures) == 2 and crow['identical_rows'] == 1, f'{copy_gene} no longer has one identical copy')
    ccounts = table(f'expression/study_counts/{cstem}', ['feature_id', copy_sample]).set_index('feature_id')[copy_sample]
    clevel = table(f'annotation/gene_counts/{cstem}', ['symbol', copy_sample]).set_index('symbol')[copy_sample]
    methods_check(float(ccounts[cfeatures[0]]) == float(ccounts[cfeatures[1]]) == float(clevel[copy_gene]), f'the count of {copy_gene} is no longer one row counted once')
    crec = methods_one(summary[(summary['accession'] == copy_study) & (summary['platform'] == copy_platform)], f'{copy_study} is no longer one gene-level table')
    features_all = pick(receipt('annotation', 3), 'features', 'annotation 03.json')
    return {
        'statuses': {k: int(status.get(k, 0)) for k in names}, 'several': several, 'hla': hla,
        'genes': int(srow['genes']), 'genes_with_several_rows': int(srow['genes_with_several_rows']),
        'matched': matched, 'renamed': dict(zip(renamed['feature_id'], renamed['symbol'])),
        'copy': {'study': copy_study, 'sample': copy_sample, 'gene': copy_gene, 'rows': cfeatures, 'count': rnd(float(clevel[copy_gene]), 2),
                 'identical': int(crec['identical_rows'])},
        'features': {k: int(v) for k, v in features_all.items()},
        'tables': {'count': len(summary), 'fewest': int(summary['genes'].min()), 'most': int(summary['genes'].max()), 'median': rnd(summary['genes'].median(), 1)},
        'rules': pick(config('annotation'), 'gene_matrices.rules', 'annotation_config.json'),
    }


def methods_apply_rules(rules, fields):
    """The class that a rule list gives a record (fields: field -> recorded wording), and the 1-based number of the rule that
    decided; the same reading as the page's simulation (case-insensitive wordings and patterns)."""
    for number, rule in enumerate(rules, 1):
        wording = fields.get(rule['from'])
        if wording is None:
            continue
        if 'patterns' in rule:
            for pattern, value in rule['patterns']:
                if re.search(pattern, wording, flags=re.I):
                    return value, number
            continue
        value = rule['values'].get(wording.strip().lower())
        if value is None or value == 'next':
            continue
        return value, number
    return 'not_recorded', 0


def methods_labels():
    """Chapter 5: the five class rules of GSE65682, every record they decide, and the record GSM1602802."""
    key, study_id, platform, sample_id = 'GSE65682|GPL13667', 'GSE65682', 'GPL13667', 'GSM1602802'
    cohort_config = config('cohort')
    rules = pick(register('cohort'), f'rules.{key}', 'curation_register.json')
    order = rules['sample_class']
    methods_check([r['from'] for r in order] == ['title', 'abdominal sepsis and controls', 'endotype class', 'pneumonia diagnoses', 'title'],
                  f'the class rules of {key} are no longer those the page walks through')
    metadata = table('cohort/sample_metadata', ['study_accession', 'platform', 'sample_accession', 'field', 'value'])
    metadata = metadata[(metadata['study_accession'] == study_id) & (metadata['platform'] == platform)]
    wanted = sorted({r['from'] for r in order} | {'gender', 'age', 'source name', 'mortality event 28days'})
    fields = metadata[metadata['field'].isin(wanted)].pivot_table(index='sample_accession', columns='field', values='value', aggfunc='first')
    sources = table('cohort/label_sources', ['study_accession', 'platform', 'sample_accession', 'variable', 'value', 'basis', 'rule'])
    classes = sources[(sources['study_accession'] == study_id) & (sources['platform'] == platform) & (sources['variable'] == 'sample_class')].set_index('sample_accession')
    methods_check(set(classes.index) == set(fields.index), f'every record of {key} needs its fields and one class')
    records, wrong = [], []
    for accession in sorted(fields.index):
        f = {k: v for k, v in fields.loc[accession].items() if isinstance(v, str)}
        value, number = methods_apply_rules(order, f)
        if (value, number) != (classes.loc[accession, 'value'], int(classes.loc[accession, 'rule'])):
            wrong.append(accession)
        records.append([accession] + [f.get(name, '') for name in ('title', 'abdominal sepsis and controls', 'endotype class', 'pneumonia diagnoses')]
                       + [classes.loc[accession, 'value'], int(classes.loc[accession, 'rule'])])
    methods_check(not wrong, f'the rules of {key} read as on the page do not give the stored class of {wrong[:3]}')
    chosen = classes.loc[sample_id]
    methods_check(chosen['rule'] == 3 and chosen['value'] == 'sepsis' and chosen['basis'] == 'sample_record_and_study_description',
                  f'{sample_id} is no longer sepsis by the third rule on its record and the study sentence')
    record = {k: v for k, v in fields.loc[sample_id].items() if isinstance(v, str)}
    labels = table('cohort/cohort_samples', ['study_accession', 'sample_accession', 'sex', 'age_years', 'age_group', 'outcome', 'outcome_definition', 'septic_shock', 'patient_id', 'time_point'])
    own = methods_one(labels[(labels['study_accession'] == study_id) & (labels['sample_accession'] == sample_id)], f'{sample_id} needs one row')
    age_groups = pick(cohort_config, 'curation.age_groups', 'cohort_config.json')
    basis = pick(receipt('cohort', 3), 'class_basis', 'cohort 03.json')
    rule_sets = pick(receipt('cohort', 3), 'rule_sets', 'cohort 03.json')
    sentences = pick(quality('cohort'), 'summary.study_sentences', 'cohort qc.json')
    return {
        'study': study_id, 'platform': platform, 'sample': sample_id,
        'rules': [{'from': r['from'], 'patterns': r.get('patterns'), 'values': r.get('values'), 'study_states': r.get('study_states', '')} for r in order],
        'records': {'__f': ['accession', 'title', 'abdominal', 'endotype', 'pneumonia', 'value', 'rule'], '__r': records},
        'record': record,
        'labels': {'sex': own['sex'], 'age': rnd(own['age_years'], 2), 'age_group': own['age_group'], 'outcome': own['outcome'],
                   'outcome_definition': own['outcome_definition'], 'septic_shock': own['septic_shock'], 'patient_id': own['patient_id'],
                   'outcome_rule': rules['outcome']},
        'age_groups': age_groups, 'basis': basis, 'rule_sets': rule_sets, 'study_sentences': sentences,
        'classes_resource': pick(receipt('cohort', 3), 'sample_class', 'cohort 03.json'),
        'infection_not_stated': pick(receipt('cohort', 3), 'infection_sepsis_not_stated', 'cohort 03.json'),
    }



def methods_repeats():
    """Chapter 6: the array of patient 11-0008 deposited in six studies, the cohort group it belongs to, the patients of
    GSE54514 with one sample marked for each, and a patient who joins two cohort groups into one patient-safe group."""
    sample_id, group_id = 'GSE4607:GSM102975', 'G024'
    expected = {'GSE4607', 'GSE8121', 'GSE9692', 'GSE13904', 'GSE26440', 'GSE66099'}
    samples = table('cohort/cohort_samples')
    records = samples[samples['sample_id'] == sample_id]
    methods_check(set(records['study_accession']) == expected and (records['linked_by'] == 'chip_type_scan_time_and_median_intensity').all(),
                  f'the records of {sample_id} are no longer the six that the scan time links')
    arrays = table('expression/affymetrix_arrays', ['study_accession', 'sample_accession', 'chip_type', 'member', 'scan_time', 'raw_pm_median'])
    found = arrays.merge(records[['study_accession', 'sample_accession']], on=['study_accession', 'sample_accession'], validate='one_to_one')
    methods_check(len(found) == len(records), f'every record of {sample_id} needs its CEL file')
    summary = table('cohort/cohort_summary', ['cohort_group', 'patient_safe_group', 'study_accession', 'platform', 'records', 'primary_records',
                                              'release_date'])
    released = summary.drop_duplicates('study_accession').set_index('study_accession')['release_date'].astype(str)
    sources = table('cohort/label_sources', ['study_accession', 'sample_accession', 'variable', 'recorded', 'study_statement'])
    sources = sources[sources['variable'] == 'sample_class'].set_index(['study_accession', 'sample_accession'])
    files = found.set_index('study_accession')
    deposits = []
    for r in records.sort_values('study_accession', key=lambda s: s.map(released)).itertuples(index=False):
        f = files.loc[r.study_accession]
        src = sources.loc[(r.study_accession, r.sample_accession)]
        deposits.append({'study': r.study_accession, 'record': r.sample_accession, 'released': released[r.study_accession][:10],
                         'file': f['member'], 'chip': f['chip_type'], 'scan_time': str(f['scan_time'])[:19], 'median': rnd(f['raw_pm_median'], 2),
                         'wording': r.authors_label, 'class': r.sample_class, 'primary': bool(r.is_primary_record), 'patient_id': r.patient_id,
                         'recorded': src['recorded'], 'statement': src['study_statement'] or ''})
    group = summary[summary['cohort_group'] == group_id]
    methods_check(set(group['study_accession']) == expected | {'GSE26378'}, f'the cohort group {group_id} has changed')
    group_rows = [{'study': r.study_accession, 'platform': r.platform, 'records': int(r.records), 'primary': int(r.primary_records),
                   'released': str(r.release_date)[:10]} for r in group.sort_values('release_date').itertuples(index=False)]
    # the patients of GSE54514: the five samples of NS_19, and the one sample marked for every patient
    part = samples[samples['study_accession'] == 'GSE54514']
    methods_check(part['patient_identified'].all() and int(part['one_per_patient'].sum()) == part['patient_key'].nunique(),
                  'GSE54514 must mark one sample for every named patient')
    own = part[part['patient_key'] == 'GSE54514:NS_19'].sort_values('time_order')
    methods_check(own['one_per_patient'].tolist() == [True] + [False] * (len(own) - 1), 'NS_19 is no longer marked on its first day')
    by_class = part.groupby('sample_class')['patient_key'].nunique()
    sizes = part.groupby('patient_key').size()
    patient = {'study': 'GSE54514', 'patient': own['patient_id'].iloc[0], 'records': len(part), 'patients': int(part['patient_key'].nunique()),
               'sepsis': int(by_class.get('sepsis', 0)), 'healthy': int(by_class.get('healthy_control', 0)),
               'most': int(sizes.max()), 'with_most': int((sizes == sizes.max()).sum()),
               'healthy_points': sorted(part.loc[part['sample_class'] == 'healthy_control', 'time_point'].unique().tolist()),
               'outcome': own['outcome'].iloc[0],
               'samples': [{'record': r.sample_accession, 'time_point': r.time_point, 'time_order': rnd(r.time_order, 2), 'unit': r.time_unit,
                            'flag': bool(r.any_flag), 'tissue': r.tissue, 'marked': bool(r.one_per_patient)}
                           for r in own.itertuples(index=False)]}
    # one patient, two cohort groups
    safe = samples[samples['patient_key'] == 'GSE13015:I001']
    methods_check(set(safe['cohort_group']) == {'G011', 'G019'} and set(safe['patient_safe_group']) == {'G011+G019'}, 'patient I001 of GSE13015 has changed')
    known = samples[samples['patient_identified']]
    spans = known.groupby('patient_key')['cohort_group'].agg(lambda v: tuple(sorted(set(v))))
    safe_rows = [{'study': r.study_accession, 'platform': r.platform, 'record': r.sample_accession, 'group': r.cohort_group, 'sample_id': r.sample_id,
                  'patient_id': r.patient_id, 'primary': bool(r.is_primary_record)}
                 for r in safe.sort_values(['cohort_group', 'study_accession']).itertuples(index=False)]
    links, cohort5 = receipt('cohort', 2), receipt('cohort', 5)
    patients = pick(cohort5, 'patients', 'cohort 05.json')
    pairs = table('validation/repeat_pairs', ['linked_by'])['linked_by'].value_counts()
    return {
        'sample_id': sample_id, 'deposits': deposits, 'group': group_id, 'group_rows': group_rows, 'patient': patient,
        'safe': {'rows': safe_rows, 'patients': int((spans == ('G011', 'G019')).sum())},
        'totals': {'samples': pick(links, 'samples', 'cohort 02.json'), 'cohorts': pick(cohort5, 'cohorts', 'cohort 05.json'),
                   'identified': patients['identified'], 'with_patient': patients['distinct_samples_with_an_identified_patient'],
                   'without_patient': patients['distinct_samples_without_a_patient_identifier'], 'one_per_patient': patients['one_per_patient']['total'],
                   'in_several_studies': patients['in_several_studies'], 'namespaces': len(patients.get('namespaces', [])),
                   'class_differs': pick(cohort5, 'distinct_samples.whose_records_differ_in_class', 'cohort 05.json')},
        'evidence': [[k, int(v)] for k, v in pairs.items()],
    }


def methods_sex():
    """Chapter 7: the score d (mean of the Y genes minus XIST) of every sample of GSE65682, the two-means split of its
    table and the comparison with the recorded sex."""
    study_id, platform = 'GSE65682', 'GPL13667'
    matrices = table('validation/sex_check_matrices')
    row = methods_one(matrices[(matrices['study_accession'] == study_id) & (matrices['platform'] == platform)], f'{study_id} needs one row of the sex check')
    methods_check(row['status'] == 'tested' and row['female_gene_read'] and not row['recorded_sex_unreliable'],
                  f'the sex check of {study_id} is no longer one split of the whole table')
    policy = pick(config('validation'), 'sex', 'validation_config.json')
    samples = table('validation/sex_check_samples')
    own = samples[(samples['study_accession'] == study_id) & (samples['platform'] == platform)].sort_values('sample_accession')
    methods_check(len(own) == row['samples'], f'the samples of the sex check of {study_id} differ from its table row')
    code = {'female': 'f', 'male': 'm', 'not_recorded': 'n', 'ambiguous': 'a', 'not_called': 'x'}
    each = [[r.sample_accession, rnd(r.male_genes_value, 4), rnd(r.female_gene_value, 4), rnd(r.difference, 4), code[r.recorded_sex], code[r.expression_sex]]
            for r in own.itertuples(index=False)]
    v2, vqc = receipt('validation', 2), quality('validation')['summary']
    stats = {}
    for k in ('samples', 'lower_group', 'higher_group', 'gap', 'separation', 'midpoint', 'called_female', 'called_male', 'ambiguous',
              'concordant', 'discordant', 'concordance'):
        stats[k] = rnd(row[k], 6) if k in ('gap', 'separation', 'midpoint', 'concordance') else int(row[k])
    return {
        'study': study_id, 'platform': platform, 'y_genes': policy['male_genes'], 'female_gene': policy['female_gene'],
        'policy': {k: policy[k] for k in ('minimum_separation', 'minimum_group', 'minimum_group_share', 'component_separation', 'ambiguous_zone',
                                          'minimum_compared', 'unreliable_below', 'minimum_samples_per_tissue')},
        'each': {'__f': ['accession', 'y', 'xist', 'd', 'recorded', 'call'], '__r': each}, 'table': stats,
        'examples': ['GSM1602801', 'GSM1602802'],
        'resource': {'matrices': pick(v2, 'matrices.total', 'validation 02.json'), 'two_sexes': pick(vqc, 'matrices_showing_two_sexes', 'validation qc.json'),
                     'samples': {k: v for k, v in pick(v2, 'samples', 'validation 02.json').items()},
                     'unreliable': pick(v2, 'unreliable_matrices', 'validation 02.json')},
    }


def methods_geometric(values):
    return math.exp(sum(math.log(v) for v in values) / len(values))


def methods_score():
    """Chapter 8: the Sepsis MetaScore of GSM1602802 from its eleven genes, and the scores of every stratum, from which the
    page counts the AUROC as the share of sepsis-control pairs in the right order."""
    study_id, platform, sample_id = 'GSE65682', 'GPL13667', 'GSM1602802'
    genes = table('validation/validation_genes')
    named = genes[genes['test'] == 'sepsis_score']
    strata = table('validation/sepsis_score_by_study')
    stratum = methods_one(strata[(strata['study_accession'] == study_id) & (strata['platform'] == platform)], f'{study_id} needs one stratum')
    methods_check(stratum['score_shift'] == 0 and stratum['score_genes_higher'] == 6 and stratum['score_genes_lower'] == 5,
                  f'the stratum of {study_id} is no longer unshifted with every gene of the score')
    matrix = table(f'annotation/gene_matrices/{methods_stem(study_id, platform)}', ['hgnc_id', sample_id]).set_index('hgnc_id')[sample_id]
    parts, means = {}, {}
    for role, short in (('higher_in_sepsis', 'up'), ('lower_in_sepsis', 'down')):
        own = named[named['role'] == role]
        values = [float(matrix[g]) for g in own['hgnc_id']]
        parts[short] = [{'published': p, 'symbol': s, 'value': rnd(v, 4)} for p, s, v in zip(own['published_symbol'], own['symbol'], values)]
        means[short] = methods_geometric(values)
    scores = table('validation/sepsis_score_samples', ['study_accession', 'platform', 'tissue', 'sample_accession', 'sample_class', 'score'])
    sample = methods_one(scores[(scores['study_accession'] == study_id) & (scores['sample_accession'] == sample_id)], f'{sample_id} needs one score')
    methods_check(abs(sample['score'] - (means['up'] - means['down'])) < 1e-9, f'the score of {sample_id} is not the difference of its geometric means')
    # every stratum: its sepsis and healthy-control scores; the pairs in the right order must give its stored AUROC
    out = []
    for r in strata[strata['score_computed']].sort_values(['study_accession', 'platform', 'tissue']).itertuples(index=False):
        own = scores[(scores['study_accession'] == r.study_accession) & (scores['platform'] == r.platform) & (scores['tissue'] == r.tissue)]
        pos = own.loc[own['sample_class'] == 'sepsis', 'score'].round(4).tolist()
        neg = own.loc[own['sample_class'] == 'healthy_control', 'score'].round(4).tolist()
        right = sum((p > q) + 0.5 * (p == q) for p in pos for q in neg)
        methods_check(len(pos) == r.sepsis and len(neg) == r.healthy_control and abs(right / (len(pos) * len(neg)) - r.auroc) < 1e-9,
                      f'the rounded scores of {r.study_accession} {r.platform} do not give its AUROC')
        out.append({'study': r.study_accession, 'platform': r.platform, 'tissue': r.tissue, 'auroc': rnd(r.auroc, 6),
                    'low': rnd(r.auroc_low, 4), 'high': rnd(r.auroc_high, 4), 'derivation': r.score_derivation or '',
                    'sepsis': sorted(pos, reverse=True), 'healthy': sorted(neg, reverse=True)})
    biology = pick(config('validation'), 'biology', 'validation_config.json')
    v4 = receipt('validation', 4)
    return {
        'study': study_id, 'platform': platform, 'sample': sample_id, 'up': parts['up'], 'down': parts['down'],
        'score': rnd(sample['score'], 6), 'minimum_genes': biology['minimum_score_genes'], 'minimum_group': biology['minimum_group'],
        'bootstrap': biology['bootstrap'], 'strata': out,
        'resource': {'strata': pick(v4, 'strata.total', 'validation 04.json'), 'score': {k: v for k, v in pick(v4, 'score', 'validation 04.json').items()
                                                                                         if not isinstance(v, (dict, list))},
                     'independent': pick(v4, 'score.independent', 'validation 04.json'),
                     'genes': {k: pick(v4, f'genes.{k}', 'validation 04.json') for k in ('comparisons', 'as_published', 'share_as_published')}},
    }


def methods_shrunk_variance(estimate, freedom):
    """The inverse-gamma prior fitted across the genes by moments, as in ComBat, and every gene's variance drawn toward it."""
    usable = estimate.notna() & (freedom >= 1)
    mean, spread = float(estimate[usable].mean()), float(estimate[usable].var(ddof=1))
    shape, rate = (2 * spread + mean ** 2) / spread, (mean * spread + mean ** 3) / spread
    return ((0.5 * freedom * estimate + rate) / (0.5 * freedom + shape - 1)).where(usable), shape, rate


def methods_shrunk_location(estimate, noise):
    usable = estimate.notna() & noise.notna()
    mean, spread = float(estimate[usable].mean()), float(estimate[usable].var(ddof=1))
    return ((spread * estimate + noise * mean) / (spread + noise)).where(usable), mean, spread


def methods_histogram(values, lo, hi, bins):
    """Counts of values in equal bins from lo to hi; values outside fall in the end bins."""
    counts = [0] * bins
    for v in values:
        counts[min(bins - 1, max(0, int((v - lo) / (hi - lo) * bins)))] += 1
    return {'lo': lo, 'hi': hi, 'counts': counts}


def methods_correction():
    """Chapter 9: CEACAM1 through the correction of GSE54514 on its healthy controls (alignment, sex effect, anchor,
    empirical Bayes location and scale, corrected value) and of E-MTAB-4421 by matched groups."""
    group, gene = 'bulk_blood', 'CEACAM1'
    study_id, platform, sepsis_id, control_id = 'GSE54514', 'GPL6947', 'GSM1317914', 'GSM1317896'
    genes = table('correction/correction_genes', ['group', 'hgnc_id', 'in_group_set'])
    order = sorted(genes.loc[(genes['group'] == group) & genes['in_group_set'], 'hgnc_id'], key=lambda g: int(g.split(':')[1]))
    model = table('correction/correction_model')
    model = model[model['group'] == group].set_index('hgnc_id').reindex(order)
    hgnc = methods_one(model[model['symbol'] == gene], f'{gene} is not a gene of the {group} model').name
    m, s, beta = (model[n] for n in ('anchor_mean', 'anchor_sd', 'sex_effect'))
    effects = table('correction/anchored_effects')
    effects = effects[(effects['group'] == group) & (effects['study_accession'] == study_id) & (effects['platform'] == platform)].set_index('hgnc_id').reindex(order)
    records = table('correction/correction_records')
    unit = records[(records['group'] == group) & (records['study_accession'] == study_id) & (records['platform'] == platform) & (records['status'] == 'in_unit')]
    unit = unit.sort_values('sample_accession', kind='stable').reset_index(drop=True)
    estimation = unit[unit['estimation']]
    controls = estimation[estimation['sample_class'] == 'healthy_control']
    reference = table('correction/alignment_reference')
    reference = reference[reference['group'] == group].sort_values('position').reset_index(drop=True)
    stem = methods_stem(study_id, platform)
    aligned = table(f'correction/aligned/{group}/{stem}').set_index('hgnc_id').reindex(order)
    matrix = table(f'annotation/gene_matrices/{stem}', ['hgnc_id', sepsis_id, control_id]).set_index('hgnc_id').reindex(order)
    corrected = table('correction/corrected/bulk_blood_expression', ['hgnc_id', *unit['sample_id']]).set_index('hgnc_id')

    def interp(x):
        p, v = reference['position'].tolist(), reference['value'].tolist()
        if x <= p[0]:
            return v[0]
        if x >= p[-1]:
            return v[-1]
        lo, hi = 0, len(p) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            lo, hi = (mid, hi) if p[mid] <= x else (lo, mid)
        return v[lo] + (v[hi] - v[lo]) * (x - p[lo]) / (p[hi] - p[lo])
    examples = []
    for role, accession in (('sepsis', sepsis_id), ('control', control_id)):
        record = methods_one(unit[unit['sample_accession'] == accession], f'{accession} is no longer a record of the unit')
        methods_check(record['female'] == 1 and record['estimation'], f'{accession} is no longer a female estimation sample')
        raw = matrix[accession].dropna()
        ranks = raw.rank(method='average')
        rank, held = float(ranks[hgnc]), int(len(raw))
        position = (rank - 0.5) / held
        value = float(aligned.loc[hgnc, accession])
        methods_check(abs(interp(position) - value) < 1e-4, f'the aligned value of {accession} is not the reference at its position')
        z = (value - beta[hgnc] - m[hgnc]) / s[hgnc]
        result = m[hgnc] + s[hgnc] * (z - effects.loc[hgnc, 'location']) / effects.loc[hgnc, 'scale'] + beta[hgnc]
        stored = float(corrected.loc[hgnc, record['sample_id']])
        methods_check(abs(result - stored) < 1e-3, f'the worked example does not give the corrected value of {accession}')
        examples.append({'role': role, 'sample': accession, 'class': record['sample_class'], 'female': float(record['female']), 'raw': rnd(raw[hgnc], 4),
                         'rank': rank, 'genes': held, 'position': rnd(position, 6), 'aligned': rnd(value, 4), 'corrected': rnd(stored, 4)})
    # the estimates before shrinkage, from the healthy controls in anchor units, and the shrinkage that gives the stored effects
    female = controls.set_index('sample_accession')['female']
    values = aligned[list(female.index)] - pd.DataFrame({a: beta * female[a] for a in female.index})
    count = values.notna().sum(axis=1).astype(float)
    standard = values.sub(m, axis=0).div(s, axis=0)
    own_location, own_variance = standard.mean(axis=1), standard.var(axis=1, ddof=1)
    variance, shape, rate = methods_shrunk_variance(own_variance, count - 1)
    noise = variance / count.clip(lower=1)
    pulled, prior_mean, prior_spread = methods_shrunk_location(own_location.where(count >= 2), noise)
    methods_check(abs(math.sqrt(variance[hgnc]) - effects.loc[hgnc, 'scale']) < 1e-6 and abs(pulled[hgnc] - effects.loc[hgnc, 'location']) < 1e-6,
                  'the shrinkage of the worked example does not give the stored effects')
    n = int(count[hgnc])
    anchored = {
        'study': study_id, 'platform': platform, 'gene': gene, 'hgnc': hgnc, 'records': len(unit), 'estimation': len(estimation),
        'controls': len(controls), 'sepsis': int((estimation['sample_class'] == 'sepsis').sum()),
        'anchor_mean': rnd(m[hgnc], 6), 'anchor_sd': rnd(s[hgnc], 6), 'sex_effect': rnd(beta[hgnc], 6), 'anchor_units': int(model.loc[hgnc, 'anchor_units']),
        'location': rnd(effects.loc[hgnc, 'location'], 6), 'scale': rnd(effects.loc[hgnc, 'scale'], 6), 'examples': examples,
        'own_location': rnd(own_location[hgnc], 6), 'own_variance': rnd(own_variance[hgnc], 6), 'n': n,
        'prior_shape': rnd(shape, 6), 'prior_rate': rnd(rate, 6), 'prior_location': rnd(prior_mean, 6), 'prior_spread': rnd(prior_spread, 6),
        'genes': len(order),
        # how the unit's own estimates spread over the genes, from which the priors are fitted
        'location_histogram': methods_histogram(own_location.dropna().tolist(), -8, 8, 64),
        'variance_histogram': methods_histogram(own_variance.dropna().tolist(), 0, 4, 64),
    }
    ref = {'positions': len(reference), 'value': [rnd(v, 4) for v in reference['value']],
           'first': rnd(reference['position'].iloc[0], 8), 'last': rnd(reference['position'].iloc[-1], 8)}
    step = (reference['position'].iloc[-1] - reference['position'].iloc[0]) / (len(reference) - 1)
    methods_check(((reference['position'].diff().dropna() - step).abs() < 1e-9).all(), 'the reference positions are no longer evenly spaced')

    # matched groups: CEACAM1 in E-MTAB-4421, adults with sepsis only
    mstudy, mplatform = 'E-MTAB-4421', 'A-MEXP-2210'
    munit = records[(records['group'] == group) & (records['study_accession'] == mstudy) & (records['platform'] == mplatform) & (records['status'] == 'in_unit')]
    mown = munit[munit['estimation']].sort_values('sample_accession', kind='stable').reset_index(drop=True)
    cells = mown.groupby(['sample_class', 'age_group']).size()
    methods_check(list(cells.index) == [('sepsis', 'adult')], 'the matched-group example expects adult sepsis samples only')
    mfemale = mown.set_index('sample_accession')['female']
    maligned = table(f'correction/aligned/{group}/{methods_stem(mstudy, mplatform)}', ['hgnc_id', *mfemale.index]).set_index('hgnc_id').reindex(order)
    mvalues = maligned[list(mfemale.index)] - pd.DataFrame({a: beta * mfemale[a] for a in mfemale.index})
    mcount = mvalues.notna().sum(axis=1).astype(float)
    mmean = mvalues.mean(axis=1)
    squares = mvalues.sub(mmean, axis=0).pow(2).sum(axis=1)
    refs = table('correction/class_references')
    refs = refs[(refs['group'] == group) & (refs['sample_class'] == 'sepsis') & (refs['age_group'] == 'adult')].set_index('hgnc_id').reindex(order)
    rho, omega = refs['mean'], refs['variance']
    floor = pick(config('correction'), 'anchored.floor_fraction', 'correction_config.json') ** 2 * float(omega.median())
    omega_used = omega.clip(lower=floor)
    z_bar = (mmean - m) / s
    estimate = (squares / s ** 2 / omega_used / (mcount - 1)).where(mcount > 1)
    mvariance, _, _ = methods_shrunk_variance(estimate, mcount - 1)
    mscale = mvariance.pow(0.5)
    raw_location = z_bar - mscale * rho
    mlocation, _, _ = methods_shrunk_location(raw_location, mvariance * omega_used / mcount)
    stored = table('correction/candidate_effects')
    stored = stored[(stored['group'] == group) & (stored['study_accession'] == mstudy) & (stored['platform'] == mplatform)].set_index('hgnc_id')
    methods_check(abs(stored.loc[hgnc, 'location'] - mlocation[hgnc]) < 1e-6 and abs(stored.loc[hgnc, 'scale'] - mscale[hgnc]) < 1e-6,
                  'the matched-group example does not give the stored effects')
    matched = {'study': mstudy, 'platform': mplatform, 'samples': int(mcount[hgnc]), 'reference_units': int(refs.loc[hgnc, 'units']),
               'mean': rnd(mmean[hgnc], 6), 'spread': rnd(squares[hgnc] / s[hgnc] ** 2 / (mcount[hgnc] - 1), 6), 'omega': rnd(omega_used[hgnc], 6),
               'rho': rnd(rho[hgnc], 6), 'variance_raw': rnd(estimate[hgnc], 6), 'scale': rnd(mscale[hgnc], 6), 'location_raw': rnd(raw_location[hgnc], 6),
               'location': rnd(mlocation[hgnc], 6)}
    r1, r4, r5, r6, r7 = (receipt('correction', i) for i in (1, 4, 5, 6, 7))
    sex = {g: pick(r4, f'groups.{g}.largest_sex_effects', 'correction 04.json') for g in ('bulk_blood', 'pbmc')}
    cfg = config('correction')
    return {
        'anchored': anchored, 'reference': ref, 'matched': matched,
        'sex_effects': {g: [[k, rnd(v, 3)] for k, v in sex[g].items()] for g in sex},
        'sex_over_one': {g: pick(r4, f'groups.{g}.genes_with_sex_effect_over_one', 'correction 04.json') for g in ('bulk_blood', 'pbmc')},
        'limits': {'minimum_unit_samples': pick(cfg, 'inputs.minimum_unit_samples', 'correction_config.json'),
                   'minimum_healthy_controls': pick(cfg, 'inputs.minimum_healthy_controls', 'correction_config.json'),
                   'minimum_gene_share': pick(cfg, 'inputs.minimum_gene_share', 'correction_config.json'),
                   'grid_points': pick(cfg, 'alignment.grid_points', 'correction_config.json')},
        'units': {g: {k: pick(r1, f'groups.{g}.{k}', 'correction 01.json') for k in ('units', 'anchored_units', 'newborn_units', 'estimation_samples', 'genes')}
                  for g in ('bulk_blood', 'pbmc')},
        'matched_units': {k: pick(r5, f'groups.bulk_blood.{k}', 'correction 05.json') for k in ('unanchored_units', 'corrected_by_matched_groups')},
        'kept': {g: pick(r6, f'groups.{g}.by_method', 'correction 06.json') for g in ('bulk_blood', 'pbmc')},
        'tables': {g: {k: pick(r7, f'groups.{g}.{k}', 'correction 07.json') for k in ('samples', 'studies', 'units', 'genes', 'by_method')}
                   for g in ('bulk_blood', 'pbmc')},
    }


def methods_heldout():
    """Chapter 10: the held-out test of every anchored unit by method (the points of Figure 6), with GSE54514 as the
    worked example."""
    rows = figure_data('Fig6', 'held_out', ['group', 'study_accession', 'platform', 'family', 'method', 'hidden_healthy', 'estimation_samples',
                                            'location_error', 'pool_dispersion', 'sampling_floor'])
    test = table('correction/held_out_test')
    methods_check(len(rows) == len(test), 'Fig6 held_out and the held-out table differ in their rows; rerun stage 08')
    methods = ['none', 'plain_combat', 'flat_genes', 'matched_groups']
    methods_check(set(rows['method']) <= set(methods), 'the held-out test names another method')
    units = []
    for (g, st, pl), part in rows.groupby(['group', 'study_accession', 'platform'], sort=True):
        by = part.set_index('method')
        first = part.iloc[0]
        units.append({'group': g, 'study': st, 'platform': pl, 'family': first['family'], 'hidden': int(first['hidden_healthy']),
                      'visible': int(first['estimation_samples']), 'floor': rnd(first['sampling_floor'], 4),
                      'error': [rnd(by.loc[k, 'location_error'], 4) if k in by.index else None for k in methods],
                      'dispersion': [rnd(by.loc[k, 'pool_dispersion'], 4) if k in by.index else None for k in methods]})
    medians = {}
    for g, part in test.groupby('group'):
        medians[g] = {'floor': rnd(part.drop_duplicates(['study_accession', 'platform'])['sampling_floor'].median(), 4),
                      'units': int(part[['study_accession', 'platform']].drop_duplicates().shape[0]),
                      'studies': int(part['study_accession'].nunique())}
        for k in methods:
            own = part[part['method'] == k]
            if len(own):
                medians[g][k] = {'error': rnd(own['location_error'].median(), 4), 'dispersion': rnd(own['pool_dispersion'].median(), 4), 'units': len(own)}
    r5 = receipt('correction', 5)
    paired = [p for p in pick(r5, 'groups.bulk_blood.paired_comparison', 'correction 05.json') if p['method'] == 'plain_combat' and p['baseline'] == 'matched_groups']
    return {'methods': methods, 'units': units, 'medians': medians, 'example': 'GSE54514',
            'flat_genes': pick(config('correction'), 'flat_genes.count', 'correction_config.json'),
            'level': pick(config('correction'), 'unanchored.comparison.level', 'correction_config.json'),
            'paired': {k: (rnd(v, 4) if isinstance(v, float) else v) for k, v in paired[0].items()} if paired else None}


def methods_bench():
    """Chapter 11: the three tasks, the reasons a record leaves a task, the split schemes, and the fold that holds out the
    patient-safe group G024 in the task of sepsis against healthy controls."""
    task, scheme, fold = 'sepsis_vs_healthy_control', 'leave_one_group_out', 'G024'
    bench = receipt('cohort', 6)
    sheet = review('cohort', 'benchmark')
    rows = sheet.set_index(['task', 'scheme', 'fold'])
    row = rows.loc[(task, scheme, fold)]
    methods_check(int(row['test_groups']) == 1 and int(row['held_out_not_scored_groups']) == 0, 'the fold of G024 no longer holds out one scored group')
    contrast = rows.loc[(task, scheme, 'G080')]
    transfer = rows.loc[(task, 'technology_transfer', 'array_to_rna_seq')]
    tasks = table('cohort/benchmark_tasks', ['study_accession', 'platform', 'task', 'patient_safe_group', 'eligible'])
    held = tasks[(tasks['task'] == task) & tasks['eligible'] & (tasks['patient_safe_group'] == fold)]
    methods_check(len(held) == int(row['test_positives']) + int(row['test_negatives']), 'the eligible samples of G024 are not those of its fold')
    by_platform = held.groupby(['study_accession', 'platform']).size()
    splits = table('cohort/benchmark_splits', ['task', 'scheme', 'fold', 'patient_safe_group', 'first_release', 'positives', 'negatives'])
    temporal = splits[(splits['task'] == task) & (splits['scheme'] == 'temporal')]
    dated = methods_one(temporal[temporal['patient_safe_group'] == fold], 'G024 needs one row of the temporal fold')
    summary = table('cohort/cohort_summary', ['patient_safe_group', 'study_accession'])
    return {
        'policy': {k: v for k, v in pick(bench, 'policy', 'cohort 06.json').items() if k != 'technologies'},
        'reasons': pick(bench, 'reasons', 'cohort 06.json'), 'schemes': pick(bench, 'schemes', 'cohort 06.json'),
        'tasks': [{'task': k, 'eligible': v['eligible'], 'positives': v['positives'], 'negatives': v['negatives'], 'groups': v['patient_safe_groups'],
                   'scored_groups': v['groups_with_the_minimum_of_each_class'], 'folds': v['folds'], 'excluded': v['excluded'], 'records': v['records']}
                  for k, v in pick(bench, 'tasks', 'cohort 06.json').items()],
        'fold': {'task': task, 'scheme': scheme, 'group': fold, 'train_groups': int(row['train_groups']), 'train_positives': int(row['train_positives']),
                 'train_negatives': int(row['train_negatives']), 'test_positives': int(row['test_positives']), 'test_negatives': int(row['test_negatives']),
                 'platforms': [[s, p, int(n)] for (s, p), n in by_platform.items()], 'first_release': str(dated['first_release'])[:10]},
        'contrast': {'group': 'G080', 'studies': sorted(summary.loc[summary['patient_safe_group'] == 'G080', 'study_accession'].unique()),
                     'positives': int(contrast['test_positives']), 'negatives': int(contrast['test_negatives'])},
        'transfer': {'test_groups': int(transfer['test_groups']), 'not_scored': int(transfer['held_out_not_scored_groups'])},
        'sepsis_only_groups': int((temporal['negatives'] == 0).sum()), 'control_only_groups': int((temporal['positives'] == 0).sum()),
    }


def methods_overview():
    """Chapter 0: the counts of the path from the searches to the corrected tables (the flow of Figure 1)."""
    flow = figure_data('Fig1', 'flow', ['step', 'count'])
    return {'flow': {r.step: int(r.count) for r in flow.itertuples(index=False)}}


def methods_vis():
    """The labels and colors of the mapping statuses, as the figures use them."""
    vis = config('visualization')
    return {'mapping_labels': pick(vis, 'mapping_labels', 'visualization_config.json'), 'mapping_colors': pick(vis, 'mapping_colors', 'visualization_config.json')}


def methods():
    """Method in detail: the inputs of every simulation and the real numbers of every worked example, one block per chapter."""
    return {'methods': {
        'overview': methods_overview(), 'search': methods_search(), 'arrays': methods_arrays(), 'counts': methods_counts(),
        'genes': methods_genes(), 'labels': methods_labels(), 'repeats': methods_repeats(), 'sex': methods_sex(), 'score': methods_score(),
        'correction': methods_correction(), 'heldout': methods_heldout(), 'bench': methods_bench(),
        'vis': methods_vis(), 'glossary': [list(entry) for entry in METHODS_GLOSSARY],
    }}


PARTS = {'core': core, 'story': story, 'resource': resource, 'studies': studies, 'methods': methods}


def main():
    parser = argparse.ArgumentParser(description='Write the data of the results page from the pipeline outputs.')
    # The part names are checked here, not with choices=: argparse before Python 3.12 rejects an empty list against choices.
    parser.add_argument('parts', nargs='*', metavar='part',
                        help=f'parts to rebuild ({", ".join(PARTS)}); all of them, then the script versions, if none is named')
    parser.add_argument('--stamp', action='store_true', help='stamp the scripts of index.html with content versions')
    args = parser.parse_args()
    unknown = [name for name in args.parts if name not in PARTS]
    if unknown:
        parser.error(f'unknown part {", ".join(unknown)}; choose from {", ".join(PARTS)}')
    for name in args.parts or PARTS:
        write_part(name, PARTS[name]())
    if args.stamp or not args.parts:
        stamp()


if __name__ == '__main__':
    main()
