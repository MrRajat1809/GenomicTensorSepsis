#!/usr/bin/env bash
# Build the resource in dependency order: acquisition, expression, annotation, cohort, validation, correction, figures, then the export of the data package.
# -E lets the ERR trap fire inside functions too, so a failed step is always reported.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="${PYTHON:-python}"
FROM="01:01"
OUTDIR="data_release"
WORKERS=2
UPDATE=0
EXPORT_ONLY=0
NO_EXPORT=0
DRY_RUN=0
CHECK_ONLY=0

usage() {
  cat <<'HELP'
Usage: bash pipeline/01_build_dataset.sh [options]

Run inside the project image, from any directory:
  docker build -t genomictensorsepsis:1.0.0 docker
  docker run --rm -v "$PWD":/workspace -w /workspace genomictensorsepsis:1.0.0 \
      bash pipeline/01_build_dataset.sh --check-only
  docker run --rm -v "$PWD":/workspace -w /workspace genomictensorsepsis:1.0.0 \
      bash pipeline/01_build_dataset.sh

The build reads the stored answers of the repositories under data/raw/: the census up to search.census_date
of the acquisition configuration, the study and sample records, the platform records, the annotation
sources and the download lock, which snapshots.tar.gz of the data package holds. It rebuilds the same
resource; missing expression files are downloaded as the lock lists them and accepted only at their SHA-256.

  --check-only        Check the environment, the configurations, the stored answers and the output paths; write nothing
  --update-search     Ask the repositories again instead of reading the stored answers, for studies released up to
                      search.census_date; a new candidate stops the build until screening_register.json decides it.
                      Move search.census_date only to extend the resource.
  --from STAGE[:NN]   Resume at a stage or at one of its scripts, e.g. 02 or 03:04 (inclusive);
                      the script must exist. A failed build prints the value to resume with.
  --export-only       Package the completed, QC-passed stages and rebuild the results page data; compute nothing
  --no-export         Run the stages and their quality control without creating the data package or the page data
  --outdir PATH       Directory of the data package; managed files are replaced (default: data_release)
  --workers N         Processes that write the per-study tables of the package (default: 2)
  --dry-run           Print the commands without running or writing anything
  -h, --help          Show this help
Order: 01 data acquisition (01-05), 02 expression processing (01-09), 03 gene annotation (01-05),
       04 cohort curation (01-07), 05 technical validation (01-05), 06 batch correction (01-08), 07 visualization (01-07),
       then the export of the data package and of the results page data (release/build_data.py).
A copy without release/ runs every stage and skips the page data.
Each stage ends with its master quality control, which must pass before the next stage starts.
Exactly these scripts run: a stage folder whose numbered scripts are not 01 to its last number, each once, stops the build.
--from runs the remaining sequence; the preceding stages must already be complete and match the
current code and configuration.

Outputs: data/raw/ (downloaded files with SHA256SUMS.txt), data/processed/{acquisition,expression,
annotation,cohort,validation,correction}/, outputs/<stage>/ with receipts and QC reports, outputs/visualization/ with the figures, and
data_release/ with the data package for GigaDB: nine tar archives of CSV.gz tables, data_dictionary.csv, cohort_flow.csv, readme.txt,
md5sum.txt and release_manifest.json; release/data.js and release/data/ for the results page, when release/ is present.
Execution logs append to outputs/logs/build_dataset.log.

A build downloads about 44 GB of expression files where they are missing and takes many hours; the
expression stage needs R with Bioconductor, which the image provides. With --update-search the searches
are repeated against the live repositories, and the records they return can differ from the stored ones.
HELP
}

need_value() { [[ $# -ge 2 && -n "$2" && "$2" != --* ]] || { echo "Missing value for $1" >&2; exit 2; }; }
while (($#)); do
  case "$1" in
    --outdir) need_value "$@"; OUTDIR="$2"; shift 2 ;;
    --from) need_value "$@"; FROM="$2"; shift 2 ;;
    --workers) need_value "$@"; WORKERS="$2"; shift 2 ;;
    --update-search) UPDATE=1; shift ;;
    --export-only) EXPORT_ONLY=1; shift ;;
    --no-export) NO_EXPORT=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    --check-only) CHECK_ONLY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ "$FROM" =~ ^0[1-7](:[0-9]{2})?$ ]] || { echo "Invalid --from: $FROM" >&2; exit 2; }
[[ "$WORKERS" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid --workers: $WORKERS" >&2; exit 2; }
if (( EXPORT_ONLY && NO_EXPORT )); then echo 'Choose either --export-only or --no-export.' >&2; exit 2; fi
if (( CHECK_ONLY && DRY_RUN )); then echo 'Choose either --check-only or --dry-run.' >&2; exit 2; fi
if (( EXPORT_ONLY )) && [[ "$FROM" != 01:01 ]]; then echo '--from cannot be used with --export-only.' >&2; exit 2; fi
FROM_STAGE="${FROM%%:*}"
FROM_SCRIPT="01"
[[ "$FROM" == *:* ]] && FROM_SCRIPT="${FROM##*:}"

STAGES=(01_data_acquisition 02_expression_processing 03_gene_annotation 04_cohort_curation 05_technical_validation 06_batch_correction 07_visualization)
# The number of scripts of each stage, the last one its master quality control; a stray or missing script stops the build before anything runs.
SCRIPTS=(5 9 5 7 5 8 7)
# The scripts that can ask a repository, each with the stored answers it reads unless --update-search is given; the download script fetches the
# files of data/raw/download_lock.csv only with --download, and lists the repositories again only with --update-search.
declare -A ANSWERS=(
  [01_data_acquisition/01_search_repositories.py]="data/raw/census/search.json data/raw/census/SHA256SUMS.txt"
  [01_data_acquisition/02_harvest_study_metadata.py]="data/raw/metadata/SHA256SUMS.txt"
  [01_data_acquisition/04_download_raw_data.py]="data/raw/download_lock.csv"
  [02_expression_processing/01_inventory_expression_files.py]="data/raw/platforms/SHA256SUMS.txt"
  [03_gene_annotation/02_collect_annotation_sources.py]="data/raw/annotation/SHA256SUMS.txt"
)

for index in "${!STAGES[@]}"; do
  stage="${STAGES[$index]}"
  position=0
  name=""
  for script in src/"$stage"/[0-9][0-9]_*.py; do
    if [[ ! -e "$script" ]]; then break; fi
    position=$((position + 1))
    name="$(basename "$script")"
    if [[ "${name%%_*}" != "$(printf '%02d' "$position")" ]]; then
      echo "src/$stage holds $name where script $(printf '%02d' "$position") belongs: remove the stray copy or renumber" >&2; exit 2
    fi
  done
  if (( position != SCRIPTS[index] )) || [[ "$name" != *_qc.py ]]; then
    echo "src/$stage must hold scripts 01 to $(printf '%02d' "${SCRIPTS[$index]}"), the last its master quality control; found $position" >&2; exit 2
  fi
  if [[ "${stage%%_*}" == "$FROM_STAGE" ]] && (( 10#$FROM_SCRIPT < 1 || 10#$FROM_SCRIPT > position )); then
    echo "Invalid --from: $FROM; src/$stage holds scripts 01 to $(printf '%02d' "$position")" >&2; exit 2
  fi
done

# Whether the build runs a script of a stage, given --from.
selected() {
  local number="${1%%_*}"
  [[ ! "$number" < "$FROM_STAGE" ]] && ! [[ "$number" == "$FROM_STAGE" && "${2%%_*}" < "$FROM_SCRIPT" ]]
}

CURRENT=""
run() {
  printf '[RUN]'; printf ' %q' "$@"; printf '\n'
  CURRENT="$*"
  if (( ! DRY_RUN )); then "$@"; fi
}
failed() {
  local resume="rerun the same command"
  if [[ "$CURRENT" =~ src/0([1-7])_[a-z_]+/([0-9][0-9])_ ]]; then
    resume="resume with --from 0${BASH_REMATCH[1]}:${BASH_REMATCH[2]}"
  elif [[ "$CURRENT" == *pipeline/export_release.py* || "$CURRENT" == *release/build_data.py* ]]; then
    resume="resume with --export-only"
  fi
  echo "[BUILD FAILED] $(date -u +%Y-%m-%dT%H:%M:%SZ) while running: ${CURRENT:-the preflight}; fix the cause and ${resume}" >&2
}
trap failed ERR

if (( ! DRY_RUN )); then
  command -v "$PY" >/dev/null 2>&1 || { echo "Python interpreter not found: $PY" >&2; exit 2; }
  "$PY" - "$EXPORT_ONLY" "$OUTDIR" <<'PY'
import importlib
import json
import shutil
import sys
from pathlib import Path

if sys.version_info < (3, 10):
    raise SystemExit('Python 3.10 or newer is required.')
export_only, release = sys.argv[1:]
root, errors = Path.cwd().resolve(), []
for folder in sorted(root.glob('src/0[1-7]_*')):
    configs = sorted(folder.glob('*_config.json'))
    if len(configs) != 1:
        errors.append(f'{folder.name}: exactly one stage configuration is required')
    for path in configs + sorted(folder.glob('*_register.json')):
        try:
            json.loads(path.read_text(encoding='utf-8'))
        except (OSError, ValueError) as error:
            errors.append(f'Required configuration {path.relative_to(root)}: {error}')
if not (root / 'pipeline/export_release.py').is_file():
    errors.append('Missing pipeline/export_release.py.')
for name in ('numpy', 'pandas', 'pyarrow') + (() if export_only == '1' else ('requests', 'lxml', 'openpyxl', 'xlrd', 'matplotlib', 'PIL')):
    try:
        importlib.import_module(name)
    except ImportError as error:
        errors.append(f'Python module {name}: {error}')
if export_only != '1' and shutil.which('Rscript') is None:
    errors.append('Rscript not found; the expression stage needs R with Bioconductor.')
target = Path(release)
target = (target if target.is_absolute() else root / target).resolve()
if any(target == root / part or (root / part) in target.parents for part in ('src', 'data', 'outputs', 'pipeline', 'docker')) or root.is_relative_to(target):
    errors.append(f'Choose a separate release destination: {target}')
if errors:
    raise SystemExit('[PREFLIGHT FAILED]\n' + '\n'.join('- ' + error for error in errors))
print('[PREFLIGHT PASS] Environment, configurations and release destination are ready.')
PY
fi
if (( ! DRY_RUN && ! UPDATE && ! EXPORT_ONLY )); then
  missing=()
  for key in $(printf '%s\n' "${!ANSWERS[@]}" | sort); do
    if selected "${key%%/*}" "${key##*/}"; then
      for path in ${ANSWERS[$key]}; do [[ -f "$path" ]] || missing+=("$path"); done
    fi
  done
  if (( ${#missing[@]} )); then
    echo '[PREFLIGHT FAILED] Stored answers of the repositories are missing:' >&2
    printf -- '- %s\n' "${missing[@]}" >&2
    echo 'Unpack snapshots.tar.gz from the GigaDB package into data/raw/, or run with --update-search to ask the repositories again (new studies then need screening decisions).' >&2
    exit 2
  fi
  echo '[PREFLIGHT PASS] The stored answers of the repositories are present; the build reads them instead of asking again.'
fi
if (( CHECK_ONLY )); then
  run "$PY" pipeline/export_release.py --outdir "$OUTDIR" --check-only || echo '[CHECK] The stages are not yet complete; a build is needed before the export.'
  exit 0
fi
if (( ! DRY_RUN )); then
  mkdir -p outputs/logs
  # The log writer runs beside the shell; at exit the shell hands back its own output and waits for the writer, so that no last line is lost (under
  # docker run the shell is the first process, and the container ends with it).
  exec 3>&1 4>&2
  exec > >(tee -a outputs/logs/build_dataset.log) 2>&1
  LOGGER=$!
  trap 'exec 1>&3 2>&4; wait "$LOGGER"' EXIT
  MODE="stored answers"
  if (( UPDATE )); then MODE="updated search"; fi
  echo "[BUILD] $(date -u +%Y-%m-%dT%H:%M:%SZ) from ${FROM}, ${MODE}, export $(( ! NO_EXPORT ))"
fi

if (( ! EXPORT_ONLY )); then
  for stage in "${STAGES[@]}"; do
    for script in src/"$stage"/[0-9][0-9]_*.py; do
      name="$(basename "$script")"
      if ! selected "$stage" "$name"; then continue; fi
      options=()
      if (( UPDATE )) && [[ -v "ANSWERS[$stage/$name]" ]]; then options+=(--update-search); fi
      if [[ "$stage/$name" == 01_data_acquisition/04_download_raw_data.py ]]; then options+=(--download); fi
      run "$PY" "$script" ${options[@]+"${options[@]}"}
    done
  done
fi
if (( ! NO_EXPORT )); then
  run "$PY" pipeline/export_release.py --outdir "$OUTDIR" --workers "$WORKERS"
  if [[ -f release/build_data.py ]]; then
    run "$PY" release/build_data.py
  else
    echo '[SKIP] No release/ folder in this copy: the results page data are not rebuilt.'
  fi
fi
echo '[BUILD] Complete.'
