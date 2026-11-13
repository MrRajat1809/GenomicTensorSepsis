# GenomicTensorSepsis: a harmonized multi-cohort blood transcriptomic resource for sepsis research

GenomicTensorSepsis (GTS) brings the public blood transcriptome studies of sepsis into one resource that can be analyzed across studies. This repository holds the pipeline that finds, reprocesses, annotates, curates, validates, corrects and benchmarks the studies, and packages the result.

## Overview

GTS holds 12,037 sample records with expression values from 96 studies in GEO (85) and ArrayExpress (11), released from 2006 to 2026, on 105 study platforms (the samples of one study on one platform). The records are 11,043 distinct samples. Every study platform is processed on its own from the lowest data level of its deposit, mapped to HGNC genes and kept on its own scale. The 73 records of the included studies that have no expression values are listed with their reason.

| Technology | Study platforms | Sample records | Distinct samples |
| --- | ---: | ---: | ---: |
| Affymetrix arrays | 32 | 3,693 | 3,068 |
| Illumina BeadChip | 21 | 3,349 | 3,196 |
| Agilent arrays | 5 | 752 | 675 |
| Other arrays | 2 | 152 | 152 |
| RNA-seq | 45 | 4,091 | 3,952 |

The distinct samples are curated from their records as sepsis (5,975), healthy control (1,189), non-septic control (737), other condition (3,139) or not recorded (3), and every label keeps the wording it rests on. Repeat deposits of one array or library, repeated samples of one patient (4,206 identified patients) and samples with quality concerns are marked, not removed. Two corrected tables place bulk blood (7,604 samples of 59 studies) and PBMC (1,013 samples of 10 studies) on a common scale. Benchmark splits for three tasks, 64 folds in all, keep patients, repeat deposits and groups of studies apart. Under a random split of the records, 34% of test records have a deposit of the same array or a sample of the same patient in training, and baseline classifiers of sepsis against non-septic controls reach a mean AUROC of 0.99, against 0.88 to 0.89 when whole patient-safe groups are held out.

![Selection of the studies and samples](Figures/Fig1.png)

The [results page](https://mrrajat1809.github.io/GenomicTensorSepsis/) explains the resource without technical background and lets you explore every study, check, fold and benchmark result.

## Using the data

The data package is deposited in Zenodo under CC0 1.0: [doi:10.5281/zenodo.24205193](https://doi.org/10.5281/zenodo.24205193). It holds the tables as gzip-compressed CSV files in tar archives, a data dictionary that defines every column, and a readme. Download every file with the Download all button of the record or with `zenodo_get 10.5281/zenodo.24205193` (`pip install zenodo-get`), or only the archives you need. Check a download with `md5sum -c md5sum.txt`; on macOS, which has no `md5sum`, run `while read -r sum file; do [ "$(md5 -q "$file")" = "$sum" ] && echo "$file: OK" || echo "$file: FAILED"; done < md5sum.txt`.

| Question | Where to look |
| --- | --- |
| Which samples are there, and what do their records say? | `cohort.tar`: `cohort_samples.csv.gz`, one row per sample record |
| What are the expression values of one study? | `gene_expression.tar`: one gene-level table per study platform |
| What are the read counts of an RNA-seq study? | `gene_counts.tar` and `probe_counts.tar` |
| Which records are the same array, and which samples share a patient? | `cohort.tar`: `cohort_sample_links.csv.gz` and the patient columns of `cohort_samples.csv.gz` |
| Where does a label come from? | `cohort.tar`: `cohort_label_sources.csv.gz` and `cohort_sample_metadata.csv.gz` |
| Which gene does a probe measure? | `gene_mapping.tar` |
| Can I compare samples of many studies on one scale? | `corrected_expression.tar`: the bulk-blood and PBMC tables |
| Which values do the benchmark baselines read? | `corrected_expression.tar`: `corrected_expression_aligned_bulk_blood_expression.csv.gz`, every bulk-blood record on one reference distribution, made without labels |
| How do I test a classifier without leakage? | `cohort.tar`: `cohort_benchmark_tasks.csv.gz` and `cohort_benchmark_splits.csv.gz` |
| How do naive and patient-safe evaluations compare? | `benchmark.tar`: linked records, baseline models, study shortcuts and published signatures |
| Which records have no expression values, and why? | `provenance.tar.gz`: `provenance_records_without_values.csv.gz` |
| What does a column mean? | `data_dictionary.csv` |

For example, to take one sample per patient of GSE65682 with its expression values:

```python
import pandas as pd

samples = pd.read_csv("cohort_samples.csv.gz")                          # from cohort.tar
expression = pd.read_csv("gene_expression_GSE65682__GPL13667.csv.gz",   # from gene_expression.tar
                         index_col="symbol").drop(columns="hgnc_id")    # genes by samples, log2 scale

chosen = samples[(samples.study_accession == "GSE65682") & (samples.one_per_patient == 1) & (samples.any_flag == 0)]
labels = chosen.set_index("sample_accession").sample_class
values = expression[labels.index.intersection(expression.columns)]
```

Use `is_primary_record` to count each sample once, `one_per_patient` to keep one sample of every patient, and `patient_safe_group` to split data so that no patient or sample appears on both sides.

## Repository Structure

```text
GenomicTensorSepsis/
├── src/
│   ├── 01_data_acquisition/       # Repository searches, screening and downloads
│   ├── 02_expression_processing/  # Processing of every study from its deposited data
│   ├── 03_gene_annotation/        # Probe-to-gene mapping and gene-level tables
│   ├── 04_cohort_curation/        # Sample labels, repeat deposits, patients, flags and benchmark splits
│   ├── 05_technical_validation/   # Recorded sex, repeat deposits and known sepsis biology
│   ├── 06_batch_correction/       # Cross-study correction of bulk blood and PBMC
│   ├── 07_benchmark/              # Split leakage, baseline models, study shortcuts and published signatures
│   └── 08_visualization/          # Main and supplementary figures
├── pipeline/
│   ├── 01_build_dataset.sh        # Run the complete workflow in dependency order
│   └── export_release.py          # Package the data for Zenodo
├── release/                       # The results page and the script that writes its data
├── .github/workflows/pages.yml    # Publishes the results page to GitHub Pages
├── docker/                        # Dockerfile and pinned package lists
├── Figures/                       # Main and supplementary figures
├── CITATION.cff
├── LICENSE
└── README.md
```

`data/`, `outputs/` and `data_release/` are local working directories excluded from Git. Each stage folder holds numbered scripts and one JSON configuration with its rules and limits.

## Workflow

1. **Acquisition:** search GEO, ArrayExpress and ENA, screen every candidate study with one recorded decision and reason, and download the expression files of the included studies.
2. **Expression processing:** process every study platform from the lowest level its files support: RMA for Affymetrix arrays, limma for other arrays, TMM-normalized log2 counts per million for RNA-seq. Values that the authors already normalized are kept where no raw data exist.
3. **Gene annotation:** map every row of every table to HGNC genes and build gene-level tables. No row is dropped; rows without a gene keep their reason.
4. **Cohort curation:** curate sample class, tissue, sex, age, time point, outcome and septic shock from the records, link repeat deposits, identify patients, flag quality concerns and define the benchmark splits.
5. **Technical validation:** compare recorded sex with expression, check that repeat deposits agree and test known sepsis biology in every study.
6. **Batch correction:** give every bulk-blood and PBMC sample one reference distribution of values, correct the studies on their healthy controls or by matched groups, and test the correction on held-out healthy controls.
7. **Benchmark:** count the test records that naive splits leave linked to a training record, compare baseline classifiers under naive and patient-safe splits, test study shortcuts and re-test published sepsis signatures within single studies.
8. **Figures:** draw the main and supplementary figures; `export_release.py` then packages the tables with a data dictionary, cohort flow and manifest.

Each stage ends with a quality-control script that must pass before the next stage starts.

## Rebuilding the Resource

```bash
git clone https://github.com/MrRajat1809/GenomicTensorSepsis.git
cd GenomicTensorSepsis
docker build -t genomictensorsepsis:1.0.0 docker
docker run --rm -it -v "$PWD":/workspace -w /workspace genomictensorsepsis:1.0.0 bash
```

The image is built for x86-64 Linux. On a Mac with Apple silicon (M1 or later), add `--platform linux/amd64` to both commands (`docker build --platform linux/amd64 ...` and `docker run --platform linux/amd64 ...`); Docker Desktop then runs the image through Rosetta. Give Docker at least 8 GB of memory (Docker Desktop: Settings, Resources).

The search is frozen: the data package stores the answers the repositories gave on the census date (11 October 2026), and the build reads them instead of searching again. Unpack `snapshots.tar.gz` from the package into `data/raw/`, then, inside the container:

```bash
bash pipeline/01_build_dataset.sh --check-only   # check the environment, configurations and stored answers
bash pipeline/01_build_dataset.sh                # build every stage and the data package
```

The build downloads the listed expression files (about 46 GB), checks each against its checksum and takes about eight hours on a desktop computer. `--update-search` asks the repositories again, limited to studies released by the census date; new studies are added only by moving that date and screening them.

| Task | Command |
| --- | --- |
| Resume at a stage or script | `bash pipeline/01_build_dataset.sh --from 02:04` |
| Redraw the figures from completed stages | `bash pipeline/01_build_dataset.sh --from 08 --no-export` |
| Package completed stages | `bash pipeline/01_build_dataset.sh --export-only` |
| List the steps without running anything (any bash, including the macOS terminal) | `bash pipeline/01_build_dataset.sh --dry-run` |
| List all options | `bash pipeline/01_build_dataset.sh --help` |

## Study Design and Interpretation

- **Labels:** a sample is sepsis only where its record, or a quoted sentence of its study, says so. Diagnoses were not re-adjudicated, and most records do not state the sepsis definition used.
- **Repeat deposits and patients:** a sample without a patient identifier counts as its own patient, so the samples marked one per patient are an upper bound on the number of patients.
- **Scales:** the per-study tables are on the scale of their processing; they are the reference for any result and are not comparable across studies without a correction.
- **Corrected tables:** studies with healthy controls are corrected on them; studies without are corrected by matched groups, which line up their class means with the other studies, so their class differences are not independent evidence. Newborn studies, plasma, swabs and sorted cells are not in the corrected tables.
- **Benchmark:** the per-study tables remain the reference for any result. The benchmark baselines read the aligned bulk-blood table, which gives every record one reference distribution without class or control information; the corrected tables used the labels and healthy controls of every study and are not inputs for the sepsis tasks. Fit any scaling on the training studies alone and score each held-out group on its own.
- **Validation:** the checks test processing and labeling against the data; they do not assess whether a deposited label is clinically correct.

## Data Availability

The data package is available in Zenodo: [doi:10.5281/zenodo.24205193](https://doi.org/10.5281/zenodo.24205193). The code of release v1.0.0, which built it, is archived in Zenodo: [doi:10.5281/zenodo.24206598](https://doi.org/10.5281/zenodo.24206598). Downloaded and generated data are not stored in this repository. The source studies remain the work of their investigators: please cite them (`cohort_studies.csv.gz` lists their accessions and addresses) together with this resource.

## How to Cite

Please cite the dataset, and the source studies whose samples you use:

Kumar P. GenomicTensorSepsis: a harmonized multi-cohort blood transcriptomic resource for sepsis research [Data set]. Version 1.0.0. Zenodo; 2026. https://doi.org/10.5281/zenodo.24205193

To cite the code: Kumar P. GenomicTensorSepsis [Software]. Version 1.0.0. Zenodo; 2026. https://doi.org/10.5281/zenodo.24206598

The article describing the resource is in preparation; its citation will be added here. `CITATION.cff` gives the citation in machine-readable form.

## License

The code is distributed under the [MIT License](LICENSE). The data package, curated metadata and documentation are dedicated to the public domain under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/). The source data remain subject to the terms of the repositories that hold them.

## Use of Generative AI

Generative AI (Claude, Anthropic) was used to assist with writing and reviewing code and documentation. The author directed the work, made the methodological decisions, checked the results, and takes responsibility for the content of this repository.

## Contact

**Author:** Priyanshu Kumar, Chandigarh University  
**Repository:** [MrRajat1809/GenomicTensorSepsis](https://github.com/MrRajat1809/GenomicTensorSepsis)
