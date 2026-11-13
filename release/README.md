# Results page

`release/` is a static page that explains GenomicTensorSepsis to readers without a technical background and lets them explore it. It shows how the resource was built, with animations of real samples, one of them of twins across a random and a patient-safe split; what it holds, chapter by chapter, with what the benchmark measures; each included study on its own, with the excluded candidates and their reasons; and every step of the method with its equation, a simulation and a worked example on real samples. It has a light and a dark theme.

The page is published at https://mrrajat1809.github.io/GenomicTensorSepsis/ by `.github/workflows/pages.yml`. Locally, open `release/index.html` in a browser or host the folder as it is; the page needs no network access. Any view or chapter can be opened directly by its address, such as `index.html#methods/m-over` or `index.html#studies/GSE65682/overview`, and `?theme=dark` or `?theme=light` sets the theme without remembering it.

All its numbers come from the page data, which `release/build_data.py` writes from the pipeline outputs: the tables in `data/processed/`, the review tables, receipts and QC reports in `outputs/`, the points of the paper's figures in `outputs/visualization/figure_data.csv` and the stage configurations. `index.html?selftest=1` loads every part, opens every view and chapter in both themes, runs the checks of each view (counts that add up to their totals, numbers recomputed from the data, worked examples recomputed in the page) and reports the result in a panel and in the page title, as `SELFTEST PASS n` or `SELFTEST FAIL k/n`.

To rebuild the page data after a pipeline run, from the repository root, on the host or in the pipeline's image (pandas and pyarrow are the only packages it uses):

```bash
python release/build_data.py               # every part, then the content versions of the scripts
python release/build_data.py methods       # one part: core, story, resource, studies or methods
```

| Location | Contents |
| --- | --- |
| `release/index.html` | The page: its styles, the header and one section per view |
| `release/common.js`, `release/main.js` | Shared helpers, addresses and the self-test; the start page |
| `release/story.js`, `resource.js`, `studies.js`, `methods.js` | One script per view: How it works, The resource, Studies, Method in detail |
| `release/data.js` | The numbers the page needs when it opens: main numbers, labels, colors, glossary and the start page |
| `release/data/` | One file per view, and one per study in `data/studies/`, loaded when the view or study is first shown |

The page is written in HTML, CSS and JavaScript, without external libraries.

The data of GenomicTensorSepsis, including the page data, are released under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/); the code, of which the page is part, under the MIT License. The expression data were deposited by the authors of the source studies; please cite the studies whose samples you use.
