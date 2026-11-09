# Base images pinned by digest as well as tag, so that a rebuild starts from the same layers. The Ubuntu snapshot
# is the one under the environment in which the resource was built.
FROM mambaorg/micromamba:2.9.0@sha256:4f29a2b62a6e65ffa3dfbccb62d8a90514e95954aae52a90afaa8e7dae347bbf AS micromamba

FROM ubuntu:jammy-20240212@sha256:f9d633ff6640178c2d0525017174a688e2c1aef28f0a0130b26bd5554491f0da

ENV PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=0 \
    OMP_NUM_THREADS=8 \
    OPENBLAS_NUM_THREADS=8 \
    MKL_NUM_THREADS=8 \
    NUMEXPR_NUM_THREADS=8 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    MAMBA_ROOT_PREFIX=/opt/micromamba

WORKDIR /workspace

# Python lives in /opt/conda, and R with its Bioconductor packages in /opt/bio. Each prefix is installed from
# an explicit list of conda packages with their checksums, so that every build installs the same builds. /opt/bio
# is appended to PATH, so that python and pip are those of /opt/conda.
COPY --from=micromamba /bin/micromamba /usr/local/bin/micromamba
COPY --from=micromamba /etc/ssl/certs/ca-certificates.crt /etc/ssl/certs/ca-certificates.crt
COPY conda_bio.txt /tmp/genomictensorsepsis-conda-bio.txt
RUN MAMBA_REMOTE_MAX_RETRIES=8 MAMBA_REMOTE_BACKOFF_FACTOR=3 \
    micromamba create -y -p /opt/bio -f /tmp/genomictensorsepsis-conda-bio.txt \
    && micromamba clean --all --yes
ENV PATH="/opt/conda/bin:${PATH}:/opt/bio/bin"
RUN Rscript -e 'for (p in c("oligo", "oligoClasses", "affxparser", "affyio", "limma", "edgeR", "AnnotationDbi", "jsonlite", "pd.hg.u133.plus.2", "pd.hg.u133a.2", "pd.hg.u219", "pd.hugene.1.0.st.v1", "pd.hugene.2.0.st", "pd.hugene.2.1.st", "pd.hta.2.0", "pd.huex.1.0.st.v2", "pd.clariom.d.human", "hgu133plus2.db", "hgu133a2.db", "hgu219.db", "hugene10sttranscriptcluster.db", "hugene20sttranscriptcluster.db", "hugene21sttranscriptcluster.db", "huex10sttranscriptcluster.db", "hta20transcriptcluster.db", "clariomdhumantranscriptcluster.db", "illuminaHumanv2.db", "illuminaHumanv3.db", "illuminaHumanv4.db", "hgug4112a.db", "HsAgilentDesign026652.db", "org.Hs.eg.db")) suppressMessages(library(p, character.only = TRUE))'

# preprocessCore is rebuilt from source without threading. The bioconda build
# stops with "return code from pthread_create() is 22" in this image, and
# oligo::rma stops with it. The source of the pinned version is fetched from
# Bioconductor 3.22, or from its archive once a later version replaces it.
RUN Rscript -e 'source <- file.path(tempdir(), "preprocessCore_1.72.0.tar.gz"); for (url in c("https://bioconductor.org/packages/3.22/bioc/src/contrib/preprocessCore_1.72.0.tar.gz", "https://bioconductor.org/packages/3.22/bioc/src/contrib/Archive/preprocessCore/preprocessCore_1.72.0.tar.gz")) if (!file.exists(source)) try(download.file(url, source, quiet = TRUE), silent = TRUE); stopifnot(file.exists(source)); install.packages(source, repos = NULL, type = "source", configure.args = "--disable-threading")' \
    && Rscript -e 'library(preprocessCore); stopifnot(packageVersion("preprocessCore") == "1.72.0", all(dim(rma.background.correct(matrix(rexp(2000, 0.01), ncol = 2))) == c(1000, 2)))'

COPY conda_python.txt /tmp/genomictensorsepsis-conda-python.txt
RUN MAMBA_REMOTE_MAX_RETRIES=8 MAMBA_REMOTE_BACKOFF_FACTOR=3 \
    micromamba create -y -p /opt/conda -f /tmp/genomictensorsepsis-conda-python.txt \
    && micromamba clean --all --yes
COPY requirements.txt /tmp/genomictensorsepsis-requirements.txt
RUN python -m pip install --no-cache-dir --no-deps --retries 10 --timeout 60 -r /tmp/genomictensorsepsis-requirements.txt \
    && python -m pip check

CMD ["bash"]
