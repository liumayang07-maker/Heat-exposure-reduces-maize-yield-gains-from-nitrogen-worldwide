# Heat exposure reduces maize yield gains from nitrogen worldwide

Data and analysis code accompanying the manuscript with this working title. This initial deposit contains **Figure 1 only**. Other analyses are not yet included. The repository is private; editor/reviewer access must be arranged separately.

Start with [Figure 1](Figure1/README.md). See [verification](VERIFICATION.md) for the exact tests performed and [sources and rights](RIGHTS.md) for attribution and unresolved permissions. No publication or archive DOI has been assigned here.

## Environment

Tested runtime: Python 3.12.5 and R 4.5.3 on Windows. Direct Python dependencies are pinned in `requirements.txt`; `requirements-upstream.txt` adds optional external-input preparation dependencies. Use a separate environment:

```text
python -m venv .venv
# Activate .venv using the command appropriate for your shell.
python -m pip install -r requirements.txt
```

In R, install required packages if unavailable:

```r
install.packages(c("terra", "sf", "ggplot2", "data.table", "nlme", "pheatmap"))
```

Verified R packages: terra 1.9.1, sf 1.1.0, ggplot2 4.0.2, data.table 1.18.2.1, nlme 3.1.168, pheatmap 1.0.13. These are recorded versions, not an R lockfile. R needs Cairo PDF support; sf/terra need their usual system geospatial libraries. Python plots use Arial; missing fonts can change appearance. Put `Rscript` on PATH or supply its path to entries accepting `--rscript`.

Scripts do not install packages automatically. Fresh environment installation and other operating systems have not been validated in this deposit.

## Reproduction levels

1. **Bundled inputs:** regenerate five panel PDFs, six GLS model rows and two source workbooks, without downloads or the original working directory.
2. **External inputs:** obtain the area HDF5 and monthly climate exports, then run preparation scripts. This requires external access and potentially substantial compute; it was not rerun for this deposit.
3. **Final layout:** the assembled manuscript figure was manually edited. Reproducing individual panels does not recreate that layout.

No open-source or blanket open-data licence has been selected. Third-party terms remain applicable; see `RIGHTS.md`.
