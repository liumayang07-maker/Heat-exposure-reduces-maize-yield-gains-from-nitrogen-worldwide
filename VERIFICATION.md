# Figure 1 deposit verification

Verified on 2026-11-08 using Python 3.12.5 and R 4.5.3 on Windows. Exact direct dependency versions are recorded in `requirements.txt`, `requirements-upstream.txt` and the repository README. Existing installed environments were used; a fresh dependency installation and other operating systems were not tested.

## Isolation and unchanged references

The 35 selected files were copied from the manuscript archive into an independent repository directory. Only documentation was revised there; scientific scripts, numerical inputs and archived outputs were not changed. A second complete Figure1 copy was used for execution, so its overwritten outputs did not replace the deposit references. The original archive is checked by SHA256 before push. The seven excluded files are the two final-layout assets and five historical audit files described in the preparation plan; no original file was deleted.

## Commands executed

The five commands in `Figure1/README.md` were executed with the actual Rscript path substituted. All returned exit code 0. Execution times were approximately 4.75 s (1a), 1.38 s (1b), 2.70 s (1b workbook), 4.94 s (1c) and 8.80 s (1d); these are local observations, not performance guarantees.

| Check | Result |
| --- | --- |
| 1b GLS model table | Six rows reproduce the archived table; numeric comparison rtol=1e-12, atol=1e-12 |
| 1b workbook | Five sheets reproduce: 240 annual rows, 240 plot rows, six GLS rows, six parameter rows, 24 covariance rows |
| 1d workbook | Five sheets reproduce: 615 annual rows, 615 heatmap rows, 15 trends, 15 parameter rows, 60 covariance rows |
| Workbook numerical tolerance | rtol=1e-12, atol=1e-12; sheet names/order and text values also compared |
| Five PDFs | One page each, equal page dimensions; full-page RGB rendering at 144 dpi has zero differing components against archived references |
| 1a raster | Single-band EPSG:4326; 599,246 valid, nonnegative pixels, including 28,070 zeros |
| 1c source table | 540,269 rows; pointwise intervals independently checked with t(0.975,38) and supplied HAC SE |
| 1c BH correction | Recomputed across all supplied raw P values; 308,069 adjusted P values below 0.05 |
| 1c rasters versus source table | All pixel coordinates mapped back to both rasters; valid-cell counts match. Slope tolerance rtol=1e-6, atol=1e-5; adjusted-P tolerance rtol=1e-6, atol=1e-7 |

Workbook/PDF metadata timestamps need not match byte-for-byte. Numerical and rendered output comparisons, rather than file hashes, establish reproduction. Input checksums in `MANIFEST.sha256` identify the deposited reference files.

## Limits

No Earth Engine task was submitted, no external climate or area data downloaded, and no complete upstream aggregation rerun. The 1a mean and 1c trends were not independently rebuilt from hourly/monthly inputs in this deposit check. Their archived inputs, plotting and source-table consistency were checked. The manually assembled final figure was not reconstructed. Source-access and redistribution questions remain documented in `RIGHTS.md`; passing numerical checks does not resolve them.
