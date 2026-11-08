# Figure 1: heat exposure across global maize systems

Scripts, archived inputs, numerical results and reference PDFs for panels 1a-d. Preserve the subdirectory structure: panels share code and geometry.

| Panel | Entry | Essential inputs | Outputs |
| --- | --- | --- | --- |
| [1a](Figure1a/README.md) | `Figure1a/plot_Figure1a.R` | mean TIFF and Natural Earth GeoPackage | PDF |
| [1b](Figure1b/README.md) | `Figure1b/plot_Figure1b.R` | continental annual CSV | PDF and GLS model CSV |
| [1c](Figure1c/README.md) | `Figure1c/plot_Figure1c.py` | slope/P TIFFs and 1a geometry | PDF |
| [1d](Figure1d/README.md) | `Figure1d/plot_Figure1d.py` | regional annual CSV and region definitions | two PDFs and source workbook |

## Run from the repository root

Archived outputs already exist. These commands explicitly replace them in your checkout. Run in a separate copy to preserve the references.

```text
Rscript Figure1/Figure1a/plot_Figure1a.R --overwrite
Rscript Figure1/Figure1b/plot_Figure1b.R --overwrite
python Figure1/build_source_data.py --rscript Rscript --overwrite
python Figure1/Figure1c/plot_Figure1c.py --data-dir Figure1/Figure1c --output Figure1/Figure1c/Figure1c.pdf --overwrite
python Figure1/Figure1d/plot_Figure1d.py --rscript Rscript --overwrite
```

`build_source_data.py` refits 1b through R and writes its source workbook. It also provides the workbook writer imported by 1d. Panel 1d imports HAC functions from `Figure1c/prepare_grid_EDH30_trends.py`; 1c shares upstream functions with `Figure1a/03_prepare_mean_EDH30.py`. The 1c map reads `Figure1a/world_boundaries.gpkg`.

Crosswalks support external-input reconstruction, preserving the actual spatial assignments. The 1c compressed source table preserves pixel-level model results beyond the two plotting TIFFs.

See [environment](../README.md), [verification](../VERIFICATION.md) and [sources/rights](../RIGHTS.md). Chinese panel notes retain detailed fields and historical checks; they do not establish that upstream reconstruction was repeated for this deposit. Individual PDFs do not recreate the manually assembled final figure.
