# Data sources, access and rights

Source pages checked on 2026-11-08. This record documents provenance and access, not legal clearance for an open release. The repository remains private. No open-source licence for author code or blanket licence for author outputs has been selected; rights are reserved pending an explicit choice. Review access does not itself grant unrestricted redistribution.

| Source | Use here | Verified terms/access | Treatment |
| --- | --- | --- | --- |
| [ERA5-Land / Earth Engine](https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY) | Hourly `temperature_2m`, converted to EDH30 | Official catalogue supplies acknowledgement and licence conditions | Derived results and code included; raw hourly/monthly files excluded |
| [Natural Earth](https://www.naturalearthdata.com/about/terms-of-use/) | Country geometry; continental and Chinese admin-1 assignments | Source states map data are public domain | Geometry-only GeoPackage and attribution included |
| [GADM](https://gadm.org/license.html) | Countries in the original regional assignment | Academic/non-commercial use permitted; redistribution of GADM data requires prior permission | No original GADM boundaries included. Derived ID-to-region table retained for private research review; unrestricted redistribution of that table is not asserted or verified |
| [Area dataset DOI](https://doi.org/10.11888/Terre.tpdc.300446), [source article](https://doi.org/10.1038/s41597-023-02526-z), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10495426/) | `Harvested_area_1961-2020.h5`, `Maize` group, mask and annual weights | DOI portal inaccessible during this check. Article CC BY 4.0 does not alone establish the HDF5 redistribution licence | HDF5 and annual area tables excluded; obtain data and applicable permission from provider |
| [Zenodo 7408852](https://zenodo.org/records/7408852) | Alternate area/fertilizer dataset record | Metadata public; files explicitly restricted on check date | Not represented as an anonymous download mirror |

ERA5-Land acknowledgement: **Generated using Copernicus Climate Change Service Information (2026).** Exposure metrics are modifications/analyses by the manuscript authors. Neither the European Commission nor ECMWF is responsible for any use made of the Copernicus information or data it contains. The year denotes this deposit, not a new climate-data version.

Made with Natural Earth. Country geometry was converted to GeoPackage without simplification or reprojection. Crosswalks contain IDs and labels, not boundary polygons or harvested-area values. They preserve actual assignments rather than substituting different boundaries.

## Before public release

Confirm area-file access and terms covering redistribution of derived spatial assignments and area-weighted results. Clarify the derived regional crosswalk's status with the applicable rights holder if needed. Select an author-approved code/output licence separately from third-party terms. Arrange editor/reviewer access while private. No provider contact or permission request has been made here.
