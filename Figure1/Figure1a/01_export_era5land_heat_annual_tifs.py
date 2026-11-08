# Figure1a upstream: ERA5-Land monthly EDH30, 1986-2025 (UTC).
# Adapted from global_maize_edh_1980_2025/scripts/gee/01_export_era5land_heat_annual_tifs.py.
# Requires earthengine-api; authenticate with: earthengine authenticate
# Run: python 01_export_era5land_heat_annual_tifs.py --project YOUR_PROJECT_ID
# Exports 40 GeoTIFFs to Google Drive, each with 12 EDH30_YYYY_MM bands.
# Only EDH30 is retained; original exports also contained HH30/HH35/EDH35.
# Annual sums, yearly maize masks and the multi-year mean are computed locally.
# Retain the exported affine transform; scale is 11132 m, not exact 0.1 degrees.
import argparse
import time


def monthly_edh30(ee, year, month):
    start = ee.Date.fromYMD(year, month, 1)
    hourly = (ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
              .filterDate(start, start.advance(1, 'month')).select('temperature_2m'))

    def excess(image):
        return (image.select('temperature_2m').subtract(273.15)
                .subtract(30.0).max(0.0).rename(f'EDH30_{year}_{month:02d}'))

    return hourly.map(excess).sum().toFloat()


def main():
    parser = argparse.ArgumentParser(description='Export monthly EDH30 for 1986-2025 to Google Drive.')
    parser.add_argument('--project', required=True, help='Earth Engine-enabled Google Cloud project ID.')
    parser.add_argument('--dry-run', action='store_true', help='Build tasks without submitting; requires EE authentication.')
    parser.add_argument('--start-year', type=int, default=1986, choices=(1984, 1985, 1986),
                        help='Use 1984 for Figure1d southern summer; default preserves Figure1a-c.')
    args = parser.parse_args()
    import ee

    ee.Initialize(project=args.project)
    region = ee.Geometry.Rectangle([-179.95, -89.95, 179.95, 89.95], None, False)
    for year in range(args.start_year, 2026):
        image = ee.Image.cat([monthly_edh30(ee, year, month) for month in range(1, 13)]).toFloat()
        name = f'era5land_heat_monthly_{year}'
        task = ee.batch.Export.image.toDrive(
            image=image, description=name, fileNamePrefix=name,
            folder='global_maize_edh_era5land_exports', region=region,
            crs='EPSG:4326', scale=11132, maxPixels=1e13,
            fileFormat='GeoTIFF', formatOptions={'cloudOptimized': True})
        if args.dry_run:
            print(f'Prepared only: {name}')
        else:
            task.start()
            print(f'Submitted: {name}; task_id={task.id}')
            time.sleep(1)


if __name__ == '__main__':
    main()
