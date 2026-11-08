from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import logging
import platform
import shutil
import sys
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / 'site_soil_configuration.csv'
SNAPSHOT = SOURCE
DATASET = 'reanalysis-era5-land-timeseries'
DATES = ['2010-01-01', '2025-12-31']
EXPECTED = pd.date_range(DATES[0], '2026-01-01', freq='h', inclusive='left')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=True), encoding='utf-8')
    tmp.replace(path)


def prepare():
    # 1. 直接复用81站点配置，不另存坐标快照，不创建小时数据目录。
    sites = pd.read_csv(
        SOURCE,
        encoding='utf-8-sig',
        usecols=['selected_site_id', 'lat', 'lon'],
    )
    if len(sites) != 81 or sites.selected_site_id.nunique() != 81:
        raise ValueError('Expected exactly 81 unique sites')
    if not sites.selected_site_id.str.fullmatch(r'APSIM_\d{3}').all():
        raise ValueError('Invalid site ID')
    if not sites.lat.between(-90, 90).all() or not sites.lon.between(-180, 180).all():
        raise ValueError('Invalid coordinates')
    return sites

def request_for(site):
    return {'variable': ['2m_temperature', 'total_precipitation'],
            'date': DATES, 'location': {'latitude': float(site.lat), 'longitude': float(site.lon)},
            'data_format': 'netcdf'}


# 3. 保留原NetCDF读取、单位转换、有效性检查与重试下载逻辑。
def read_raw(folder, site):
    payload = folder / 'response.download'
    files = []
    if zipfile.is_zipfile(payload):
        extracted = folder / 'netcdf'
        extracted.mkdir(exist_ok=True)
        with zipfile.ZipFile(payload) as archive:
            for member in archive.infolist():
                if member.is_dir() or not member.filename.lower().endswith('.nc'):
                    continue
                target = extracted / Path(member.filename).name
                if target in files:
                    raise ValueError('Duplicate NetCDF names in archive')
                with archive.open(member) as src, target.open('wb') as dst:
                    shutil.copyfileobj(src, dst)
                files.append(target)
    else:
        files = [payload]
    arrays = {}
    coords = []
    units = {}
    for path in files:
        with path.open('rb') as stream, xr.open_dataset(stream, engine='h5netcdf') as ds:
            for short, candidates in [('t2m', ['t2m', '2m_temperature']),
                                      ('tp', ['tp', 'total_precipitation'])]:
                names = [v for v in candidates if v in ds.data_vars]
                if not names:
                    continue
                if short in arrays:
                    raise ValueError(f'Duplicate variable: {short}')
                da = ds[names[0]]
                coord = []
                for choices in [('latitude', 'lat'), ('longitude', 'lon')]:
                    name = next((n for n in choices if n in ds), None)
                    if name is None or ds[name].size != 1:
                        raise ValueError('Expected one explicit latitude/longitude')
                    coord.append(float(ds[name].values.reshape(-1)[0]))
                coord[1] = (coord[1] + 180) % 360 - 180
                coords.append(coord)
                tname = next((n for n in ['time', 'valid_time'] if n in da.dims), None)
                if tname is None:
                    raise ValueError('Missing hourly time dimension')
                for dim in list(da.dims):
                    if dim != tname:
                        if da.sizes[dim] != 1:
                            raise ValueError(f'Unexpected dimension: {dim}')
                        da = da.isel({dim: 0})
                index = pd.DatetimeIndex(da[tname].values)
                if not index.equals(EXPECTED):
                    raise ValueError(f'{short}: hourly time axis differs from 2010-2025 UTC')
                values = np.asarray(da.values, dtype=float)
                if not np.isfinite(values).all():
                    raise ValueError(f'{short}: missing or non-finite values')
                units[short] = str(da.attrs.get('units', ''))
                arrays[short] = values
    if set(arrays) != {'t2m', 'tp'}:
        raise ValueError(f'Missing expected variables; found {list(arrays)}')
    if not all(np.allclose(c, coords[0], atol=1e-6, rtol=0) for c in coords):
        raise ValueError('Temperature/precipitation grid coordinates differ')
    lat, lon = coords[0]
    if abs(lat - float(site.lat)) > 0.05001 or abs(lon - float(site.lon)) > 0.05001:
        raise ValueError('Returned grid is not the nearest 0.1-degree grid point')
    if units['t2m'] not in ['K', 'kelvin'] or units['tp'] not in ['m', 'm of water equivalent']:
        raise ValueError(f'Unexpected units: {units}')
    frame = pd.DataFrame({'time_utc': EXPECTED, 'selected_site_id': site.selected_site_id,
                          'site_lat': site.lat, 'site_lon': site.lon,
                          'grid_lat': lat, 'grid_lon': lon,
                          't2m_C': arrays['t2m'] - 273.15, 'tp_mm': arrays['tp'] * 1000})
    # ARCO precipitation is already de-accumulated; retain negative values for QC.
    qc = {'selected_site_id': site.selected_site_id, 'rows': len(frame),
          'start_utc': str(EXPECTED[0]), 'end_utc': str(EXPECTED[-1]),
          'grid_lat': lat, 'grid_lon': lon, 'missing_values': 0,
          'temperature_min_C': float(frame.t2m_C.min()), 'temperature_max_C': float(frame.t2m_C.max()),
          'precipitation_min_mm': float(frame.tp_mm.min()), 'precipitation_max_mm': float(frame.tp_mm.max()),
          'negative_precipitation_hours': int((frame.tp_mm < 0).sum()),
          'source_units': units}
    return frame, qc


def validate_csv(path, expected_frame):
    saved = pd.read_csv(path)
    saved['time_utc'] = pd.to_datetime(saved.time_utc, utc=True).dt.tz_localize(None)
    pd.testing.assert_frame_equal(saved, expected_frame, check_dtype=False, check_exact=False,
                                  rtol=1e-9, atol=1e-8)


def transfer_payload(receipt, part):
    offset = part.stat().st_size if part.exists() else 0
    total = int(receipt.get('content_length') or 0)
    if total and offset == total:
        return
    if total and offset > total:
        raise ValueError('Partial file exceeds CDS content length')
    started = time.monotonic()
    with requests.get(receipt['download_url'], stream=True, timeout=(15, 45),
                      headers={'Range': f'bytes={offset}-'}) as response:
        response.raise_for_status()
        if response.status_code == 206:
            if not response.headers.get('Content-Range', '').startswith(f'bytes {offset}-'):
                raise ValueError('Unexpected Content-Range in resumed download')
            mode = 'ab' if offset else 'wb'
        else:
            mode = 'wb'
        with part.open(mode) as stream:
            for block in response.iter_content(16 * 1024):
                stream.write(block)
                if time.monotonic() - started > 180:
                    raise TimeoutError('Transfer exceeded 180 seconds; saved bytes will be resumed')
    if total and part.stat().st_size != total:
        raise ValueError('Downloaded byte count differs from CDS content length')
    if part.stat().st_size == 0:
        raise ValueError('Empty download')


def process(site, client, verify_only=False):
    sid = site.selected_site_id
    folder = ROOT / 'raw' / sid
    folder.mkdir(parents=True, exist_ok=True)
    request = request_for(site)
    reqpath = folder / 'request.json'
    if reqpath.exists() and json.loads(reqpath.read_text()) != request:
        raise ValueError(f'{sid}: existing request differs')
    write_json(reqpath, request)
    payload = folder / 'response.download'
    if not payload.exists():
        if verify_only:
            raise FileNotFoundError(f'{sid}: raw download missing')
        part = folder / 'response.download.part'
        logging.info('%s requesting full 2010-2025 series', sid)
        receipt_path = folder / 'cds_receipt.json'
        receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
        if not receipt.get('download_url'):
            if receipt.get('request_id'):
                remote = client.client.get_remote(receipt['request_id'])
            else:
                remote = client.retrieve(DATASET, request)
                receipt = {'request_id': remote.request_id,
                           'submitted_utc': datetime.now(timezone.utc).isoformat()}
                write_json(receipt_path, receipt)
            result = remote.get_results()
            receipt.update({'download_url': result.location, 'content_length': result.content_length})
            write_json(receipt_path, receipt)
        transfer_payload(receipt, part)
        part.replace(payload)
    frame, qc = read_raw(folder, site)
    csv = ROOT / 'hourly_csv' / f'{sid}_2010_2025_hourly_utc.csv'
    if csv.exists():
        validate_csv(csv, frame)
    elif verify_only:
        raise FileNotFoundError(f'{sid}: converted CSV missing')
    else:
        tmp = csv.with_suffix('.csv.part')
        frame.to_csv(tmp, index=False, date_format='%Y-%m-%dT%H:%M:%SZ')
        # Read back as UTC-naive values for comparison with source netCDF time.
        saved = pd.read_csv(tmp)
        saved['time_utc'] = pd.to_datetime(saved.time_utc, utc=True).dt.tz_localize(None)
        pd.testing.assert_frame_equal(saved, frame, check_dtype=False, rtol=1e-9, atol=1e-8)
        tmp.replace(csv)
    qc.update({'status': 'passed_with_flags' if qc['negative_precipitation_hours'] else 'passed',
               'raw_sha256': sha256(payload), 'csv_sha256': sha256(csv),
               'input_sha256': getattr(site, 'input_sha256', sha256(SNAPSHOT)),
               'checked_utc': datetime.now(timezone.utc).isoformat()})
    if hasattr(site, 'input_snapshot'):
        qc['input_snapshot'] = site.input_snapshot
    write_json(ROOT / 'qc' / f'{sid}.json', qc)
    logging.info('%s validated: %d hours; negative precipitation hours=%d', sid, len(frame), qc['negative_precipitation_hours'])
    return qc


def run_site(site, verify_only):
    client = None
    if not verify_only and not (ROOT / 'raw' / site.selected_site_id / 'response.download').exists():
        import cdsapi
        client = cdsapi.Client(timeout=60, retry_max=3, sleep_max=10, quiet=True,
                               progress=False, wait_until_complete=False)
    for attempt in range(1, 4):
        try:
            return process(site, client, verify_only)
        except Exception as exc:
            message = str(exc)
            logging.error('%s attempt %d: %s', site.selected_site_id, attempt, message)
            write_json(ROOT / 'metadata' / f'{site.selected_site_id}_last_error.json',
                       {'site': site.selected_site_id, 'attempt': attempt, 'error': message,
                        'utc': datetime.now(timezone.utc).isoformat()})
            fatal = any(s in message.lower() for s in ['403', '401', 'licence', 'license', 'terms', '400', 'invalid request'])
            if fatal or verify_only or attempt == 3:
                raise
            time.sleep(10 * attempt)


def main():
    global ROOT
    # 2. 预检只显示请求；实际下载或核查必须显式指定外部数据目录。
    parser = argparse.ArgumentParser(description='Download ERA5-Land hourly data for the 81 Figure4 sites.')
    parser.add_argument('--site-id', action='append')
    parser.add_argument('--output-dir', type=Path, help='External directory for raw, hourly CSV and QC files.')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--workers', type=int, default=1, choices=[1, 2, 3, 4])
    args = parser.parse_args()
    sites = prepare()
    if args.site_id:
        unknown = set(args.site_id) - set(sites.selected_site_id)
        if unknown:
            raise ValueError(f'Unknown sites: {unknown}')
        sites = sites[sites.selected_site_id.isin(args.site_id)]
    if args.dry_run:
        print(json.dumps({
            'dataset': DATASET,
            'site_count': len(sites),
            'expected_hours_per_site': len(EXPECTED),
            'coordinates_source': str(SOURCE),
            'requests': {row.selected_site_id: request_for(row) for row in sites.itertuples()},
        }, indent=2))
        return
    if args.output_dir is None:
        parser.error('--output-dir is required unless --dry-run is used')
    ROOT = args.output_dir.resolve()
    figure_directory = SOURCE.parent.parent.resolve()
    if ROOT == figure_directory or figure_directory in ROOT.parents:
        parser.error('--output-dir must be outside the Figure4 archive')
    for name in ['raw', 'hourly_csv', 'metadata', 'qc', 'logs']:
        (ROOT / name).mkdir(parents=True, exist_ok=True)
    provenance = ROOT / 'metadata/input_provenance.json'
    if provenance.exists():
        previous = json.loads(provenance.read_text(encoding='utf-8'))
        if previous['sha256'] != sha256(SOURCE):
            raise ValueError('Site configuration changed since download initialization')
    write_json(provenance, {
        'source': str(SOURCE),
        'sha256': sha256(SOURCE),
    })
    all_sites = sites.copy()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.FileHandler(ROOT / 'logs' / 'download.log', encoding='utf-8'), logging.StreamHandler()])
    write_json(ROOT / 'metadata' / 'environment.json', {'python': platform.python_version(),
               'packages': {p: importlib.metadata.version(p) for p in ['cdsapi', 'pandas', 'numpy', 'xarray', 'h5netcdf', 'h5py']}})
    write_json(ROOT / 'metadata' / 'download_plan.json', {'dataset': DATASET, 'expected_hours_per_site': len(EXPECTED),
               'requests': {row.selected_site_id: request_for(row) for row in all_sites.itertuples()}})
    outcomes = []
    pool = ThreadPoolExecutor(max_workers=args.workers)
    futures = [pool.submit(run_site, site, args.verify_only) for site in sites.itertuples()]
    try:
        for future in as_completed(futures):
            outcomes.append(future.result())
            pd.DataFrame([{k: v for k, v in q.items() if k != 'source_units'} for q in outcomes]).sort_values(
                'selected_site_id').to_csv(ROOT / 'qc' / 'latest_run_summary.csv', index=False)
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
    logging.info('Completed and validated %d requested sites', len(outcomes))


if __name__ == '__main__':
    main()
