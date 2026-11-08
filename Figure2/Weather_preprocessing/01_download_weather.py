"""Download planned ERA5-Land point time series with visible progress."""

from __future__ import annotations

import argparse
import os
import random
import sys
import tempfile
import threading
import time
import zipfile
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


# 1. 使用已核验的最终任务表；大型下载文件和状态日志必须写到外部目录。
MODULE_DIR = Path(__file__).resolve().parent
TASKS_CSV = MODULE_DIR / "download_tasks.csv"
STATUS_DIR: Path
DATASET = "reanalysis-era5-land-timeseries"
VARIABLES = ["2m_temperature", "total_precipitation"]
PRINT_LOCK = threading.Lock()
STATUS_WARNING_SHOWN = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download the prepared 0.1-degree ERA5-Land point time-series tasks."
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-workers", type=int, default=2)
    parser.add_argument("--retries", type=int, default=6)
    parser.add_argument("--test-n-tasks", type=int)
    parser.add_argument(
        "--task-ids",
        help="Optional comma-separated task IDs to download or validate.",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate the plan and existing files without contacting CDS.",
    )
    return parser.parse_args()


def validate_csv(
    path: Path,
    request_start: pd.Timestamp,
    request_end: pd.Timestamp,
    expected_lat: float,
    expected_lon: float,
) -> tuple[bool, str]:
    if not path.exists() or path.stat().st_size < 100:
        return False, "missing_or_empty"
    try:
        data = pd.read_csv(path)
        required = {"valid_time", "t2m", "tp", "latitude", "longitude"}
        missing = required - set(data.columns)
        if missing:
            return False, "missing_columns:" + ",".join(sorted(missing))
        time_utc = pd.to_datetime(data["valid_time"], errors="raise", utc=True)
        expected = pd.date_range(
            request_start,
            request_end + pd.Timedelta(days=1),
            freq="h",
            inclusive="left",
            tz="UTC",
        )
        if len(data) != len(expected):
            return False, f"row_count:{len(data)}_expected:{len(expected)}"
        if not time_utc.is_unique or set(time_utc) != set(expected):
            return False, "timestamps_incomplete_or_duplicated"
        if not np.isfinite(data[["t2m", "tp", "latitude", "longitude"]].apply(pd.to_numeric, errors="coerce").to_numpy()).all():
            return False, "missing_weather_values"
        lat = pd.to_numeric(data["latitude"], errors="coerce")
        lon = ((pd.to_numeric(data["longitude"], errors="coerce") + 180) % 360) - 180
        if lat.sub(expected_lat).abs().max() > 1e-6 or lon.sub(expected_lon).abs().max() > 1e-6:
            return False, "returned_coordinates_do_not_match_download_grid"
        return True, "valid"
    except Exception as exc:
        return False, f"read_error:{type(exc).__name__}:{exc}"


def normalize_download_table(data: pd.DataFrame) -> pd.DataFrame:
    """Normalize one ERA5-Land CSV group to the final audit schema."""
    aliases = {
        "time": "valid_time",
        "2m_temperature": "t2m",
        "total_precipitation": "tp",
    }
    data = data.rename(columns={key: value for key, value in aliases.items() if key in data.columns})
    keys = ["valid_time", "latitude", "longitude"]
    missing_keys = set(keys) - set(data.columns)
    if missing_keys:
        raise ValueError(f"Missing key columns: {sorted(missing_keys)}")
    variables = [name for name in ("t2m", "tp") if name in data.columns]
    if not variables:
        raise ValueError("CSV group contains neither t2m nor tp")
    out = data[keys + variables].copy()
    if out.duplicated(keys).any():
        raise ValueError("CSV group has duplicated time-coordinate rows")
    return out


def merge_download_tables(tables: list[pd.DataFrame]) -> pd.DataFrame:
    """Merge variable-group CSV files returned by the ERA5-Land service."""
    if not tables:
        raise RuntimeError("No CSV tables were returned")
    keys = ["valid_time", "latitude", "longitude"]
    merged = normalize_download_table(tables[0])
    for table in tables[1:]:
        current = normalize_download_table(table)
        duplicate_variables = (set(merged.columns) & set(current.columns)) - set(keys)
        if duplicate_variables:
            current = current.drop(columns=sorted(duplicate_variables))
        if len(current.columns) > len(keys):
            merged = merged.merge(current, on=keys, how="outer", validate="one_to_one")
    return merged.sort_values("valid_time", kind="stable").reset_index(drop=True)


def read_download(path: Path) -> pd.DataFrame:
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            bad = archive.testzip()
            if bad is not None:
                raise RuntimeError(f"Corrupted archive entry: {bad}")
            csv_names = [name for name in archive.namelist() if name.lower().endswith(".csv")]
            if not csv_names:
                raise RuntimeError("No CSV files found in downloaded archive")
            tables = []
            for name in csv_names:
                with archive.open(name) as handle:
                    tables.append(pd.read_csv(handle))
            return merge_download_tables(tables)
    return merge_download_tables([pd.read_csv(path)])


def validate_candidate(data: pd.DataFrame, row: dict) -> None:
    required = {"valid_time", "t2m", "tp", "latitude", "longitude"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    time_utc = pd.to_datetime(data["valid_time"], errors="raise", utc=True)
    expected = pd.date_range(
        row["request_start_utc"],
        row["request_end_utc"] + pd.Timedelta(days=1),
        freq="h",
        inclusive="left",
        tz="UTC",
    )
    if len(data) != len(expected) or not time_utc.is_unique or set(time_utc) != set(expected):
        raise ValueError(
            f"Incomplete UTC hours: rows={len(data)}, expected={len(expected)}, "
            f"unique={time_utc.nunique()}"
        )
    weather = data[["t2m", "tp", "latitude", "longitude"]].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(weather.to_numpy()).all():
        raise ValueError("Missing temperature or precipitation values")
    lat = pd.to_numeric(data["latitude"], errors="coerce")
    lon = ((pd.to_numeric(data["longitude"], errors="coerce") + 180) % 360) - 180
    if (
        lat.sub(row["download_gp_lat"]).abs().max() > 1e-6
        or lon.sub(row["download_gp_lon"]).abs().max() > 1e-6
    ):
        raise ValueError("Returned coordinates do not match the requested ERA5 grid point")


def download_task(row: dict, raw_dir: Path, retries: int) -> dict:
    destination = raw_dir / row["output_filename"]
    valid, reason = validate_csv(
        destination,
        row["request_start_utc"],
        row["request_end_utc"],
        row["download_gp_lat"],
        row["download_gp_lon"],
    )
    if valid:
        return {"task_id": row["task_id"], "status": "skipped_valid", "path": str(destination), "attempts": 0, "message": reason}
    if destination.exists():
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        destination.replace(destination.with_suffix(destination.suffix + f".invalid_{stamp}"))

    request = {
        "variable": VARIABLES,
        "location": {
            "latitude": float(row["download_gp_lat"]),
            "longitude": float(row["download_gp_lon"]),
        },
        "date": [f"{row['request_start_utc']:%Y-%m-%d}/{row['request_end_utc']:%Y-%m-%d}"],
        "data_format": "csv",
    }
    messages = []
    for attempt in range(1, retries + 1):
        part = destination.with_suffix(destination.suffix + ".part")
        try:
            if part.exists():
                part.unlink()
            # 只在真正下载时加载 CDS 客户端；--check-only 不需凭据且不联网。
            import cdsapi
            client = cdsapi.Client(quiet=True)
            client.retrieve(DATASET, request).download(str(part))
            data = read_download(part)
            validate_candidate(data, row)
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".csv", delete=False, dir=raw_dir,
                encoding="utf-8", newline=""
            ) as handle:
                temporary_csv = Path(handle.name)
                data.to_csv(handle, index=False)
            os.replace(temporary_csv, destination)
            part.unlink(missing_ok=True)
            return {"task_id": row["task_id"], "status": "downloaded", "path": str(destination), "attempts": attempt, "message": "validated"}
        except Exception as exc:
            messages.append(f"attempt {attempt}: {type(exc).__name__}: {exc}")
            part.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(min(180, 5 * (2 ** (attempt - 1))) + random.uniform(0, 3))
    return {"task_id": row["task_id"], "status": "failed", "path": str(destination), "attempts": retries, "message": " | ".join(messages)}


def progress(completed: int, total: int, counts: dict[str, int], elapsed: float) -> None:
    width = 36
    fraction = completed / total if total else 1.0
    filled = min(width, int(width * fraction))
    bar = "#" * filled + "-" * (width - filled)
    message = (
        f"\r[{bar}] {completed:>4}/{total:<4} {fraction:6.1%} | "
        f"downloaded {counts.get('downloaded', 0)} | "
        f"skipped {counts.get('skipped_valid', 0)} | "
        f"failed {counts.get('failed', 0)} | elapsed {elapsed/60:.1f} min"
    )
    with PRINT_LOCK:
        print(message, end="", flush=True)


def write_status(results: list[dict]) -> Path:
    global STATUS_WARNING_SHOWN
    target = STATUS_DIR / "download_status.csv"
    table = pd.DataFrame(results)
    if target.exists():
        previous = pd.read_csv(target)
        previous = previous.loc[~previous["task_id"].isin(table["task_id"])]
        table = pd.concat([previous, table], ignore_index=True)
    table = table.sort_values("task_id", kind="stable")
    temporary = STATUS_DIR / f"download_status.{os.getpid()}.tmp"
    table.to_csv(temporary, index=False, encoding="utf-8-sig")
    for attempt in range(1, 11):
        try:
            os.replace(temporary, target)
            return target
        except PermissionError:
            if attempt < 10:
                time.sleep(0.1 * attempt)

    checkpoint = STATUS_DIR / f"download_status_checkpoint_{os.getpid()}.csv"
    table.to_csv(checkpoint, index=False, encoding="utf-8-sig")
    temporary.unlink(missing_ok=True)
    if not STATUS_WARNING_SHOWN:
        with PRINT_LOCK:
            print(
                f"\nWarning: {target.name} is locked; status is being saved to "
                f"{checkpoint.name}. Downloads will continue.",
                flush=True,
            )
        STATUS_WARNING_SHOWN = True
    return checkpoint


def main() -> None:
    global STATUS_DIR
    args = parse_args()
    if args.max_workers < 1 or args.retries < 1:
        raise ValueError("--max-workers and --retries must be positive integers")
    # 2. 最终任务表已包含实际陆地网格和跨年补齐范围，不再依赖主表快照。
    tasks = pd.read_csv(
        TASKS_CSV,
        parse_dates=["window_start", "window_end", "request_start_utc", "request_end_utc"],
    )
    if tasks.task_id.duplicated().any() or tasks.output_filename.duplicated().any():
        raise ValueError("Duplicate task IDs or output filenames in the task table")
    if not np.isfinite(tasks[["download_gp_lat", "download_gp_lon"]].to_numpy()).all():
        raise ValueError("Missing ERA5-Land download coordinates in the task table")
    if tasks[["window_start", "window_end", "request_start_utc", "request_end_utc"]].isna().any().any():
        raise ValueError("Missing dates in the task table")
    if any(Path(name).name != name for name in tasks.output_filename):
        raise ValueError("Task filenames must not contain directory components")
    if len(tasks) != 775:
        raise ValueError("Expected 775 final download tasks")
    if args.task_ids:
        requested_ids = {item.strip() for item in args.task_ids.split(",") if item.strip()}
        unknown_ids = requested_ids - set(tasks["task_id"])
        if unknown_ids:
            raise ValueError(f"Unknown --task-ids: {sorted(unknown_ids)}")
        tasks = tasks.loc[tasks["task_id"].isin(requested_ids)].copy()
    if args.test_n_tasks is not None:
        tasks = tasks.head(args.test_n_tasks).copy()

    raw_dir = args.output_dir / "raw_utc"
    if args.output_dir.resolve().is_relative_to(MODULE_DIR.resolve()):
        raise ValueError("Store large weather downloads outside Weather_preprocessing")
    print("Formal sowing-to-harvest environments: 856")
    print(f"ERA5-Land 0.1-degree point tasks selected: {len(tasks)}")
    print(f"Output: {raw_dir}")
    print(f"Workers: {args.max_workers}; retries: {args.retries}")
    print("Variables: 2m_temperature, total_precipitation")
    fallback_count = int((
        (tasks["download_gp_lat"].sub(tasks["gp_lat"]).abs() > 1e-9)
        | (tasks["download_gp_lon"].sub(tasks["gp_lon"]).abs() > 1e-9)
    ).sum())
    print(f"Nearest-valid-land fallback tasks: {fallback_count}")

    if args.check_only:
        states = []
        for row in tasks.to_dict("records"):
            valid, reason = validate_csv(
                raw_dir / row["output_filename"],
                row["request_start_utc"], row["request_end_utc"],
                row["download_gp_lat"], row["download_gp_lon"],
            )
            states.append({"task_id": row["task_id"], "status": "valid" if valid else "pending", "message": reason})
        status = pd.DataFrame(states)
        print(status.status.value_counts().to_string())
        print("Check only: CDS was not contacted.")
        pending = status[status.status != "valid"]
        if not pending.empty:
            print(pending.head(10).to_string(index=False))
            sys.exit(1)
        return

    # 3. 下载可断点重跑；既有有效文件跳过，状态日志只保存在外部下载目录。
    raw_dir.mkdir(parents=True, exist_ok=True)
    STATUS_DIR = args.output_dir / "logs"
    STATUS_DIR.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    results: list[dict] = []
    counts: dict[str, int] = {}
    status_path = STATUS_DIR / "download_status.csv"
    records = tasks.to_dict("records")
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        pending = {executor.submit(download_task, row, raw_dir, args.retries) for row in records}
        progress(0, len(records), counts, 0)
        while pending:
            completed_now, pending = wait(pending, timeout=1, return_when=FIRST_COMPLETED)
            for future in completed_now:
                result = future.result()
                results.append(result)
                counts[result["status"]] = counts.get(result["status"], 0) + 1
                status_path = write_status(results)
            progress(len(results), len(records), counts, time.monotonic() - started)
    print()
    failed = counts.get("failed", 0)
    print(f"Finished: {len(results) - failed}/{len(results)} valid; failures: {failed}")
    print(f"Status table: {status_path}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
