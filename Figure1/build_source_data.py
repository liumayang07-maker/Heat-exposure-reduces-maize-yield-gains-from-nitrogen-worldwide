"""Build panel-specific source workbooks from the actual plotting scripts."""

import argparse
from pathlib import Path
import subprocess
import tempfile

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill


ROOT = Path(__file__).resolve().parent
def write_workbook(path, sheets):
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=name, index=False)
            sheet = writer.sheets[name]
            sheet.freeze_panes = "A2"
            sheet.auto_filter.ref = sheet.dimensions
            for cell in sheet[1]:
                cell.font = Font(name="Arial", bold=True, color="FFFFFF")
                cell.fill = PatternFill("solid", fgColor="245A64")
            for row in sheet.iter_rows(min_row=2):
                for cell in row:
                    cell.font = Font(name="Arial", size=10)
                    cell.alignment = Alignment(vertical="top", wrap_text=True)
                    if isinstance(cell.value, float):
                        field = sheet.cell(1, cell.column).value
                        cell.number_format = ("0.000000000000E+00" if field in
                                              ("p_value", "p_BH_six_continents", "p_raw", "p_BH")
                                              else "0.###############")
            for column in sheet.columns:
                width = min(85, max(16, max(len(str(c.value or "")) for c in column) + 2))
                sheet.column_dimensions[column[0].column_letter].width = width
    # Check the actual serialized workbook, not only the in-memory frames.
    loaded = pd.read_excel(path, sheet_name=None, keep_default_na=False)
    assert list(loaded) == list(sheets)
    for name, frame in sheets.items():
        pd.testing.assert_frame_equal(loaded[name], frame, check_dtype=False,
                                      check_exact=False, rtol=1e-12, atol=1e-12)


def panel_b(folder, export):
    names = ["Annual_data", "Plot_data", "GLS_models", "Model_parameters",
             "Coefficient_covariance"]
    tables = {name: pd.read_csv(export / (name + ".csv")) for name in names}
    for sheet, file in [("Annual_data", "continent_annual_EDH30_1986_2025.csv"),
                        ("GLS_models", "GLS_models.csv")]:
        pd.testing.assert_frame_equal(tables[sheet], pd.read_csv(folder / file),
                                      check_exact=False, rtol=1e-12, atol=1e-12)
    assert len(tables["Plot_data"]) == 240 and len(tables["Model_parameters"]) == 6
    return tables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    targets = [ROOT / "Figure1b" / "Figure1b_Source_Data.xlsx"]
    if not args.overwrite and any(path.exists() for path in targets):
        parser.error("Source workbooks exist; use --overwrite explicitly.")
    with tempfile.TemporaryDirectory(prefix="figure1_source_") as temporary:
        temporary = Path(temporary)
        for target in targets:
            panel = "Figure1b"
            folder = ROOT / panel
            export = temporary / panel
            export.mkdir()
            subprocess.run([args.rscript, str(folder / ("plot_" + panel + ".R")),
                            "--source-data-dir=" + str(export)], check=True)
            write_workbook(temporary / target.name, panel_b(folder, export))
        for target in targets:
            target.write_bytes((temporary / target.name).read_bytes())
            print(f"Verified: {target.name} ({target.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
