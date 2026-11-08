"""从数据重新绘制左右图，合并矢量输出；不使用 Figure4.png。"""
from pathlib import Path
import argparse
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

import pymupdf
from PIL import Image
import matplotlib.pyplot as plt


# 1. 原始 Python 数据绘图负责左图，原生 R 脚本负责右图。
base = Path(__file__).resolve().parent
root = base.parent
resolution = 300
left_width_mm = 297
right_width_mm = 99
height_mm = 210
# 通过参数或 PATH 定位 Rscript，不保存整理者的本地安装路径。
parser = argparse.ArgumentParser(description='Draw the combined Figure4 from archived inputs.')
parser.add_argument(
    '--rscript',
    default='Rscript',
    help='Rscript executable name or path; default: Rscript on PATH.',
)
args = parser.parse_args()
rscript = shutil.which(args.rscript)
if rscript is None:
    parser.error('Rscript was not found; add it to PATH or supply --rscript.')
import plot_Figure4 as main

subprocess.run(
    [rscript, str(base / 'plot_warming_management_comparison.R')],
    check=True,
)
right_stem = base / 'Supplementary_N2O_warming_management'
output = root / 'Figure4_combined'

with tempfile.TemporaryDirectory(prefix='figure4_vector_') as folder:
    temporary = Path(folder)
    figure = main.draw_figure('totals', return_figure=True)
    for extension in ['pdf', 'svg']:
        figure.savefig(temporary / f'left.{extension}', facecolor='white')
    plt.close(figure)

    # 2. PDF 页面对象保留两侧路径、文本和字体，保持 3:1 的物理宽度。
    factor = 72 / 25.4
    left_width = left_width_mm * factor
    right_width = right_width_mm * factor
    height = height_mm * factor
    document = pymupdf.open()
    page = document.new_page(width=left_width + right_width, height=height)
    with pymupdf.open(temporary / 'left.pdf') as left, pymupdf.open(right_stem.with_suffix('.pdf')) as right:
        page.show_pdf_page(pymupdf.Rect(0, 0, left_width, height), left, 0)
        page.show_pdf_page(pymupdf.Rect(left_width, 0, left_width + right_width, height), right, 0)
    assert not page.get_images(full=True), 'Composite PDF unexpectedly contains raster images'
    assert len(page.get_text().strip()) > 100
    document.save(output.with_suffix('.pdf'), garbage=4, deflate=True)

    # 3. 从同一矢量页面渲染 300 dpi PNG/TIFF。
    pixels = page.get_pixmap(dpi=resolution, alpha=False)
    raster = Image.frombytes('RGB', (pixels.width, pixels.height), pixels.samples)
    raster.save(output.with_suffix('.png'), dpi=(resolution, resolution))
    raster.save(output.with_suffix('.tiff'), dpi=(resolution, resolution), compression='tiff_lzw')
    document.close()

    # 4. SVG 使用嵌套矢量画布；前缀隔离 ID，避免裁切路径重名。
    namespace = 'http://www.w3.org/2000/svg'
    ET.register_namespace('', namespace)
    canvas = ET.Element(f'{{{namespace}}}svg', {
        'width': '396mm', 'height': '210mm',
        'viewBox': f'0 0 {left_width + right_width} {height}',
    })
    for prefix, path, x, width in [
        ('left_', temporary / 'left.svg', 0, left_width),
        ('right_', right_stem.with_suffix('.svg'), left_width, right_width),
    ]:
        text = path.read_text(encoding='utf-8')
        tree = ET.fromstring(text)
        identifiers = [node.attrib['id'] for node in tree.iter() if 'id' in node.attrib]
        for identifier in sorted(identifiers, key=len, reverse=True):
            text = text.replace(f'id="{identifier}"', f'id="{prefix}{identifier}"')
            text = text.replace(f'#{identifier}"', f'#{prefix}{identifier}"')
            text = text.replace(f'#{identifier})', f'#{prefix}{identifier})')
        tree = ET.fromstring(text)
        tree.set('x', str(x))
        tree.set('y', '0')
        tree.set('width', str(width))
        tree.set('height', str(height))
        canvas.append(tree)
    assert not list(canvas.iter(f'{{{namespace}}}image'))
    ET.ElementTree(canvas).write(output.with_suffix('.svg'), encoding='utf-8', xml_declaration=True)

print('PASS: both panels redrawn from data; vector PDF/SVG; 300 dpi PNG/TIFF; width ratio 3:1.')
