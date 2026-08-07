#!/usr/bin/env python3
"""Extract and plot candidate-gene PHx time courses from GSE95135."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


MATRIX_PATTERN = "GSE95135_Rib_et_al.RPKM_log2*.csv.gz"
SOFT_NAME = "GSE95135_family.soft.gz"
GENE_SPECS = {
    "Adcy7": {"aliases": {"Adcy7"}, "ensembl": "ENSMUSG00000031659"},
    "Piezo1": {"aliases": {"Piezo1", "Fam38a"}, "ensembl": "ENSMUSG00000014444"},
    "Cxcl2": {"aliases": {"Cxcl2"}, "ensembl": "ENSMUSG00000058427"},
    "Ifrd1": {"aliases": {"Ifrd1"}, "ensembl": "ENSMUSG00000001627"},
}


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def find_matrix(raw_dir: Path) -> Path:
    candidates = sorted(raw_dir.glob(MATRIX_PATTERN), key=lambda p: p.stat().st_size, reverse=True)
    if not candidates:
        raise FileNotFoundError(f"No matrix matching {MATRIX_PATTERN!r} in {raw_dir}")
    return candidates[0]


def extract_gene_rows(matrix_path: Path, requested: list[str]) -> tuple[list[str], dict[str, list[str]]]:
    aliases = {
        gene: {alias.casefold() for alias in GENE_SPECS[gene]["aliases"]}
        for gene in requested
    }
    found: dict[str, list[str]] = {}
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        for row in reader:
            if len(row) < 3:
                continue
            matrix_symbol = row[1].casefold()
            for current_symbol, accepted_aliases in aliases.items():
                expected_ensembl = GENE_SPECS[current_symbol]["ensembl"]
                if matrix_symbol in accepted_aliases or row[0] == expected_ensembl:
                    if len(row) != len(header):
                        raise ValueError(
                            f"The {current_symbol} row has {len(row)} columns; expected {len(header)}."
                        )
                    found[current_symbol] = row
            if len(found) == len(requested):
                break
    missing = sorted(set(requested) - set(found))
    if missing:
        raise ValueError(f"Genes not found in matrix: {', '.join(missing)}")
    return header, found


def parse_matrix_sample(sample_code: str) -> dict[str, object]:
    match = re.fullmatch(r"([CSX])(\d+)([HW])(\d+)", sample_code)
    if match is None:
        raise ValueError(f"Unrecognized matrix sample code: {sample_code}")
    prefix, time_value, unit, replicate = match.groups()
    group = {"C": "Control", "S": "Sham", "X": "Post-PH"}[prefix]
    time_value_int = int(time_value)
    time_hours = time_value_int if unit == "H" else time_value_int * 168
    time_label = f"{time_value_int} h" if unit == "H" else f"{time_value_int} wk"
    return {
        "sample_code": sample_code,
        "group": group,
        "time_value": time_value_int,
        "time_unit": "hour" if unit == "H" else "week",
        "time_label": time_label,
        "time_hours": time_hours,
        "replicate": int(replicate),
    }


def build_expression_table(header: list[str], rows: dict[str, list[str]]) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for current_symbol, row in rows.items():
        for sample_code, value in zip(header[3:], row[3:]):
            record = parse_matrix_sample(sample_code)
            record.update(
                {
                    "gene_id": row[0],
                    "gene_symbol_current": current_symbol,
                    "gene_symbol_matrix": row[1],
                    "gene_type": row[2],
                    "log2_rpkm": float(value),
                }
            )
            records.append(record)
    return pd.DataFrame.from_records(records)


def parse_geo_catalog(soft_path: Path) -> pd.DataFrame:
    if not soft_path.exists():
        return pd.DataFrame(columns=["gsm", "title"])
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)
                current = {"gsm": line.split("=", 1)[1].strip()}
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
        if current is not None:
            records.append(current)
    return pd.DataFrame.from_records(records)


def summarize(expression: pd.DataFrame) -> pd.DataFrame:
    summary = (
        expression.groupby(
            ["gene_symbol_current", "gene_symbol_matrix", "gene_id", "group", "time_hours", "time_label"],
            as_index=False,
        )
        .agg(n=("log2_rpkm", "size"), mean_log2_rpkm=("log2_rpkm", "mean"), sd_log2_rpkm=("log2_rpkm", "std"))
        .sort_values(["gene_symbol_current", "group", "time_hours"])
    )
    summary["sem_log2_rpkm"] = summary["sd_log2_rpkm"] / np.sqrt(summary["n"])
    baseline = (
        summary[(summary["group"] == "Control") & (summary["time_hours"] == 0)]
        .set_index("gene_symbol_current")["mean_log2_rpkm"]
    )
    summary["delta_from_control_0h"] = summary.apply(
        lambda row: row["mean_log2_rpkm"] - baseline[row["gene_symbol_current"]], axis=1
    )
    return summary


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
        Path(
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
            if bold
            else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
        ),
    ]
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size=size)
    return ImageFont.load_default()


def save_png_pdf(image: Image.Image, output_png: Path, output_pdf: Path) -> None:
    image.convert("RGB").save(output_png, dpi=(300, 300), optimize=True)
    pdf = canvas.Canvas(str(output_pdf), pagesize=(9.2 * 72, 5.6 * 72))
    pdf.drawImage(ImageReader(image.convert("RGB")), 0, 0, width=9.2 * 72, height=5.6 * 72)
    pdf.showPage()
    pdf.save()


def new_plot_canvas(title: str, y_min: float, y_max: float):
    width, height = 2760, 1680
    left, right, top, bottom = 250, 110, 190, 330
    plot_left, plot_right = left, width - right
    plot_top, plot_bottom = top, height - bottom
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(56, bold=True)
    axis_font = load_font(36)
    tick_font = load_font(30)
    note_font = load_font(26)
    legend_font = load_font(29)
    draw.text((plot_left, 62), title, fill=(20, 28, 39), font=title_font)

    def x_px(value: float) -> float:
        return plot_left + (value / 72.0) * (plot_right - plot_left)

    def y_px(value: float) -> float:
        return plot_bottom - ((value - y_min) / (y_max - y_min)) * (plot_bottom - plot_top)

    y_step = 1.0 if y_max - y_min <= 9 else 2.0
    y_ticks = np.arange(math.ceil(y_min / y_step) * y_step, y_max + 0.001, y_step)
    for y_tick in y_ticks:
        y_coord = y_px(float(y_tick))
        draw.line((plot_left, y_coord, plot_right, y_coord), fill=(221, 225, 230), width=2)
        label = f"{y_tick:.0f}"
        box = draw.textbbox((0, 0), label, font=tick_font)
        draw.text((plot_left - 28 - (box[2] - box[0]), y_coord - 17), label, fill=(75, 85, 99), font=tick_font)
    draw.line((plot_left, plot_top, plot_left, plot_bottom), fill=(45, 55, 72), width=4)
    draw.line((plot_left, plot_bottom, plot_right, plot_bottom), fill=(45, 55, 72), width=4)
    ticks = [0, 1, 4, 10, 20, 28, 36, 44, 48, 72]
    for tick in ticks:
        x_coord = x_px(tick)
        draw.line((x_coord, plot_bottom, x_coord, plot_bottom + 12), fill=(45, 55, 72), width=3)
        label = str(tick)
        box = draw.textbbox((0, 0), label, font=tick_font)
        draw.text((x_coord - (box[2] - box[0]) / 2, plot_bottom + 22), label, fill=(55, 65, 81), font=tick_font)
    return image, draw, x_px, y_px, (width, height, plot_left, plot_right, plot_top, plot_bottom), (axis_font, note_font, legend_font)


def finish_axes(image, draw, geometry, fonts, y_label: str, note: str) -> None:
    width, height, plot_left, _plot_right, _plot_top, plot_bottom = geometry
    axis_font, note_font, _legend_font = fonts
    x_label = "Time after surgery (hours)"
    x_box = draw.textbbox((0, 0), x_label, font=axis_font)
    draw.text(((width - (x_box[2] - x_box[0])) / 2, plot_bottom + 78), x_label, fill=(35, 45, 60), font=axis_font)
    label_layer = Image.new("RGBA", (1000, 80), (255, 255, 255, 0))
    ImageDraw.Draw(label_layer).text((0, 10), y_label, fill=(35, 45, 60), font=axis_font)
    rotated = label_layer.rotate(90, expand=True)
    image.alpha_composite(rotated, (32, int((height - rotated.height) / 2 - 50)))
    draw = ImageDraw.Draw(image)
    draw.multiline_text((plot_left, height - 205), note, fill=(75, 85, 99), font=note_font, spacing=10)


def plot_timecourse(expression: pd.DataFrame, gene: str, output_png: Path, output_pdf: Path) -> None:
    rng = np.random.default_rng(7)
    colors = {"Post-PH": (209, 73, 91), "Sham": (39, 125, 161), "Control": (107, 114, 128)}
    baseline = expression[(expression["group"] == "Control") & (expression["time_hours"] == 0)]
    baseline_mean, baseline_sem = baseline["log2_rpkm"].mean(), baseline["log2_rpkm"].sem()
    focus = expression[
        ((expression["group"] == "Control") & (expression["time_hours"] == 0))
        | (expression["group"].isin(["Post-PH", "Sham"]) & expression["time_hours"].between(1, 72))
    ]
    y_min = math.floor((focus["log2_rpkm"].min() - 0.35) * 2) / 2
    y_max = math.ceil((focus["log2_rpkm"].max() + 0.35) * 2) / 2
    image, draw, x_px, y_px, geometry, fonts = new_plot_canvas(
        f"{gene} expression during mouse liver regeneration (GSE95135)", y_min, y_max
    )
    overlay = Image.new("RGBA", image.size, (255, 255, 255, 0))
    overlay_draw = ImageDraw.Draw(overlay)

    def point(x_value, y_value, color, radius=12):
        x_coord, y_coord = x_px(float(x_value)), y_px(float(y_value))
        draw.ellipse((x_coord-radius, y_coord-radius, x_coord+radius, y_coord+radius), fill=color, outline="white", width=2)

    for jitter, value in zip(rng.uniform(-0.22, 0.22, len(baseline)), baseline["log2_rpkm"]):
        point(jitter, value, colors["Control"])
    for group in ["Post-PH", "Sham"]:
        subset = expression[(expression["group"] == group) & expression["time_hours"].between(1, 72)].copy()
        for time_hours, time_frame in subset.groupby("time_hours"):
            for offset, value in zip(rng.uniform(-0.28, 0.28, len(time_frame)), time_frame["log2_rpkm"]):
                point(time_hours + offset, value, colors[group])
        group_summary = subset.groupby("time_hours", as_index=False).agg(mean=("log2_rpkm", "mean"), sem=("log2_rpkm", "sem")).sort_values("time_hours")
        line_x = np.concatenate(([0.0], group_summary["time_hours"].to_numpy(float)))
        line_y = np.concatenate(([baseline_mean], group_summary["mean"].to_numpy(float)))
        line_sem = np.concatenate(([baseline_sem], group_summary["sem"].to_numpy(float)))
        upper = [(x_px(x), y_px(y+e)) for x, y, e in zip(line_x, line_y, line_sem)]
        lower = [(x_px(x), y_px(y-e)) for x, y, e in zip(line_x[::-1], line_y[::-1], line_sem[::-1])]
        overlay_draw.polygon(upper + lower, fill=colors[group] + (34,))
        coords = [(x_px(x), y_px(y)) for x, y in zip(line_x, line_y)]
        draw.line(coords, fill=colors[group], width=8, joint="curve")
        for x_coord, y_coord in coords:
            draw.ellipse((x_coord-14, y_coord-14, x_coord+14, y_coord+14), fill=colors[group], outline="white", width=3)
    image = Image.alpha_composite(image.convert("RGBA"), overlay)
    draw = ImageDraw.Draw(image)
    legend_font = fonts[2]
    legend_x, legend_y = geometry[3] - 930, 132
    for label, color in [("Control 0 h samples", colors["Control"]), ("Post-PH mean", colors["Post-PH"]), ("Sham mean", colors["Sham"])]:
        draw.line((legend_x, legend_y+14, legend_x+52, legend_y+14), fill=color, width=7)
        draw.ellipse((legend_x+18, legend_y, legend_x+46, legend_y+28), fill=color, outline="white", width=2)
        draw.text((legend_x+68, legend_y-4), label, fill=(55, 65, 81), font=legend_font)
        legend_x += 310
    note = (
        "Points are individual samples; lines are means; shaded bands are SEM.  Untreated Control 0 h is the shared baseline.\n"
        "Sham was measured at 1, 4, 10, 20 and 48 h; missing Sham time points were not imputed."
    )
    finish_axes(image, draw, geometry, fonts, f"{gene} expression (log2 RPKM)", note)
    save_png_pdf(image, output_png, output_pdf)


def plot_comparison(summary: pd.DataFrame, output_png: Path, output_pdf: Path) -> None:
    colors = {
        "Adcy7": (110, 84, 148),
        "Piezo1": (42, 157, 143),
        "Cxcl2": (242, 142, 43),
        "Ifrd1": (209, 73, 91),
    }
    focus = summary[(summary["group"] == "Post-PH") & summary["time_hours"].between(1, 72)].copy()
    y_min = math.floor(min(-0.5, focus["delta_from_control_0h"].min() - 0.4))
    y_max = math.ceil(max(0.5, focus["delta_from_control_0h"].max() + 0.4))
    image, draw, x_px, y_px, geometry, fonts = new_plot_canvas(
        "Candidate-gene responses after partial hepatectomy (GSE95135)", y_min, y_max
    )
    if y_min < 0 < y_max:
        draw.line((geometry[2], y_px(0), geometry[3], y_px(0)), fill=(90, 100, 115), width=4)
    legend_font = fonts[2]
    legend_x, legend_y = geometry[3] - 915, 132
    for gene in ["Adcy7", "Piezo1", "Cxcl2", "Ifrd1"]:
        gene_frame = focus[focus["gene_symbol_current"] == gene].sort_values("time_hours")
        line_x = np.concatenate(([0.0], gene_frame["time_hours"].to_numpy(float)))
        line_y = np.concatenate(([0.0], gene_frame["delta_from_control_0h"].to_numpy(float)))
        coords = [(x_px(x), y_px(y)) for x, y in zip(line_x, line_y)]
        draw.line(coords, fill=colors[gene], width=9, joint="curve")
        for x_coord, y_coord in coords:
            draw.ellipse((x_coord-14, y_coord-14, x_coord+14, y_coord+14), fill=colors[gene], outline="white", width=3)
        draw.line((legend_x, legend_y+14, legend_x+52, legend_y+14), fill=colors[gene], width=8)
        draw.text((legend_x+65, legend_y-4), gene, fill=(55, 65, 81), font=legend_font)
        legend_x += 220
    note = (
        "Each curve is the Post-PH group mean minus that gene's untreated Control 0 h mean.\n"
        "Positive values indicate higher expression than baseline; absolute expression levels differ among genes."
    )
    image = image.convert("RGBA")
    finish_axes(image, ImageDraw.Draw(image), geometry, fonts, "Change from Control 0 h (log2 RPKM)", note)
    save_png_pdf(image, output_png, output_pdf)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=project_root())
    parser.add_argument("--genes", nargs="+", default=list(GENE_SPECS))
    args = parser.parse_args()
    unknown = sorted(set(args.genes) - set(GENE_SPECS))
    if unknown:
        raise ValueError(f"Unknown configured genes: {', '.join(unknown)}")

    root = args.project_root.resolve()
    raw_dir = root / "datasets" / "GSE95135" / "raw"
    processed_dir = root / "datasets" / "GSE95135" / "processed"
    metadata_dir = root / "datasets" / "GSE95135" / "metadata"
    figure_dir = root / "figures" / "GSE95135"
    summary_dir = root / "results" / "summary" / "GSE95135"
    for directory in [processed_dir, metadata_dir, figure_dir, summary_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    matrix_path = find_matrix(raw_dir)
    header, rows = extract_gene_rows(matrix_path, args.genes)
    expression = build_expression_table(header, rows)
    summary = summarize(expression)
    expression.to_csv(processed_dir / "candidate_genes_log2RPKM_long.tsv", sep="\t", index=False)
    expression.drop(columns=["log2_rpkm", "gene_id", "gene_symbol_current", "gene_symbol_matrix", "gene_type"]).drop_duplicates().to_csv(
        metadata_dir / "matrix_sample_metadata.tsv", sep="\t", index=False
    )
    parse_geo_catalog(raw_dir / SOFT_NAME).to_csv(metadata_dir / "geo_sample_catalog.tsv", sep="\t", index=False)
    summary.to_csv(summary_dir / "candidate_genes_timecourse_summary.tsv", sep="\t", index=False)

    for gene in args.genes:
        gene_expression = expression[expression["gene_symbol_current"] == gene].copy()
        gene_expression.to_csv(processed_dir / f"{gene}_log2RPKM_long.tsv", sep="\t", index=False)
        summary[summary["gene_symbol_current"] == gene].to_csv(
            summary_dir / f"{gene}_timecourse_summary.tsv", sep="\t", index=False
        )
        plot_timecourse(
            gene_expression,
            gene,
            figure_dir / f"{gene}_PHx_vs_sham_0_72h.png",
            figure_dir / f"{gene}_PHx_vs_sham_0_72h.pdf",
        )

    plot_comparison(
        summary,
        figure_dir / "candidate_genes_Post-PH_delta_from_baseline_0_72h.png",
        figure_dir / "candidate_genes_Post-PH_delta_from_baseline_0_72h.pdf",
    )
    print(f"Matrix: {matrix_path}")
    for gene in args.genes:
        row = rows[gene]
        print(f"{gene}: matrix symbol {row[1]}, {row[0]}")
    print(f"Expression rows: {len(expression)}")
    print(f"Summary tables: {summary_dir}")
    print(f"Figures: {figure_dir}")


if __name__ == "__main__":
    main()
