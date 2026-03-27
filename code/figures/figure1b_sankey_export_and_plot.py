#!/usr/bin/env python3
# Public release script for the climate-AMR manuscript

"""Build Figure 1B Sankey outputs from the public source-data table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

TOTAL_STUDIES = 185
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_CSV = PROJECT_ROOT / "figure_source_data" / "source_data_figure1b_sankey_assignments.csv"
DEFAULT_OUTPUT_CSV = PROJECT_ROOT / "generated_figures" / "Figure1B_Sankey_assignments_supplement_ready.csv"
DEFAULT_OUTPUT_XLSX = PROJECT_ROOT / "generated_figures" / "Figure1B_Sankey_assignments_supplement_ready.xlsx"
DEFAULT_OUTPUT_HTML = PROJECT_ROOT / "generated_figures" / "Figure1B_sankey_diagram_185studies.html"
DEFAULT_OUTPUT_PDF = PROJECT_ROOT / "generated_figures" / "Figure1B_sankey_diagram_185studies.pdf"

DESIGN_ORDER = [
    "Cross-sectional study",
    "Time-series study",
    "Ecological study",
    "Modeling study",
    "Laboratory experiment",
    "Cohort study",
]

CLIMATE_ORDER = [
    "Ambient temperature",
    "Seasonality",
    "Precipitation",
    "Extreme weather events",
    "Relative humidity",
    "Interaction/Composite exposure",
]

SAMPLE_ORDER = [
    "Clinical isolates",
    "Environmental water",
    "Wastewater",
    "Soil/agricultural sources",
    "Laboratory isolates",
    "Animal sources",
    "Air/aerosol",
]

COLORS = {
    "Cross-sectional study": "#91CC75",
    "Time-series study": "#5470C6",
    "Ecological study": "#FAC858",
    "Modeling study": "#73C0DE",
    "Laboratory experiment": "#EE6666",
    "Cohort study": "#9A60B4",
    "Ambient temperature": "#F28E2B",
    "Seasonality": "#4E79A7",
    "Precipitation": "#59A14F",
    "Extreme weather events": "#E15759",
    "Relative humidity": "#76B7B2",
    "Interaction/Composite exposure": "#D62728",
    "Clinical isolates": "#F28E8B",
    "Environmental water": "#5DA5DA",
    "Wastewater": "#7488D8",
    "Soil/agricultural sources": "#9C755F",
    "Laboratory isolates": "#F1CE63",
    "Animal sources": "#8E63B6",
    "Air/aerosol": "#B07AA1",
}

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT_CSV)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT_CSV)
    parser.add_argument("--output-xlsx", type=Path, default=DEFAULT_OUTPUT_XLSX)
    parser.add_argument("--output-html", type=Path, default=DEFAULT_OUTPUT_HTML)
    parser.add_argument("--output-pdf", type=Path, default=DEFAULT_OUTPUT_PDF)
    return parser.parse_args()

def percent_string(count: int, total: int = TOTAL_STUDIES) -> str:
    return f"{(count / total) * 100:.2f}%"

def split_climate_factors(value: str) -> list[str]:
    if pd.isna(value) or not str(value).strip():
        return []
    return [part.strip() for part in str(value).split(";") if part.strip()]

def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group_name, order, column in [
        ("Study Design", DESIGN_ORDER, "study_design"),
        ("Sample Source", SAMPLE_ORDER, "sample_source"),
    ]:
        for item in order:
            count = int((df[column] == item).sum())
            rows.append(
                {
                    "category_group": group_name,
                    "label": item,
                    "count_n": count,
                    "percentage": round((count / TOTAL_STUDIES) * 100, 2),
                    "figure_1b_label": f"{item}\nn={count} ({percent_string(count)})",
                }
            )

    expanded = df.assign(climate_factor=df["climate_factors"].map(split_climate_factors)).explode("climate_factor")
    for item in CLIMATE_ORDER:
        count = int((expanded["climate_factor"] == item).sum())
        rows.append(
            {
                "category_group": "Climate Factor",
                "label": item,
                "count_n": count,
                "percentage": round((count / TOTAL_STUDIES) * 100, 2),
                "figure_1b_label": f"{item}\nn={count} ({percent_string(count)})",
            }
        )

    return pd.DataFrame(rows)

def build_links(df: pd.DataFrame) -> tuple[list[str], list[int], list[int], list[int]]:
    expanded = df.assign(climate_factor=df["climate_factors"].map(split_climate_factors)).explode("climate_factor")

    node_order = DESIGN_ORDER + CLIMATE_ORDER + SAMPLE_ORDER
    node_index = {name: idx for idx, name in enumerate(node_order)}
    sources: list[int] = []
    targets: list[int] = []
    values: list[int] = []

    left_counts = (
        expanded.groupby(["study_design", "climate_factor"]).size().reset_index(name="n")
    )
    right_counts = (
        expanded.groupby(["climate_factor", "sample_source"]).size().reset_index(name="n")
    )

    for row in left_counts.itertuples(index=False):
        sources.append(node_index[row.study_design])
        targets.append(node_index[row.climate_factor])
        values.append(int(row.n))

    for row in right_counts.itertuples(index=False):
        sources.append(node_index[row.climate_factor])
        targets.append(node_index[row.sample_source])
        values.append(int(row.n))

    return node_order, sources, targets, values

def build_node_labels(summary: pd.DataFrame, node_order: list[str]) -> list[str]:
    summary_map = summary.set_index("label")[["count_n", "percentage"]].to_dict("index")
    labels = []
    for name in node_order:
        entry = summary_map[name]
        labels.append(f"{name}<br>n={entry['count_n']} ({entry['percentage']:.2f}%)")
    return labels

def write_outputs(df: pd.DataFrame, summary: pd.DataFrame, args: argparse.Namespace) -> None:
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_xlsx.parent.mkdir(parents=True, exist_ok=True)
    args.output_html.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    with pd.ExcelWriter(args.output_xlsx, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="Supplement_Ready_Table", index=False)
        summary.to_excel(writer, sheet_name="Category_Summary", index=False)

def make_figure(df: pd.DataFrame, summary: pd.DataFrame, args: argparse.Namespace) -> None:
    node_order, sources, targets, values = build_links(df)
    node_labels = build_node_labels(summary, node_order)
    node_colors = [COLORS[name] for name in node_order]

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="fixed",
                node=dict(label=node_labels, color=node_colors, pad=18, thickness=20),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=[COLORS[node_order[s]] for s in sources],
                ),
            )
        ]
    )
    fig.update_layout(width=1400, height=900, font=dict(family="Times New Roman", size=16))
    fig.write_html(args.output_html)
    try:
        fig.write_image(args.output_pdf)
    except Exception:
        pass

def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    df = df.sort_values("study_id").reset_index(drop=True)
    summary = build_summary(df)
    write_outputs(df, summary, args)
    make_figure(df, summary, args)

if __name__ == "__main__":
    main()