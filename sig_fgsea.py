import argparse
import csv
import math
import os
import sys
from pathlib import Path


REQUIRED_COLUMNS = ("pathway", "padj", "NES")
MANIFEST_FIELDS = ("contrast_type", "object_name", "group_by", "ident1", "ident2")
EXTRA_FIELDS = (
    "contrast_id",
    "source_layout",
    "source_fgsea_csv",
    "contrast_type",
    "object_name",
    "group_by",
    "ident1",
    "ident2",
    "direction"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Export significant FGSEA pathways from existing hallmark_fgsea.csv files "
            "without rerunning the full analysis pipeline."
        )
    )
    parser.add_argument(
        "--input-root",
        default="results/gsea",
        help="Directory containing current and legacy FGSEA CSV outputs."
    )
    parser.add_argument(
        "--padj-threshold",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff for significance (strictly less than this value)."
    )
    parser.add_argument(
        "--summary-output",
        default="results/gsea/_summary/all_significant_pathways.csv",
        help="Path for the combined significant-pathways CSV."
    )
    return parser.parse_args()


def load_manifest(manifest_path):
    if not manifest_path.exists():
        return {}

    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = {}
        for row in reader:
            contrast_id = row.get("contrast_id", "")
            if contrast_id:
                rows[contrast_id] = row
        return rows


def display_path(path):
    return os.path.relpath(path, Path.cwd())


def discover_sources(input_root):
    sources = []

    for csv_path in sorted(input_root.glob("*/hallmark_fgsea.csv")):
        contrast_id = csv_path.parent.name
        output_path = csv_path.parent/"significant_pathways.csv"
        sources.append({"contrast_id": contrast_id,
                        "source_layout": "current",
                        "input_path": csv_path,
                        "output_path": output_path})

    for csv_path in sorted(input_root.glob("*_hallmark_fgsea.csv")):
        contrast_id = csv_path.name[:-len("_hallmark_fgsea.csv")]
        output_path = csv_path.with_name(f"{contrast_id}_significant_pathways.csv")
        sources.append({"contrast_id": contrast_id,
                        "source_layout": "legacy",
                        "input_path": csv_path,
                        "output_path": output_path})

    return sources


def parse_float(value, field_name, csv_path):
    cleaned = "" if value is None else str(value).strip()
    if cleaned == "" or cleaned.upper() == "NA":
        return math.nan

    return float(cleaned)

def validate_columns(fieldnames, csv_path):
    if fieldnames is None:
        raise SystemExit(f"{csv_path}: missing CSV header row")

    missing = []
    for column in REQUIRED_COLUMNS:
        if column not in fieldnames:
            missing.append(column)
    if missing:
        missing_str = ", ".join(missing)
        raise SystemExit(f"{csv_path}: missing required columns: {missing_str}")

    return fieldnames


def direction_for_row(nes, metadata):
    ident1 = metadata.get("ident1", "")
    ident2= metadata.get("ident2", "")
    if not ident1 or not ident2 or math.isnan(nes):
        return ""
    return ident1 if nes > 0 else ident2


def sort_key(row):
    padj = parse_float(row.get("padj", ""), "padj", Path(row.get("source_fgsea_csv", "<row>")))
    nes = parse_float(row.get("NES", ""), "NES", Path(row.get("source_fgsea_csv", "<row>")))
    contrast_id = row.get("contrast_id", "")
    pathway = row.get("pathway", "")
    return (padj, -abs(nes), contrast_id, pathway)


def build_output_row(base_row, contrast_id, source_layout, source_fgsea_csv, metadata):
    nes = parse_float(base_row.get("NES", ""), "NES", Path(source_fgsea_csv))
    manifest_values = {}
    for field in MANIFEST_FIELDS:
        manifest_values[field] = metadata.get(field, "")

    output_row = {
        "contrast_id": contrast_id,
        "source_layout": source_layout,
        "source_fgsea_csv": source_fgsea_csv,
        **manifest_values,
        "direction": direction_for_row(nes, metadata)
    }
    output_row.update(base_row)
    return output_row


def write_rows(output_path, fieldnames, rows):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def process_source(source, manifest_rows, padj_threshold):
    csv_path = Path(source["input_path"])
    contrast_id = str(source["contrast_id"])
    source_layout = str(source["source_layout"])
    output_path = Path(source["output_path"])
    metadata = manifest_rows.get(contrast_id, {})

    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        original_fieldnames = validate_columns(reader.fieldnames, csv_path)

        total_rows = 0
        significant_rows = []
        source_fgsea_csv = display_path(csv_path)

        for row in reader:
            total_rows += 1
            padj = parse_float(row.get("padj", ""), "padj", csv_path)
            if math.isnan(padj) or padj >= padj_threshold:
                continue
            significant_rows.append(
                build_output_row(
                    base_row=row,
                    contrast_id=contrast_id,
                    source_layout=source_layout,
                    source_fgsea_csv=source_fgsea_csv,
                    metadata=metadata
                )
            )

    ordered_rows = sorted(significant_rows, key=sort_key)
    fieldnames = list(EXTRA_FIELDS) + original_fieldnames
    write_rows(output_path, fieldnames, ordered_rows)

    print(
        f"{source_fgsea_csv}: total_rows={total_rows} "
        f"significant_rows={len(ordered_rows)} "
        f"output={display_path(output_path)}"
    )
    return ordered_rows, total_rows, len(ordered_rows)


def main():
    args = parse_args()
    input_root = Path(args.input_root).resolve()
    summary_output = Path(args.summary_output).resolve()

    if not input_root.exists():
        raise SystemExit(f"Input root does not exist: {input_root}")
    if not input_root.is_dir():
        raise SystemExit(f"Input root is not a directory: {input_root}")
    if args.padj_threshold <= 0:
        raise SystemExit("--padj-threshold must be greater than 0")

    manifest_path = input_root/"_summary"/"run_manifest.csv"
    manifest_rows = load_manifest(manifest_path)
    sources = discover_sources(input_root)

    all_rows = []
    combined_fieldnames = None

    for source in sources:
        rows, _, _ = process_source(
            source=source,
            manifest_rows=manifest_rows,
            padj_threshold=args.padj_threshold
        )
        if combined_fieldnames is None:
            input_path = Path(source["input_path"])
            with input_path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                combined_fieldnames = list(EXTRA_FIELDS) + validate_columns(reader.fieldnames, input_path)
        all_rows.extend(rows)

    if combined_fieldnames is None:
        combined_fieldnames = list(EXTRA_FIELDS)

    summary_rows = sorted(all_rows, key=sort_key)
    write_rows(summary_output, combined_fieldnames, summary_rows)
    print(f"{display_path(summary_output)}: combined_significant_rows={len(summary_rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
