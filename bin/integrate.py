#!/usr/bin/env python3
"""
Integrate RNA-seq expression/count information with eggNOG annotation and optional
metabolic model information.

Expected main output:
    final_complete_multiomics.csv
"""

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd


def read_table_auto(path: str) -> pd.DataFrame:
    """Read CSV/TSV-like tables with simple delimiter inference."""
    path_obj = Path(path)
    if not path_obj.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    if path_obj.suffix.lower() in {".tsv", ".txt", ".annotations"}:
        return pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    return pd.read_csv(path, low_memory=False)


def load_counts(counts_file: str) -> pd.DataFrame:
    print(f"Loading count/expression file: {counts_file}")
    df = read_table_auto(counts_file)

    if df.empty:
        raise ValueError("Counts/expression file is empty.")

    df = df.dropna(subset=[df.columns[0]]).copy()

    rename_map = {
        "Geneid": "RNA_locus_tag",
        "geneid": "RNA_locus_tag",
        "locus_tag": "RNA_locus_tag",
        "Count": "expression",
        "count": "expression",
        "counts": "expression",
        "name": "gene_symbol",
        "desc": "description",
        "KEGG": "KEGG_ID",
        "pathway": "KEGG_Pathway_id",
    }
    df = df.rename(columns={c: rename_map.get(c, c) for c in df.columns})

    if "RNA_locus_tag" not in df.columns:
        df = df.rename(columns={df.columns[0]: "RNA_locus_tag"})

    if "Met_locus_tag" not in df.columns:
        df["Met_locus_tag"] = df["RNA_locus_tag"]

    return df


def load_eggnog(eggnog_file: str) -> pd.DataFrame:
    """Load eggNOG-mapper annotation table and normalise key column names."""
    print(f"Loading eggNOG annotation file: {eggnog_file}")

    path_obj = Path(eggnog_file)
    if not path_obj.exists():
        raise FileNotFoundError(f"eggNOG annotation file not found: {eggnog_file}")

    # eggNOG files often have comment lines before the real header.
    header_line = None
    with open(eggnog_file, "r", encoding="utf-8", errors="replace") as handle:
        for i, line in enumerate(handle):
            if line.startswith("#query") or line.startswith("query"):
                header_line = i
                break

    if header_line is not None:
        egg = pd.read_csv(eggnog_file, sep="\t", skiprows=header_line, low_memory=False)
        egg.columns = [str(c).lstrip("#") for c in egg.columns]
    else:
        egg = pd.read_csv(eggnog_file, sep="\t", comment="#", low_memory=False)

    if egg.empty:
        print("WARNING: eggNOG annotation file is empty. Continuing without eggNOG merge.")
        return pd.DataFrame()

    egg = egg.rename(
        columns={
            "query": "RNA_locus_tag",
            "seed_ortholog": "seed_ortholog",
            "evalue": "eggnog_evalue",
            "score": "eggnog_score",
            "eggNOG_OGs": "eggNOG_OGs",
            "max_annot_lvl": "max_annot_lvl",
            "COG_category": "COG",
            "Description": "description",
            "Preferred_name": "gene_symbol",
            "GOs": "GO",
            "KEGG_ko": "KEGG_ID",
            "KEGG_Pathway": "KEGG_Pathway_id",
        }
    )

    keep_cols = [
        "RNA_locus_tag",
        "gene_symbol",
        "description",
        "COG",
        "GO",
        "KEGG_ID",
        "KEGG_Pathway_id",
    ]
    keep_cols = [c for c in keep_cols if c in egg.columns]

    if "RNA_locus_tag" not in keep_cols:
        print("WARNING: no query/RNA_locus_tag column found in eggNOG file. Skipping eggNOG merge.")
        return pd.DataFrame()

    return egg[keep_cols].drop_duplicates(subset=["RNA_locus_tag"])


def merge_eggnog(counts_df: pd.DataFrame, eggnog_df: pd.DataFrame) -> pd.DataFrame:
    if eggnog_df.empty:
        return counts_df

    print("Merging expression table with eggNOG annotations...")
    merged = counts_df.merge(eggnog_df, on="RNA_locus_tag", how="left", suffixes=("", "_eggnog"))

    for col in ["gene_symbol", "description", "COG", "GO", "KEGG_ID", "KEGG_Pathway_id"]:
        egg_col = f"{col}_eggnog"
        if egg_col in merged.columns:
            if col in merged.columns:
                merged[col] = merged[col].fillna(merged[egg_col])
                merged = merged.drop(columns=[egg_col])
            else:
                merged = merged.rename(columns={egg_col: col})

    return merged


def load_metabolic_model(metabolic_model: str, locus_prefix: str) -> pd.DataFrame:
    if not metabolic_model:
        print("No metabolic model supplied. Skipping metabolic integration.")
        return pd.DataFrame()

    metabolic_path = Path(os.path.expanduser(metabolic_model))
    if not metabolic_path.exists():
        print(f"WARNING: metabolic model not found: {metabolic_model}. Skipping metabolic integration.")
        return pd.DataFrame()

    print(f"Loading metabolic model: {metabolic_path}")
    model = pd.read_csv(metabolic_path, header=2, low_memory=False)
    model.columns = [str(c).strip() for c in model.columns]

    # Fallback positional names based on the source metabolic model layout.
    positional_names = {
        0: "Metabolite_ID",
        1: "Metabolite_name",
        7: "Genes",
        8: "Subsystem",
    }
    for pos, name in positional_names.items():
        if pos < len(model.columns) and name not in model.columns:
            model = model.rename(columns={model.columns[pos]: name})

    required = {"Metabolite_name", "Genes", "Subsystem"}
    missing = required.difference(model.columns)
    if missing:
        print(f"WARNING: metabolic model missing columns {sorted(missing)}. Skipping metabolic integration.")
        return pd.DataFrame()

    model = model[["Metabolite_name", "Genes", "Subsystem"]].dropna(subset=["Genes"])
    model = model[model["Genes"].astype(str).str.strip() != ""].copy()

    model = model.assign(locus_tag=model["Genes"].astype(str).str.split()).explode("locus_tag")
    model["locus_tag"] = model["locus_tag"].astype(str).str.strip()

    if locus_prefix:
        model = model[model["locus_tag"].str.startswith(locus_prefix)]

    def extract_pathway(subsystem):
        if pd.isna(subsystem):
            return pd.NA, pd.NA
        first = str(subsystem).split(";")[0].strip()
        match = re.match(r"(rn\d+|map\d+)\s+(.*)", first)
        if match:
            return match.group(1).replace("rn", "map"), match.group(2).strip()
        return pd.NA, first

    extracted = model["Subsystem"].apply(extract_pathway)
    model["metabolite_pathway_id"] = [x[0] for x in extracted]
    model["metabolite_pathway_name"] = [x[1] for x in extracted]

    model = model.rename(columns={"Metabolite_name": "metabolite"})
    model = model[["locus_tag", "metabolite", "metabolite_pathway_id", "metabolite_pathway_name"]]

    print(f"Metabolic pairs: {len(model)}")
    print(f"Unique metabolic genes: {model['locus_tag'].nunique()}")
    return model


def merge_metabolic_data(df: pd.DataFrame, metabolic_df: pd.DataFrame) -> pd.DataFrame:
    if metabolic_df.empty:
        return df

    print("Merging expression table with metabolic model information...")
    final = df.merge(metabolic_df, left_on="RNA_locus_tag", right_on="locus_tag", how="left")
    return final.drop(columns=["locus_tag"], errors="ignore")


def main() -> int:
    parser = argparse.ArgumentParser(description="Integrate RNA expression, eggNOG annotation, and metabolic model data.")
    parser.add_argument("--counts", required=True, help="Input count/expression table, for example important_genes.csv")
    parser.add_argument("--eggnog", required=True, help="eggNOG-mapper .emapper.annotations file")
    parser.add_argument("--metabolic-model", default="", help="Optional metabolic model CSV file")
    parser.add_argument("--output", default="final_complete_multiomics.csv", help="Output CSV filename")
    parser.add_argument("--locus-prefix", default="DJ41_", help="Optional locus-tag prefix for metabolic gene filtering")
    args = parser.parse_args()

    counts_df = load_counts(args.counts)
    eggnog_df = load_eggnog(args.eggnog)
    merged_df = merge_eggnog(counts_df, eggnog_df)

    metabolic_df = load_metabolic_model(args.metabolic_model, args.locus_prefix)
    final = merge_metabolic_data(merged_df, metabolic_df)

    preferred_order = [
        "RNA_locus_tag",
        "Met_locus_tag",
        "expression",
        "gene_symbol",
        "description",
        "COG",
        "GO",
        "KEGG_ID",
        "KEGG_Pathway_id",
        "metabolite",
        "metabolite_pathway_id",
        "metabolite_pathway_name",
    ]

    ordered_cols = [c for c in preferred_order if c in final.columns]
    remaining_cols = [c for c in final.columns if c not in ordered_cols]
    final = final[ordered_cols + remaining_cols]

    final.to_csv(args.output, index=False)

    print(f"DONE: {args.output}")
    print(f"Rows: {len(final)}")
    print(f"Columns: {len(final.columns)}")
    if "gene_symbol" in final.columns:
        print(f"Annotated gene symbols: {final['gene_symbol'].notna().sum()}")
    if "metabolite" in final.columns:
        print(f"Rows with metabolite information: {final['metabolite'].notna().sum()}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
