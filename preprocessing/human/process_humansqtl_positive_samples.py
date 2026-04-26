import glob
import os
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from splice_junction_utils import distance_bin, load_splice_sites_from_gtf, nearest_site_distance

SQTL_PIP_DIR = "./human_raw_data/sQTL_pip"
OUTPUT_DIR = "./processed_data"
OUTPUT_POSITIVE_FILE = os.path.join(OUTPUT_DIR, "human_positive_samples_with_tissue_pipvcf.tsv")
OUTPUT_SUMMARY_FILE = os.path.join(OUTPUT_DIR, "human_pip_variant_summary_pipvcf.tsv")
GTF_FILE = "./annotation/gencode.v38.annotation.gtf.gz"

POSITIVE_PIP_THRESHOLD = 0.9
SUPPORT_PIP_THRESHOLD = 0.1
MAX_SPLICE_DISTANCE_BP = 2000

GTEX_TISSUE_ID_MAP = {
    "Adipose_Subcutaneous": 0,
    "Adipose_Visceral_Omentum": 1,
    "Adrenal_Gland": 2,
    "Artery_Aorta": 3,
    "Artery_Coronary": 4,
    "Artery_Tibial": 5,
    "Brain_Cortex": 6,
    "Brain_Frontal_Cortex_BA9": 7,
    "Brain_Putamen_basal_ganglia": 8,
    "Colon_Sigmoid": 9,
    "Heart_Left_Ventricle": 10,
    "Liver": 11,
    "Lung": 12,
    "Muscle_Skeletal": 13,
    "Whole_Blood": 14,
}
DEFAULT_TISSUE_ID = 15


def parse_gtex_variant_id(variant_id: str) -> Optional[Tuple[str, int, str, str]]:
    parts = str(variant_id).split("_")
    if len(parts) < 4:
        return None
    chrom = parts[0]
    try:
        pos = int(parts[1])
    except ValueError:
        return None
    ref = parts[2].upper()
    alt = parts[3].upper()
    if len(ref) != 1 or len(alt) != 1:
        return None
    if ref not in "ACGT" or alt not in "ACGT":
        return None
    return chrom, pos, ref, alt


def get_tissue_name_from_pip_filename(path: str) -> str:
    base = os.path.basename(path)
    if ".v8." in base:
        return base.split(".v8.")[0]
    if ".variants_pip" in base:
        return base.split(".variants_pip")[0].replace(".v8.sqtl_signifpairs.txt.gz", "")
    return base.split(".")[0]


def get_tissue_id(tissue_name: str) -> int:
    return GTEX_TISSUE_ID_MAP.get(tissue_name, DEFAULT_TISSUE_ID)


def select_available_column(df: pd.DataFrame, candidates) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def aggregate_pip_variants(pip_files, splice_sites) -> pd.DataFrame:
    agg: Dict[Tuple[str, int, str, str], dict] = {}

    for path in pip_files:
        tissue = get_tissue_name_from_pip_filename(path)
        tissue_id = get_tissue_id(tissue)

        try:
            df = pd.read_csv(path, sep="\t", compression="infer")
        except Exception:
            continue

        variant_col = select_available_column(df, ["variant_id", "variant", "snp", "vid"])
        pip_col = select_available_column(df, ["pip", "PIP", "posterior_prob", "posterior_inclusion_prob"])
        rank_col = select_available_column(df, ["rank", "Rank", "posterior_rank"])
        if variant_col is None or pip_col is None:
            logger.warning(
                "  Skipping %s (required columns missing). Found columns: %s",
                os.path.basename(path), ",".join(df.columns[:20].tolist()),
            )
            continue

        keep_cols = [variant_col, pip_col]
        if rank_col:
            keep_cols.append(rank_col)
        df = df[keep_cols].copy()
        df = df.rename(columns={variant_col: "variant_id", pip_col: "pip"})
        if rank_col:
            df = df.rename(columns={rank_col: "rank"})
        else:
            df["rank"] = np.nan

        df["pip"] = pd.to_numeric(df["pip"], errors="coerce")
        df["rank"] = pd.to_numeric(df["rank"], errors="coerce")
        df = df.dropna(subset=["variant_id", "pip"])

        df = df.groupby("variant_id", as_index=False).agg({"pip": "max", "rank": "min"})

        for row in df.itertuples(index=False):
            parsed = parse_gtex_variant_id(row.variant_id)
            if parsed is None:
                continue
            chrom, pos, ref, alt = parsed
            key = (chrom, pos, ref, alt)
            pip_val = float(row.pip)
            rank_val = float(row.rank) if pd.notna(row.rank) else np.nan

            if key not in agg:
                agg[key] = {
                    "CHR": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "best_tissue": tissue,
                    "best_tissue_id": tissue_id,
                    "max_pip": pip_val,
                    "min_rank": rank_val,
                    "tissue_support_any": 1 if pip_val > 0 else 0,
                    "tissue_support_high": 1 if pip_val >= SUPPORT_PIP_THRESHOLD else 0,
                }
            else:
                rec = agg[key]
                if pip_val > rec["max_pip"]:
                    rec["max_pip"] = pip_val
                    rec["best_tissue"] = tissue
                    rec["best_tissue_id"] = tissue_id
                if pd.notna(rank_val):
                    if pd.isna(rec["min_rank"]) or rank_val < rec["min_rank"]:
                        rec["min_rank"] = rank_val
                if pip_val > 0:
                    rec["tissue_support_any"] += 1
                if pip_val >= SUPPORT_PIP_THRESHOLD:
                    rec["tissue_support_high"] += 1

    if not agg:
        raise RuntimeError("No valid variants aggregated from PIP files.")

    summary_df = pd.DataFrame(list(agg.values()))
    summary_df["nearest_splice_dist"] = summary_df.apply(
        lambda r: nearest_site_distance(str(r["CHR"]), int(r["POS"]), splice_sites),
        axis=1,
    )
    summary_df = summary_df[summary_df["nearest_splice_dist"].notna()].copy()
    summary_df["nearest_splice_dist"] = summary_df["nearest_splice_dist"].astype(int)
    summary_df["distance_bin"] = summary_df["nearest_splice_dist"].apply(distance_bin)
    summary_df["variant_id"] = summary_df.apply(
        lambda r: f"{r['CHR']}:{int(r['POS'])}_{r['REF']}/{r['ALT']}",
        axis=1,
    )
    return summary_df


def process_positive_samples():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    pip_files = sorted(glob.glob(os.path.join(SQTL_PIP_DIR, "*.variants_pip.txt*")))
    if not pip_files:
        raise FileNotFoundError(f"No PIP files found in {SQTL_PIP_DIR}")

    if not os.path.exists(GTF_FILE):
        raise FileNotFoundError(f"GTF not found: {GTF_FILE}")
    splice_sites = load_splice_sites_from_gtf(GTF_FILE)

    summary_df = aggregate_pip_variants(pip_files, splice_sites)
    summary_df = summary_df.sort_values(["max_pip", "tissue_support_high"], ascending=[False, False]).reset_index(drop=True)
    summary_df.to_csv(OUTPUT_SUMMARY_FILE, sep="\t", index=False)

    positive_df = summary_df[
        (summary_df["max_pip"] >= POSITIVE_PIP_THRESHOLD)
        & (summary_df["nearest_splice_dist"] <= MAX_SPLICE_DISTANCE_BP)
    ].copy()

    positive_df = positive_df.rename(
        columns={
            "best_tissue": "Tissue",
            "best_tissue_id": "tissue_id",
        }
    )

    positive_out = positive_df[
        [
            "CHR",
            "POS",
            "REF",
            "ALT",
            "Tissue",
            "tissue_id",
            "max_pip",
            "tissue_support_any",
            "tissue_support_high",
            "min_rank",
            "nearest_splice_dist",
            "distance_bin",
        ]
    ].copy()

    positive_out.to_csv(OUTPUT_POSITIVE_FILE, sep="\t", index=False)


if __name__ == "__main__":
    process_positive_samples()
