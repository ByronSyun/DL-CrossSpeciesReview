import os
import math
import csv
import pandas as pd
from pyfaidx import Fasta
from tqdm import tqdm

POSITIVE_FILE  = './processed_data/human_positive_samples_with_tissue_pipvcf.tsv'
FASTA_FILE     = './reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
OUTPUT_FILE    = './processed_data/human_positive_samples_sequences_pipvcf.tsv'

SEQ_LENGTH     = 8192
CENTER_IDX     = SEQ_LENGTH // 2 - 1   # = 4095 (0-based), same as ratGTEx standard
TARGET_N       = 20_000
RANDOM_SEED    = 42


def fetch_slice(fasta: Fasta, chrom: str, start_1based: int, end_1based: int):
    try:
        return str(fasta[chrom][start_1based - 1:end_1based]).upper()
    except (KeyError, IndexError):
        return None


def get_ref_alt_sequences(fasta: Fasta, chrom: str, pos: int, ref: str, alt: str):
    start = pos - (SEQ_LENGTH // 2) + 1
    end   = pos + (SEQ_LENGTH // 2)

    if start < 1:
        return None, None

    for key in [chrom,
                chrom.replace('chr', '') if chrom.startswith('chr') else f'chr{chrom}']:
        seq = fetch_slice(fasta, key, start, end)
        if seq is None:
            continue
        if len(seq) != SEQ_LENGTH:
            continue
        # Verify reference allele at center position
        if seq[CENTER_IDX] != ref.upper():
            return None, None
        alt_seq = seq[:CENTER_IDX] + alt.upper() + seq[CENTER_IDX + 1:]
        return seq, alt_seq

    return None, None


def stratified_sample(df: pd.DataFrame, target_n: int, seed: int) -> pd.DataFrame:
    tissues = df['Tissue'].unique()
    n_tissues = len(tissues)
    per_tissue = max(1, math.floor(target_n / n_tissues))

    sampled_parts = []
    used_indices = set()

    for tissue in tissues:
        tissue_df = df[df['Tissue'] == tissue]
        take_n = min(len(tissue_df), per_tissue)
        sampled = tissue_df.sample(n=take_n, random_state=seed)
        sampled_parts.append(sampled)
        used_indices.update(sampled.index)

    sampled_df = pd.concat(sampled_parts, ignore_index=False)
    n_so_far = len(sampled_df)

    remaining_need = target_n - n_so_far
    if remaining_need > 0:
        remaining_pool = df[~df.index.isin(used_indices)]
        if len(remaining_pool) > 0:
            top_up = remaining_pool.sample(
                n=min(remaining_need, len(remaining_pool)),
                random_state=seed
            )
            sampled_df = pd.concat([sampled_df, top_up], ignore_index=False)

    return sampled_df.sample(frac=1, random_state=seed).reset_index(drop=True)


def generate_sequences():
    os.makedirs(os.path.dirname(OUTPUT_FILE) if os.path.dirname(OUTPUT_FILE) else '.', exist_ok=True)

    pos_df = pd.read_csv(POSITIVE_FILE, sep='\t')
    sample_df = stratified_sample(pos_df, TARGET_N, RANDOM_SEED)

    if not os.path.exists(FASTA_FILE):
        raise FileNotFoundError(f"FASTA not found: {FASTA_FILE}")
    fasta = Fasta(
        FASTA_FILE,
        key_function=lambda k: k.split(' ')[0],
        sequence_always_upper=True
    )
    records = []
    stats = dict(total=0, boundary=0, ref_mismatch=0, success=0)

    for _, row in tqdm(sample_df.iterrows(), total=len(sample_df), desc="Fetching sequences"):
        stats['total'] += 1
        chrom     = str(row['CHR'])
        pos       = int(row['POS'])
        ref       = str(row['REF']).upper()
        alt       = str(row['ALT']).upper()
        tissue_id = int(row.get('tissue_id', 15))

        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom

        ref_seq, alt_seq = get_ref_alt_sequences(fasta, chrom, pos, ref, alt)

        if ref_seq is None:
            stats['boundary'] += 1
            continue

        if ref_seq[CENTER_IDX] != ref:
            stats['ref_mismatch'] += 1
            continue

        variant_id = f"{chrom}:{pos}_{ref}/{alt}"
        records.append([variant_id, ref_seq, alt_seq, 1, tissue_id])
        stats['success'] += 1

    if not records:
        raise RuntimeError("No records generated.")

    out_df = pd.DataFrame(
        records,
        columns=['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    )
    out_df.to_csv(OUTPUT_FILE, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    generate_sequences()
