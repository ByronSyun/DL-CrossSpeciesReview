import os
import gzip
import csv
import random
import pandas as pd
from pyfaidx import Fasta
from tqdm import tqdm

POSITIVE_FILE      = './processed_data/human_positive_samples_with_tissue_pipvcf.tsv'
VCF_FILE           = './human_raw_data/1kg_grch38_vcf/1kGP.chrALL.GRCh38.vcf.gz'  # ALL chromosomes
FASTA_FILE         = './reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
OUTPUT_FILE        = './processed_data/human_negative_pool_pipvcf.tsv'

SEQ_LENGTH         = 8192
CENTER_IDX         = SEQ_LENGTH // 2 - 1   # = 4095 (0-based) — ratGTEx standard
UNKNOWN_TISSUE_ID  = 15     # human, tissue unspecified
MAX_NEGATIVE_SAMPLES = 30_000   # pool size; balanced set takes 20k matched to positives

def get_ref_alt_sequences(chrom: str, pos: int, ref: str, alt: str, fasta_handle: Fasta):
    start = pos - (SEQ_LENGTH // 2) + 1
    end   = pos + (SEQ_LENGTH // 2)
    if start < 1:
        return None, None

    for key in [chrom,
                chrom.replace('chr', '') if chrom.startswith('chr') else f'chr{chrom}']:
        try:
            if key not in fasta_handle.keys():
                continue
            seq = fasta_handle[key][start - 1:end].seq.upper()
            if len(seq) != SEQ_LENGTH:
                continue
            if seq[CENTER_IDX] != ref.upper():
                return None, None   # ref mismatch
            alt_seq = seq[:CENTER_IDX] + alt.upper() + seq[CENTER_IDX + 1:]
            return seq, alt_seq
        except (KeyError, IndexError):
            continue
    return None, None


def create_negative_pool():
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    pos_df = pd.read_csv(POSITIVE_FILE, sep='\t')
    positive_locations = set()
    for _, row in pos_df.iterrows():
        try:
            positive_locations.add((str(row['CHR']), int(row['POS'])))
        except Exception:
            pass

    if not os.path.exists(FASTA_FILE):
        raise FileNotFoundError(f"FASTA not found: {FASTA_FILE}")
    fasta_handle = Fasta(
        FASTA_FILE,
        key_function=lambda k: k.split(' ')[0],
        sequence_always_upper=True
    )

    if not os.path.exists(VCF_FILE):
        raise FileNotFoundError(f"VCF not found: {VCF_FILE}")

    random.seed(42)

    negative_samples = []
    stats = dict(lines=0, snps=0, excluded_pos=0, boundary=0, ref_mismatch=0, n_content=0, valid_candidates=0)

    with gzip.open(VCF_FILE, 'rt') as fh:
        for line in tqdm(fh, desc="Scanning VCF"):
            if line.startswith('#'):
                continue
            stats['lines'] += 1
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue

            chrom, pos_str, _, ref, alt_field = parts[0], parts[1], parts[2], parts[3], parts[4]

            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom

            try:
                pos = int(pos_str)
            except ValueError:
                continue

            alt = alt_field.split(',')[0]   # take first alt allele if multi-allelic

            if not (len(ref) == 1 and len(alt) == 1):
                continue
            if ref not in 'ACGT' or alt not in 'ACGT':
                continue
            stats['snps'] += 1

            if (chrom, pos) in positive_locations:
                stats['excluded_pos'] += 1
                continue

            ref_seq, alt_seq = get_ref_alt_sequences(chrom, pos, ref, alt, fasta_handle)
            if ref_seq is None:
                stats['boundary'] += 1
                continue

            if ref_seq.count('N') > SEQ_LENGTH * 0.1:
                stats['n_content'] += 1
                continue

            variant_id = f"{chrom}:{pos}_{ref}/{alt}"
            candidate = {
                'variant_id':   variant_id,
                'ref_sequence': ref_seq,
                'alt_sequence': alt_seq,
                'label':        0,
                'tissue_id':    UNKNOWN_TISSUE_ID,
            }
            stats['valid_candidates'] += 1

            # Reservoir sampling to ensure UNIFORM distribution across all chromosomes
            if len(negative_samples) < MAX_NEGATIVE_SAMPLES:
                negative_samples.append(candidate)
            else:
                j = random.randint(0, stats['valid_candidates'] - 1)
                if j < MAX_NEGATIVE_SAMPLES:
                    negative_samples[j] = candidate

    if not negative_samples:
        raise RuntimeError("No negative samples generated.")

    neg_df = pd.DataFrame(negative_samples)
    cols   = ['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    neg_df[cols].to_csv(OUTPUT_FILE, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    create_negative_pool()
