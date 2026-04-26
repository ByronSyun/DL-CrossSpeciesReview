import os
import pandas as pd

POSITIVE_FILE = './processed_data/human_positive_samples_sequences_pipvcf.tsv'
NEGATIVE_FILE = './processed_data/human_negative_pool_pipvcf.tsv'
OUTPUT_FILE   = './processed_data/humansqtl_silver_benchmark_balanced_pipvcf.tsv'
RANDOM_SEED   = 42


def create_balanced():
    pos_df = pd.read_csv(POSITIVE_FILE, sep='\t', header=None)
    n_pos  = len(pos_df)

    neg_df = pd.read_csv(NEGATIVE_FILE, sep='\t', header=None)

    if len(neg_df) < n_pos:
        raise RuntimeError(
            f"Not enough negatives ({len(neg_df):,}) to match positives ({n_pos:,})."
        )

    neg_sampled = neg_df.sample(n=n_pos, random_state=RANDOM_SEED)

    combined = pd.concat([pos_df, neg_sampled], ignore_index=True)
    combined = combined.sample(frac=1, random_state=RANDOM_SEED).reset_index(drop=True)

    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    combined.to_csv(OUTPUT_FILE, sep='\t', index=False, header=False)


if __name__ == '__main__':
    create_balanced()
