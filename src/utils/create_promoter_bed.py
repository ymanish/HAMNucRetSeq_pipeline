#/src/utils/create_promoter_bed.py
# Created on 2025-04-06

import pandas as pd
from src.config import path
from itertools import combinations


gtf_path = path.DATA_DIR / 'processed/protein_coding_transcripts.gtf'
gtf = pd.read_csv(
    gtf_path, sep='\t', comment='#', header=None,
    names=['seqname','source','feature','start','end','score','strand','frame','attributes']
)

# Keep only transcript lines
tx = gtf[gtf['feature']=='transcript'].copy()

tx['transcript_id'] = tx['attributes'].str.extract(r'transcript_id "([^\"]+)"')
tx['gene_id']       = tx['attributes'].str.extract(r'gene_id "([^\"]+)"')


up, down = 1000, 1000
# Convert to 0-based BED coordinates
# BED: [start, end) where start is 0-based, end is 1-based
# So [TSS - 1000, TSS + 1000] becomes [TSS - 1000 - 1, TSS + 1000)

tx['TSS'] = tx.apply(lambda r: r.start if r.strand=='+' else r.end, axis=1)
# Define promoter bounds (0-based half-open)
tx['prom_start'] = (tx['TSS'] - up - 1).clip(lower=0)
tx['prom_end']   = tx['TSS'] + down

# Keep only needed columns
prom = tx[['seqname','prom_start','prom_end','gene_id','transcript_id']].copy()



def find_redundant_pairs(intervals):
    """
    intervals: list of (tx_id, start, end)
    returns: list of (txA, txB, overlap_len, small_len)
    """
    pairs = []
    for (txA, sA, eA), (txB, sB, eB) in combinations(intervals, 2):
        # compute overlap
        common = max(0, min(eA, eB) - max(sA, sB))
        small = min(eA-sA, eB-sB)
        if common > 0.5 * small:
            pairs.append((txA, txB, common, small))
    return pairs



def compute_overlap_coverage(intervals):
    """
    Given a list of (start, end) intervals, compute the percent of the union
    covered by overlapping (depth >=2) regions.
    """
    pts = sorted({p for iv in intervals for p in iv})
    total_union = 0
    overlap_union = 0
    for x, y in zip(pts, pts[1:]):
        seg_len = y - x
        depth = sum(1 for s, e in intervals if s <= x and e >= y)
        if depth >= 1:
            total_union += seg_len
        if depth >= 2:
            overlap_union += seg_len
    return (overlap_union / total_union * 100) if total_union else 0

def compute_gene_metrics(intervals, tx_pairs):
    """
    Compute metrics for a gene given its promoter intervals and redundant pairs.
    intervals: list of (start, end) for all transcripts
    tx_pairs: list of redundant (txA, txB, overlap, small)
    Returns dict with num_pairs, num_txs_involved, overlap_coverage_pct
    """
    num_pairs = len(tx_pairs)
    txs_involved = {a for a, *_ in tx_pairs} | {b for _, b, *_ in tx_pairs}
    coverage_pct = compute_overlap_coverage(intervals)
    return {
        'num_redundant_pairs': num_pairs,
        'num_txs_involved': len(txs_involved),
        'overlap_coverage_pct': coverage_pct
    }


# Collect metrics and redundant transcripts per gene
metrics = []
cnt = 0
keep_tx = set(prom.transcript_id)
for gene, grp in prom.groupby('gene_id'):
    # Build list of (tx_id, start, end)
    ints = list(grp[['transcript_id','prom_start','prom_end']].itertuples(index=False, name=None))
    # Find redundant transcript pairs
    redundant = find_redundant_pairs(ints)
    # Compute gene metrics
    intervals = [(s, e) for _, s, e in ints]
    gene_metrics = compute_gene_metrics(intervals, redundant)
    gene_metrics['gene_id'] = gene
    metrics.append(gene_metrics)
    # Discard redundant transcripts (keep first of each pair)
    for a, b, *_ in redundant:
        keep_tx.discard(b)
    
  
# # Save metrics to TSV
# metrics_df = pd.DataFrame(metrics)
# metrics_df.to_csv(path.DATA_DIR / 'processed/promoter_overlap_metrics.tsv', sep='\t', index=False)

# --- Step 4: Write final promoter BED ---
final = prom[prom.transcript_id.isin(keep_tx)][['seqname','prom_start','prom_end', 'transcript_id']]
# Remove mitochondrial chromosome (chrM)
final = final[final['seqname'] != 'chrM']
final = final.drop_duplicates()
final.to_csv(path.DATA_DIR / 'processed/promoter_regions.bed', sep='\t', header=False, index=False)

print(f"Wrote {len(final)} promoter regions (post-filter).")
