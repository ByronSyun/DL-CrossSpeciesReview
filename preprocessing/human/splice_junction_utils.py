#!/usr/bin/env python3
"""
Utilities for building splice-junction-aware human sQTL benchmarks.
"""

from __future__ import annotations

import bisect
import gzip
from collections import defaultdict
from typing import Dict, List, Tuple


def open_text_auto(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def load_splice_sites_from_gtf(gtf_path: str) -> Dict[str, List[int]]:
    """
    Parse exon records from a GTF and collect exon boundary coordinates as
    splice-site proxies. Returned coordinates are 1-based genomic positions.
    """
    sites = defaultdict(set)

    with open_text_auto(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue

            chrom = parts[0]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue

            sites[chrom].add(start)
            sites[chrom].add(end)

    return {chrom: sorted(values) for chrom, values in sites.items()}


def nearest_site_distance(chrom: str, pos: int, splice_sites: Dict[str, List[int]]) -> int | None:
    """
    Return absolute distance from a genomic position to the nearest splice site.
    """
    keys_to_try = [chrom]
    if chrom.startswith("chr"):
        keys_to_try.append(chrom.replace("chr", ""))
    else:
        keys_to_try.append(f"chr{chrom}")

    sites = None
    for key in keys_to_try:
        if key in splice_sites:
            sites = splice_sites[key]
            break

    if not sites:
        return None

    idx = bisect.bisect_left(sites, pos)
    candidates = []
    if idx < len(sites):
        candidates.append(abs(sites[idx] - pos))
    if idx > 0:
        candidates.append(abs(sites[idx - 1] - pos))
    return min(candidates) if candidates else None


def distance_bin(distance: int | None) -> str:
    if distance is None:
        return "unannotated"
    if distance <= 2:
        return "bin_000_002"
    if distance <= 10:
        return "bin_003_010"
    if distance <= 50:
        return "bin_011_050"
    if distance <= 250:
        return "bin_051_250"
    if distance <= 1000:
        return "bin_251_1000"
    if distance <= 10000:
        return "bin_1001_10000"
    return "bin_gt_10000"
