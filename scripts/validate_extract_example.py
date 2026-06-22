#!/usr/bin/env python3
"""Regression check: our Tier-A candidate universe must recover SComatic's PASS calls.

Compares the candidates emitted by `extract_tier_a.py` on SComatic's bundled example
BAM against SComatic's own final PASS SNVs (Step4 calling.step2.tsv). The extractor is
*permissive* by design, so the requirement is: our candidates ⊇ SComatic PASS calls.

Run in the `clonemut` env (needs pyarrow), after extract_tier_a has been run:
    conda run -n clonemut python scripts/validate_extract_example.py
"""
import sys
import pyarrow.parquet as pq

CAND = "results/example_extract/Example.candidates.parquet"
SCOMATIC = ("external/SComatic/example_data/results/"
            "Step4_VariantCalling/Example.calling.step2.tsv")


def scomatic_pass():
    out = set()
    with open(SCOMATIC) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            chrom, start, _end, ref, alt, filt = c[0], c[1], c[2], c[3], c[4], c[5]
            if filt == "PASS":
                out.add((chrom, int(start), ref, alt))
    return out


def ours():
    t = pq.read_table(CAND).to_pydict()
    return {(t["chrom"][i], int(t["pos"][i]), t["ref"][i], t["alt"][i])
            for i in range(len(t["chrom"]))}


def main():
    truth, got = scomatic_pass(), ours()
    missing = truth - got
    extra = got - truth
    for s in sorted(truth):
        print(("  OK " if s in got else " MISS") + f" {s[0]}:{s[1]} {s[2]}>{s[3]}")
    print(f"\nSComatic PASS={len(truth)}  ours={len(got)}  "
          f"recovered={len(truth & got)}  missing={len(missing)}  extra={len(extra)}")
    if extra:
        print("extra (permissive candidates beyond SComatic PASS):",
              sorted(f"{c}:{p} {r}>{a}" for c, p, r, a in extra))
    if missing:
        sys.exit(f"FAIL: {len(missing)} SComatic PASS call(s) not in our candidates")
    print("PASS: candidate universe recovers all SComatic PASS calls")


if __name__ == "__main__":
    main()
