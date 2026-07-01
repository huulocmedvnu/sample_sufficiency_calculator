#!/usr/bin/env python
"""Pass 0 (study-design layout) — derive the full experimental layout of Tahoe-100M from
metadata/obs_metadata.parquet (one row per cell, 100.6M rows, INCLUDING pre-filter cells).

This is independent of the expression stream (Pass 2). It emits the physical-well / plate layout,
the raw (pre-filter) vs retained (post-filter) cell-count distributions, per-plate variation, and the
per-(drug x cell_line) cell counts used as the quota baseline N0.

A Tahoe "well" (`sample`) is one drug on one plate, with all cell lines pooled and computationally
demultiplexed; `pass_filter=="full"` are the cells present in the expression release.

Outputs -> $OUT/layout_summary.json, $OUT/per_plate.csv, $OUT/condition_counts.csv
           (drug, cell_line, cells_pre, cells_post).
"""
import os, json, numpy as np, pandas as pd, pyarrow.parquet as pq

OUT = os.environ.get("OUT", "/mnt/hdd2/loc-tran/tahoe_work/out")
OBS = "/mnt/hdd2/loc-tran/tahoe_work/meta/metadata/obs_metadata.parquet"
FULL = "full"


def dist(s):
    return dict(min=int(s.min()), median=int(s.median()), mean=round(float(s.mean()), 1),
                max=int(s.max()), total=int(s.sum()), n=int(s.size))


def main():
    df = pq.read_table(OBS, columns=["sample", "plate", "drug", "cell_line", "cell_name",
                                     "pass_filter"]).to_pandas()
    full = df["pass_filter"] == FULL
    wells = df.groupby("sample", observed=True)
    pre = wells.size()
    post = df[full].groupby("sample", observed=True).size()

    per_plate = df.groupby("plate", observed=True).agg(
        wells=("sample", "nunique"), cells_pre=("sample", "size"))
    per_plate["cells_post"] = df[full].groupby("plate", observed=True).size()
    per_plate["pct_filtered"] = (100 * (1 - per_plate.cells_post / per_plate.cells_pre)).round(2)
    per_plate.to_csv(os.path.join(OUT, "per_plate.csv"))

    cc = df[full].groupby(["drug", "cell_line"], observed=True).size().rename("cells_post").reset_index()
    ccpre = df.groupby(["drug", "cell_line"], observed=True).size().rename("cells_pre").reset_index()
    cc = cc.merge(ccpre, on=["drug", "cell_line"])
    cc.to_csv(os.path.join(OUT, "condition_counts.csv"), index=False)

    summary = dict(
        total_cells_prefilter=int(len(df)),
        pass_filter_counts={str(k): int(v) for k, v in df.pass_filter.value_counts().items()},
        pct_prefilter_excess=round(100 * (1 - int(full.sum()) / len(df)), 2),
        n_wells=int(df["sample"].nunique()),
        n_plates=int(df.plate.nunique()),
        wells_per_plate=int(per_plate.wells.iloc[0]) if per_plate.wells.nunique() == 1 else None,
        n_drugs=int(df.drug.nunique()), n_cell_lines=int(df.cell_line.nunique()),
        n_conditions=int(len(cc)),
        cells_per_well_prefilter=dist(pre),
        cells_per_well_postfilter=dist(post),
        cells_per_condition_postfilter=dist(cc["cells_post"]),
        per_plate={str(pl): dict(wells=int(r.wells), cells_pre=int(r.cells_pre),
                                 cells_post=int(r.cells_post), pct_filtered=float(r.pct_filtered))
                   for pl, r in per_plate.iterrows()},
    )
    json.dump(summary, open(os.path.join(OUT, "layout_summary.json"), "w"), indent=2)
    print(json.dumps(summary, indent=2))
    print("\n=== per plate ===\n", per_plate.to_string())
    print(f"\n[wrote] layout_summary.json, per_plate.csv, condition_counts.csv in {OUT}")


if __name__ == "__main__":
    main()
