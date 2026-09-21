# Calculation pipeline

The intended contract is coordinates/topology + explicit analysis settings → intermediate
observables → statistical tables → figures. Cached research results are not required for
the implemented numerical analyses. No script silently fetches or publishes research data.

## Configuration

Use the generated synthetic `study.json` as an executable full-schema example. For private
inputs create `inputs/study.json`; paths resolve relative to that JSON. Each system requires:

```json
{
  "dataset": "A",
  "system": "WT",
  "state": "Ag",
  "topology": "WT_A_Ag.psf",
  "trajectories": ["WT_A_Ag_TM_aligned.xtc"],
  "protein_selection": "segid PROA",
  "first": 0,
  "last": 9999,
  "stride": 1,
  "coordinates_rotated": true,
  "contact_chain": "A",
  "contacts": {
    "selection": "protein",
    "solvent": "resname WAT",
    "first": 0,
    "last": 8999,
    "stride": 1
  },
  "ligand": {
    "reference": "WT_A_Ag_start.pdb",
    "selection": "resname P0G and not name H*",
    "cation_selection": "resname P0G and name N1",
    "periodic_mode": "whole"
  }
}
```

These names are examples, **not verified settings for every package**. Supply the matching
IAg entry and WT entries for both states in every dataset. Complete regional plots need all
eight systems in each dataset. Ligand atom names, reference topology, protonation, solvent
name, chain, numbering and selections must be checked against each topology. `protein_selection`
and ligand selections use MDAnalysis syntax; `contacts.selection` uses VMD syntax.

At the top level, `systems` contains these entries; defaults are:

- `block_ns`: 100; `microswitch_block_ns`: same; `pair_threshold_A`: 0.25.
- `response_threshold_A`: 0.1 for the separately defined WT-referenced microswitch response.
- `bootstrap_seed`: 20260903; `bootstrap_resamples`: 10000 for microswitch block statistics.
- `contact_profile`: `pooled_S`. AL: S≥1 Å; IL: S≤−1 Å; WT/WT-like excluded.
- `require_ligand`: false; set true for a complete ligand analysis.
- `synthetic`: false; true only for artificial test inputs.

Dataset A P183H requires `source_revision: "replacement_2"`, an explicit source-selection
attestation. The general JSON reader cannot verify that claim from an arbitrary filename;
inspect the preparation manifest. The separate package extractor checks replacement paths.

## Trajectory requirements

All receptor CA residues must be ordered 1–282, with Leu47 and Leu215. The protein must be
whole across periodic images. The calculation script does **not** autoimage, fit or rewrite
trajectories. Atom topology must match the trajectory exactly. Pass numbered files in
chronological order; no files are skipped. Times must be uniformly spaced and increasing.
Reset segment timestamps, discontinuities, incomplete frames and wrong atom mappings must
be resolved in the preparation stage rather than concealed by this analysis.

Coordinate windows are zero-based and inclusive. Use an even number of frames for the
half-trajectory comparison and exact multiples of the block duration for kink statistics.
At least two complete blocks are required. Each half uses its **corresponding WT half** when
calculating S. Framewise SD describes fluctuations, not uncertainty in an independent mean.

Kinks use the acute angle between centered SVD axes of six CA atoms before and six after
the site, excluding the central residue. This is the quantitative Figure 2c definition;
it is **not a replacement for the Bendix rendering algorithm** used in Figure 2a.

Regional distances use arithmetic three-CA centroids at the documented helix ends/middle.
Microswitch distances use the minimum over all atom pairs selected with `not name H*`,
matching the archived convention. Verify that the topology names all hydrogens accordingly.

Ligand RMSD uses equal-weight TM CA Kabsch fitting, then heavy-atom ligand RMSD without
symmetry reassignment. `whole` requires an already whole protein–ligand complex. For an
unrotated periodic trajectory, `minimum_image` translates the entire ligand to the nearest
image around Asp85 before measurement. It does not reconstruct a split ligand internally.
Never apply original unrotated box vectors to fitted/rotated coordinates.

## GetContacts stage

Provide top-level `getcontacts` settings: `directory` (external checkout), `python` (an
interpreter with vmd-python and engine dependencies), `cores`, and preferably `expected_commit`.
The reviewed engine must support `hp` as well as the other requested types. The installed
engine's own README is authoritative for its environment. This repository does not vendor it.

```sh
python calculate.py inputs/study.json --out generated/study --stages coordinates
python calculate.py inputs/study.json --out generated/study --stages contacts
python calculate.py inputs/study.json --out generated/study --stages contacts --execute-contacts
python calculate.py inputs/study.json --out generated/study --stages score
```

Without `--execute-contacts`, only commands/provenance are prepared. All Python engine files
are hashed; a supplied commit must match and tracked modifications cause rejection. Use one
explicit detector trajectory per job. If coordinate analysis uses numbered segments, supply
an independently validated `contacts.trajectory` with the same atom order and mapped window.
The wrapper deliberately does not silently concatenate or transform coordinates.

The `contacts` window can differ from the coordinate window, e.g. archived 9000 contact
frames versus 10000 marker frames. The difference remains explicit in the configuration.
The event header count and every event frame are checked. The archived detector's `--end 0`
behavior is unsafe and rejected. Engine frame-range behavior must be tested before production.

For existing event files, supply `contact_events: {"strict": "...", "default": "...",
"lenient": "..."}` in each job. Default frame labels are absolute source indices; specify
`contact_frame_labels: "ordinal"` only if events were explicitly renumbered.

Each pair/type is counted once per frame, including first and last events. Empty frames
remain in the denominator. Water-bridge endpoint atoms are the first two fields; mediator
water atoms are not interpreted as extra receptor contacts. Extended/two-water bridges are
not among the ten archived types. Missing contacts mean zero only in successfully parsed
complete files; missing states/files cause errors. Numeric residue ordering avoids swapped
keys across mutations. Edge names currently follow the first observed naming; review any
edge involving the substituted residue before assigning a WT label.

## Scores and statistical conventions

CP = contact-containing frames / selected frames. ΔCP = CP_Ag − CP_IAg.
For each class, discard |ΔCP|≤10⁻⁶, choose the majority sign, and take the signed geometric
mean of |ΔCP|+10⁻⁶ among observations supporting that sign. Ties give G=0. F=G_AL−G_IL.
The small epsilon is preserved from the archived calculation, not silently removed.
Edges require ≥3 supporting observations in each class, nonadjacent residues, both endpoints
in the defined TM/H8 regions, and |F|≥0.3. S6 further excludes VDW. Class sizes are data-dependent;
the historical five/five split and 42/34 counts are **not enforced as expected answers**.

`legacy_A` is an explicitly different historical profile: six AL observations (including WT)
versus P183A/P183F, with 4/2 agreement. Do not mix it with pooled S-classified outputs.
Sensitivity changes twelve geometric floating-point criteria by ±5%; sequence-separation
exclusions remain fixed. Calling this a change to *every numerical parameter* is inaccurate.

For each microswitch, shared 60-bin histograms yield peak shifts; medians are also reported.
Reference orientation signs are explicit in `analysis/statistics.py`. R subtracts the
same-dataset WT oriented shift. Correlations exclude WT; BH adjustments are within each
dataset and correlation method among finite tests. Constant inputs produce undefined values.
Bootstrap CIs refer to the **difference of mean block medians**, not the pooled median
difference: states are independently resampled and a final ≥half-sized block is retained.
The seed/resample count are configurable; historical row ordering may change exact bootstrap
Monte Carlo endpoints. Ligand plotting separately uses 500-frame means and 20000 resamples.
Use the original 0.1-ns spacing for 50-ns ligand blocks; short synthetic examples are tests,
not useful uncertainty estimates. None of these comparisons establishes convergence.

## Outputs and evidence

`calculation_config.json`, `coordinate_provenance.json`, and engine command/hash records
identify inputs and methods. Input topology is hashed; large trajectory paths/sizes/mtime
are recorded, not rehashed. Keep preparation checksums for stronger identity verification.
These metadata may contain private paths; all generated outputs are ignored by Git.

`data/current/cache` stores marker/kink arrays; `data/historical/regions` stores regional
arrays; `data/microswitches` stores all 21 pair traces and statistical summaries;
`full_and_half_trajectory_scores.csv` and `structural_shift_scores.csv` report S;
`contact_probabilities` stores per-type CP; contact score tables include per-observation ΔCP;
ligand arrays are separated by dataset. The plotter consumes these calculated outputs.

Tests cover synthetic coordinates and artificial contact events. They do not certify the
manuscript numbers. Full source-consistent scientific validation remains required.
