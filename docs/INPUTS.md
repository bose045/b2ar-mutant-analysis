# Input formats

`--input-root` selects a local directory with the following layout. Only the files required
by your chosen mode need to exist. Paths inside this directory should remain relative.
The example generator provides executable examples of every statistical input schema.
`calculate.py` now creates the coordinate, regional, microswitch, ligand and contact inputs
from trajectories/events. See [CALCULATIONS.md](CALCULATIONS.md) for upstream configuration;
this page describes the plotting interchange format, not a requirement to precompute results elsewhere.

```text
inputs/
  EXAMPLE_DATA.json                  # include ONLY for synthetic data; enables figure watermark
  source_manifest.json              # verify mode: [{destination: relative path, sha256: ...}]
  regional_scores.json               # {A: {P60F: S, ...}, B: {...}}; all seven mutants per dataset
  assets/marker_structure.png        # optional for current mode
  data/
    current/
      statistics_and_provenance.json
      cache/<system>_<dataset>_<state>.npz
    jones_EC100.json
    historical/
      regions/dataset_A_framewise_interhelix_regions.npz
      regions/dataset_B_framewise_interhelix_regions.npz
      contacts/combined_AB_nonVDW_edges.csv
      contacts/functional_shift_scores_{strict,default,lenient}.csv
      contacts/display_edges_F0.3_{strict,default,lenient}.csv
      ligand/<system>_<Active|Inactive>_ligand_pose.csv
```

## Current marker/kink plots

Systems: `WT, P60F, P60M, P183A, P183F, P183H, P228A, P228M`.
Dataset labels: `A, B`. State labels: `Ag, IAg`.

Each NPZ contains equally sized `time_ns`, `distance_A`, and applicable `kink_60_deg`,
`kink_183_deg`, `kink_228_deg` arrays. Angles at all three sites are used for WT; mutant
panels use the corresponding mutated site. The fixed plot axes show 0–1000 ns and 8–25 Å.

Summary JSON has `results` and `missing` lists. Each result contains `label`, `system`,
`dataset`, `state`, `distance_mean_A`, and `angles`. Each applicable `angles` key is a
string residue number mapping to `mean` and `block_sd` (plus optional statistics).
`missing` lists unavailable system labels. The example includes 32 fictional systems.
Use the same frames to calculate JSON statistics and NPZ arrays; the plotter does not
automatically verify their agreement. S is not calculated for a dataset without both WT states.

## Regional plots

Regional NPZ keys are `<system>_<Active|Inactive>_<region>_TM<i>_TM<j>`, where i < j,
region is `extracellular`, `middle` or `intracellular`, and each value is a framewise
distance array in Å. All eight systems, two states and 21 helix pairs are required per dataset.
The arrays should describe distances between three-Cα regional centroids. See the residue
windows in `analysis/geometry.py`; `calculate.py` reads coordinates, while the plotter reads its outputs.

`regional_scores.json` supplies S values separately per dataset. WT is excluded from
correlations and functional classes. AL uses S >= 1 Å; IL uses S <= -1 Å. A resolved regional
shift has magnitude >= 0.25 Å. Consensus requires at least three observations, intended for
five observations in each class. Review this fixed rule before using a different study design.

## Contacts

Keys: `res1_num,res2_num,contact_type`. Network rows also require `res1_name,res2_name`,
`functional_shift_msgm` (F), `abs_functional_shift_msgm`, `msgm_act,msgm_inact`,
`maj_sign_act,maj_sign_inact`. Cross-criterion robust tables also require
`n_agree_act,n_agree_inact`. CSVs must have unique keys. Names are three-letter residue codes.

Supported contact codes: hbbb, hbsb, hbss, hp, wb, sb, pc, ps, ts, vdw. S6 excludes vdw.
Residue numbering and helix boundaries are study-specific. `analysis/contacts.py` computes
these scores and class membership; `calculate.py --stages score` writes the required tables.
S9 needs only the three keys and F. S10 uses prefiltered edges passing |F| >= 0.3 under each criterion.

## Ligand and assay examples

Ligand table columns: `ligand_rmsd_A`, `asp113_ionic_distance_A`, `salt_bridge_le_4A`.
New calculations store these under `ligand/dataset_A/` and `ligand/dataset_B/`; choose
`--ligand-dataset A` or `B` when plotting. Legacy direct-folder tables remain supported.
The legacy distance column is displayed as ligand N–simulation residue Asp85; verify this
mapping for your own input. Occupancy must be 0 or 1. Histidine mutant filename is `P183HSP`.
Bootstrap uses contiguous 500-frame means and 20,000 resamples (seed 20260820).
The current implementation includes a final partial block if present; supply equally sized
blocks for the study protocol and review uncertainty estimation before changing this convention.

Assay JSON: `{"rows": [{"AA": "F", "Pos": 88, "Norm": 1.2}, ...]}`.
Positions are 88,168,211,288,323. Supply only the intended 0.625 μM condition.
`Norm` is the already normalized activity; this code does not normalize experimental data.
Synthetic assay values are fictional and are not Jones experimental values.

## Structural/snapshot modes (no demo inputs supplied)

S12 requires `structural/publication_matrices.txt`, `publication_coordinates.tsv`,
`mapped_pairs.csv`, `publication_switch_map_highres.tga` and `publication_switch_map_highres.dat`.
The compositor expects the original VMD/Tachyon orthographic convention and table schema.
It requires tkinter/Tcl, generally provided by the Python distribution. It is a legacy study
compositor, not a generic molecular renderer.

Snapshot requires `figure_inventory.json`: a mapping of document identifiers to objects
with a `figures` list; each figure has `figure`, `caption` and `assets` (relative image paths).
Only supply material you have permission to distribute. Neither mode recreates molecular structures.

## Completed-package extraction

The standalone extractor requires `dataset_A/<label>` and `dataset_B/<label>`, each with
`manifest.json`, `COMPLETE.json`, and `upload/` files. Manifest fields: `label`, `frame_count`,
`source_files.trajectory`, `gpcrmd_files.topology`, `gpcrmd_files.trajectory`.
`COMPLETE.json` must contain `checksums` mapping filenames to SHA-256 hashes. Topology is
checked; the full trajectory hash is not recomputed. The 282 Cα atoms must be selected by
`segid PROA and name CA`, with consecutive simulation residue numbers 1–282.
Only validated complete packages should be used; this tool is not a replacement for upstream
atom-identity, frame-continuity and processing-provenance validation.
