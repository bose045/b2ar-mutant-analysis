# Figure-to-code mapping

| Figure | Launcher mode | Required input |
|---|---|---|
| 2c, 3 | `current` | Distance/kink NPZ arrays plus summary JSON; optional structural schematic |
| 4 | `historical --figures 4` | Regional framewise NPZ and dataset-specific S JSON |
| S1 | `historical --figures S1` | Assay matrix JSON; actual Jones data not included |
| S3–S5 | `historical --figures S3 S4 S5` | Regional framewise NPZ and S JSON |
| S6 | `historical --figures S6` | Precomputed non-VDW edge table with F and class values |
| S7 | `historical --figures S7` | Precomputed ligand metric tables |
| S9–S10 | `historical --figures S9 S10` | Precomputed strict/default/lenient score and edge tables |
| S12 | `structural` | VMD image, projection matrices, coordinates and residue-pair table |
| All supplied embedded figures | `snapshot` | Figure inventory and corresponding images; exports only |

## Upstream calculation coverage

| Result | Code that calculates it | Starting input |
|---|---|---|
| Figure 2c kink means/block SD | `analysis/geometry.py:kink,summarize` | Whole receptor coordinates |
| Figure 3 distances/means/SD/distributions; S | `calculate.py:coordinates,structural_scores` | Coordinates for matched Ag/IAg and WT |
| Reviewer full/half comparisons | `analysis/statistics.py:half_summaries` | Calculated marker traces |
| Figure 4 and S3–S5 | `analysis/geometry.py:regional_distances`; regional vendor helpers | Coordinates → 63 regional distances/frame → shifts, correlations, consensus |
| 21 established microswitches (tables/text) | `analysis/geometry.py:minimum_distance,pair_summary`; `analysis/statistics.py` | Heavy-atom coordinates → minima, peaks, medians, block uncertainty, R, correlations |
| Figure S6 | `calculate.py:contacts,contact_scores`; `analysis/contacts.py` | External GetContacts detector → events → CP → ΔCP → G → F → graph |
| Figure S7 | `analysis/geometry.py:ligand_metrics`; `reproduce.py:ligand` | Coordinates + explicit reference → aligned RMSD, ionic distance, occupancy and bootstrap |
| Figures S9–S10 | Same detector/scoring at three criteria; `reproduce.py:contacts` | Recalculated strict/default/lenient F and eligible edge sets |

The word `historical` is a retained plot-mode name. It can now read newly calculated inputs;
it does not mean the inputs are necessarily archived results. Archive reproduction and corrected
source-consistent recomputation must still be labeled distinctly.

Figures 1, 2a–b and S2 have no numerical or structural regeneration pipeline here.
Figure 2a specifically needs the original Bendix version/settings and selected structures.
Figure 2b needs the original MD/REMD/AlphaFold3 models and alignment/rendering specification.
S2's old caption refers to a Figure 5 absent from the current main inventory: its source
network/profile must be resolved before claiming exact reproduction. S1 consumes the published
normalized assay values; it does not reproduce the external laboratory/sequencing workflow.
`structural` does not launch VMD and is not exercised by the synthetic example.
The original 42-edge and 34-sign-flip counts are research results, not general input requirements.
Synthetic network examples intentionally contain a different number of edges.
