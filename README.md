# β₂AR mutant analysis

Python postprocessing and plotting code for a β₂-adrenergic receptor mutation study.
This repository contains analysis code and synthetic examples.

## Quick start

Python 3.11 or newer. Commands below work from the repository root on Windows, macOS or Linux.

### Calculate from coordinates

```sh
python -m pip install -r requirements-trajectories.txt
python examples/make_trajectory_demo.py
python run_study.py generated/trajectory_demo/study.json --out generated/end_to_end --contact-mode existing-events
```

This reads **32 synthetic coordinate trajectories**, calculates the observables, counts
contact-events, derives CP → ΔCP → G → F, and makes the plots.
Use a fresh output directory for each run. See [calculation configuration](docs/CALCULATIONS.md)
to run your own trajectories, including `--contact-mode run-detector` with an external GetContacts
installation. All scientific inputs remain local and ignored by Git.

### Plot-only example

```sh
python -m pip install -r requirements.txt
python examples/make_demo.py
python reproduce.py verify --input-root examples/demo_inputs
python reproduce.py current --input-root examples/demo_inputs
python reproduce.py historical --input-root examples/demo_inputs --figures 4 S1 S6 S7 S9 S10
python -m unittest discover -s tests -v
```

Outputs appear in `generated/current` and `generated/historical` (PNG, SVG and PDF).
To demonstrate the larger regional histogram panels separately:

```sh
python reproduce.py historical --input-root examples/demo_inputs --figures S3 S4 S5
```

The example generator refuses to overwrite a nonempty input folder. Reuse existing examples,
or pass a new `--out` directory. Plotting commands overwrite their own generated figures.
Use `--output-root` to retain separate runs. `--figures` selects panels only in `historical` mode.

## Code map

| File | Responsibility |
|---|---|
| `run_study.py` | Sequential coordinate → contact → score → figure launcher |
| `calculate.py` | Explicit trajectory inputs, GetContacts commands, intermediate outputs |
| `analysis/geometry.py` | Cα distances, kink axes, regional centroids, heavy-atom minima, ligand alignment/metrics |
| `analysis/statistics.py` | Full/half S, microswitch block CIs, reference-oriented/WT-referenced shifts, correlations and BH q |
| `analysis/contacts.py` | Deduplicated contact events → CP → ΔCP → class consensus → F |
| `reproduce.py` | CLI, file loading, S calculation, assay heatmap, sensitivity comparison, ligand bootstrap |
| `extract_observables.py` | Cα distances and six-Cα SVD kink angles |
| `vendor/current_replot.py` | Kink and activation-marker plots |
| `vendor/analyze_interhelix_regions.py` | Regional summary statistics and Spearman correlations |
| `vendor/plot_regional_histogram_matrices.py` | Regional distributions and class-consensus arrows |
| `vendor/network_plotting_common.py` | Contact networks and across-criterion intersection |
| `vendor/make_datasetB_figure5.py` | Connected-component diagrams |
| `structural/compose_publication_map.py` | Annotates an externally supplied VMD rendering; does not run VMD |
| `example_label.py` | Adds synthetic-example labels to figures |

`vendor` contains extracted study plotting helpers.
Dependencies are installed separately.

Place your inputs in `inputs/` (ignored by Git), following [the input schema](docs/INPUTS.md).

```sh
python reproduce.py current --input-root inputs --output-root generated/my_run
python reproduce.py historical --input-root inputs --figures 4 S6
```

Optional trajectory reading requires MDAnalysis and the study's completed GPCRmd package layout:

```sh
python -m pip install -r requirements-trajectories.txt
python extract_observables.py /path/to/processed_trajectories --check-only
python extract_observables.py /path/to/processed_trajectories --out inputs/data/current
python reproduce.py current --input-root inputs
```

The file extractor is study-specific:
282 receptor residues, 10,000 frames, 0.1-ns spacing, and ten 100-ns blocks. 

## Scientific scope and limitations

- Structural shift: `S = (mean dAg − mean dIAg)mutant − (mean dAg − mean dIAg)WT`, separately per dataset.
- Kink: acute angle between SVD axes through six Cα atoms on each side of the specified site, excluding the central residue.
- Kink error bars: SD of ten block means; distance error bars: framewise sample SD.
- `calculate.py` derives dataset-specific S and intermediate arrays/tables directly from coordinates.
- GetContacts detection is an external dependency; the repo supplies commands, criteria, event counting and all subsequent contact-score calculations.
- Ligand RMSD is calculated after TM Cα alignment to an explicit reference, with ligand–Asp85 distance and ≤4 Å occupancy.
- Microswitch outputs include 21 minimum-heavy-atom distances, shared-bin peak shifts, median shifts, block-median difference CIs and WT-referenced responses.

See [figure coverage](docs/FIGURE_MAP.md).
