# Signal MC studies

This directory groups studies performed on scouting signal MC samples to understand acceptance, reconstruction, category migrations, and signal shapes.

## Layout

### `GenLevel-studies/`

Generator-level sanity checks when new signals are (privately) generated and comparisons across formats.

Main scripts:

- `mcAna_genLevelStudies.py`
- `mcAna_genLevelStudies_wAccCuts.py`
- `mcAna_miniAOD.py`
- `plot_mcAna_genLevel.py`
- `plot_mcAna_massComp.py`
- `plot_mcAna_miniAOD.py`
- `compPlot_mcAna_genLevel.py`

This area is useful for validating signal samples and preliminarily checking the signal topology before moving to full reco-level studies.

### `VBF-CAT/`

Signal studies for the VBF category on ScoutNano VBF and non-VBF signal MC.
This area is currently including:

- GEN-level VBF signal studies;
- reco-level central-widejet plus forward-pair reconstruction;
- working-point scans balancing VBF efficiency against non-VBF contamination;
- non-VBF contamination summaries and plots;
- VBF purity summaries based on LHE matching;
- standalone signal-shape plotting and fitting.

See [VBF-CAT/README.md](./VBF-CAT/README.md) for the detailed workflow.

### `FullyRES-CAT/`

Signal studies for the fully resolved pp-like category.

Main scripts:

- `fullyResCategory_dijetMass_ZPrimeToQQ.py`: builds dijet-mass outputs for the
  resolved category;
- `plot_allSignals_ZPrimeToQQ.py`: overlays signal shapes;
- `fit_signalPeak_ZPrimeToQQ.py`: performs simple signal-peak fitting.

