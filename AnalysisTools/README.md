# AnalysisTools

This directory collects a set of standalone analysis utilities to perform different physics tasks.

## Structure

The tools in this directory assume the local Run 3 scouting analysis environment used in this repository, including ROOT, NanoAODTools, and `correctionlib` when needed.

### `TriggerEfficiencies/`

Tools to measure and compare trigger efficiencies in scouting data and MC.
This area contains:

- the JSON inputs used by the trigger studies, including the HLT JEC payload and the golden JSON;
- scripts to compute inclusive trigger efficiencies in data and MC;
- scripts to compare categories such as VBF-like, boosted+ISR, resolved+ISR, and generic pp topologies;
- Condor submission helpers for large-scale processing.

See [TriggerEfficiencies/README.md](./TriggerEfficiencies/README.md) for the current workflow and submission notes.

### `SigMC-studies/`

Signal-MC validation and optimization studies for the scouting dijet analysis.
This area is split into:

- `GenLevel-studies/`: quick generator-level sanity checks and acceptance studies;
- `FullyRES-CAT/`: fully resolved pp-category signal-shape studies;
- `VBF-CAT/`: VBF-category reconstruction, optimization, contamination, purity, and signal-modeling studies.

See [SigMC-studies/README.md](./SigMC-studies/README.md) for a map of the available studies.