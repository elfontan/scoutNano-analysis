# Run 3 scouting dijet analysis: VBF-category signal studies

This directory contains the signal-side studies used to define and validate the VBF category in the scouting dijet analysis.

The workflow is organized in layers:

1. GEN-level validation of the VBF topology;
2. reco-level reconstruction of a central dijet plus forward VBF pair;
3. optimization of the VBF working point;
4. checks against non-VBF signal contamination;
5. purity studies based on LHE matching;
6. standalone signal-shape visualization and fitting.

## Directory content

### `GENLevel/`

Generator-level studies for VBF signals.

Main scripts:

- `studyVBFSignals_miniaod.py`: sanity checks on MINIAODSIM when ScoutNano is
  not yet available;
- `studyVBFSignals_scoutNano.py`: GEN-level validation directly on ScoutNano;
- `plot_genVBFHiggs_scoutNano.py`: plotting helper for the GEN-level outputs.

## Notes

- The scripts assume the local scouting analysis environment, including ROOT, NanoAODTools, and `correctionlib`.
- Some outputs are stored locally while most plots are exported to EOS/web areas. Adapt commands to your working environment and preferredoutput locations, where needed

### Main reco-level scripts

- `studyCentralWideJetVBF_scoutNano.py`
  Builds the central-widejet plus forward-pair reco study output tree and
  summary histograms.
- `plot_centralWideJetVBF_scoutNano.py`
  Plots the content of the reco-level study output.
- `scan_jointVBF_vs_nonVBF_500_scoutNano.py`
  Scans forward cuts on a VBF and a non-VBF 500 GeV sample to identify useful
  working points.
- `plot_scan_jointVBF_vs_nonVBF_500_scoutNano.py`
  Visualization of the working-point scan.
- `summarize_nonVBFSignalContamination_scoutNano.py`
  Measures how often non-VBF signal samples pass the VBF selection.
- `plot_nonVBFSignalContamination_scoutNano.py`
  Plots non-VBF contamination either by working point or by signal family.
- `summarize_centralWideJetVBF_vbfMatching.py`
  Summarizes Central-to-VBFTag purity with respect to LHE VBF quarks.
- `plot_centralWideJetVBF_vbfPuritySummary.py`
  Plots the purity summary versus signal mass.
- `plot_savedWideJetMassesAll.py`
  Overlays saved signal peaks from several masses.
- `fit_signalPeak_VBF.py`
  Performs simple signal-peak fits for the VBF category.

### `OUTDATED/`

Older scan and reco-study scripts kept only for reference. They should not be used as the default workflow unless a historical comparison is needed.

## Naming convention

Current reco-level ROOT outputs typically follow a pattern like:

```text
CentralVBFHTo2B_M500_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root
```

The encoded working point is:

- central pair: `centPt30-Deta1p3`
- forward VBF selection: `fwdPt24_VBF-Deta4-Mjj500`

In the current interpretation this means:

- central pair built from jets with `pT > 30 GeV`;
- central pair constrained to `|DeltaEta| < 1.3`;
- forward pair required to satisfy `pT > 24 GeV`;
- VBF separation `|DeltaEta| > 4.0`;
- additional forward-pair mass requirement `mjj > 500 GeV`.

## Typical workflow

### 1. GEN-level study

MINIAOD-based sanity check:

```bash
python3 GENLevel/studyVBFSignals_miniaod.py \
  --input "PATH-To/RunIII2024Summer24MiniAOD*root" \
  --out VBFHTo2B_M1000_genHiggs.root
```

ScoutNano-based GEN study:

```bash
python3 GENLevel/studyVBFSignals_scoutNano.py \
  --decay 2B \
  --mass 500 \
  --vbfTagMode maxDEta
```

Plotting:

```bash
python3 GENLevel/plot_genVBFHiggs_scoutNano.py \
  --input VBFHTo2B_M500_genHiggs_scoutNano_VBFTagMode-maxDEta.root \
  --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/VBFCategory-GENLevel/M500 \
  --label "VBFHTo2B M500 (2024)" \
  --vbfTagMode maxDEta
```

### 2. Reco-level central-widejet plus VBF study

Example:

```bash
python3 studyCentralWideJetVBF_scoutNano.py \
  --decay 2B \
  --mass 500 \
  --requireTriggerBaseline \
  --forwardPairAlgo maxDEta \
  --widejetPairDEtaMax 1.3 \
  --forwardJetPtMin 24 \
  --forwardPairDEtaMin 4 \
  --forwardPairExtraMjjMin 500 \
  --output CentralVBFHTo2B_M500_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root
```

This step produces the event tree used by the downstream contamination and
purity summaries.

Plotting:

```bash
python3 plot_centralWideJetVBF_scoutNano.py \
  --input CentralVBFHTo2B_M500_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root \
  --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/CentralSamples-VBFCategory-FinalSELECTION/M500 \
  --label "VBFHTo2B M500 (2024)"
```

### 3. Joint VBF versus non-VBF optimization

The current scan tool optimizes the Central-to-VBFTag selection on one VBF and
one non-VBF 500 GeV sample:

```bash
python3 scan_jointVBF_vs_nonVBF_500_scoutNano.py \
  --outputCsv scan_jointVBF_vs_nonVBF_500.csv \
  --outputParetoCsv scan_jointVBF_vs_nonVBF_500_pareto.csv
```

Then plot the scan:

```bash
python3 plot_scan_jointVBF_vs_nonVBF_500_scoutNano.py \
  --csv scan_jointVBF_vs_nonVBF_500.csv
```

### 4. Non-VBF signal contamination

The contamination study checks how often non-VBF signals pass the VBF-category
selection. The currently used presets are:

- `fwdPt24_VBF-Deta4-Mjj500`
- `fwdPt21_VBF-Deta4p5-Mjj450`
- `fwdPt20_VBF-Deta5`

Examples:

```bash
python3 summarize_nonVBFSignalContamination_scoutNano.py \
  --sampleSet zprime_qq \
  --workingPointPreset fwdPt24_VBF-Deta4-Mjj500 fwdPt21_VBF-Deta4p5-Mjj450 fwdPt20_VBF-Deta5 \
  --outputCsv summary_nonVBFSignalContamination_zprime_3wps.csv \
  --outputRoot summary_nonVBFSignalContamination_zprime_3wps.root
```

```bash
python3 summarize_nonVBFSignalContamination_scoutNano.py \
  --sampleSet glugluspin0_2b \
  --workingPointPreset fwdPt24_VBF-Deta4-Mjj500 fwdPt21_VBF-Deta4p5-Mjj450 fwdPt20_VBF-Deta5 \
  --outputCsv summary_nonVBFSignalContamination_gluglu2b_3wps.csv \
  --outputRoot summary_nonVBFSignalContamination_gluglu2b_3wps.root
```

Plotting the combined outputs:

```bash
python3 plot_nonVBFSignalContamination_scoutNano.py \
  --csv summary_nonVBFSignalContamination_gluglu2b_3wps.csv \
  --csv summary_nonVBFSignalContamination_rsgrav2glu_3wps.csv \
  --csv summary_nonVBFSignalContamination_rsgrav2q_3wps.csv \
  --csv summary_nonVBFSignalContamination_zprime_3wps.csv \
  --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/nonVBFSignalContamination/ \
  --label "2024 (13.6 TeV)"
```

The current plotter provides both views automatically:

- one canvas per VBF working point with all signal families overlaid;
- one canvas per signal family with the different VBF selections overlaid.

### 5. VBF purity

Purity is evaluated only for the Central-to-VBFTag strategy in the current
setup. The summary and plot scripts can now auto-expand the input ROOT files
from a single working-point tag:

```bash
python3 summarize_centralWideJetVBF_vbfMatching.py \
  --selection centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500 \
  --writeCSV \
  --csv /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/CentralSamples-VBFCategory-FinalSELECTION/summary_centralWideJetVBF_vbfMatching_fwdPt24_VBF-Deta4-Mjj500.csv
```

```bash
python3 plot_centralWideJetVBF_vbfPuritySummary.py \
  --selection centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500 \
  --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/CentralSamples-VBFCategory-FinalSELECTION \
  --label "VBFHTo2B" \
  --outname summary_vbfPurity_vs_mass_fwdPt24_VBF-Deta4-Mjj500 \
  --csv /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/CentralSamples-VBFCategory-FinalSELECTION/summary_vbfPurity_vs_mass_fwdPt24_VBF-Deta4-Mjj500.csv
```

This is usually more convenient than listing all mass files explicitly by hand.

### 6. Standalone signal-shape visualization and fitting

Overlay the saved mass peaks:

```bash
python3 plot_savedWideJetMassesAll.py \
  --workingPoint fwdPt24_VBF-Deta4-Mjj500 \
  --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/CentralSamples-VBFCategory-FinalSELECTION \
  --label "2024 (13.6 TeV)"
```

Fit the signal peaks:

```bash
python3 fit_signalPeak_VBF.py \
  --indir ./ \
  --algo fwdPt24_VBF-Deta4-Mjj500 \
  --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/CentralSamples-VBFCategory-FinalSELECTION/
```

