# Run 3 scouting dijet analysis: signal studies 

## Reco-level study (with and without gen-matching): ScoutingNano 

To perform reco-level studies (with and without gen-matching), a script is provided to study the forward jets from the VBF topology and to optimize the reconstruction of the central dijet pair.

How to run:
``` 
python3 studyVBFReco_scoutNano.py --decay 2B --mass 500  --requireTriggerBaseline --centralPairAlgo leading  --centralEtaMax 2.5  --output VBFHTo2B_M500_recoVBF_scoutNano_leading_massDepSel.root
```

The option `requireTriggerBaseline` applies the trigger selection. There is also the possibility to simply save the tribber bit decision and apply it offline..

The algorithm for pairing the central dijet system can be chosen `--centralPairAlgo`:
* leading
* minDR
* minDPhi
* massClosest

The option `centralPairSelection` allows to choose between a newly tested `massDependent` approach and a more standard angular approach to apply a selection on the central pair:
* angular mode uses `dEta < args.centralPairDEtaMax and dPhi > args.centralPairDPhiMin`;
* massDependent mode uses `dR > dr_min and apt < apt_max`.

By default, a dR and leading/subleading jet compatibility with a dependency on the mass is applied:
* `DeltaR > dR_min(mass)`
* `A_pt < A_pt_max(mass)`, with A_pt = |pT1 - pT2| / (pT1 + pT2)

with reference values of:
* 300 GeV -> (1.5, 0.60)
* 500 GeV -> (2.0, 0.45)
* 1000 GeV -> (2.5, 0.35)

This behaviour can be overwritten:
* `--centralPairDRMin`
* `--centralPairPtImbalanceMax`


A plotter to examine the output is provided:
``` 
python3 plot_recoVBF_scoutNano.py --input VBFHTo2B_M500_recoVBF_scoutNano.root --outdir  /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B//M500/RECOLevel --label "VBFHTo2B M500 (2024)"
``` 

Signal mass histograms (inclusive and matched only) are saved as output of the plotting step, e.g.: `VBF-dijetMass-Histos_ForFIT/vbf-m500-leading.root`, with a 5 GeV bining:
* `h_signal_peak_selected_5GeV`	
* `h_signal_peak_matched_5GeV`

