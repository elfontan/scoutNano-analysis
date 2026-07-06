# Run 3 scouting dijet analysis: signal studies 


## Gen-level study: MINIAODSIM

If ScoutingNano are not available yet, a script to perform sanity checks on VBF signals on MINIAODSIM is provided.

How to run:
``` 
python3  GENLevel/studyVBFSignals_miniaod.py --input "PATH-To/RunIII2024Summer24MiniAOD*root" --out VBFHTo2B_M1000_genHiggs.root
```


## Gen-level study: ScoutingNano 

The same type of checks from the script above is implemented to run on ScoutingNano.

Some options can be specified from command line:
* mass;
* flavour can be specified: 2B, 2C, 2Q, 2Glu;
* a specific algorithm for taggin VBF, e.g. maxDEta or maxMjj

How to run:
``` 
python3 GENLevel/studyVBFSignals_scoutNano.py --decay 2B --mass 500 --vbfTagMode maxDEta
``` 

A plotter to examine the output is provided:
``` 
python3 GENLevel/plot_genVBFHiggs_scoutNano.py --input VBFHTo2B_M500_genHiggs_scoutNano_VBFTagMode-maxDEta.root --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/VBFCategory-GENLevel/M500 --label "VBFHTo2B M500 (2024)" --vbfTagMode maxDEta
``` 


## Reco-level study (with and without gen-matching): ScoutingNano 

Moving to reco-level studies (with and without gen-matching), a script is provided to study the forward jets from the VBF topology and to optimize the reconstruction of the central dijet pair.

How to run:
``` 
python3 studyCentralWideJetVBF_scoutNano.py --decay 2B --mass 500 --requireTriggerBaseline --forwardPairAlgo maxDEta --centralPairDEtaMax 1.3 --forwardJetPtMin 24 --forwardPairDEtaMin 4 --forwardPairExtraMjjMin 500 --output CentralVBFHTo2B_M500_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root
```

The option `requireTriggerBaseline` applies the trigger selection. There is also the possibility to simply save the trigger bit decision and apply it offline.

Various parameters can be tuned to optimize either the central pair selection:

* forwardPairAlgo (maxDEta is assumed as default)
* centralPairDEtaMax (1.3 is assumed as default)

or the VBF selection:

* forwardJetPtMin
* forwardPairDEtaMin
* forwardPairExtraMjjMin

A plotter to examine the output is provided:
``` 
python3 plot_centralWideJetVBF_scoutNano.py --input CentralVBFHTo2B_M500_centralWideJetVBF_scoutNano_centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.root --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/CentralSamples-VBFCategory-FinalSELECTION/M500  --label "VBFHTo2B M500 (2024)"
``` 



## Standalone signal modeling

From the output of the previous step, signal peaks at different mass hypotheses can be plotted together:

```
python3 plot_savedWideJetMassesAll.py --workingPoint fwdPt24_VBF-Deta4-Mjj500 --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/CentralSamples-VBFCategory-FinalSELECTION --label "2024 (13.6 TeV)" (--plotMatched)
```
and a basic signal modelling based on a double-sided Crystal Ball can be done for every signal mass under consideration:
```
python3 fit_signalPeak_VBF.py --indir ./ --algo fwdPt24_VBF-Deta4-Mjj500 --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/CentralSamples-VBFCategory-FinalSELECTION/
```