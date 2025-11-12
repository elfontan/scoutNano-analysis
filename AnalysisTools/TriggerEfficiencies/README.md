# Dijet analysis with scouting data: trigger efficiencies

## Trigger efficiency code

Setup the proxy:
```
voms-proxy-init --voms cms --valid 168:00 
cp /tmp/x509up_u101050 ~/
export X509_USER_PROXY=~/x509up_u101050
```
or:
```
voms-proxy-init --voms cms --valid 168:00 -out $HOME/private/.proxy
export X509_USER_PROXY=$HOME/private/.proxy
```

A simple script `getEffs.py` and the corresponding `plotDijetHTEff.py` are provided as an example.

The code at the moment compute trigger efficiencies for the JetHT scouting path (`PFScouting_JetHT`), using the scouting single-muon trigger (`PFScouting_SingleMuon`) as a reference, for different categories of events:
* VBF-like topology;
* boosted + ISR;
* resolved + ISR;
* rest of pp events;

as a function of:
* m_jj;
* H_T.

The trigger efficiency as a function of H_T is also studied inclusively.

## Condor submission

A setup for condor submission to study trigger efficiency on partial or full dataset is provided in the repository `condor_sub`.
The ingredients (that can be customised according to the needs) are:
* A `data2024_triggEff.submit` for condor submission. Note: err/out/log files are disabled as they occupy a lot of space and the fill up the disk available in afs quickly.
* A `run_data2024_triggEff.sh` as executable running the python script to compute trigger efficiencies for every file in the chosen list (see below).
* The python script to compute trigger efficiencies `data2024_triggEff.py`: note that the only difference with respect to the script provided in the main directory is the parsing of some options to allow to run on every single file separately.
* A list of files to consider. Two examples are provided: a reduced list (`red_list_data2024I.txt`) and the full list for 2024I data (`list_data2024I.txt`). Note that the prefix `root://eoscms.cern.ch//eos/cms` is needed.
