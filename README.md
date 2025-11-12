# ScoutingNano analysis framework

This framework allows to run on Run 3 ScoutingNano data exploiting RDF within the cmgrdf-prototype (included as back-end).

## Setup of the working area

Clone this repository with the --recursive option, so cmgrdf-prototype is also cloned.
Then, follow the setup instructions as described in [that repository](https://gitlab.cern.ch/cms-new-cmgtools/cmgrdf-prototype).

After the first installation, the two following commands are needed to properly configure the environment:s
```
voms-proxy-init --voms cms --valid 168:00
source /cvmfs/sft.cern.ch/lcg/views/LCG_106a_cuda/x86_64-el9-gcc11-opt/setup.sh
cd cmgrdf-prototype/
eval $(make env)
```

Copy the script `run_simple_scoutNano_data.py` from the main repository to `examples`:
```
cp ../run_simple_scoutNano_data.py examples/.
cd examples
python3 run_simple_scoutNano_data.py
```

A similar example to run on signal MC events is provided: `run_mcSig_scoutNano.py`.
It can be tested in the same way:
```
cp ../run_mcSig_scoutNano.py examples/.
cd examples
python3 run_mcSig_scoutNano.py
```

## AnalysisTools directory

Some script to study signal MC samples and trigger efficiencies are collected in the `AnalysisTools` directory. There are two separate repositories:

* `SigMC-studies`
* `TriggerEfficiencies`