# ScoutingNano analysis framework

This framework allows to run on Run 3 ScoutingNano data exploiting RDF within the cmgrdf-prototype (included as back-end).

## Setup of the working area

Clone this repository with the --recursive option, so cmgrdf-prototype is also cloned.
Then, follow the setup instructions as described in [that repository](https://gitlab.cern.ch/cms-new-cmgtools/cmgrdf-prototype).

After the first installation, the two following commands are needed to properly configure the environment:s
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_106a_cuda/x86_64-el9-gcc11-opt/setup.sh
eval $(make env)
```

Copy the script `run_simple_scoutNano_data.py`:
```
cp run_simple_scoutNano_data.py cmgrdf-prototype/examples/.
cd cmgrdf-prototype/examples
python3 run_simple_scoutNano_data.py
```