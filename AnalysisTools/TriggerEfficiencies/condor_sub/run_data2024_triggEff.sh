#!/bin/bash
# run_data2024_triggEff.sh
# Arguments: inputFile outputFile

# Set up CMS environment                                                                                                                      
# ----------------------
cd /afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/TRIGGER_EFF/CMSSW_14_0_12/src/
source /cvmfs/cms.cern.ch/cmsset_default.sh                                                                                              
eval `scramv1 runtime -sh`

INPUTFILE=$1
BASENAME=$(basename ${INPUTFILE} .root)
OUTPUTFILE="histos_${BASENAME}.root"

echo "[INFO] Running data2024_triggEff.py"
echo "       Input:  ${INPUTFILE}"
echo "       Output: ${OUTPUTFILE}"

# Run triggEff code
cd /afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/TRIGGER_EFF/condor_sub/
python3 /afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/TRIGGER_EFF/condor_sub/data2024_triggEff.py inputFile=$INPUTFILE outputFile=file:$OUTPUTFILE

# Stage output
OUTDIR=/afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/TRIGGER_EFF/condor_sub/outputHistos
#cp -f ${OUTPUTFILE} ${OUTDIR}/
mv ${OUTPUTFILE} ${OUTDIR}/
echo "[INFO] Job done. Output copied to ${OUTDIR}/${OUTPUTFILE}"

