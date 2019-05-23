#!/bin/tcsh
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.csh  ## if a bash script, use .sh instead of .csh
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/CMSSW.tgz CMSSW.tgz
tar -xf CMSSW.tgz
rm CMSSW.tgz
setenv SCRAM_ARCH slc6_amd64_gcc530
cd CMSSW_8_0_30/src/
scramv1 b ProjectRename
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
cd BoostedDiTau/BoostedDiTauReco/SubmitMiniAODSIM
echo "Arguments passed to this script are: for 1: $1"
cmsRun ${1}
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_8_0_30
