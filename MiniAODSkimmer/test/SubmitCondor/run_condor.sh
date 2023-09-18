#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a bash script, use .sh instead of .csh
echo ">>> copying code package from EOS to local..."
xrdcp root://cmseos.fnal.gov//store/user/zhangj/CMSSW.tgz CMSSW.tgz

echo ">>> unzipping the package locally..."
tar -xf CMSSW.tgz
echo ">>> finish unzipping the package!"
rm CMSSW.tgz

export SCRAM_ARCH=slc7_amd64_gcc900
cd CMSSW_12_1_1/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
cd BoostedDiTau/MiniAODSkimmer/test/SubmitCondor/

echo "Arguments passed to this script are: "
echo "  for 1 (input): $1"
echo "  for 2 (output): $2"

cmsRun rerunTauRecoOnMiniAOD_WithClean_Custom.py file:root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/UL2017/TCP_m_50_w_1_htj_400toInf_slc6_amd64_gcc630_MINIAOD/${1}

xrdcp -f mttNtuple.root root://cmseos.fnal.gov//eos/uscms/store/group/lpcsusyhiggs/Jingyu/svfit/2017/TCP_m_50_w_1_htj_400toInf/${2}

cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_12_1_1
