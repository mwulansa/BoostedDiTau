#!/bin/tcsh
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.csh  ## if a bash script, use .sh instead of .csh
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/TCPAnalysis/CMSSW.tgz CMSSW.tgz
tar -xf CMSSW.tgz
rm CMSSW.tgz
setenv SCRAM_ARCH slc7_amd64_gcc630
cd CMSSW_9_4_15/src/
scramv1 b ProjectRename
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
cd BoostedDiTau/BoostedDiTauReco/SubmitMiniAODSIM
ls ./filelists
ls ./filelists/Debug
ls ./filelists/Debug/debug
echo "Arguments passed to this script are: "
echo "  for 1: $1"
echo "  for 2: $2"
echo "  for 3: $3"
python ${1} ${2} ${3}
echo /eos/uscms/store/user/zhangj/TCPAnalysis/Plots94X/h_analyzeMC_${4}
#eosrm /eos/uscms/store/user/zhangj/TCPAnalysis/Plots/h_${4}
#xrdcp h_${4} root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/TCPAnalysis/Plots/h_${4}
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_9_4_15
