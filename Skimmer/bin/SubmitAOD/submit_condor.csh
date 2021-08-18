#!/bin/tcsh

setenv CMSSW_BASE /uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_10_6_12

cd $CMSSW_BASE/src

tar -zcvf ../../CMSSW.tgz ../../CMSSW_10_6_12/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"

eosrm /eos/uscms/store/user/zhangj/events/ALP/CMSSW.tgz

xrdcp ../../CMSSW.tgz root://cmseos.fnal.gov//store/user/zhangj/events/ALP/CMSSW.tgz

cd $CMSSW_BASE/src/BoostedDiTau/Skimmer/bin/SubmitAOD/

set cfgDir="./configs/"
foreach MASS (30 50)
    foreach JOB (`seq 1 100`)
	setenv CFG ${cfgDir}ALP_m${MASS}_w1_htjmin400_RunIISummer19UL17RECO_MINIAODSIM_${JOB}.py
	#setenv CFG ${cfgDir}UpsilonToTauTau_pthatmin400_RunIISummer19UL17RECO_MINIAODSIM_${JOB}.py
	setenv JOBNUMBER m${MASS}_j${JOB}
	#setenv JOBNUMBER upsilon_j${JOB}
	echo $CFG
	condor_submit condor.jdl
    end
end
