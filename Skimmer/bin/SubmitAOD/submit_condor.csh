#!/bin/tcsh

setenv CMSSW_BASE /uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_9_4_15

cd $CMSSW_BASE/src

tar -zcvf ../../CMSSW.tgz ../../CMSSW_9_4_15/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"

eosrm /eos/uscms/store/user/zhangj/events/ALP/CMSSW.tgz

xrdcp ../../CMSSW.tgz root://cmseos.fnal.gov//store/user/zhangj/events/ALP/CMSSW.tgz

cd $CMSSW_BASE/src/BoostedDiTau/Skimmer/bin/SubmitAOD/

set cfgDir="./configs/"
foreach MASS (10 30 50)
#foreach MASS (30)
    foreach JOB (`seq 1 100`)
    #foreach JOB (2 68 92 93 94 95 96 97 98 99 100)
    #foreach JOB (31)
	setenv CFG ${cfgDir}ALP_m${MASS}_w1_htjmin400_RunIISummer17DR94Premix_MINIAODSIM_${JOB}.py
	setenv JOBNUMBER m${MASS}_j${JOB}
	echo $CFG
	condor_submit condor.jdl
    end
end
