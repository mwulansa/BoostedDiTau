#!/bin/tcsh

setenv CMSSW_BASE /uscms_data/d3/jingyu/TCP/Generator/CMSSW_9_4_9

cd $CMSSW_BASE/src

#tar -zcvf ../../CMSSW.tgz ../../CMSSW_9_4_9/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"

eosrm /eos/uscms/store/user/zhangj/events/ALP/CMSSW.tgz

xrdcp ../../CMSSW.tgz root://cmseos.fnal.gov//store/user/zhangj/events/ALP/CMSSW.tgz

cd $CMSSW_BASE/src/BoostedDiTau/SignalGeneration/SubmitAODSIM

set cfgDir="./configs/"
#foreach MASS (10 30 50)
foreach MASS (10)
    #foreach JOB (`seq 1 100`)
    #foreach JOB (41 49)
    #foreach JOB (19 23 49 50 51 52 56 57 80)
    #foreach JOB (4 30 38 39 3 6 82 83 84 86 88 9)
    foreach JOB (41)
	setenv CFG ${cfgDir}ALP_m${MASS}_w1_htjmin400_RunIISummer17DR94Premix_AODSIM_${JOB}.py
	setenv JOBNUMBER m${MASS}_j${JOB}
	echo $CFG
	condor_submit condor.jdl
    end
end
