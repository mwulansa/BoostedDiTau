#!/bin/tcsh

setenv CMSSW_BASE /uscms_data/d3/jingyu/TCP/Generator/CMSSW_8_0_25

cd $CMSSW_BASE/src

tar -zcvf ../../CMSSW.tgz ../../CMSSW_8_0_25/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"

eosrm /eos/uscms/store/user/zhangj/events/ALP/CMSSW.tgz

xrdcp ../../CMSSW.tgz root://cmseos.fnal.gov//store/user/zhangj/events/ALP/CMSSW.tgz

cd $CMSSW_BASE/src/BoostedDiTau/SignalGeneration/SubmitAODSIM

set cfgDir="./configs/"
#foreach MASS (30 50)
foreach MASS (10)
    #foreach JOB (`seq 1 100`)
    #foreach JOB (2 68 92 93 94 95 96 97 98 99 100)
    foreach JOB (31)
	setenv CFG ${cfgDir}ALP_m${MASS}_w1_htjmin400_RunIISummer16DR80Premix_AODSIM_${JOB}.py
	setenv JOBNUMBER m${MASS}_j${JOB}
	echo $CFG
	condor_submit condor.jdl
    end
end
