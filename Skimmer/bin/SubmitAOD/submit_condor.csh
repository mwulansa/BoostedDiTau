#!/bin/tcsh


setenv MAKETARBALL 1
setenv SIGNAL 0
setenv BKG 1


#setenv SAMPLE DYJetsToLL
setenv SAMPLE QCD


setenv CMSSW_BASE /uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_9_4_15

cd $CMSSW_BASE/src

tar -zcvf ../../CMSSW.tgz ../../CMSSW_9_4_15/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"


setenv CMSSW_BASE /uscms/home/mwulansa/nobackup/TCP/boostedDiTauReco/CMSSW_10_6_16

if ($MAKETARBALL == 1) then

    cd $CMSSW_BASE/src


    tar -zcvf ../../CMSSW.tgz ../../CMSSW_10_6_16/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout" --exclude="*.txt"

    eosrm /eos/uscms/store/user/mwulansa/DIS/TCPAnalysis/CMSSW.tgz

    xrdcp ../../CMSSW.tgz root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCPAnalysis/CMSSW.tgz

    cd $CMSSW_BASE/src/BoostedDiTau/Skimmer/bin/SubmitAOD/

endif

if ($SIGNAL == 1 ) then
    set cfgDir="./configs/"
    foreach MASS (10)
    #foreach MASS (10 30 50)
	foreach JOB (`seq 1 100`)
	#foreach JOB (2 68 92 93 94 95 96 97 98 99 100)
	    setenv CFG ${cfgDir}ALP_m${MASS}_w1_htjmin400_RunIISummer17DR94Premix_miniAODSIM_${JOB}.py
	    eosrm store/user/mwulansa/DIS/TCPAnalysis/RunIISummer17DR94Premix/ALP_m${MASS}_w1_htjmin400_RunIISummer17DR94Premix_MINIAODSIM_Slimmed_${JOB}.root
	    setenv JOBNUMBER m${MASS}_j${JOB}
	    echo $CFG
	    condor_submit condor_signal.jdl
	end
    end
endif


if ($BKG == 1) then
    foreach Mass (`ls filelists/$SAMPLE`)
	setenv MASS $Mass
	setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
	echo $NQueue
	echo ./configs/${SAMPLE}_${MASS}_Process.py

	condor_submit condor.jdl
    end
endif
