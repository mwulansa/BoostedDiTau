#!/bin/tcsh

setenv ISGEN 0
setenv MAKETARBALL 1

#setenv SAMPLE ST
#setenv SAMPLE DYJetsToLLNLO
#setenv SAMPLE TTJets

#setenv SAMPLE TCP
#setenv SAMPLE DYJetsToLL
#setenv SAMPLE QCD_Pt
#setenv SAMPLE Diboson
setenv SAMPLE WJetsToLNu

#setenv OutputPrefix ./
setenv OutputPrefix root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/TCPAnalysis/Plots94X/

#setenv CMSSW_BASE /uscms_data/d3/jingyu/TCP/boostedDiTauReco/CMSSW_8_0_30
setenv CMSSW_BASE /uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_9_4_15

if ($MAKETARBALL == 1) then 
    cd $CMSSW_BASE/src

    tar -zcvf ../../CMSSW.tgz ../../CMSSW_9_4_15/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"

    eosrm /eos/uscms/store/user/zhangj/TCPAnalysis/CMSSW.tgz

    xrdcp ../../CMSSW.tgz root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/TCPAnalysis/CMSSW.tgz

    cd $CMSSW_BASE/src/BoostedDiTau/BoostedDiTauReco/SubmitMiniAODSIM
endif

foreach Mass (`ls filelists/$SAMPLE`)
#foreach Mass (WW)
    setenv MASS $Mass
    setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
    echo Number of Jobs: $NQueue
    echo ./filelists/$SAMPLE/$MASS/${SAMPLE}_${MASS}_Process.txt $OutputPrefix
    if ($ISGEN == 1) then
	condor_submit condorGen.jdl
    else
	condor_submit condor.jdl
    endif
end
