#!/bin/tcsh

setenv ISGEN 0
setenv MAKETARBALL 0

#setenv SAMPLE TCP
#setenv SAMPLE ST
#setenv SAMPLE DYJetsToLL
#setenv SAMPLE DYJetsToLLNLO
setenv SAMPLE TTJets
#setenv SAMPLE Diboson

#setenv OutputPrefix ./
setenv OutputPrefix root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/TCPAnalysis/Plots/

setenv CMSSW_BASE /uscms_data/d3/jingyu/TCP/boostedDiTauReco/CMSSW_8_0_30

if ($MAKETARBALL == 1) then 
    cd $CMSSW_BASE/src

    tar -zcvf ../../CMSSW.tgz ../../CMSSW_8_0_30/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"

    eosrm /eos/uscms/store/user/zhangj/TCPAnalysis/CMSSW.tgz

    xrdcp ../../CMSSW.tgz root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/TCPAnalysis/CMSSW.tgz

    cd $CMSSW_BASE/src/BoostedDiTau/BoostedDiTauReco/SubmitMiniAODSIM
endif

foreach Mass (`ls filelists/$SAMPLE`)
#foreach Mass (WW)
    setenv MASS $Mass
    setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
    echo $NQueue
    echo ./filelists/$SAMPLE/$MASS/${SAMPLE}_${MASS}_Process.txt $OutputPrefix
    if ($ISGEN == 1) then
	condor_submit condorGen.jdl
    else
	condor_submit condor.jdl
    endif
end
