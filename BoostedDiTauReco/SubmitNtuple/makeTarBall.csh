#!/bin/tcsh

setenv MAKETARBALL 1

setenv CMSSW_BASE /uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_12_1_1

if ($MAKETARBALL == 1) then

    cd $CMSSW_BASE/src

    tar -zcvf ../../CMSSW_12X.tgz ../../CMSSW_12_1_1/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.Log" --exclude="*stderr" --exclude="*stdout" --exclude="*.log" --exclude="../../CMSSW_12_1_1/src/BoostedDiTau/MiniAODSkimmer/test/crabConfig"

    #eosrm /store/user/mwulansa/DIS/TCPAnalysis/CMSSW_12X.tgz

    xrdcp -f ../../CMSSW_12X.tgz root://cmseos.fnal.gov//store/user/zhangj/CMSSW_12X.tgz

    cd $CMSSW_BASE/src/BoostedDiTau/BoostedDiTauReco/SubmitNtuple

endif
