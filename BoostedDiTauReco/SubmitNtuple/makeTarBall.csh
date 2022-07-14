#!/bin/tcsh

setenv MAKETARBALL 1

setenv CMSSW_BASE /uscms/home/mwulansa/nobackup/TCPNtuple/CMSSW_12_1_0_pre3

if ($MAKETARBALL == 1) then

    cd $CMSSW_BASE/src

    tar -zcvf ../../CMSSW_12X.tgz ../../CMSSW_12_1_0_pre3/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.Log" --exclude="*stderr" --exclude="*stdout" --exclude="*.log" --exclude="../../CMSSW_12_1_0_pre3/src/BoostedDiTau/MiniAODSkimmer/test/crabConfig"

    eosrm /store/user/mwulansa/DIS/TCPAnalysis/CMSSW_12X.tgz

    xrdcp ../../CMSSW_12X.tgz root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCPAnalysis/CMSSW_12X.tgz

    cd $CMSSW_BASE/src/BoostedDiTau/BoostedDiTauReco/SubmitNtuple

endif
