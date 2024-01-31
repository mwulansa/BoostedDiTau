#!/bin/tcsh

setenv MAKETARBALL 1

setenv CMSSW_BASE /uscms/home/mwulansa/nobackup/TCPNtuple/CMSSW_12_1_0_pre3

if ($MAKETARBALL == 1) then

    cp EMu_OS_BTag_Efficiency.root $CMSSW_BASE/src

    cd $CMSSW_BASE/src

    tar -zcvf ../../CMSSW_12X.tgz ../../CMSSW_12_1_0_pre3/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.Log" --exclude="*stderr" --exclude="*stdout" --exclude="*.log" --exclude="../../CMSSW_12_1_0_pre3/src/BoostedDiTau/MiniAODSkimmer/test/crabConfig" --exclude="../../CMSSW_12_1_0_pre3/src/BoostedDiTau/BoostedDiTauReco/SubmitNtuple/filelists_bkp" --exclude="../../CMSSW_12_1_0_pre3/src/BoostedDiTau/BoostedDiTauReco/SubmitNtuple/filelists_old"

    eosrm /store/user/mwulansa/DIS/TCPAnalysis/CMSSW_12X.tgz

    cp ../../CMSSW_12X.tgz ../../tarToUpdate

    rm ../../CMSSW_12X.tgz

    cp EMu_OS_BTag_Efficiency.root ../../tarToUpdate

    cd ../../tarToUpdate
    
    ls -l

    tar -zxvf CMSSW_12X.tgz

    cp EMu_OS_BTag_Efficiency.root  CMSSW_12_1_0_pre3/src/BoostedDiTau/BoostedDiTauReco/SubmitNtuple/

    tar -zcvf ../CMSSW_12X.tgz CMSSW_12_1_0_pre3/

    xrdcp ../CMSSW_12X.tgz root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCPAnalysis/CMSSW_12X.tgz

    cd $CMSSW_BASE/src/BoostedDiTau/BoostedDiTauReco/SubmitNtuple

endif
