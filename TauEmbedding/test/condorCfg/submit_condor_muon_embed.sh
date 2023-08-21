#!/bin/bash

export CMSSW_BASE=/uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_12_1_1

cd $CMSSW_BASE/src

tar -zcvf ../../CMSSW.tgz ../../CMSSW_12_1_1/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"

xrdcp -f ../../CMSSW.tgz root://cmseos.fnal.gov//store/user/zhangj/CMSSW.tgz

cd $CMSSW_BASE/src/BoostedDiTau/TauEmbedding/test/condorCfg/

for i in `seq 1 1 500`
#for i in `seq 398 1 398`
do
    export INPUT=YMuMu_pth400_${i}.root
    echo Arguments: ${INPUT}
    condor_submit condor_muon_embed.jdl
done
