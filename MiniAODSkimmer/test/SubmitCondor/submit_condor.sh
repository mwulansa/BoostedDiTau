#!/bin/bash

export CMSSW_BASE=/uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_12_1_1

cd $CMSSW_BASE/src

tar -zcvf ../../CMSSW.tgz ../../CMSSW_12_1_1/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"

xrdcp -f ../../CMSSW.tgz root://cmseos.fnal.gov//store/user/zhangj/CMSSW.tgz

cd $CMSSW_BASE/src/BoostedDiTau/MiniAODSkimmer/test/SubmitCondor/

for i in `seq 1 1 500`
do
    export INPUT=${i}
    echo Arguments: ${INPUT}
    condor_submit condor.jdl
done
