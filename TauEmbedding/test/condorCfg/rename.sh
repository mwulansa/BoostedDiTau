#!/bin/bash

for i in `seq 1 1 500`
do
    eos root://cmseos.fnal.gov mv /eos/uscms/store/user/zhangj/events/UpsilonTauTau/UL2017/YTauTau_${i}.root /eos/uscms/store/user/zhangj/events/UpsilonTauTau/UL2017/YTauTau_pth400_${i}.root
done
