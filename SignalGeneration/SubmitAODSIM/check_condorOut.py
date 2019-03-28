import os, sys, ROOT
import numpy as np

masses=[10, 30, 50]

jobs=np.linspace(100,1,100)

outputPrefix="root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer16DR80Premix/"

for mass in masses:
    for job in jobs:
        filename="ALP_m"+str(mass)+"_w1_htjmin400_RunIISummer16DR80Premix_AODSIM_"+str(int(job))+".root"

        f=ROOT.TFile(outputPrefix+filename)

        if f.IsZombie():
            print filename

        
