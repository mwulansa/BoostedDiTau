import subprocess
import sys,string,math,os
import ConfigParser
import glob
import numpy as np
from makeFileLists import *

isCopy = True
isGen = False

if isGen:
    plotDir="./plots/"+version+"Gen"
else:
    plotDir="./plots/"+version

checkAndMakeDir(plotDir)

if isGen:
    outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/zhangj/TCPAnalysis/Plots94X/ | grep h_Gen_"+Sample).read().split()
else:
    outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/zhangj/TCPAnalysis/Plots94X/ | grep h_analyzeMC_"+Sample).read().split()

print outputfiles
    
if isCopy:
    for fil in outputfiles:
        print fil
        os.system("xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/TCPAnalysis/Plots94X/"+fil+" "+plotDir+"/"+fil)
    
#sys.exit()
    
for mass in masses:
    if isGen:
        searchString='h_Gen_'+Sample+'_'+mass+'_*'
    else:
        searchString='h_analyzeMC_'+Sample+'_'+mass+'_*'
    print 'ls '+plotDir+'/'+searchString
    os.system('ls '+plotDir+'/'+searchString)
    if os.path.exists(searchString.replace("*",version+".root")):
        os.remove(searchString.replace("*",version+".root"))
    os.system('hadd '+searchString.replace("*",version+".root")+' '+plotDir+'/'+searchString)
    
