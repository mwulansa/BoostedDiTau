import subprocess
import sys,string,math,os
import ConfigParser
import glob
import numpy as np
from sampleAndMasses import *

filesPerList=50


def checkAndMakeDir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def clearDir(dir):
    for fil in glob.glob(dir+"/*"):
        os.remove(fil)

if __name__ == "__main__":
    if Sample == 'TCP':
        for mass in masses:
            fileListDir="./filelists/TCP/"+mass+"/"
            checkAndMakeDir("./filelists/TCP/")
            checkAndMakeDir(fileListDir)
            clearDir(fileListDir)
            for job in range(1, 101): 
                filelistIdx=int((job-1)/filesPerList)
                if job%filesPerList==1:
                    out=open(fileListDir+"TCP_"+mass+"_"+str(filelistIdx)+".txt","w")
                out.write(prefix+"ALP_"+mass+"_w1_htjmin400_RunIISummer17DR94Premix_MINIAODSIM_Cleaned_"+str(int(job))+".root\n")
                
    else:
        for TMass in masses:
            fileListDir="./filelists/"+Sample+"/"+TMass+"/"
            checkAndMakeDir("./filelists/"+Sample)
            checkAndMakeDir(fileListDir)
            clearDir(fileListDir)
            #print mass
            searchString=preSearchString.replace("REPLACEME",TMass)
            os.system('dasgoclient --query "'+searchString+'"')
            query = 'dasgoclient --query "file dataset='+searchString+' instance=prod/phys03"'
            #print query
            files=os.popen(query).read().split()
            
            for nf in range(1, len(files)+1):
                filelistIdx=int((nf-1)/filesPerList)
                if nf%filesPerList==1:
                    out=open(fileListDir+Sample+"_"+TMass+"_"+str(filelistIdx)+".txt","w")
                out.write(prefix+files[nf-1]+"\n")
            
