import sys,string,math,os,glob
import numpy as np

prefix = "root://xrootd.unl.edu/"
filesPerList=50

def checkAndMakeDir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def clearDir(dir):
    for fil in glob.glob(dir+"/*"):
        os.remove(fil)


masses = ['30']

bins = ['0to100','100to400','400toInf']

for mass in masses:
    fileListDir="./filelists/TCP/"+mass+"/"
    checkAndMakeDir("./filelists/TCP")
    checkAndMakeDir(fileListDir)
    clearDir(fileListDir)
    for b in bins:
        searchString = '/store/user/nbower/Events/TCP_m_{}_w_1_{}_slc6_amd64_gcc630_MINIAOD/'.format(mass, b)
        query = 'eos root://cmseos.fnal.gov ls '+searchString
        files = os.popen(query).read().split()

        for nf in range(1, len(files)+1):
            filelistIdx=int((nf-1)/filesPerList)
            if nf%filesPerList==1:
                out=open(fileListDir+'TCP_m_{}_w_1_htj_{}_{}.txt'.format(mass, b, str(filelistIdx)), 'w')
            out.write(prefix+searchString+files[nf-1]+"\n")

