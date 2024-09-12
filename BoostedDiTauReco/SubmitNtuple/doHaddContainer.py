import subprocess
import sys,string,math,os
import glob
import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To copy output histograms from eos space in preparation for hadd")
    parser.add_argument("--era", type=str, help="Era. e.g. 2016preVFP, 2016postVFP, 2017, 2018")
    parser.add_argument("-v", "--version", type=str, help="version. User defined")
    parser.add_argument("--sample", nargs="+", type=str, help="sample. e.g. DYJetsToLL_M-50, TTTo2L2Nu, etc. Able take more than one")
    parser.add_argument("--fname", type=str, help="output name. e.g. plotBoostedTauTau")
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    version = args.version
    fname = args.fname
    plotDir = "./output/"+version        

    for s in args.sample:

        smpl = s
        era = args.era
        if era == "2016preVFP" or era == "2016postVFP":
            hist = "h_"+fname+"_"+era+"_"+smpl
        else:
            hist = "h_"+fname+"_"+era+"_Ntuple_"+smpl
        
        searchString = hist+"*"

        print('hadd '+searchString.replace("*","_"+version+".root")+' '+plotDir+'/'+searchString)                   
        os.system('hadd '+searchString.replace("*","_"+version+".root")+' '+plotDir+'/'+searchString)   

        # parts = os.listdir('filelists/'+s+'/UL'+era+'/')
        # print(parts)

        # for part in parts:
        #     p = part.split("_")[-3]
        #     print(p)
        #     print('hadd '+searchString.replace("*","_"+p+"_"+version+".root")+' '+plotDir+'/'+searchString+p+"*")
        #     os.system('hadd '+searchString.replace("*","_"+p+"_"+version+".root")+' '+plotDir+'/'+searchString+p+"*")


        
