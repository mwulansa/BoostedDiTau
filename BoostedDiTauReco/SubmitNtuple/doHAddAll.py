import sys,os

#dataset = sys.argv[1]
bkgSample = ['DYJetsToLL_M-50','DYJetsToLL_M-4to50','QCD','WJetsToLNu']

for sample in bkgSample:
    os.system("./runDoHAdd.sh "+sample+"_HTBinned")

os.system("./runDoHAddFlat.sh")
#os.system("python3 doHadd.py "+dataset)
