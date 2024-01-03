import sys,os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Do HAdd. Example python3 doHAddAll --version v1 --filename studySVFit")
    parser.add_argument("-v", "--version", type=str, help="Version for doHadd")
    parser.add_argument("-f", "--filename", type=str, help="filename from condor output")
    parser.add_argument("-s", "--sampleName", type=str, help="To do HAdd to one sample. Specify the sample name")
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    print("version", args.version)

    #bkgSample = ['DYJetsToLL_M-50','DYJetsToLL_M-4to50','QCD','WJetsToLNu']
    #bkgSample = ['DYJetsToLL_M-50','DYJetsToLL_M-4to50','WJetsToLNu','QCD']

    bkgSample = []

    if args.sampleName is not None:
        bkgSample+=[args.sampleName]
        print(bkgSample)
    else:
        bkgSample = ['DYJetsToLL_M-50','DYJetsToLL_M-4to50','WJetsToLNu','YMuMu', 'QCD', 'SingleMuon']
        print(bkgSample)

    for sample in bkgSample:
        print("./runDoHAdd.sh "+sample+" "+args.version+" "+args.filename)
        os.system("./runDoHAdd.sh "+sample+" "+args.version+" "+args.filename)

    os.system("./runDoHAddFlat.sh "+args.version+" "+args.filename)
    #os.system("python3 doHadd.py "+dataset)
