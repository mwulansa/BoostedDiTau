import ROOT,sys,os

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

prefix="root://cmseos.fnal.gov//eos/uscms/store/user/mwulansa/DIS/TCPAnalysis/Backgrounds/RunIIUL17/"

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    inputFileName = prefix+inputFileName.replace("\n","")

    f = ROOT.TFile.Open(inputFileName)

    print inputFileName

