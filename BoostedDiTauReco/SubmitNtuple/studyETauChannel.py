import ROOT, sys, os
import numpy as np
import time

start_time = time.time()

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

outputTitle = "h_studyETauChannel"

if "-b" in opts:
    isData = 0

ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/JetInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/MuonInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/ElectronInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/TauInfoDS.h"')

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

if len(sys.argv)>2 and not sys.argv[2].startswith("-"):
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./output/"

out=ROOT.TFile.Open(outputFileName,'recreate')
print(outputFileName)

fchain = ROOT.TChain('tcpNtuples/analysisTree')
chain2 = ROOT.TChain('tcpTrigNtuples/triggerTree')
if isData == 0:
    chain3 = ROOT.TChain('lumiSummary/lumiTree')
    chain4 = ROOT.TChain('tcpGenNtuples/genTree')

h = {}

event_cut = {
    'jetPt': 70,
    'dRl': 0.4,
    'dRltau': 0.05,
    'dRlj': 0.8,
    'metcut': 100.0,
    'mtcut': 50.0,
    'dPhiml': 1,
    'dPhimj': 2,
    'mass' : 5
}

