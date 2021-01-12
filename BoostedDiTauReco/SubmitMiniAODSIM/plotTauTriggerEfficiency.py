

import ROOT,sys,os
from DataFormats.FWLite import Events, Handle
from looseElectron import *
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

if len(sys.argv)>2:
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./plots/"


outputFileName = outputFileDir+"h_TauTriggerEfficiency_ETau.root"
print outputFileName
    

out=ROOT.TFile.Open(outputFileName,'recreate')

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleBoostedTaus = Handle ('vector<pat::Tau>')
labelBoostedTaus = ('slimmedTausBoosted')

handleElectronCleanedTaus = Handle ('vector<pat::Tau>')
labelElectronCleanedTaus = ('slimmedTausElectronCleaned')

handleMuonCleanedTaus = Handle ('vector<pat::Tau>')
labelMuonCleanedTaus = ('slimmedTausMuonCleaned')

handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlineSlimmedPrimaryVertices')

handleBs = Handle ('reco::BeamSpot')
labelBs = ("offlineBeamSpot")

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults','','HLT')

handleConv = Handle ('vector<reco::Conversion>')
labelConv = ('reducedEgamma', 'reducedConversions')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('slimmedJets')

handleGenJet = Handle ('vector<reco::GenJet>')
labelGenJet = ('slimmedGenJets')

#<Histograms>

h = {}

h['hNEvent'] = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)




h['hETaudR_M'] = ROOT.TH1F ("hETau_dR_M", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETaudR_ePt'] = ROOT.TH1F ("hETau_dR_ePt", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETaudR_tauPt'] = ROOT.TH1F ("hETau_dR_tauPt", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETaudR_jPt'] = ROOT.TH1F ("hETau_dR_jPt", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETaudR_dR'] = ROOT.TH1F ("hETau_dR_dR", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETaudR_dRlj'] = ROOT.TH1F ("hETau_dR_dRlj", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2)

h['hETauTrig_M_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_M_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_dR_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2)

h['hETauTrig_M_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_M_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_dR_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1)

h['hETauTrig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1)

h['hETauTrig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1)

h['hETauTrig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1)

h['hETauTrig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5)

h['hETauTrig_M_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hETau_Trig_M_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hETau_Trig_dR_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7)

h['hETauTrig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hETau_Trig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hETau_Trig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7)

h['hETauTrig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1)

h['hETauTrig_M_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hETau_Trig_M_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hETau_Trig_dR_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

h['hETauTrig_M_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_M_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dR_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

h['hETauTrig_M_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_M_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dR_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

h['hETauTrig_M_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_M_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dR_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

h['hETauTrig_M_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_M_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dR_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5)

h['hETauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hETau_Trig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hETau_Trig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5)

h['hETauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hETau_Trig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hETau_Trig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5)

h['hETauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hETau_Trig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hETau_Trig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5)

h['hETauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hETau_Trig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hETau_Trig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_v7)

h['hETauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_PFTau120_eta2p1_v5)

h['hETauTrig_M_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_M_HLT_PFTau120_eta2p1_v5", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_PFTau120_eta2p1_v5", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_PFTau120_eta2p1_v5", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_PFTau120_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_dR_HLT_PFTau120_eta2p1_v5", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_PFTau120_eta2p1_v5", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_PFTau140_eta2p1_v5)

h['hETauTrig_M_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_M_HLT_PFTau140_eta2p1_v5", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_PFTau140_eta2p1_v5", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_PFTau140_eta2p1_v5", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_PFTau140_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_dR_HLT_PFTau140_eta2p1_v5", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_PFTau140_eta2p1_v5", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5)

h['hETauTrig_M_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_M_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_dR_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

        #-----Passing Trigger (HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5)

h['hETauTrig_M_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_M_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_ePt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_tauPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_jPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_dR_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hETau_Trig_dRlj_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)
OB
    

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    print inputFileName.replace("
","")
    inputFileName=inputFileName.replace("
","")
    f=ROOT.TFile.Open(inputFileName)

    if not f.IsZombie():
        events=Events(inputFileName)
    else:
        print "Can't Open File: "+inputFileName
        continue

    for event in events:

        event.getByLabel(labelVertex, handleVertex)
        vertex=handleVertex.product()

        event.getByLabel(labelMuon, handleMuon)
        muons=handleMuon.product()

        event.getByLabel(labelElectron, handleElectron)
        electrons=handleElectron.product()

        event.getByLabel(labelJet, handleJet)
        jets=handleJet.product()

        event.getByLabel(labelHLT, handleHLT)
        triggerResults=handleHLT.product()
        names = event.object().triggerNames(triggerResults)

        event.getByLabel(labelBs, handleBs)
        bs=handleBs.product()

        event.getByLabel(labelConv, handleConv)
        convs=handleConv.product()

        event.getByLabel(labelRho, handleRho)
        rho=handleRho.product()[0]

        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()

        event.getByLabel(labelElectronCleanedTaus, handleElectronCleanedTaus)
        etaus=handleElectronCleanedTaus.product()

        event.getByLabel(labelMuonCleanedTaus, handleMuonCleanedTaus)
        mutaus=handleMuonCleanedTaus.product()

        event.getByLabel(labelBoostedTaus, handleBoostedTaus)
        btaus=handleBoostedTaus.product()
        
        h['hNEvent'].Fill(0.5, 1)
        h['hNEvent'].Fill(1.5, genweight)
        
#<muonSelection>

        selected_muons=[]

        for muon in muons:
            if not muon.isLooseMuon(): continue
            if muon.pt()>3 or muon.eta()<2.4:
                selected_muons+=[muon]

        selected_muons.sort(key=lambda x: x.pt(), reverse=True)

#</muonSelection>

#<electronSelection>

        selected_electrons=[]
        for electron in electrons:
            if electron.pt()<7: continue
            if abs(electron.eta())>2.5: continue
            if electron.isEB():
                if electron.full5x5_sigmaIetaIeta()<0.011 and electron.hadronicOverEm()<0.298 and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.222 and GsfEleEInverseMinusPInverse(electron)<0.241 and abs(dEtaInSeed(electron))<0.00477 and GsfEleMissingHitsCut(electron)<=1 and electron.passConversionVeto():
                    selected_electrons+=[electron]
            if electron.isEE():
                if electron.full5x5_sigmaIetaIeta()<0.0314 and electron.hadronicOverEm()<0.101 and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.213 and GsfEleEInverseMinusPInverse(electron)<0.14 and abs(dEtaInSeed(electron))<0.00868 and GsfEleMissingHitsCut(electron)<=1 and electron.passConversionVeto():
                    selected_electrons+=[electron]

        selected_electrons.sort(key=lambda x: x.pt(), reverse=True)

#</electronSelection>

#<tauSelection>


        selected_etaus=[]

        for tau in etaus:
            if tau.pt()<20 or abs(tau.eta())>2.3: continue
            if not tau.tauID("decayModeFinding"): continue
            if not tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"): continue
            selected_etaus+=[tau]

        selected_etaus.sort(key=lambda x: x.pt(), reverse=True)



#</tauSelection>

#<jetSelection>

        selected_jets=[]
        for jet in jets:
            if jet.pt()<20 or abs(jet.eta())>2.5: continue
            NHF  = jet.neutralHadronEnergyFraction()
            NEMF = jet.neutralEmEnergyFraction()
            CHF  = jet.chargedHadronEnergyFraction()
            MUF  = jet.muonEnergyFraction()
            CEMF = jet.chargedEmEnergyFraction()
            NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity()
            NumNeutralParticles =jet.neutralMultiplicity()
            CHM      = jet.chargedMultiplicity()
            if (NHF<0.99 and NEMF<0.99 and NumConst>1) and ((abs(jet.eta())<=2.4 and CHF>0 and CHM>0 and CEMF<0.99) or abs(jet.eta())>2.4) and abs(jet.eta())<=2.7:
                selected_jets+=[jet]
        selected_jets.sort(key=lambda x: x.pt(), reverse=True)

#</jetSelection>



#----------ETau----------

        if len(selected_etaus)>=1 and len(selected_electrons)>=1 and selected_etaus[0].charge()*selected_electrons[0].charge()<0 and len(selected_jets)>=1:

            etau = ROOT.TLorentzVector()
            etau.SetPtEtaPhiM(selected_etaus[0].pt(), selected_etaus[0].eta(), selected_etaus[0].phi(), selected_etaus[0].mass())

            e = ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())

            if len(selected_jets)==1:
                j0=ROOT.TLorentzVector()
                j0.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j=j0
            else:
                j1=ROOT.TLorentzVector()
                j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j2=ROOT.TLorentzVector()
                j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())
                if etau.DeltaR(j1)>0.3:
                    j=j1
                else:
                    j=j2

            if etau.DeltaR(e)<0.4 and etau.DeltaR(j)>0.8 and e.DeltaR(j)>0.8:

                h['hETaudR_M'].Fill((etau+e).M(), genweight)
                h['hETaudR_ePt'].Fill(e.Pt(), genweight) 
                h['hETaudR_tauPt'].Fill(etau.Pt(), genweight)
                h['hETaudR_jPt'].Fill(j.Pt(), genweight)
                h['hETaudR_dR'].Fill(e.DeltaR(etau), genweight)
                h['hETaudR_dRlj'].Fill((e+etau).DeltaR(j), genweight)

    

    #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2)

                if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2')):

                    h['hETauTrig_M_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2)

                if triggerResults.accept(names.triggerIndex('HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2')):

                    h['hETauTrig_M_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1)

                if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1')):

                    h['hETauTrig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1)

                if triggerResults.accept(names.triggerIndex('HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1')):

                    h['hETauTrig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1)

                if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1')):

                    h['hETauTrig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1)

                if triggerResults.accept(names.triggerIndex('HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1')):

                    h['hETauTrig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5)

                if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5')):

                    h['hETauTrig_M_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7)

                if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7')):

                    h['hETauTrig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7)

                if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7')):

                    h['hETauTrig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1)

                if triggerResults.accept(names.triggerIndex('HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1')):

                    h['hETauTrig_M_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

                if triggerResults.accept(names.triggerIndex('HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1')):

                    h['hETauTrig_M_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

                if triggerResults.accept(names.triggerIndex('HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1')):

                    h['hETauTrig_M_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

                if triggerResults.accept(names.triggerIndex('HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1')):

                    h['hETauTrig_M_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

                if triggerResults.accept(names.triggerIndex('HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1')):

                    h['hETauTrig_M_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5)

                if triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5')):

                    h['hETauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5)

                if triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5')):

                    h['hETauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5)

                if triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5')):

                    h['hETauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5)

                if triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5')):

                    h['hETauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_v7)

                if triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_v7')):

                    h['hETauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_PFTau120_eta2p1_v5)

                if triggerResults.accept(names.triggerIndex('HLT_PFTau120_eta2p1_v5')):

                    h['hETauTrig_M_HLT_PFTau120_eta2p1_v5'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_PFTau120_eta2p1_v5'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_PFTau120_eta2p1_v5'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_PFTau120_eta2p1_v5'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_PFTau120_eta2p1_v5'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_PFTau120_eta2p1_v5'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_PFTau140_eta2p1_v5)

                if triggerResults.accept(names.triggerIndex('HLT_PFTau140_eta2p1_v5')):

                    h['hETauTrig_M_HLT_PFTau140_eta2p1_v5'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_PFTau140_eta2p1_v5'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_PFTau140_eta2p1_v5'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_PFTau140_eta2p1_v5'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_PFTau140_eta2p1_v5'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_PFTau140_eta2p1_v5'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5)

                if triggerResults.accept(names.triggerIndex('HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5')):

                    h['hETauTrig_M_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill((e+etau).DeltaR(j), genweight)


    #-----Passing Trigger (HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5)

                if triggerResults.accept(names.triggerIndex('HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5')):

                    h['hETauTrig_M_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill((e+etau).DeltaR(j), genweight)


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

