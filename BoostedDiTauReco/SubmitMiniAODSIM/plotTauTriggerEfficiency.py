

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


outputFileName = outputFileDir+"h_TauTriggerEfficiency_BTau_m10.root"
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

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handlePatMETs = Handle("vector<pat::MET>")
labelPatMETs = ('slimmedMETs')

#<Histograms>

h = {}

h['hNEvent'] = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)




h['hBTaudR_M'] = ROOT.TH1F ("hBTau_dR_M", "#tau - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hBTaudR_tau1Pt'] = ROOT.TH1F ("hBTau_dR_tau1Pt", "tau1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTaudR_tau2Pt'] = ROOT.TH1F ("hBTau_dR_tau2Pt", "tau2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTaudR_jPt'] = ROOT.TH1F ("hBTau_dR_jPt", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTaudR_dR'] = ROOT.TH1F ("hBTau_dR_dR", "#tau1 #tau2 #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTaudR_dRlj'] = ROOT.TH1F ("hBTau_dR_dRlj", "#tau#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauGen_M'] = ROOT.TH1F ("hBTau_Gen_M", "#tau - #tau mass;M_{#mu#tau};N_{events}", 1000, 0, 200)
h['hBTauGen_tau1Pt'] = ROOT.TH1F ("hBTau_Gen_tau1Pt", "tau1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauGen_tau2Pt'] = ROOT.TH1F ("hBTau_Gen_tau2Pt", "tau2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauGen_jPt'] = ROOT.TH1F ("hBTau_Gen_jPt", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauGen_dR'] = ROOT.TH1F ("hBTau_Gen_dR", "#tau1 #tau2 #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauGen_dRlj'] = ROOT.TH1F ("hBTau_Gen_dRlj", "#tau#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigJet_M'] = ROOT.TH1F ("hBTau_TrigJet_M", "#tau - #tau mass;M_{#mu#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigJet_tau1Pt'] = ROOT.TH1F ("hBTau_TrigJet_tau1Pt", "tau1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigJet_tau2Pt'] = ROOT.TH1F ("hBTau_TrigJet_tau2Pt", "tau2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigJet_jPt'] = ROOT.TH1F ("hBTau_TrigJet_jPt", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigJet_dR'] = ROOT.TH1F ("hBTau_TrigJet_dR", "#tau1 #tau2 #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigJet_dRlj'] = ROOT.TH1F ("hBTau_TrigJet_dRlj", "#tau#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2)

h['hBTauTrig_M_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_M_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2)

h['hBTauTrig_M_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_M_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1)

h['hBTauTrig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1)

h['hBTauTrig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1)

h['hBTauTrig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1)

h['hBTauTrig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5)

h['hBTauTrig_M_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_Trig_M_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7)

h['hBTauTrig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_Trig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7)

h['hBTauTrig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1)

h['hBTauTrig_M_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_Trig_M_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

h['hBTauTrig_M_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_M_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

h['hBTauTrig_M_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_M_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

h['hBTauTrig_M_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_M_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

h['hBTauTrig_M_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_M_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5)

h['hBTauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_Trig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5)

h['hBTauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_Trig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5)

h['hBTauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_Trig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5)

h['hBTauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_Trig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_v7)

h['hBTauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_PFTau120_eta2p1_v5)

h['hBTauTrig_M_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_M_HLT_PFTau120_eta2p1_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_PFTau120_eta2p1_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_PFTau120_eta2p1_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_PFTau120_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_PFTau120_eta2p1_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_PFTau120_eta2p1_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_PFTau120_eta2p1_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_PFTau120_eta2p1_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_PFTau120_eta2p1_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_PFTau120_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_PFTau120_eta2p1_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_PFTau120_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_PFTau120_eta2p1_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_PFTau140_eta2p1_v5)

h['hBTauTrig_M_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_M_HLT_PFTau140_eta2p1_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_PFTau140_eta2p1_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_PFTau140_eta2p1_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_PFTau140_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_PFTau140_eta2p1_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_PFTau140_eta2p1_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_PFTau140_eta2p1_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_PFTau140_eta2p1_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_PFTau140_eta2p1_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_PFTau140_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_PFTau140_eta2p1_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_PFTau140_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_PFTau140_eta2p1_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5)

h['hBTauTrig_M_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_M_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

        #-----Passing Trigger (HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5)

h['hBTauTrig_M_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_M_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_jPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_dR_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_Trig_dRlj_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_M_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    print inputFileName.replace("\n","")
    inputFileName=inputFileName.replace("\n","")
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

        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()

        h['hNEvent'].Fill(0.5, 1)
        h['hNEvent'].Fill(1.5, genweight)

#<genInfo>

        genMuons=[]
        genTaus=[]
        genElectrons=[]
        genNutaus=[]

        for particle in particles:
            if abs(particle.pdgId())==15 and particle.mother().pdgId()==9999:
                genTaus+=[particle]
            if abs(particle.pdgId())==11 and particle.isDirectHardProcessTauDecayProductFinalState():
                genElectrons+=[particle]
            if abs(particle.pdgId())==13 and particle.isDirectHardProcessTauDecayProductFinalState():
                genMuons+=[particle]
            if abs(particle.pdgId())==16 and particle.isDirectHardProcessTauDecayProductFinalState():
                genNutaus+=[particle]

#</genInfo>
        
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


        selected_btaus=[]

        for tau in btaus:
            if tau.pt()<20 or abs(tau.eta())>2.3: continue
            if not tau.tauID("decayModeFinding"): continue
            if not tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"): continue
            selected_btaus+=[tau]

        selected_btaus.sort(key=lambda x: x.pt(), reverse=True)



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



#----------BTau----------

        if len(selected_btaus)>=2 and selected_btaus[0].charge()*selected_btaus[1].charge()<0 and len(selected_jets)>=1:

            tau2 = ROOT.TLorentzVector()
            tau2.SetPtEtaPhiM(selected_btaus[1].pt(), selected_btaus[1].eta(), selected_btaus[1].phi(), selected_btaus[1].mass())

            tau1 = ROOT.TLorentzVector()
            tau1.SetPtEtaPhiM(selected_btaus[0].pt(), selected_btaus[0].eta(), selected_btaus[0].phi(), selected_btaus[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
            
            if tau2.DeltaR(tau1)<0.4 and tau2.DeltaR(j)>0.8 and tau1.DeltaR(j)>0.8:

                h['hBTaudR_M'].Fill((tau1+tau2).M(), genweight)
                h['hBTaudR_tau2Pt'].Fill(tau2.Pt(), genweight) 
                h['hBTaudR_tau1Pt'].Fill(tau1.Pt(), genweight)
                h['hBTaudR_jPt'].Fill(j.Pt(), genweight)
                h['hBTaudR_dR'].Fill(tau2.DeltaR(tau1), genweight)
                h['hBTaudR_dRlj'].Fill((tau2+tau1).DeltaR(j), genweight)

                if len(genTaus) == 2 and len(genMuons) == 0 and len(genElectrons) == 0:
                    btau1 = genTaus[0]
                    btau2 = genTaus[1]

                    genTau1 = ROOT.TLorentzVector()
                    genTau2 = ROOT.TLorentzVector()
                    genNt1 = ROOT.TLorentzVector()
                    genNt2 = ROOT.TLorentzVector()

                    genTau1.SetPtEtaPhiM(btau1.pt(), btau1.eta(), btau1.phi(), btau1.mass()) if btau1.pdgId()==15 else genTau1.SetPtEtaPhiM(btau2.pt(), btau2.eta(), btau2.phi(), btau2.mass())
                    genNt1.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId==16 else genNt1.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                    genTau2.SetPtEtaPhiM(btau1.pt(), btau1.eta(), btau1.phi(), btau1.mass()) if btau1.pdgId()==-15 else genTau2.SetPtEtaPhiM(btau2.pt(), btau2.eta(), btau2.phi(), btau2.mass())
                    genNt2.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId==-16 else genNt2.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
 
                    if (tau1.DeltaR(genTau1-genNt1)<0.4 or tau1.DeltaR(genTau2-genNt2)<0.4) and (tau2.DeltaR(genTau1-genNt1)<0.4 or tau2.DeltaR(genTau2-genNt2)<0.4):                       
                        h['hBTauGen_M'].Fill((tau1+tau2).M(), genweight)
                        h['hBTauGen_tau2Pt'].Fill(tau2.Pt(), genweight) 
                        h['hBTauGen_tau1Pt'].Fill(tau1.Pt(), genweight)
                        h['hBTauGen_jPt'].Fill(j.Pt(), genweight)
                        h['hBTauGen_dR'].Fill(tau2.DeltaR(tau1), genweight)
                        h['hBTauGen_dRlj'].Fill((tau2+tau1).DeltaR(j), genweight)

                        if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))):
                            h['hBTauTrigJet_M'].Fill((tau1+tau2).M(), genweight)
                            h['hBTauTrigJet_tau2Pt'].Fill(tau2.Pt(), genweight) 
                            h['hBTauTrigJet_tau1Pt'].Fill(tau1.Pt(), genweight)
                            h['hBTauTrigJet_jPt'].Fill(j.Pt(), genweight)
                            h['hBTauTrigJet_dR'].Fill(tau2.DeltaR(tau1), genweight)
                            h['hBTauTrigJet_dRlj'].Fill((tau2+tau1).DeltaR(j), genweight)                    

    

    #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2)

                            if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2')):
                                h['hBTauTrig_M_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'))):
                                h['hBTauTrigTauJet_M_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2)

                            if triggerResults.accept(names.triggerIndex('HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2')):
                                h['hBTauTrig_M_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'))):
                                h['hBTauTrigTauJet_M_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_v2'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1)

                            if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1')):
                                h['hBTauTrig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'))):
                                h['hBTauTrigTauJet_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1)

                            if triggerResults.accept(names.triggerIndex('HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1')):
                                h['hBTauTrig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'))):
                                h['hBTauTrigTauJet_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1)

                            if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1')):
                                h['hBTauTrig_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'))):
                                h['hBTauTrigTauJet_M_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1)

                            if triggerResults.accept(names.triggerIndex('HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1')):
                                h['hBTauTrig_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'))):
                                h['hBTauTrigTauJet_M_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_v1'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5)

                            if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5')):
                                h['hBTauTrig_M_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'))):
                                h['hBTauTrigTauJet_M_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v5'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7)

                            if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7')):
                                h['hBTauTrig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'))):
                                h['hBTauTrigTauJet_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v7'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7)

                            if triggerResults.accept(names.triggerIndex('HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7')):
                                h['hBTauTrig_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'))):
                                h['hBTauTrigTauJet_M_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_v7'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1)

                            if triggerResults.accept(names.triggerIndex('HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1')):
                                h['hBTauTrig_M_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'))):
                                h['hBTauTrigTauJet_M_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_v1'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

                            if triggerResults.accept(names.triggerIndex('HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1')):
                                h['hBTauTrig_M_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'))):
                                h['hBTauTrigTauJet_M_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

                            if triggerResults.accept(names.triggerIndex('HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1')):
                                h['hBTauTrig_M_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'))):
                                h['hBTauTrigTauJet_M_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

                            if triggerResults.accept(names.triggerIndex('HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1')):
                                h['hBTauTrig_M_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'))):
                                h['hBTauTrigTauJet_M_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1)

                            if triggerResults.accept(names.triggerIndex('HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1')):
                                h['hBTauTrig_M_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'))):
                                h['hBTauTrigTauJet_M_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_v1'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5)

                            if triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5')):
                                h['hBTauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'))):
                                h['hBTauTrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v5'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5)

                            if triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5')):
                                h['hBTauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'))):
                                h['hBTauTrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5)

                            if triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5')):
                                h['hBTauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'))):
                                h['hBTauTrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v5'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5)

                            if triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5')):
                                h['hBTauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'))):
                                h['hBTauTrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v5'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_LooseIsoPFTau50_Trk30_eta2p1_v7)

                            if triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_v7')):
                                h['hBTauTrig_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'))):
                                h['hBTauTrigTauJet_M_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_LooseIsoPFTau50_Trk30_eta2p1_v7'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_PFTau120_eta2p1_v5)

                            if triggerResults.accept(names.triggerIndex('HLT_PFTau120_eta2p1_v5')):
                                h['hBTauTrig_M_HLT_PFTau120_eta2p1_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_PFTau120_eta2p1_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_PFTau120_eta2p1_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_PFTau120_eta2p1_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_PFTau120_eta2p1_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_PFTau120_eta2p1_v5'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_PFTau120_eta2p1_v5'))):
                                h['hBTauTrigTauJet_M_HLT_PFTau120_eta2p1_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_PFTau120_eta2p1_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_PFTau120_eta2p1_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_PFTau120_eta2p1_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_PFTau120_eta2p1_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_PFTau120_eta2p1_v5'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_PFTau140_eta2p1_v5)

                            if triggerResults.accept(names.triggerIndex('HLT_PFTau140_eta2p1_v5')):
                                h['hBTauTrig_M_HLT_PFTau140_eta2p1_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_PFTau140_eta2p1_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_PFTau140_eta2p1_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_PFTau140_eta2p1_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_PFTau140_eta2p1_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_PFTau140_eta2p1_v5'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_PFTau140_eta2p1_v5'))):
                                h['hBTauTrigTauJet_M_HLT_PFTau140_eta2p1_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_PFTau140_eta2p1_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_PFTau140_eta2p1_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_PFTau140_eta2p1_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_PFTau140_eta2p1_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_PFTau140_eta2p1_v5'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5)

                            if triggerResults.accept(names.triggerIndex('HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5')):
                                h['hBTauTrig_M_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'))):
                                h['hBTauTrigTauJet_M_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_VLooseIsoPFTau120_Trk50_eta2p1_v5'].Fill((tau2+tau1).DeltaR(j), genweight)



    #-----Passing Trigger (HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5)

                            if triggerResults.accept(names.triggerIndex('HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5')):
                                h['hBTauTrig_M_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrig_tau2Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrig_tau1Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrig_jPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrig_dR_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrig_dRlj_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill((tau2+tau1).DeltaR(j), genweight)

                            if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'))):
                                h['hBTauTrigTauJet_M_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill((tau1+tau2).M(), genweight)
                                h['hBTauTrigTauJet_tau2Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(tau2.Pt(), genweight)
                                h['hBTauTrigTauJet_tau1Pt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(tau1.Pt(), genweight)
                                h['hBTauTrigTauJet_jPt_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(j.Pt(), genweight)
                                h['hBTauTrigTauJet_dR_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill(tau2.DeltaR(tau1), genweight)
                                h['hBTauTrigTauJet_dRlj_HLT_VLooseIsoPFTau140_Trk50_eta2p1_v5'].Fill((tau2+tau1).DeltaR(j), genweight)



out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

