import os, sys

ch = "BTau"

script = open ("plotTauTriggerEfficiency.py", "w")
script.writelines("""

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

""")

script.close()

script = open ("plotTauTriggerEfficiency.py", "a")

if ch == "ETau":
    script.writelines("""
outputFileName = outputFileDir+"h_TauTriggerEfficiency_ETau.root"
print outputFileName
    """)

if ch == "BTau":
    script.writelines("""
outputFileName = outputFileDir+"h_TauTriggerEfficiency_BTau.root"
print outputFileName
    """)

if ch == "MuTau":
    script.writelines("""
outputFileName = outputFileDir+"h_TauTriggerEfficiency_MuTau.root"
print outputFileName
    """)

script.close()
script = open ("plotTauTriggerEfficiency.py", "a")
script.writelines("""

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

""")

script.close()
script = open ("plotTauTriggerEfficiency.py", "a")

if ch == "ETau":
    script.writelines("""


h['hETaudR_M'] = ROOT.TH1F ("hETau_dR_M", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETaudR_ePt'] = ROOT.TH1F ("hETau_dR_ePt", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETaudR_tauPt'] = ROOT.TH1F ("hETau_dR_tauPt", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETaudR_jPt'] = ROOT.TH1F ("hETau_dR_jPt", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETaudR_dR'] = ROOT.TH1F ("hETau_dR_dR", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETaudR_dRlj'] = ROOT.TH1F ("hETau_dR_dRlj", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    """)

if ch == "MuTau":
    script.writelines("""


h['hMuTaudR_M'] = ROOT.TH1F ("hMuTau_dR_M", "#mu - #tau mass;M_{#mu#tau};N_{events}", 1000, 0, 200)
h['hMuTaudR_muPt'] = ROOT.TH1F ("hMuTau_dR_muPt", "muon P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hMuTaudR_tauPt'] = ROOT.TH1F ("hMuTau_dR_tauPt", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hMuTaudR_jPt'] = ROOT.TH1F ("hMuTau_dR_jPt", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hMuTaudR_dR'] = ROOT.TH1F ("hMuTau_dR_dR", "#mu #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hMuTaudR_dRlj'] = ROOT.TH1F ("hMuTau_dR_dRlj", "#e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    """)

if ch == "BTau":
    script.writelines("""


h['hBTaudR_M'] = ROOT.TH1F ("hBTau_dR_M", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hBTaudR_tau1Pt'] = ROOT.TH1F ("hBTau_dR_tau1Pt", "tau1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTaudR_tau2Pt'] = ROOT.TH1F ("hBTau_dR_tau2Pt", "tau2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTaudR_jPt'] = ROOT.TH1F ("hBTau_dR_jPt", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTaudR_dR'] = ROOT.TH1F ("hBTau_dR_dR", "#tau1 #tau2 #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTaudR_dRlj'] = ROOT.TH1F ("hBTau_dR_dRlj", "#tau#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    """)


script.close()

if ch == "ETau":

    triggerNames  = open("triggerTauList.txt", "r")
    for triggerName in triggerNames:
        script = open ("plotTauTriggerEfficiency.py", "a")
        script.writelines("""

        #-----Passing Trigger ("""+triggerName.replace("\n","")+""")

h['hETauTrig_M_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_Trig_M_"""+triggerName.replace("\n","")+"""", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrig_ePt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_Trig_ePt_"""+triggerName.replace("\n","")+"""", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_tauPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_Trig_tauPt_"""+triggerName.replace("\n","")+"""", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrig_jPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_Trig_jPt_"""+triggerName.replace("\n","")+"""", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrig_dR_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_Trig_dR_"""+triggerName.replace("\n","")+"""", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrig_dRlj_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_Trig_dRlj_"""+triggerName.replace("\n","")+"""", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hETauTrigJet_M_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigJet_M_"""+triggerName.replace("\n","")+"""", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrigJet_ePt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigJet_ePt_"""+triggerName.replace("\n","")+"""", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrigJet_tauPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigJet_tauPt_"""+triggerName.replace("\n","")+"""", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrigJet_jPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigJet_jPt_"""+triggerName.replace("\n","")+"""", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrigJet_dR_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigJet_dR_"""+triggerName.replace("\n","")+"""", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrigJet_dRlj_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigJet_dRlj_"""+triggerName.replace("\n","")+"""", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hETauTrigTauJet_M_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigTauJet_M_"""+triggerName.replace("\n","")+"""", "e - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hETauTrigTauJet_ePt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigTauJet_ePt_"""+triggerName.replace("\n","")+"""", "electron P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrigTauJet_tauPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigTauJet_tauPt_"""+triggerName.replace("\n","")+"""", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hETauTrigTauJet_jPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigTauJet_jPt_"""+triggerName.replace("\n","")+"""", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hETauTrigTauJet_dR_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigTauJet_dR_"""+triggerName.replace("\n","")+"""", "e #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hETauTrigTauJet_dRlj_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hETau_TrigTauJet_dRlj_"""+triggerName.replace("\n","")+"""", "e#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    """)

    script.close()

if ch == "BTau":

    triggerNames  = open("triggerTauList.txt", "r")
    for triggerName in triggerNames:
        script = open ("plotTauTriggerEfficiency.py", "a")
        script.writelines("""

        #-----Passing Trigger ("""+triggerName.replace("\n","")+""")

h['hBTauTrig_M_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_Trig_M_"""+triggerName.replace("\n","")+"""", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrig_tau1Pt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_Trig_tau1Pt_"""+triggerName.replace("\n","")+"""", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_tau2Pt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_Trig_tau2Pt_"""+triggerName.replace("\n","")+"""", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrig_jPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_Trig_jPt_"""+triggerName.replace("\n","")+"""", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrig_dR_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_Trig_dR_"""+triggerName.replace("\n","")+"""", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrig_dRlj_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_Trig_dRlj_"""+triggerName.replace("\n","")+"""", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigJet_M_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigJet_M_"""+triggerName.replace("\n","")+"""", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigJet_tau1Pt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigJet_tau1Pt_"""+triggerName.replace("\n","")+"""", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigJet_tau2Pt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigJet_tau2Pt_"""+triggerName.replace("\n","")+"""", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigJet_jPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigJet_jPt_"""+triggerName.replace("\n","")+"""", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigJet_dR_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigJet_dR_"""+triggerName.replace("\n","")+"""", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigJet_dRlj_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigJet_dRlj_"""+triggerName.replace("\n","")+"""", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hBTauTrigTauJet_M_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigTauJet_M_"""+triggerName.replace("\n","")+"""", "#tau mass;M_{#tau};N_{events}", 1000, 0, 200)
h['hBTauTrigTauJet_tau1Pt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigTauJet_tau1Pt_"""+triggerName.replace("\n","")+"""", "tau 1 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_tau2Pt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigTauJet_tau2Pt_"""+triggerName.replace("\n","")+"""", "tau 2 P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hBTauTrigTauJet_jPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigTauJet_jPt_"""+triggerName.replace("\n","")+"""", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hBTauTrigTauJet_dR_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigTauJet_dR_"""+triggerName.replace("\n","")+"""", "#tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hBTauTrigTauJet_dRlj_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hBTau_TrigTauJet_dRlj_"""+triggerName.replace("\n","")+"""", "#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    """)

    script.close()

if ch == "MuTau":

    triggerNames  = open("triggerTauList.txt", "r")
    for triggerName in triggerNames:
        script = open ("plotTauTriggerEfficiency.py", "a")
        script.writelines("""

        #-----Passing Trigger ("""+triggerName.replace("\n","")+""")

h['hMuTauTrig_M_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_Trig_M_"""+triggerName.replace("\n","")+"""", "#mu - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hMuTauTrig_muPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_Trig_muPt_"""+triggerName.replace("\n","")+"""", "#mu P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hMuTauTrig_tauPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_Trig_tauPt_"""+triggerName.replace("\n","")+"""", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hMuTauTrig_jPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_Trig_jPt_"""+triggerName.replace("\n","")+"""", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hMuTauTrig_dR_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_Trig_dR_"""+triggerName.replace("\n","")+"""", "#mu #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hMuTauTrig_dRlj_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_Trig_dRlj_"""+triggerName.replace("\n","")+"""", "#mu#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hMuTauTrigJet_M_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigJet_M_"""+triggerName.replace("\n","")+"""", "#mu - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hMuTauTrigJet_muPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigJet_muPt_"""+triggerName.replace("\n","")+"""", "#mu P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hMuTauTrigJet_tauPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigJet_tauPt_"""+triggerName.replace("\n","")+"""", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hMuTauTrigJet_jPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigJet_jPt_"""+triggerName.replace("\n","")+"""", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hMuTauTrigJet_dR_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigJet_dR_"""+triggerName.replace("\n","")+"""", "#mu #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hMuTauTrigJet_dRlj_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigJet_dRlj_"""+triggerName.replace("\n","")+"""", "#mu#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

h['hMuTauTrigTauJet_M_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigTauJet_M_"""+triggerName.replace("\n","")+"""", "#mu - #tau mass;M_{e#tau};N_{events}", 1000, 0, 200)
h['hMuTauTrigTauJet_muPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigTauJet_muPt_"""+triggerName.replace("\n","")+"""", "#mu P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hMuTauTrigTauJet_tauPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigTauJet_tauPt_"""+triggerName.replace("\n","")+"""", "tau P_{t} ; P_{t} ; N_{events}", 500, 0, 500)
h['hMuTauTrigTauJet_jPt_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigTauJet_jPt_"""+triggerName.replace("\n","")+"""", "jet P_{t}; P_{t};N_{events}", 2000, 0, 2000)
h['hMuTauTrigTauJet_dR_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigTauJet_dR_"""+triggerName.replace("\n","")+"""", "#mu #tau #Delta R;#Delta R;N_{events}", 100, 0, 5)
h['hMuTauTrigTauJet_dRlj_"""+triggerName.replace("\n","")+"""'] = ROOT.TH1F ("hMuTau_TrigTauJet_dRlj_"""+triggerName.replace("\n","")+"""", "#mu#tau and jet #Delta R;#Delta R;N_{events}", 100, 0, 5)

    """)

    script.close()

script = open ("plotTauTriggerEfficiency.py", "a")
script.writelines("""

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

        event.getByLabel(labelPatMETs, handlePatMETs)
        met=handlePatMETs.product().front()
        mets=[]
        mets+=[met]
        mets.sort(key=lambda x: x.pt(), reverse=True)

        m=ROOT.TLorentzVector()
        m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

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

""")

script.close()

script = open ("plotTauTriggerEfficiency.py", "a")

if ch == "BTau":
    script.writelines("""
        selected_btaus=[]

        for tau in btaus:
            if tau.pt()<20 or abs(tau.eta())>2.3: continue
            if not tau.tauID("decayModeFinding"): continue
            if not tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"): continue
            selected_btaus+=[tau]

        selected_btaus.sort(key=lambda x: x.pt(), reverse=True)
""")

if ch == "ETau":
    script.writelines("""
        selected_etaus=[]

        for tau in etaus:
            if tau.pt()<20 or abs(tau.eta())>2.3: continue
            if not tau.tauID("decayModeFinding"): continue
            if not tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"): continue
            selected_etaus+=[tau]

        selected_etaus.sort(key=lambda x: x.pt(), reverse=True)
""")

if ch == "MuTau":
    script.writelines("""
        selected_mutaus=[]

        for tau in mutaus:
            if tau.pt()<20 or abs(tau.eta())>2.3: continue
            if not tau.tauID("decayModeFinding"): continue
            if not tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"): continue
            selected_mutaus+=[tau]

        selected_mutaus.sort(key=lambda x: x.pt(), reverse=True)
""")

script.close()

script = open ("plotTauTriggerEfficiency.py", "a")
script.writelines("""


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

""")

script.close()
script = open ("plotTauTriggerEfficiency.py", "a")
if ch == "ETau":
    script.writelines("""

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



    """)

    script.close()

    triggerNames  = open("triggerTauList.txt", "r")
    for triggerName in triggerNames:
        print triggerName
        script = open ("plotTauTriggerEfficiency.py", "a")
        script.writelines("""

    #-----Passing Trigger ("""+triggerName.replace("\n","")+""")

                if triggerResults.accept(names.triggerIndex('"""+triggerName.replace("\n","")+"""')):

                    h['hETauTrig_M_"""+triggerName.replace("\n","")+"""'].Fill((etau+e).M(), genweight)
                    h['hETauTrig_ePt_"""+triggerName.replace("\n","")+"""'].Fill(e.Pt(), genweight)
                    h['hETauTrig_tauPt_"""+triggerName.replace("\n","")+"""'].Fill(etau.Pt(), genweight)
                    h['hETauTrig_jPt_"""+triggerName.replace("\n","")+"""'].Fill(j.Pt(), genweight)
                    h['hETauTrig_dR_"""+triggerName.replace("\n","")+"""'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrig_dRlj_"""+triggerName.replace("\n","")+"""'].Fill((e+etau).DeltaR(j), genweight)

                    if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('"""+triggerName.replace("\n","")+"""'))):
                    h['hETauTrigTauJet_M_"""+triggerName.replace("\n","")+"""'].Fill((etau+e).M(), genweight)
                    h['hETauTrigTauJet_ePt_"""+triggerName.replace("\n","")+"""'].Fill(e.Pt(), genweight)
                    h['hETauTrigTauJet_tauPt_"""+triggerName.replace("\n","")+"""'].Fill(etau.Pt(), genweight)
                    h['hETauTrigTauJet_jPt_"""+triggerName.replace("\n","")+"""'].Fill(j.Pt(), genweight)
                    h['hETauTrigTauJet_dR_"""+triggerName.replace("\n","")+"""'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrigTauJet_dRlj_"""+triggerName.replace("\n","")+"""'].Fill((e+etau).DeltaR(j), genweight)

                if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))):
                    h['hETauTrigJet_M_"""+triggerName.replace("\n","")+"""'].Fill((etau+e).M(), genweight)
                    h['hETauTrigJet_ePt_"""+triggerName.replace("\n","")+"""'].Fill(e.Pt(), genweight)
                    h['hETauTrigJet_tauPt_"""+triggerName.replace("\n","")+"""'].Fill(etau.Pt(), genweight)
                    h['hETauTrigJet_jPt_"""+triggerName.replace("\n","")+"""'].Fill(j.Pt(), genweight)
                    h['hETauTrigJet_dR_"""+triggerName.replace("\n","")+"""'].Fill(e.DeltaR(etau), genweight)
                    h['hETauTrigJet_dRlj_"""+triggerName.replace("\n","")+"""'].Fill((e+etau).DeltaR(j), genweight)
""")

        script.close()

if ch == "BTau":
    script.writelines("""

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

    """)

    script.close()

    triggerNames  = open("triggerTauList.txt", "r")
    for triggerName in triggerNames:
        print triggerName
        script = open ("plotTauTriggerEfficiency.py", "a")
        script.writelines("""

    #-----Passing Trigger ("""+triggerName.replace("\n","")+""")

                if triggerResults.accept(names.triggerIndex('"""+triggerName.replace("\n","")+"""')):

                    h['hBTauTrig_M_"""+triggerName.replace("\n","")+"""'].Fill((tau1+tau2).M(), genweight)
                    h['hBTauTrig_tau2Pt_"""+triggerName.replace("\n","")+"""'].Fill(tau2.Pt(), genweight)
                    h['hBTauTrig_tau1Pt_"""+triggerName.replace("\n","")+"""'].Fill(tau1.Pt(), genweight)
                    h['hBTauTrig_jPt_"""+triggerName.replace("\n","")+"""'].Fill(j.Pt(), genweight)
                    h['hBTauTrig_dR_"""+triggerName.replace("\n","")+"""'].Fill(tau2.DeltaR(tau1), genweight)
                    h['hBTauTrig_dRlj_"""+triggerName.replace("\n","")+"""'].Fill((tau2+tau1).DeltaR(j), genweight)

                if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('"""+triggerName.replace("\n","")+"""'))):
                    h['hBTauTrigTauJet_M_"""+triggerName.replace("\n","")+"""'].Fill((tau1+tau2).M(), genweight)
                    h['hBTauTrigTauJet_tau2Pt_"""+triggerName.replace("\n","")+"""'].Fill(tau2.Pt(), genweight)
                    h['hBTauTrigTauJet_tau1Pt_"""+triggerName.replace("\n","")+"""'].Fill(tau1.Pt(), genweight)
                    h['hBTauTrigTauJet_jPt_"""+triggerName.replace("\n","")+"""'].Fill(j.Pt(), genweight)
                    h['hBTauTrigTauJet_dR_"""+triggerName.replace("\n","")+"""'].Fill(tau2.DeltaR(tau1), genweight)
                    h['hBTauTrigTauJet_dRlj_"""+triggerName.replace("\n","")+"""'].Fill((tau2+tau1).DeltaR(j), genweight)

                if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))):
                    h['hBTauTrigJet_M_"""+triggerName.replace("\n","")+"""'].Fill((tau1+tau2).M(), genweight)
                    h['hBTauTrigJet_tau2Pt_"""+triggerName.replace("\n","")+"""'].Fill(tau2.Pt(), genweight)
                    h['hBTauTrigJet_tau1Pt_"""+triggerName.replace("\n","")+"""'].Fill(tau1.Pt(), genweight)
                    h['hBTauTrigJet_jPt_"""+triggerName.replace("\n","")+"""'].Fill(j.Pt(), genweight)
                    h['hBTauTrigJet_dR_"""+triggerName.replace("\n","")+"""'].Fill(tau2.DeltaR(tau1), genweight)
                    h['hBTauTrigJet_dRlj_"""+triggerName.replace("\n","")+"""'].Fill((tau2+tau1).DeltaR(j), genweight)
""")

        script.close()

if ch == "MuTau":
    script.writelines("""

#----------MuTau----------

        if len(selected_mutaus)>=1 and len(selected_muons)>=1 and selected_mutaus[0].charge()*selected_muons[0].charge()<0 and len(selected_jets)>=1:

            mutau = ROOT.TLorentzVector()
            mutau.SetPtEtaPhiM(selected_mutaus[0].pt(), selected_mutaus[0].eta(), selected_mutaus[0].phi(), selected_mutaus[0].mass())

            mu = ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())

            if len(selected_jets)==1:
                j0=ROOT.TLorentzVector()
                j0.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j=j0
            else:
                j1=ROOT.TLorentzVector()
                j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j2=ROOT.TLorentzVector()
                j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())
                if mutau.DeltaR(j1)>0.3:
                    j=j1
                else:
                    j=j2

            if mutau.DeltaR(mu)<0.4 and mutau.DeltaR(j)>0.8 and mu.DeltaR(j)>0.8:

                h['hMuTaudR_M'].Fill((mutau+mu).M(), genweight)
                h['hMuTaudR_muPt'].Fill(mu.Pt(), genweight) 
                h['hMuTaudR_tauPt'].Fill(mutau.Pt(), genweight)
                h['hMuTaudR_jPt'].Fill(j.Pt(), genweight)
                h['hMuTaudR_dR'].Fill(mu.DeltaR(mutau), genweight)
                h['hMuTaudR_dRlj'].Fill((mu+mutau).DeltaR(j), genweight)

    """)

    script.close()

    triggerNames  = open("triggerTauList.txt", "r")
    for triggerName in triggerNames:
        print triggerName
        script = open ("plotTauTriggerEfficiency.py", "a")
        script.writelines("""

    #-----Passing Trigger ("""+triggerName.replace("\n","")+""")

                if triggerResults.accept(names.triggerIndex('"""+triggerName.replace("\n","")+"""')):

                    h['hMuTauTrig_M_"""+triggerName.replace("\n","")+"""'].Fill((mutau+mu).M(), genweight)
                    h['hMuTauTrig_muPt_"""+triggerName.replace("\n","")+"""'].Fill(mu.Pt(), genweight)
                    h['hMuTauTrig_tauPt_"""+triggerName.replace("\n","")+"""'].Fill(mutau.Pt(), genweight)
                    h['hMuTauTrig_jPt_"""+triggerName.replace("\n","")+"""'].Fill(j.Pt(), genweight)
                    h['hMuTauTrig_dR_"""+triggerName.replace("\n","")+"""'].Fill(mu.DeltaR(mutau), genweight)
                    h['hMuTauTrig_dRlj_"""+triggerName.replace("\n","")+"""'].Fill((mu+mutau).DeltaR(j), genweight)

                if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (triggerResults.accept(names.triggerIndex('"""+triggerName.replace("\n","")+"""'))) or (mu.Pt()>26 and (triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")))) or (mu.Pt()>50 and (triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")))):
                    h['hMuTauTrigTauJet_M_"""+triggerName.replace("\n","")+"""'].Fill((mutau+mu).M(), genweight)
                    h['hMuTauTrigTauJet_muPt_"""+triggerName.replace("\n","")+"""'].Fill(mu.Pt(), genweight)
                    h['hMuTauTrigTauJet_tauPt_"""+triggerName.replace("\n","")+"""'].Fill(mutau.Pt(), genweight)
                    h['hMuTauTrigTauJet_jPt_"""+triggerName.replace("\n","")+"""'].Fill(j.Pt(), genweight)
                    h['hMuTauTrigTauJet_dR_"""+triggerName.replace("\n","")+"""'].Fill(mu.DeltaR(mutau), genweight)
                    h['hMuTauTrigTauJet_dRlj_"""+triggerName.replace("\n","")+"""'].Fill((mu+mutau).DeltaR(j), genweight)

                if (j.Pt() > 500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))) or (mu.Pt()>26 and (triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")))) or (mu.Pt()>50 and (triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")))):
                    h['hMuTauTrigJet_M_"""+triggerName.replace("\n","")+"""'].Fill((mutau+mu).M(), genweight)
                    h['hMuTauTrigJet_muPt_"""+triggerName.replace("\n","")+"""'].Fill(mu.Pt(), genweight)
                    h['hMuTauTrigJet_tauPt_"""+triggerName.replace("\n","")+"""'].Fill(mutau.Pt(), genweight)
                    h['hMuTauTrigJet_jPt_"""+triggerName.replace("\n","")+"""'].Fill(j.Pt(), genweight)
                    h['hMuTauTrigJet_dR_"""+triggerName.replace("\n","")+"""'].Fill(mu.DeltaR(mutau), genweight)
                    h['hMuTauTrigJet_dRlj_"""+triggerName.replace("\n","")+"""'].Fill((mu+mutau).DeltaR(j), genweight)

""")

        script.close()

script = open ("plotTauTriggerEfficiency.py", "a")
script.writelines("""

out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

""")

script.close()
