import ROOT, sys, math, gc
from DataFormats.FWLite import Events, Handle
import numpy as np

from looseElectron import *

sample = "QCD"

handleTaus = Handle ('vector<pat::Tau>')
labelTaus = ('slimmedTaus')

handleTausNewID = Handle ('vector<pat::Tau>')
labelTausNewID = ('slimmedTausNewID')

handleTausElectronCleaned = Handle ('vector<pat::Tau>')
labelTausElectronCleaned = ('slimmedTausElectronCleaned', '', 'PAT')

handleTausNewIDElectronCleaned = Handle ('vector<pat::Tau>')
labelTausNewIDElectronCleaned = ('slimmedTausNewIDElectronCleaned', '', 'PAT')

handleTausMuonCleaned = Handle ('vector<pat::Tau>')
labelTausMuonCleaned = ('slimmedTausMuonCleaned', '', 'PAT')

handleTausNewIDMuonCleaned = Handle ('vector<pat::Tau>')
labelTausNewIDMuonCleaned = ('slimmedTausNewIDMuonCleaned', '', 'PAT')

handleTausBoosted = Handle ('vector<pat::Tau>')
labelTausBoosted = ('slimmedTausBoosted', '', 'PAT')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

prefix="root://cmseos.fnal.gov//eos/uscms/store/user/mwulansa/DIS/TCPAnalysis/Backgrounds/RunIIUL17/"

jobs=np.linspace(1, 0, 2)

outputFileName = "h_"+sample+"_cleaningValidation.root"
print outputFileName

out=ROOT.TFile.Open(outputFileName,'recreate')

h = {}

h['hNEvent'] = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)
h['hNTaus'] = ROOT.TH1F ("hNTaus", "Number of Taus;;N_{events}", 7, 0, 7)

h['hTauPt'] = ROOT.TH1F ("hTau_Pt", "#tau P_t; P_t;", 500, 0, 500)
h['heTauPt'] = ROOT.TH1F ("heTau_Pt", "#tau P_t; P_t;", 500, 0, 500)
h['hmTauPt'] = ROOT.TH1F ("hmTau_Pt", "#tau P_t; P_t;", 500, 0, 500)

h['hTauNewIDPt'] = ROOT.TH1F ("hTau_NewIDPt", "#tau P_t; P_t;", 500, 0, 500)
h['heTauNewIDPt'] = ROOT.TH1F ("heTau_NewIDPt", "#tau P_t; P_t;", 500, 0, 500)
h['hmTauNewIDPt'] = ROOT.TH1F ("hmTau_NewIDPt", "#tau P_t; P_t;", 500, 0, 500)

h['hbTauPt'] = ROOT.TH1F ("hbTau_Pt", "#tau P_t; P_t;", 500, 0, 500)

for job in jobs:

    filename = "cleaned_"+sample+"_Pt-15to7000_"+str(int(job))+".root"

    inputFileName = prefix+filename

    print inputFileName

    f = ROOT.TFile.Open(inputFileName)

    if not f.IsZombie():
        events=Events(inputFileName)
    else:
        print "Can't Open File: "+inputFileName
        continue

    for event in events:
    
        event.getByLabel(labelTausElectronCleaned, handleTausElectronCleaned)
        tausElectronCleaned=handleTausElectronCleaned.product()

        event.getByLabel(labelTausMuonCleaned, handleTausMuonCleaned)
        tausMuonCleaned=handleTausMuonCleaned.product()

        event.getByLabel(labelTaus, handleTaus)
        taus=handleTaus.product()

        event.getByLabel(labelTausNewIDElectronCleaned, handleTausNewIDElectronCleaned)
        tausNewIDElectronCleaned=handleTausNewIDElectronCleaned.product()

        event.getByLabel(labelTausNewIDMuonCleaned, handleTausNewIDMuonCleaned)
        tausNewIDMuonCleaned=handleTausNewIDMuonCleaned.product()

        event.getByLabel(labelTausNewID, handleTausNewID)
        tausNewID=handleTausNewID.product()

        event.getByLabel(labelTausBoosted, handleTausBoosted)
        tausBoosted = handleTausBoosted.product()

        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()

        h['hNEvent'].Fill(0.5, 1)
        h['hNEvent'].Fill(1.5, genweight)

        selected_taus=[]
        if taus.size()>0:
            for tau in taus:
                if tau.pt()>20 and abs(tau.eta())<2.3:
                    selected_taus+=[tau]

        if len(selected_taus)>0:
            selected_taus.sort(key=lambda x: x.pt(), reverse=True)
            h['hNTaus'].Fill(0.5 , 1)
            h['hTauPt'].Fill(selected_taus[0].pt())

        selected_etaus=[]
        if tausElectronCleaned.size()>0:
            for tau in tausElectronCleaned:
                if tau.pt()>20 and abs(tau.eta())<2.3:
                    selected_etaus+=[tau]
                    
        if len(selected_etaus)>0:
            selected_etaus.sort(key=lambda x: x.pt(), reverse=True)
            h['hNTaus'].Fill(1.5 , 1)
            h['heTauPt'].Fill(selected_etaus[0].pt())

        selected_mtaus=[]
        if tausMuonCleaned.size()>0:
            for tau in tausMuonCleaned:
                if tau.pt()>20 and abs(tau.eta())<2.3:
                    selected_mtaus+=[tau]
            
        if len(selected_mtaus)>0:
            selected_mtaus.sort(key=lambda x: x.pt(), reverse=True)
            h['hNTaus'].Fill(2.5 , 1)
            h['hmTauPt'].Fill(selected_mtaus[0].pt())

        selected_taus_newID=[]
        if tausNewID.size()>0:
            for tau in tausNewID:
                if tau.pt()>20 and abs(tau.eta())<2.3:
                    selected_taus_newID+=[tau]
        
        if len(selected_taus_newID)>0:            
            selected_taus_newID.sort(key=lambda x: x.pt(), reverse=True)
            h['hNTaus'].Fill(3.5 , 1)
            h['hTauNewIDPt'].Fill(selected_taus_newID[0].pt())

        selected_etaus_newID=[]
        if tausNewIDElectronCleaned.size()>0:
            for tau in tausNewIDElectronCleaned:
                if tau.pt()>20 and abs(tau.eta())<2.3:
                    selected_etaus_newID+=[tau]

        if len(selected_etaus_newID)>0:
            selected_etaus_newID.sort(key=lambda x: x.pt(), reverse=True)
            h['hNTaus'].Fill(4.5 , 1)
            h['heTauNewIDPt'].Fill(selected_etaus_newID[0].pt())

        selected_mtaus_newID=[]
        if tausNewIDMuonCleaned.size()>0:
            for tau in tausNewIDMuonCleaned:
                if tau.pt()>20 and abs(tau.eta())<2.3:
                    selected_mtaus_newID+=[tau]

        if len(selected_mtaus_newID)>0:
            selected_mtaus_newID.sort(key=lambda x: x.pt(), reverse=True)
            h['hNTaus'].Fill(5.5 , 1)
            h['hmTauNewIDPt'].Fill(selected_mtaus_newID[0].pt())

        selected_btaus = []
        if tausBoosted.size>0:
            for tau in tausBoosted:
                if tau.pt()>20 and abs(tau.eta())<2.3:
                    selected_btaus+=[tau]

        if len(selected_btaus)>0:
            selected_btaus.sort(key=lambda x: x.pt(), reverse=True)
            h['hbTauPt'].Fill(selected_btaus[0].pt())
            h['hNTaus'].Fill(6.5 , 1)


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
