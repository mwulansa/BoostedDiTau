import ROOT,sys
from DataFormats.FWLite import Events, Handle
from looseElectron import *
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

jobs=np.linspace(100, 1, 100)
masses=[10, 30, 50]
#jobs=[1]
#masses=[30]

handleMuon = Handle ('vector<reco::Muon>')
labelMuon = ('muons')

handleElectron = Handle ('vector<reco::GsfElectron>')
labelElectron = ('gedGsfElectrons')

handleJet = Handle ('vector<reco::PFJet>')
labelJet = ('ak4PFJetsCHS')

handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlinePrimaryVerticesWithBS')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('genParticles')

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults')

handleBs = Handle ('reco::BeamSpot')
labelBs = ("offlineBeamSpot")

handleConv = Handle ('vector<reco::Conversion>')
labelConv = ("particleFlowEGamma")

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

prefix="root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer16DR80Premix/"

for mass in masses:

    out=ROOT.TFile("h_plotSignalReco_m"+str(int(mass))+".root",'recreate')

    hMuMuBaseline_M = ROOT.TH1F ("hMuMuBaseline_M", "#mu - #mu inv. mass;M_{#mu#mu};N_{events}", 100, 0, 100)
    hMuMu_M = ROOT.TH1F ("hMuMu_M", "#mu - #mu inv. mass;M_{#mu#mu};N_{events}", 100, 0, 100)
    hMuMuTrig_M = ROOT.TH1F ("hMuMuTrig_M", "#mu - #mu inv. mass;M_{#mu#mu};N_{events}", 100, 0, 100)
    
    hEMuBaseline_M = ROOT.TH1F ("hEMuBaseline_M", "e - #mu inv. mass;M_{e#mu};N_{events}", 100, 0, 100)
    hEMu_M = ROOT.TH1F ("hEMu_M", "e - #mu inv. mass;M_{e#mu};N_{events}", 100, 0, 100)
    hEMuTrig_M = ROOT.TH1F ("hEMuTrig_M", "e - #mu inv. mass;M_{e#mu};N_{events}", 100, 0, 100)
    
    hMuMuBaseline_lPt = ROOT.TH1F ("hMuMuBaseline_lPt", "leading muon pt;P_{t};N_{events}", 50, 0, 500)
    hMuMu_lPt = ROOT.TH1F ("hMuMu_lPt", "leading muon pt;P_{t};N_{events}", 50, 0, 500)
    hMuMuTrig_lPt = ROOT.TH1F ("hMuMuTrig_lPt", "leading muon pt;P_{t};N_{events}", 50, 0, 500)
    
    hEMuBaseline_muPt = ROOT.TH1F ("hEMuBaseline_muPt", "muon pt;P_{t};N_{events}", 50, 0, 500)
    hEMu_muPt = ROOT.TH1F ("hEMu_muPt", "muon pt;P_{t};N_{events}", 50, 0, 500)
    hEMuTrig_muPt = ROOT.TH1F ("hEMuTrig_muPt", "muon pt;P_{t};N_{events}", 50, 0, 500)
    
    hMuMuBaseline_tPt = ROOT.TH1F ("hMuMuBaseline_tPt", "sub-leading muon pt;P_{t};N_{events}", 50, 0, 500)
    hMuMu_tPt = ROOT.TH1F ("hMuMu_tPt", "sub-leading muon pt;P_{t};N_{events}", 50, 0, 500)
    hMuMuTrig_tPt = ROOT.TH1F ("hMuMuTrig_tPt", "sub-leading muon pt;P_{t};N_{events}", 50, 0, 500)
    
    hEMuBaseline_ePt = ROOT.TH1F ("hEMuBaseline_ePt", "electron pt;P_{t};N_{events}", 50, 0, 500)
    hEMu_ePt = ROOT.TH1F ("hEMu_ePt", "electron pt;P_{t};N_{events}", 50, 0, 500)
    hEMuTrig_ePt = ROOT.TH1F ("hEMuTrig_ePt", "electron pt;P_{t};N_{events}", 50, 0, 500)
    
    hMu1Pt=ROOT.TH1F ("hMu1Pt", "muon pt;P_{t};N_{events}", 50, 0, 500)
    hMu1PtIsoMuTrig= ROOT.TH1F ("hMu1PtIsoMuTrig", "muon pt;P_{t};N_{events}", 50, 0, 500)
    hMu1PtMuTrig= ROOT.TH1F ("hMu1PtMuTrig", "muon pt;P_{t};N_{events}", 50, 0, 500)
    
    hJet1Pt=ROOT.TH1F ("hJet1Pt", "jet pt;P_{t};N_{events}", 200, 0, 2000)
    hJet2Pt=ROOT.TH1F ("hJet2Pt", "jet pt;P_{t};N_{events}", 200, 0, 2000)
    hJet1PtTrig=ROOT.TH1F ("hJetTrigPt", "triggered jet pt;P_{t};N_{events}", 200, 0, 2000)

    for job in jobs:
    
        events=Events(prefix+"ALP_m"+str(int(mass))+"_w1_htjmin400_RunIISummer16DR80Premix_AODSIM_"+str(int(job))+".root")

        print events.size()
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
            #print rho[1]
    
            selected_muons=[]
            for muon in muons:
                if not muon.isTrackerMuon() or not muon.isGlobalMuon(): continue
                if not muon.isPFMuon(): continue
                if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue
                if muon.pt()<3 or muon.eta()>2.4: continue
                if (muon.pfIsolationR04().sumChargedHadronPt+max(0,muon.pfIsolationR04().sumPhotonEt+muon.pfIsolationR04().sumNeutralHadronEt-0.5*muon.pfIsolationR04().sumPUPt))/muon.pt()<0.25:
                    selected_muons+=[muon]
    
            selected_muons.sort(key=lambda x: x.pt(), reverse=True)

            if len(selected_muons)>0:
                hMu1Pt.Fill(selected_muons[0].pt())
                if triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")):
                    hMu1PtMuTrig.Fill(selected_muons[0].pt())
                if triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")):
                    hMu1PtIsoMuTrig.Fill(selected_muons[0].pt())
    
            selected_electrons=[]
            for electron in electrons:
                if electron.pt()<7: continue
                if electron.eta()>2.5: continue
                if electron.isEB():
                    if not electron.full5x5_sigmaIetaIeta()<0.011 or not electron.hadronicOverEm()<0.289 or not abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.222 or not GsfEleEInverseMinusPInverse(electron)<0.241 or not abs(dEtaInSeed(electron))<0.00477 or not GsfEleMissingHitsCut(electron)<=1 or not GsfEleConversionVetoCut(electron, bs, convs) or not GsfEleEffAreaPFIsoCut(electron, rho)<0.0994 or not dZ(electron, vertex[0])<0.1 or not dXY(electron, vertex[0])<0.05:
                        continue
                if electron.isEE():
                    if not electron.full5x5_sigmaIetaIeta()<0.03 or not electron.hadronicOverEm()<0.101 or not abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.213 or not GsfEleEInverseMinusPInverse(electron)<0.14 or not abs(dEtaInSeed(electron))<0.00868 or not GsfEleMissingHitsCut(electron)<=1 or not GsfEleConversionVetoCut(electron, bs, convs) or not GsfEleEffAreaPFIsoCut(electron, rho)<0.107 or not dZ(electron, vertex[0])<0.2 or not dXY(electron, vertex[0])<0.1:
                        continue
                selected_electrons+=[electron]
    
            selected_electrons.sort(key=lambda x: x.pt(), reverse=True)
    
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
    
            if len(selected_jets)>=1:
                hJet1Pt.Fill(selected_jets[0].pt())
                if len(selected_jets)>=2:
                    hJet2Pt.Fill(selected_jets[1].pt())
                if triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")):
                    hJet1PtTrig.Fill(selected_jets[0].pt())
            
            if len(selected_muons)>=2 and len(selected_jets)>=1:
    
                mu1=ROOT.TLorentzVector()
                mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())
    
                mu2=ROOT.TLorentzVector()
                mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass())
    
                j=ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
    
                hMuMuBaseline_M.Fill((mu1+mu2).M())
                hMuMuBaseline_lPt.Fill(mu1.Pt())
                hMuMuBaseline_tPt.Fill(mu2.Pt())
    
                #if (selected_muons[0].pfIsolationR04().sumChargedHadronPt+max(0,selected_muons[0].pfIsolationR04().sumPhotonEt+selected_muons[0].pfIsolationR04().sumNeutralHadronEt-0.5*selected_muons[0].pfIsolationR04().sumPUPt))/selected_muons[0].pt()<0.25 and (selected_muons[1].pfIsolationR04().sumChargedHadronPt+max(0,selected_muons[1].pfIsolationR04().sumPhotonEt+selected_muons[1].pfIsolationR04().sumNeutralHadronEt-0.5*selected_muons[1].pfIsolationR04().sumPUPt))/selected_muons[1].pt()<0.25:
                    #hMuMuIsolation_M.Fill((mu1+mu2).M())
                    #hMuMuIsolation_lPt.Fill(mu1.Pt())
                    #hMuMuIsolation_tPt.Fill(mu2.Pt())
    
                if mu1.DeltaR(mu2)<0.4 and mu1.DeltaR(j)>0.8 and mu2.DeltaR(j)>0.8:
                    hMuMu_M.Fill((mu1+mu2).M())
                    hMuMu_lPt.Fill(mu1.Pt())
                    hMuMu_tPt.Fill(mu2.Pt())
                    
                    if triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")):
                        hMuMuTrig_M.Fill((mu1+mu2).M())
                        hMuMuTrig_lPt.Fill(mu1.Pt())
                        hMuMuTrig_tPt.Fill(mu2.Pt())
    
            if len(selected_muons)>=1 and len(selected_electrons)>=1 and len(selected_jets)>=1:
                mu=ROOT.TLorentzVector()
                mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())
    
                e=ROOT.TLorentzVector()
                e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())
    
                j=ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
            
                hEMuBaseline_M.Fill((mu+e).M())
                hEMuBaseline_muPt.Fill(mu.Pt())
                hEMuBaseline_ePt.Fill(e.Pt())
             
                #if (selected_muons[0].pfIsolationR04().sumChargedHadronPt+max(0,selected_muons[0].pfIsolationR04().sumPhotonEt+selected_muons[0].pfIsolationR04().sumNeutralHadronEt-0.5*selected_muons[0].pfIsolationR04().sumPUPt))/selected_muons[0].pt()<0.25 and mu.DeltaR(e)<0.4 and mu.DeltaR(j)>0.8 and e.DeltaR(j)>0.8:
                    
                if mu.DeltaR(e)<0.4 and mu.DeltaR(j)>0.8 and e.DeltaR(j)>0.8:
                    hEMu_M.Fill((mu+e).M())
                    hEMu_muPt.Fill(mu.Pt())
                    hEMu_ePt.Fill(e.Pt())
                    if triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")):
                        hEMuTrig_M.Fill((mu+e).M())
                        hEMuTrig_muPt.Fill(mu.Pt())
                        hEMuTrig_ePt.Fill(e.Pt())
    
    out.cd()
    hMuMuBaseline_M.Write()
    hMuMuBaseline_lPt.Write()
    hMuMuBaseline_tPt.Write()
    #hMuMuIsolation_M.Write()
    #hMuMuIsolation_lPt.Write()
    #hMuMuIsolation_tPt.Write()
    hMuMu_M.Write()
    hMuMu_lPt.Write()
    hMuMu_tPt.Write()
    hMuMuTrig_M.Write()
    hMuMuTrig_lPt.Write()
    hMuMuTrig_tPt.Write()
    
    
    hEMuBaseline_M.Write()
    hEMuBaseline_muPt.Write()
    hEMuBaseline_ePt.Write()
    hEMu_M.Write()
    hEMu_muPt.Write()
    hEMu_ePt.Write()
    hEMuTrig_M.Write()
    hEMuTrig_muPt.Write()
    hEMuTrig_ePt.Write()
    
    hJet1Pt.Write()
    hJet2Pt.Write()
    hJet1PtTrig.Write()

    hMu1Pt.Write()
    hMu1PtIsoMuTrig.Write()
    hMu1PtMuTrig.Write()

    out.Close()
    
        

    
