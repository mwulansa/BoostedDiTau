import ROOT,sys
from DataFormats.FWLite import Events, Handle
from looseElectron import *
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

jobs=np.linspace(100, 1, 100)
#masses=[10, 30, 50]
#jobs=[1]
masses=[30]

## handleMuon = Handle ('vector<reco::Muon>')
## labelMuon = ('muons')
## 
## handleElectron = Handle ('vector<reco::GsfElectron>')
## labelElectron = ('gedGsfElectrons')
## 
## handleJet = Handle ('vector<reco::PFJet>')
## labelJet = ('ak4PFJetsCHS')
## 
## handleVertex = Handle ('vector<reco::Vertex>')
## labelVertex = ('offlinePrimaryVerticesWithBS')
## 
## handleGenParticle = Handle ('vector<reco::GenParticle>')
## labelGenParticle = ('genParticles')
## 
## handleHLT = Handle ('edm::TriggerResults')
## labelHLT = ('TriggerResults')
## 
## handleBs = Handle ('reco::BeamSpot')
## labelBs = ("offlineBeamSpot")
## 
## handleConv = Handle ('vector<reco::Conversion>')
## labelConv = ('conversions')
## 
## handleRho = Handle ('double')
## labelRho = ('fixedGridRhoFastjetAll')

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('slimmedJets')

handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlineSlimmedPrimaryVertices')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('genParticles')

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults','','HLT')

handleBs = Handle ('reco::BeamSpot')
labelBs = ("offlineBeamSpot")

handleConv = Handle ('vector<reco::Conversion>')
labelConv = ('reducedEgamma', 'reducedConversions')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handleGenJet = Handle ('vector<reco::GenJet>')
labelGenJet = ('slimmedGenJets')

prefix="root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer16DR80Premix/"

for mass in masses:

    out=ROOT.TFile("h_plotBoostedDiLepReco_m"+str(int(mass))+"_v1.root",'recreate')

    hMu1Pt=ROOT.TH1F ("hMu1Pt", "muon pt;P_{t};N_{events}", 50, 0, 500)
    hMu1PtIsoMuTrig= ROOT.TH1F ("hMu1PtIsoMuTrig", "muon pt;P_{t};N_{events}", 50, 0, 500)
    hMu1PtMuTrig= ROOT.TH1F ("hMu1PtMuTrig", "muon pt;P_{t};N_{events}", 50, 0, 500)
    
    hJet1Pt=ROOT.TH1F ("hJet1Pt", "jet pt;P_{t};N_{events}", 200, 0, 2000)
    hJet2Pt=ROOT.TH1F ("hJet2Pt", "jet pt;P_{t};N_{events}", 200, 0, 2000)
    hJet1PtTrig=ROOT.TH1F ("hJetTrigPt", "triggered jet pt;P_{t};N_{events}", 200, 0, 2000)

    hJet1Pt_genMatched=ROOT.TH1F ("hJet1Pt_genMatched", "jet pt;P_{t};N_{events}", 200, 0, 2000)

    hMuMuBaseline_M = ROOT.TH1F ("hMuMuBaseline_M", "#mu - #mu inv. mass;M_{#mu#mu};N_{events}", 100, 0, 100)
    hMuMu_M = ROOT.TH1F ("hMuMu_M", "#mu - #mu inv. mass;M_{#mu#mu};N_{events}", 100, 0, 100)
    hMuMuTrig_M = ROOT.TH1F ("hMuMuTrig_M", "#mu - #mu inv. mass;M_{#mu#mu};N_{events}", 100, 0, 100)
    hMuMu_M_genMatched = ROOT.TH1F ("hMuMu_M_genMatched", "#mu - #mu inv. mass;M_{#mu#mu};N_{events}", 100, 0, 100)

    hMuMuBaseline_lPt = ROOT.TH1F ("hMuMuBaseline_Mu1Pt", "leading muon pt;P_{t};N_{events}", 50, 0, 500)
    hMuMu_lPt = ROOT.TH1F ("hMuMu_Mu1Pt", "leading muon pt;P_{t};N_{events}", 50, 0, 500)
    hMuMu_lPt_genMatched = ROOT.TH1F ("hMuMu_Mu1Pt_genMatched", "leading muon pt;P_{t};N_{events}", 50, 0, 500)
    hMuMuTrig_lPt = ROOT.TH1F ("hMuMuTrig_Mu1Pt", "leading muon pt;P_{t};N_{events}", 50, 0, 500)

    hMuMuBaseline_tPt = ROOT.TH1F ("hMuMuBaseline_Mu2Pt", "sub-leading muon pt;P_{t};N_{events}", 50, 0, 500)
    hMuMu_tPt = ROOT.TH1F ("hMuMu_Mu2Pt", "sub-leading muon pt;P_{t};N_{events}", 50, 0, 500)
    hMuMu_tPt_genMatched = ROOT.TH1F ("hMuMu_Mu2Pt_genMatched", "sub-leading muon pt;P_{t};N_{events}", 50, 0, 500)
    hMuMuTrig_tPt = ROOT.TH1F ("hMuMuTrig_Mu2Pt", "sub-leading muon pt;P_{t};N_{events}", 50, 0, 500)
    
    hEMuBaseline_M = ROOT.TH1F ("hEMuBaseline_M", "e - #mu inv. mass;M_{e#mu};N_{events}", 100, 0, 100)
    hEMu_M = ROOT.TH1F ("hEMu_M", "e - #mu inv. mass;M_{e#mu};N_{events}", 100, 0, 100)
    hEMuTrig_M = ROOT.TH1F ("hEMuTrig_M", "e - #mu inv. mass;M_{e#mu};N_{events}", 100, 0, 100)
    hEMu_M_genMatched = ROOT.TH1F ("hEMu_M_genMatched", "e - #mu inv. mass;M_{e#mu};N_{events}", 100, 0, 100)
    
    hEMuBaseline_muPt = ROOT.TH1F ("hEMuBaseline_MuPt", "muon pt;P_{t};N_{events}", 50, 0, 500)
    hEMu_muPt = ROOT.TH1F ("hEMu_MuPt", "muon pt;P_{t};N_{events}", 50, 0, 500)
    hEMu_muPt_genMatched = ROOT.TH1F ("hEMu_MuPt_genMatched", "muon pt;P_{t};N_{events}", 50, 0, 500)
    hEMuTrig_muPt = ROOT.TH1F ("hEMuTrig_MuPt", "muon pt;P_{t};N_{events}", 50, 0, 500)
    
    hEMuBaseline_ePt = ROOT.TH1F ("hEMuBaseline_EPt", "electron pt;P_{t};N_{events}", 50, 0, 500)
    hEMu_ePt = ROOT.TH1F ("hEMu_EPt", "electron pt;P_{t};N_{events}", 50, 0, 500)
    hEMu_ePt_genMatched = ROOT.TH1F ("hEMu_EPt_genMatched", "electron pt;P_{t};N_{events}", 50, 0, 500)
    hEMuTrig_ePt = ROOT.TH1F ("hEMuTrig_EPt", "electron pt;P_{t};N_{events}", 50, 0, 500)

    hMuMu_dR = ROOT.TH1F ("hMuMu_dR", "#mu - #mu delta R;#delta R;N_{events}", 20, 0, 1)
    hEMu_dR = ROOT.TH1F ("hEMu_dR", "e - #mu delta R;#delta R;N_{events}", 20, 0, 1)

    hMuMu_dR_genMatched = ROOT.TH1F ("hMuMu_dR_genMatched", "#mu - #mu delta R;#delta R;N_{events}", 20, 0, 1)
    hEMu_dR_genMatched = ROOT.TH1F ("hEMu_dR_genMatched", "e - #mu delta R;#delta R;N_{events}", 20, 0, 1)
    

    for job in jobs:
        print prefix+"ALP_m"+str(int(mass))+"_w1_htjmin400_RunIISummer16DR80Premix_miniAODSIM_"+str(int(job))+".root"
        events=Events(prefix+"ALP_m"+str(int(mass))+"_w1_htjmin400_RunIISummer16DR80Premix_miniAODSIM_"+str(int(job))+".root")

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

            event.getByLabel(labelGenParticle, handleGenParticle)
            particles=handleGenParticle.product()

            # gen information
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

            event.getByLabel(labelGenJet, handleGenJet)
            #jets=handleGenJet.product()
        
            genJets=[]
            for jet in handleGenJet.product():
                genJets+=[jet]
                
            genJets.sort(key=lambda x: x.pt(), reverse=True)
        
            genJet1=ROOT.TLorentzVector(genJets[0].px(), genJets[0].py(), genJets[0].pz(), genJets[0].energy())

            # muon selection
            selected_muons=[]
            for muon in muons:
                if not muon.isLooseMuon(): continue
                if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue
                if muon.pt()<3 or muon.eta()>2.4: continue
                if (muon.pfIsolationR04().sumChargedHadronPt+max(0,muon.pfIsolationR04().sumPhotonEt+muon.pfIsolationR04().sumNeutralHadronEt-0.5*muon.pfIsolationR04().sumPUPt))/muon.pt()<0.4:
                    selected_muons+=[muon]
    
            selected_muons.sort(key=lambda x: x.pt(), reverse=True)
            
            if len(selected_muons)>0:
                hMu1Pt.Fill(selected_muons[0].pt())
                if triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")):
                    hMu1PtMuTrig.Fill(selected_muons[0].pt())
                if triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")):
                    hMu1PtIsoMuTrig.Fill(selected_muons[0].pt())

            # electron selection
            selected_electrons=[]
            for electron in electrons:
                if electron.pt()<7: continue
                if abs(electron.eta())>2.5: continue
                if electron.isEB():
                    if not electron.full5x5_sigmaIetaIeta()<0.011 or not electron.hadronicOverEm()<0.289 or not abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.222 or not GsfEleEInverseMinusPInverse(electron)<0.241 or not abs(dEtaInSeed(electron))<0.00477 or not GsfEleMissingHitsCut(electron)<=1 or not electron.passConversionVeto() or not GsfEleEffAreaPFIsoCut(electron, rho)<0.0994 or not abs(electron.gsfTrack().dz(vertex[0].position()))<0.1 or not abs(electron.gsfTrack().dxy(vertex[0].position()))<0.05:
                        continue
                if electron.isEE():
                    if not electron.full5x5_sigmaIetaIeta()<0.03 or not electron.hadronicOverEm()<0.101 or not abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.213 or not GsfEleEInverseMinusPInverse(electron)<0.14 or not abs(dEtaInSeed(electron))<0.00868 or not GsfEleMissingHitsCut(electron)<=1 or not electron.passConversionVeto() or not GsfEleEffAreaPFIsoCut(electron, rho)<0.107 or not abs(electron.gsfTrack().dz(vertex[0].position()))<0.2 or not abs(electron.gsfTrack().dxy(vertex[0].position()))<0.1:
                        continue
                selected_electrons+=[electron]
    
            selected_electrons.sort(key=lambda x: x.pt(), reverse=True)

            # jet selection
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
                j=ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                if genJet1.DeltaR(j)<0.4:
                    hJet1Pt_genMatched.Fill(selected_jets[0].pt())
                if len(selected_jets)>=2:
                    hJet2Pt.Fill(selected_jets[1].pt())
                if triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")):
                    hJet1PtTrig.Fill(selected_jets[0].pt())

            # tau_mu tau_mu selection
            if len(selected_muons)>=2 and len(selected_jets)>=1:
                mu1=ROOT.TLorentzVector()
                mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())
    
                mu2=ROOT.TLorentzVector()
                mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass())

                hMuMuBaseline_M.Fill((mu1+mu2).M())
                hMuMuBaseline_lPt.Fill(mu1.Pt())
                hMuMuBaseline_tPt.Fill(mu2.Pt())

                j=ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
    
                #if mu1.DeltaR(mu2)<0.4 and mu1.DeltaR(j)>0.8 and mu2.DeltaR(j)>0.8:
                if mu1.DeltaR(mu2)<0.4 and mu1.DeltaR(j)>0.8 and mu2.DeltaR(j)>0.8 and (selected_muons[0].pfIsolationR04().sumChargedHadronPt+max(0,selected_muons[0].pfIsolationR04().sumPhotonEt+selected_muons[0].pfIsolationR04().sumNeutralHadronEt-0.5*selected_muons[0].pfIsolationR04().sumPUPt))/selected_muons[0].pt()<0.25:
                    hMuMu_dR.Fill(mu1.DeltaR(mu2))
                    hMuMu_M.Fill((mu1+mu2).M())
                    hMuMu_lPt.Fill(mu1.Pt())
                    hMuMu_tPt.Fill(mu2.Pt())

                    # gen match
                    evtMatched=False
                    if len(genMuons)==2 and len(genJets)>0:
                        if not len(genElectrons)==0: print "Debug 1: Something is Wrong..."
                        genMu1=ROOT.TLorentzVector()
                        genMu1.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

                        genMu2=ROOT.TLorentzVector()
                        genMu2.SetPtEtaPhiM(genMuons[1].pt(), genMuons[1].eta(), genMuons[1].phi(), genMuons[1].mass())

                        if (mu1.DeltaR(genMu1)<0.3 or mu1.DeltaR(genMu2)<0.3) and (mu2.DeltaR(genMu1)<0.3 or mu2.DeltaR(genMu2)<0.3):
                            hMuMu_M_genMatched.Fill((mu1+mu2).M())
                            hMuMu_dR_genMatched.Fill(mu1.DeltaR(mu2))
                            hMuMu_lPt_genMatched.Fill(mu1.Pt())
                            hMuMu_tPt_genMatched.Fill(mu2.Pt())
                        
                    # trig
                    if triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")):
                        hMuMuTrig_M.Fill((mu1+mu2).M())
                        hMuMuTrig_lPt.Fill(mu1.Pt())
                        hMuMuTrig_tPt.Fill(mu2.Pt())

                        
            # tau_mu tau_e selection
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
                
                if mu.DeltaR(e)<0.4 and mu.DeltaR(j)>0.8 and e.DeltaR(j)>0.8:
                    hEMu_M.Fill((mu+e).M())
                    hEMu_muPt.Fill(mu.Pt())
                    hEMu_ePt.Fill(e.Pt())
                    hEMu_dR.Fill(mu.DeltaR(e))

                    # gen match
                    if len(genMuons)==1 and len(genElectrons)==1:
                        genMu=ROOT.TLorentzVector()
                        genMu.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

                        genEle=ROOT.TLorentzVector()
                        genEle.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())

                        if e.DeltaR(genEle)<0.3 and mu.DeltaR(genMu)<0.3:
                            hEMu_M_genMatched.Fill((e+mu).M())
                            hEMu_dR_genMatched.Fill(mu.DeltaR(e))
                            hEMu_muPt_genMatched.Fill(mu.Pt())
                            hEMu_ePt_genMatched.Fill(e.Pt())
                            
                    # trig
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
    hMuMu_M_genMatched.Write()
    hMuMu_lPt.Write()
    hMuMu_tPt.Write()
    hMuMu_lPt_genMatched.Write()
    hMuMu_tPt_genMatched.Write()
    hMuMuTrig_M.Write()
    hMuMuTrig_lPt.Write()
    hMuMuTrig_tPt.Write()
    
    hEMuBaseline_M.Write()
    hEMuBaseline_muPt.Write()
    hEMuBaseline_ePt.Write()
    hEMu_M.Write()
    hEMu_M_genMatched.Write()
    hEMu_muPt.Write()
    hEMu_ePt.Write()
    hEMu_muPt_genMatched.Write()
    hEMu_ePt_genMatched.Write()
    hEMuTrig_M.Write()
    hEMuTrig_muPt.Write()
    hEMuTrig_ePt.Write()
    
    hJet1Pt.Write()
    hJet1Pt_genMatched.Write()
    hJet2Pt.Write()
    hJet1PtTrig.Write()

    hMu1Pt.Write()
    hMu1PtIsoMuTrig.Write()
    hMu1PtMuTrig.Write()

    hMuMu_dR.Write()
    hEMu_dR.Write()

    hMuMu_dR_genMatched.Write()
    hEMu_dR_genMatched.Write()

    out.Close()
    
        

    
