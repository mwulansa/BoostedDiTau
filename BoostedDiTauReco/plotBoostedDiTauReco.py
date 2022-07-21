import ROOT,sys
from DataFormats.FWLite import Events, Handle
import numpy as np

from looseElectron import *

jobs=np.linspace(100, 1, 100)
#jobs=[1]

if len(sys.argv)>1:
    mass=sys.argv[1]
else:
    mass='10'
    print "Use default signal mass: 10 GeV."

prefix="root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer19UL17RECO/"

out=ROOT.TFile("h_plotBoostedDiTauReco_m"+mass+"_v1.root",'recreate')

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
labelHLT = ('TriggerResults')

handleBs = Handle ('reco::BeamSpot')
labelBs = ("offlineBeamSpot")

handleConv = Handle ('vector<reco::Conversion>')
labelConv = ('reducedEgamma', 'reducedConversions')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleTaus = Handle ('vector<pat::Tau>')
labelTaus = ('slimmedTaus')

handleElectronCleanedTaus = Handle ('vector<pat::Tau>')
labelElectronCleanedTaus = ('slimmedTausElectronCleaned')

handleMuonCleanedTaus = Handle ('vector<pat::Tau>')
labelMuonCleanedTaus = ('slimmedTausMuonCleaned')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

# book histograms
hETau_M = ROOT.TH1F ("hETau_M", "e - #tau mass;M_{e#tau};N_{events}", 100, 0, 100)
hMuTau_M = ROOT.TH1F ("hMuTau_M", "mu - #tau mass;M_{#mu#tau};N_{events}", 100, 0, 100)

hETau_M_eCleaned = ROOT.TH1F ("hETau_ECleaned_M", "e - #tau mass;M_{e#tau};N_{events}", 100, 0, 100)
hETau_M_eCleaned_genMatched = ROOT.TH1F ("hETau_ECleaned_M_genMatched", "e - #tau mass;M_{e#tau};N_{events}", 100, 0, 100)

hMuTau_M_muCleaned = ROOT.TH1F ("hMuTau_MuCleaned_M", "#mu - #tau mass;M_{#mu#tau};N_{events}", 100, 0, 100)
hMuTau_M_muCleaned_genMatched = ROOT.TH1F ("hMuTau_MuCleaned_M_genMatched", "#mu - #tau mass;M_{#mu#tau};N_{events}", 100, 0, 100)

hETau_dR_eCleaned = ROOT.TH1F ("hETau_ECleaned_dR", "e - #Delta R;M_{e#tau};N_{events}", 20, 0, 1)
hMuTau_dR_muCleaned = ROOT.TH1F ("hMuTau_MuCleaned_dR", "#Delta R;M_{#mu#tau};N_{events}", 20, 0, 1)

hETau_dR_eCleaned_genMatched = ROOT.TH1F ("hETau_ECleaned_dR_genMatched", "e - #Delta R;M_{e#tau};N_{events}", 20, 0, 1)
hMuTau_dR_muCleaned_genMatched = ROOT.TH1F ("hMuTau_MuCleaned_dR_genMatched", "#Delta R;M_{#mu#tau};N_{events}", 20, 0, 1)

hETau_tauPt_eCleaned = ROOT.TH1F ("hETau_ECleaned_TauPt", "Tau Pt;P_{t,#tau};N_{events}", 50, 0, 500)
hMuTau_tauPt_muCleaned = ROOT.TH1F ("hMuTau_MuCleaned_TauPt", "Tau Pt;P_{t,#tau};N_{events}", 50, 0, 500)

hETau_ePt_eCleaned = ROOT.TH1F ("hETau_ECleaned_EPt", "Electron Pt;P_{t,e};N_{events}", 50, 0, 500)
hMuTau_muPt_muCleaned = ROOT.TH1F ("hMuTau_MuCleaned_MuPt", "Muon Pt;P_{t,#mu};N_{events}", 50, 0, 500)

hETau_tauPt_eCleaned_genMatched = ROOT.TH1F ("hETau_ECleaned_TauPt_genMatched", "Tau Pt;P_{t,#tau};N_{events}", 50, 0, 500)
hMuTau_tauPt_muCleaned_genMatched = ROOT.TH1F ("hMuTau_MuCleaned_TauPt_genMatched", "Tau Pt;P_{t,#tau};N_{events}", 50, 0, 500)

hETau_ePt_eCleaned_genMatched = ROOT.TH1F ("hETau_ECleaned_EPt_genMatched", "Electron Pt;P_{t,e};N_{events}", 50, 0, 500)
hMuTau_muPt_muCleaned_genMatched = ROOT.TH1F ("hMuTau_MuCleaned_MuPt_genMatched", "Muon Pt;P_{t,#mu};N_{events}", 50, 0, 500)

hNEvent = ROOT.TH1F ("hNEvent","Number of Events;;N{event}",1,0,2)

# loop over events
for job in jobs:
    print prefix+"ALP_m"+mass+"_w1_htjmin400_RunIISummer16DR80Premix_miniAODSIM_"+str(int(job))+".root"
    events=Events(prefix+"ALP_m"+mass+"_w1_htjmin400_RunIISummer16DR80Premix_miniAODSIM_"+str(int(job))+".root")

    ntot=0 
    nmatch1=0
    nmatch2=0
    
    for event in events:
        ntot+=1

        hNEvent.Fill(1)

        #if ntot>4: continue

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

        event.getByLabel(labelTaus, handleTaus)
        taus=handleTaus.product()

        event.getByLabel(labelElectronCleanedTaus, handleElectronCleanedTaus)
        etaus=handleElectronCleanedTaus.product()

        event.getByLabel(labelMuonCleanedTaus, handleMuonCleanedTaus)
        mutaus=handleMuonCleanedTaus.product()

        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()

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

    #Muon selection
        selected_muons=[]
        for muon in muons:
            if not muon.isLooseMuon(): continue
            if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue 
            if muon.pt()<3 or muon.eta()>2.4: continue
            selected_muons+=[muon]
        selected_muons.sort(key=lambda x: x.pt(), reverse=True) 
    
    #Electron selection
        selected_electrons=[]
        for electron in electrons:
            if electron.pt()<7: continue
            if abs(electron.eta())>2.5: continue
            if electron.isEB(): 
                if not electron.full5x5_sigmaIetaIeta()<0.011 or not electron.hadronicOverEm()<0.298 or not abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.222 or not GsfEleEInverseMinusPInverse(electron)<0.241 or not abs(dEtaInSeed(electron))<0.00477 or not GsfEleMissingHitsCut(electron)<=1 or not electron.passConversionVeto() or not abs(electron.gsfTrack().dz(vertex[0].position()))<0.1 or not abs(electron.gsfTrack().dxy(vertex[0].position()))<0.05:
                    continue
            if electron.isEE(): 
                if not electron.full5x5_sigmaIetaIeta()<0.0314 or not electron.hadronicOverEm()<0.101 or not abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.213 or not GsfEleEInverseMinusPInverse(electron)<0.14 or not abs(dEtaInSeed(electron))<0.00868 or not GsfEleMissingHitsCut(electron)<=1 or not electron.passConversionVeto() or not abs(electron.gsfTrack().dz(vertex[0].position()))<0.2 or not abs(electron.gsfTrack().dxy(vertex[0].position()))<0.1:
                    continue
            selected_electrons+=[electron]
        selected_electrons.sort(key=lambda x: x.pt(), reverse=True) 

    #jet selection
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

    #tau selection        
        selected_taus=[]
        for tau in taus:
            if tau.pt()<10 or abs(tau.eta())>2.3 or tau.leadChargedHadrCand().get().dxy(vertex[0].position())>0.2 or tau.leadChargedHadrCand().get().dz(vertex[0].position())>0.5 or not tau.tauID("decayModeFinding") or tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") < -0.5 or not tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
                continue
            selected_taus+=[tau]
        selected_taus.sort(key=lambda x: x.pt(), reverse=True)
            
        selected_etaus=[] 
        for tau in etaus:
            if tau.pt()<10 or abs(tau.eta())>2.3 or tau.leadChargedHadrCand().get().dxy(vertex[0].position())>0.2 or tau.leadChargedHadrCand().get().dz(vertex[0].position())>0.5 or not tau.tauID("decayModeFinding") or tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") < -0.5 or not tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
                continue
            selected_etaus+=[tau]
        selected_etaus.sort(key=lambda x: x.pt(), reverse=True)

        selected_mutaus=[] 
        for tau in mutaus:
            if tau.pt()<10 or abs(tau.eta())>2.3 or abs(tau.leadChargedHadrCand().get().dxy(vertex[0].position()))>0.2 or abs(tau.leadChargedHadrCand().get().dz(vertex[0].position()))>0.5 or not tau.tauID("decayModeFinding") or tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") < -0.5 or not tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
                continue
            selected_mutaus+=[tau]
        selected_mutaus.sort(key=lambda x: x.pt(), reverse=True)
        
    # tau_mu tau_had
        if len(selected_muons)>0 and len(selected_taus)>0 and len(selected_jets)>0: 
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            if mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                hMuTau_M.Fill((mu+tau).M()) 

    # tau_e tau_had
        if len(selected_electrons)>0 and len(selected_taus)>0 and len(selected_jets)>0: 
            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            if e.DeltaR(tau)<0.4 and e.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                hETau_M.Fill((e+tau).M()) 

    # muon cleaned tau_mu tau_had
        if len(selected_muons)>0 and len(selected_mutaus)>0 and len(selected_jets)>0:
        #if len(selected_muons)>0 and len(selected_mutaus)>0 and len(selected_jets)>0 and len(selected_electrons)==0:
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_mutaus[0].pt(), selected_mutaus[0].eta(), selected_mutaus[0].phi(), selected_mutaus[0].mass()) 
            
            if len(selected_jets)==1:
                j0=ROOT.TLorentzVector()
                j0.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j=j0
            else: 
                j1=ROOT.TLorentzVector()
                j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j2=ROOT.TLorentzVector()
                j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())
                if tau.DeltaR(j1)>0.3: 
                    j=j1
                else:
                    j=j2 

            if mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                hMuTau_M_muCleaned.Fill((mu+tau).M()) 
                hMuTau_dR_muCleaned.Fill(tau.DeltaR(mu))
                hMuTau_tauPt_muCleaned.Fill(tau.Pt())
                hMuTau_muPt_muCleaned.Fill(mu.Pt())

                # gen match
                #evtMatched=False
                if len(genElectrons)==0 and len(genMuons)==1: 
                #    evtMatched=True
                    nmatch1+=1
                    genMu=ROOT.TLorentzVector()
                    genMu.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())
                    genTau=ROOT.TLorentzVector()
                    genNt=ROOT.TLorentzVector()
                    if genMuons[0].pdgId()==13:
                        genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                        genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass()) 
                    else:
                        genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                        genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                        
                    if mu.DeltaR(genMu)<0.3 and tau.DeltaR(genTau-genNt)<0.4:
                        nmatch2+=1
                        hMuTau_M_muCleaned_genMatched.Fill((mu+tau).M())
                        hMuTau_dR_muCleaned_genMatched.Fill(tau.DeltaR(mu))
                        hMuTau_tauPt_muCleaned_genMatched.Fill(tau.Pt())
                        hMuTau_muPt_muCleaned_genMatched.Fill(mu.Pt())
                        
                #if not evtMatched:
                #    print "tauMu_tauHad", len(genElectrons), len(genMuons), len(genTaus), len(genNutaus)

        # electron cleaned tau_e tau_had
        if len(selected_electrons)>0 and len(selected_etaus)>0 and len(selected_jets)>0: 
        #if len(selected_electrons)>0 and len(selected_etaus)>0 and len(selected_jets)>0 and len(selected_muons)==0:
            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_etaus[0].pt(), selected_etaus[0].eta(), selected_etaus[0].phi(), selected_etaus[0].mass()) 

            if len(selected_jets)==1:
                j0=ROOT.TLorentzVector()
                j0.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j=j0
            else:
                j1=ROOT.TLorentzVector()
                j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j2=ROOT.TLorentzVector()
                j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())
                if tau.DeltaR(j1)>0.3:
                    j=j1
                else:
                    j=j2

            if e.DeltaR(tau)<0.4 and e.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8:
                hETau_M_eCleaned.Fill((e+tau).M())
                hETau_dR_eCleaned.Fill(tau.DeltaR(e))
                hETau_tauPt_eCleaned.Fill(tau.Pt())
                hETau_ePt_eCleaned.Fill(e.Pt())

                # gen match
                #evtMatched=False
                if len(genElectrons)==1 and len(genMuons)==0:
                    #evtMatched=True
                    genEle=ROOT.TLorentzVector()
                    genEle.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                    genTau=ROOT.TLorentzVector()
                    genNt=ROOT.TLorentzVector()
                    if genElectrons[0].pdgId()==11:
                        genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                        genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                    else:
                        genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                        genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                        
                    if e.DeltaR(genEle)<0.3 and tau.DeltaR(genTau-genNt)<0.4:
                        hETau_M_eCleaned_genMatched.Fill((e+tau).M())
                        hETau_dR_eCleaned_genMatched.Fill(tau.DeltaR(e))
                        hETau_tauPt_eCleaned_genMatched.Fill(tau.Pt())
                        hETau_ePt_eCleaned_genMatched.Fill(e.Pt())
                #if not evtMatched:
                #    print "tauE_tauHad", len(genElectrons), len(genMuons), len(genTaus), len(genNutaus)

out.cd()
hETau_M.Write()
hMuTau_M.Write()

hETau_M_eCleaned.Write()
hETau_dR_eCleaned.Write()
hETau_tauPt_eCleaned.Write()
hETau_ePt_eCleaned.Write()
hETau_M_eCleaned_genMatched.Write()
hETau_dR_eCleaned_genMatched.Write()
hETau_tauPt_eCleaned_genMatched.Write()
hETau_ePt_eCleaned_genMatched.Write()

hMuTau_M_muCleaned.Write()
hMuTau_dR_muCleaned.Write()
hMuTau_tauPt_muCleaned.Write()
hMuTau_muPt_muCleaned.Write()
hMuTau_M_muCleaned_genMatched.Write()
hMuTau_dR_muCleaned_genMatched.Write()
hMuTau_tauPt_muCleaned_genMatched.Write()
hMuTau_muPt_muCleaned_genMatched.Write()

hNEvent.Write()
