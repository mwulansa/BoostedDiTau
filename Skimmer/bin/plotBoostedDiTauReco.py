import ROOT, sys, math, gc
from DataFormats.FWLite import Events, Handle
import numpy as np

from looseElectron import *

# book histograms
hETau_M = ROOT.TH1F ("hETau_M", "e - #tau mass;M_{e#tau};N_{events}", 100, 0, 100)
hMuTau_M = ROOT.TH1F ("hMuTau_M", "mu - #tau mass;M_{#mu#tau};N_{events}", 100, 0, 100)

h_tauPt = ROOT.TH1F ("h_tauPt", "", 500, 0, 500)
h_etauPt = ROOT.TH1F ("h_etauPt", "", 500, 0, 500)
h_mtauPt = ROOT.TH1F ("h_mtauPt", "", 500, 0, 500)

h_isoEPt = ROOT.TH1F ("h_isoEPt", "", 500, 0, 500)

hETau_M_eCleaned = ROOT.TH1F ("hETau_ECleaned_M", "e - #tau mass;M_{e#tau};N_{events}", 100, 0, 100)
hETau_M_eCleaned_genMatched = ROOT.TH1F ("hETau_ECleaned_M_genMatched", "e - #tau mass;M_{e#tau};N_{events}", 100, 0, 100)

hMuTau_M_muCleaned = ROOT.TH1F ("hMuTau_MuCleaned_M", "#mu - #tau mass;M_{#mu#tau};N_{events}", 100, 0, 100)
hMuTau_M_muCleaned_genMatched = ROOT.TH1F ("hMuTau_MuCleaned_M_genMatched", "#mu - #tau mass;M_{#mu#tau};N_{events}", 100, 0, 100)

hETau_dR_eCleaned = ROOT.TH1F ("hETau_ECleaned_dR", "e - #Delta R;M_{e#tau};N_{events}", 20, 0, 1)
hMuTau_dR_muCleaned = ROOT.TH1F ("hMuTau_MuCleaned_dR", "#Delta R;M_{#mu#tau};N_{events}", 20, 0, 1)

hETau_dR_eCleaned_genMatched = ROOT.TH1F ("hETau_ECleaned_dR_genMatched", "e - #Delta R;M_{e#tau};N_{events}", 20, 0, 1)
hMuTau_dR_muCleaned_genMatched = ROOT.TH1F ("hMuTau_MuCleaned_dR_genMatched", "#Delta R;M_{#mu#tau};N_{events}", 20, 0, 1)

hETau_tauPt_eCleaned = ROOT.TH1F ("hETau_ECleaned_TauPt", "Tau Pt;P_{t,#tau};N_{events}", 500, 0, 500)
hMuTau_tauPt_muCleaned = ROOT.TH1F ("hMuTau_MuCleaned_TauPt", "Tau Pt;P_{t,#tau};N_{events}", 500, 0, 500)

hETau_ePt_eCleaned = ROOT.TH1F ("hETau_ECleaned_EPt", "Electron Pt;P_{t,e};N_{events}", 500, 0, 500)
hMuTau_muPt_muCleaned = ROOT.TH1F ("hMuTau_MuCleaned_MuPt", "Muon Pt;P_{t,#mu};N_{events}", 500, 0, 500)

hETau_tauPt_eCleaned_genMatched = ROOT.TH1F ("hETau_ECleaned_TauPt_genMatched", "Tau Pt;P_{t,#tau};N_{events}", 500, 0, 500)
hMuTau_tauPt_muCleaned_genMatched = ROOT.TH1F ("hMuTau_MuCleaned_TauPt_genMatched", "Tau Pt;P_{t,#tau};N_{events}", 500, 0, 500)

hETau_ePt_eCleaned_genMatched = ROOT.TH1F ("hETau_ECleaned_EPt_genMatched", "Electron Pt;P_{t,e};N_{events}", 500, 0, 500)
hMuTau_muPt_muCleaned_genMatched = ROOT.TH1F ("hMuTau_MuCleaned_MuPt_genMatched", "Muon Pt;P_{t,#mu};N_{events}", 500, 0, 500)

hNEvent = ROOT.TH1F ("hNEvent","Number of Events;;N{event}",1,0,2)
hNElectrons = ROOT.TH1F ("hNElectrons","Number of Events;;N{event}",4,0,4)
hNIsoElectrons = ROOT.TH1F ("hNIsoElectrons","",4,0,4)

hNMuons = ROOT.TH1F ("hNMuons","",4,0,4)
hNTaus = ROOT.TH1F ("hNTaus","",4,0,4)
hNeTaus = ROOT.TH1F ("hNeTaus","",4,0,4)
hNmTaus = ROOT.TH1F ("hNmTaus","",4,0,4)


handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('slimmedJets')

handleTaus = Handle ('vector<pat::Tau>')
labelTaus = ('slimmedTaus')

handleElectronCleanedTaus = Handle ('vector<pat::Tau>')
labelElectronCleanedTaus = ('slimmedTausElectronCleaned', '', 'PAT')

handleMuonCleanedTaus = Handle ('vector<pat::Tau>')
labelMuonCleanedTaus = ('slimmedTausMuonCleaned', '', 'PAT')

handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlineSlimmedPrimaryVertices')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('genParticles')

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults')

#handleBs = Handle ('reco::BeamSpot')
#labelBs = ("offlineBeamSpot")

#handleConv = Handle ('vector<reco::Conversion>')
#labelConv = ('reducedEgamma', 'reducedConversions')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')


#prefix=""
prefix="root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer17DR94Premix/"

mass="50"

out=ROOT.TFile("h_plotBoostedDiTauReco_m"+mass+"_v4_94X_3.root",'recreate')

#jobs=np.linspace(100, 61, 40)
#jobs=np.linspace(50, 1, 50)
jobs=np.linspace(60, 51, 10)
print jobs
#sys.exit()
#jobs=[16]

#out=ROOT.TFile(filename.replace("MINIAODSIM","hists"),'recreate')

for job in jobs:
    filename = "ALP_m"+mass+"_w1_htjmin400_RunIISummer17DR94Premix_MINIAODSIM_Cleaned_"+str(int(job))+".root"
    #filename = "ALP_m"+mass+"_w1_htjmin400_RunIISummer17DR94Premix_MINIAODSIM_test.root"
    #filename = "mumutautau_MultipleID_1.root"
    print(prefix+filename)

    #if int(job)==16: continue
    
    events=Events(prefix+filename)

    print events, sys.getsizeof(events)
    
    ntot=0
    
    for event in events:
    
        ntot+=1
        #if not ntot==879: continue
        #print "Begin processing Event", ntot
        hNEvent.Fill(1)
    
        event.getByLabel(labelJet, handleJet)
        jets=handleJet.product()
            
        event.getByLabel(labelVertex, handleVertex)
        vertex=handleVertex.product()
        pv = vertex[0].position()
    
        #slimmedMuons=None
        event.getByLabel(labelMuon, handleMuon)
        muons=handleMuon.product()
        #print "muons", muons
    
        hNMuons.Fill(0.5, muons.size())
    
        event.getByLabel(labelTaus, handleTaus)
        taus=handleTaus.product()
        #print ("taus", taus)
    
        hNTaus.Fill(0.5, taus.size())
    
        event.getByLabel(labelElectronCleanedTaus, handleElectronCleanedTaus)
        etaus=handleElectronCleanedTaus.product()
    
        hNeTaus.Fill(0.5, etaus.size())
    
        event.getByLabel(labelMuonCleanedTaus, handleMuonCleanedTaus)
        mtaus=handleMuonCleanedTaus.product()
    
        hNmTaus.Fill(0.5, mtaus.size())
    
        event.getByLabel(labelElectron, handleElectron)
        electrons=handleElectron.product()
        #print "electrons", electrons
    
        hNElectrons.Fill(0.5, electrons.size())
        hNIsoElectrons.Fill(0.5, electrons.size())
    
        #event.getByLabel(labelBs, handleBs)
        #bs=handleBs.product()
        #print ("beamspots", bs)
        
        #event.getByLabel(labelConv, handleConv)
        #convs=handleConv.product()
    
        #print "convs:", convs
    
        event.getByLabel(labelRho, handleRho)
        if len(handleRho.product())>0:
            rho = handleRho.product()[0]
        else:
            rho = 0
        #print "rho", rho
    
        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()
    
        genMuons=[]
        genTaus=[]
        genElectrons=[]
        genNutaus=[]
        for particle in particles: 
            if abs(particle.pdgId())==15 and particle.isHardProcess(): 
                genTaus+=[particle]     
            if abs(particle.pdgId())==11 and particle.isDirectHardProcessTauDecayProductFinalState(): 
                genElectrons+=[particle]
            if abs(particle.pdgId())==13 and particle.isDirectHardProcessTauDecayProductFinalState(): 
                genMuons+=[particle]
            if abs(particle.pdgId())==16 and particle.isDirectHardProcessTauDecayProductFinalState(): 
                genNutaus+=[particle]
    
    
        #Muon Selection
        selected_muons=[]
        if muons.size()>0:
            for muon in muons:
                #print "muon", muon
                if muon.pt()>3 and abs(muon.eta())<2.4:
                    hNMuons.Fill(1.5)
                    if muon.isLooseMuon():
                        hNMuons.Fill(2.5)
                        dxy = abs(muon.innerTrack().dxy(pv))
                        dz = abs(muon.innerTrack().dz(pv))
                        if dxy<0.2 and dz<0.5:
                            selected_muons+=[muon]
                            hNMuons.Fill(3.5)
        selected_muons.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_muons)>0:
            hMuTau_muPt_muCleaned.Fill(selected_muons[0].pt())
    
        #Electron Selection
        selected_electrons=[]
        if electrons.size()>0:
            for electron in electrons:
                #print("electron", electron)
                if electron.pt()>7 and abs(electron.eta())<2.5:
                    hNElectrons.Fill(1.5)
                    E_c = electron.superCluster().energy()
                    if electron.isEB():
                        if electron.full5x5_sigmaIetaIeta()<0.0112 \
                        and GsfEleEInverseMinusPInverse(electron)<0.193 \
                        and abs(dEtaInSeed(electron))<0.00377     \
                        and GsfEleMissingHitsCut(electron)<=1 \
                        and electron.passConversionVeto() \
                        and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.0884 \
                        and electron.hadronicOverEm()< (0.05 + 1.16/E_c + 0.0324*rho/E_c):
                            #and electron.hadronicOverEm() < 0.236:
                            hNElectrons.Fill(2.5)
                            if abs(electron.gsfTrack().dz(pv))<0.1 \
                            and abs(electron.gsfTrack().dxy(pv))<0.05:
                                hNElectrons.Fill(3.5)
                                selected_electrons+=[electron]
                                
                    if electron.isEE():
                        if electron.full5x5_sigmaIetaIeta() < 0.0425 \
                        and GsfEleEInverseMinusPInverse(electron) < 0.111 \
                        and abs(dEtaInSeed(electron)) < 0.00674 \
                        and GsfEleMissingHitsCut(electron)<=1 \
                        and electron.passConversionVeto() \
                        and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.169 \
                        and electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c):
                            #and electron.hadronicOverEm() < 0.0801:
                        #and GsfEleConversionVetoCut(electron, bs, convs):
                            hNElectrons.Fill(2.5)
                            if abs(electron.gsfTrack().dz(pv)) < 0.2 \
                            and abs(electron.gsfTrack().dxy(pv)) < 0.1:
                                hNElectrons.Fill(3.5)
                                selected_electrons+=[electron]
        selected_electrons.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_electrons) > 0:
            hETau_ePt_eCleaned.Fill(selected_electrons[0].pt())

        #Isolated electron Selection
        selected_electrons_iso=[]
        if electrons.size()>0:
            for electron in electrons:
                #print("electron", electron)
                if electron.pt()>7 and abs(electron.eta())<2.5:
                    hNIsoElectrons.Fill(1.5)
                    E_c = electron.superCluster().energy()
                    if electron.isEB():
                        if electron.full5x5_sigmaIetaIeta()<0.0112 \
                        and GsfEleEInverseMinusPInverse(electron)<0.193 \
                        and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.0884 \
                        and abs(dEtaInSeed(electron))<0.00377     \
                        and GsfEleMissingHitsCut(electron)<=1 \
                        and GsfEleEffAreaPFIsoCut(electron, rho) < 0.112+0.506/electron.pt() \
                        and electron.passConversionVeto() :
                        #
                        #and electron.hadronicOverEm()< (0.05 + 1.16/E_c + 0.0324*rho/E_c) \
                            hNIsoElectrons.Fill(2.5)
                            if abs(electron.gsfTrack().dz(pv))<0.1 \
                            and abs(electron.gsfTrack().dxy(pv))<0.05:
                                hNIsoElectrons.Fill(3.5)
                                selected_electrons_iso+=[electron]
                                
                    if electron.isEE():
                        if electron.full5x5_sigmaIetaIeta() < 0.0425 \
                        and GsfEleEInverseMinusPInverse(electron) < 0.111 \
                        and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.169 \
                        and abs(dEtaInSeed(electron)) < 0.00674 \
                        and GsfEleMissingHitsCut(electron)<=1 \
                        and GsfEleEffAreaPFIsoCut(electron, rho) < 0.108+0.963/electron.pt() \
                        and electron.passConversionVeto() :
                        #                     
                        #and electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c) \
                        #and GsfEleConversionVetoCut(electron, bs, convs):
                            hNIsoElectrons.Fill(2.5)
                            if abs(electron.gsfTrack().dz(pv)) < 0.2 \
                            and abs(electron.gsfTrack().dxy(pv)) < 0.1:
                                hNIsoElectrons.Fill(3.5)
                                selected_electrons_iso+=[electron]
        selected_electrons.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_electrons_iso) > 0:
            h_isoEPt.Fill(selected_electrons_iso[0].pt())
    
    
        #Jet Selection
        selected_jets=[]
        if jets.size() > 0:
            for jet in jets:
                if jet.pt()>20 and abs(jet.eta())<2.5: 
                    NHF  = jet.neutralHadronEnergyFraction()
                    NEMF = jet.neutralEmEnergyFraction()
                    CHF  = jet.chargedHadronEnergyFraction()
                    MUF  = jet.muonEnergyFraction()
                    CEMF = jet.chargedEmEnergyFraction()
                    NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity()
                    NumNeutralParticles =jet.neutralMultiplicity()
                    CHM      = jet.chargedMultiplicity()
                    if (NHF<0.9 and NEMF<0.9 and NumConst>1 and MUF<0.8) and ((abs(jet.eta())<=2.4 and CHF>0 and CHM>0 and CEMF<0.8) or abs(jet.eta())>2.4) and abs(jet.eta())<=2.7:
                        selected_jets+=[jet]     
        selected_jets.sort(key=lambda x: x.pt(), reverse=True)
    
        selected_taus=[]
        if taus.size()>0:
            for tau in taus:
                if tau.pt()>10 and abs(tau.eta())<2.3:
                    hNTaus.Fill(1.5)
                    if tau.tauID("decayModeFinding")>0.5 \
                    and tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") > -0.5 \
                    and tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
                        hNTaus.Fill(2.5)
                        selected_taus+=[tau]
                        dxy = abs(tau.leadChargedHadrCand().get().dxy(pv))
                        dz = abs(tau.leadChargedHadrCand().get().dz(pv))
                        if dxy <0.2 and dz < 0.5:
                            hNTaus.Fill(3.5)
        selected_taus.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_taus) > 0:
            h_tauPt.Fill(selected_taus[0].pt())
    
        selected_etaus=[]
        if etaus.size()>0:
            for tau in etaus:
                if tau.pt()>10 and abs(tau.eta())<2.3:
                    hNeTaus.Fill(1.5)
                    if tau.tauID("decayModeFinding")>0.5 \
                    and tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") > -0.5 \
                    and tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
                        hNeTaus.Fill(2.5)
                        selected_etaus+=[tau]
                        dxy = abs(tau.leadChargedHadrCand().get().dxy(pv))
                        dz = abs(tau.leadChargedHadrCand().get().dz(pv))
                        if dxy <0.2 and dz < 0.5:
                            hNeTaus.Fill(3.5)
        selected_etaus.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_etaus) > 0:
            h_etauPt.Fill(selected_etaus[0].pt())
            hETau_tauPt_eCleaned.Fill(selected_etaus[0].pt())
            
    
        selected_mtaus=[]
        if mtaus.size()>0:
            for tau in mtaus:
                if tau.pt()>10 and abs(tau.eta())<2.3:
                    hNmTaus.Fill(1.5)
                    if tau.tauID("decayModeFinding")>0.5 \
                    and tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") > -0.5 \
                    and tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
                        hNmTaus.Fill(2.5)
                        selected_mtaus+=[tau]
                        dxy = abs(tau.leadChargedHadrCand().get().dxy(pv))
                        dz = abs(tau.leadChargedHadrCand().get().dz(pv))
                        if dxy <0.2 and dz < 0.5:
                            hNmTaus.Fill(3.5)
        selected_mtaus.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_mtaus) > 0:
            h_mtauPt.Fill(selected_mtaus[0].pt())
            
    
    
        # tau_mu tau_had
        if len(selected_muons)>0 and len(selected_taus)>0 and len(selected_jets)>0: 
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())
             
            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass())
             
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
                hMuTau_M.Fill((mu+tau).M())
                hMuTau_tauPt_muCleaned.Fill(tau.Pt())
    
        # tau_e tau_had
        if len(selected_electrons)>0 and len(selected_taus)>0 and len(selected_jets)>0: 
            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())
            
            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass())
            
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
                hETau_M.Fill((e+tau).M())
                hETau_ePt_eCleaned.Fill(tau.Pt())
                
    
        # muon cleaned tau_mu tau_had
        if len(selected_muons)>0 and len(selected_mtaus)>0 and len(selected_jets)>0:
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())
            
            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_mtaus[0].pt(), selected_mtaus[0].eta(), selected_mtaus[0].phi(), selected_mtaus[0].mass()) 
            
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
    
                # gen match
                #evtMatched=False
                if len(genElectrons)==0 and len(genMuons)==1: 
                    #    evtMatched=True
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
                        
                    if mu.DeltaR(genMu)<0.2 and tau.DeltaR(genTau-genNt)<0.3:
                        hMuTau_M_muCleaned_genMatched.Fill((mu+tau).M())
                        hMuTau_dR_muCleaned_genMatched.Fill(tau.DeltaR(mu))
                        hMuTau_tauPt_muCleaned_genMatched.Fill(tau.Pt())
                        hMuTau_muPt_muCleaned_genMatched.Fill(mu.Pt())
    
                
    
        # electron cleaned tau_e tau_had
        if len(selected_electrons)>0 and len(selected_etaus)>0 and len(selected_jets)>0: 
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
                         
                     if e.DeltaR(genEle)<0.2 and tau.DeltaR(genTau-genNt)<0.3:
                         hETau_M_eCleaned_genMatched.Fill((e+tau).M())
                         hETau_dR_eCleaned_genMatched.Fill(tau.DeltaR(e))
                         hETau_tauPt_eCleaned_genMatched.Fill(tau.Pt())
                         hETau_ePt_eCleaned_genMatched.Fill(e.Pt())
        
    
        #print ("------------")

out.cd()
hETau_M.Write()
hMuTau_M.Write()

h_tauPt.Write()
h_etauPt.Write()
h_mtauPt.Write()

h_isoEPt.Write()

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
hNElectrons.Write()
hNIsoElectrons.Write()

hNMuons.Write()
hNTaus.Write()
hNeTaus.Write()
hNmTaus.Write()

out.Close()
