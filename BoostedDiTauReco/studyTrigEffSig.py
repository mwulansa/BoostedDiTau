import ROOT,sys, os
from DataFormats.FWLite import Events, Handle
import numpy as np
from looseElectron import *

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()
ROOT.TH1.SetDefaultSumw2()

if len(sys.argv)<2:
    fin = 'filelist.txt'
else:
    fin = sys.argv[1]


#prefix = "root://xrootd.unl.edu/"

fileList = open(fin, 'r')

out=ROOT.TFile('h_trigStudy_'+fin.split('/')[-1].replace('.txt','.root'),'recreate')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handleGenMet = Handle ('vector<reco::GenMET>')
labelGenMet = ('genMetTrue')

handleGenJet = Handle ('vector<reco::GenJet>')
labelGenJet = ('slimmedGenJets')

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults', '', 'HLT')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('slimmedJets')

handleTaus = Handle ('vector<pat::Tau>')
labelTaus = ('slimmedTaus')

handleBoostedTau = Handle ('vector<pat::Tau>')
labelBoostedTau = ('slimmedTausBoosted')

handleMet = Handle ('vector<pat::MET>')
labelMet = ('slimmedMETs')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

h={}

# book histograms

histname = 'MuTau_M_SingleMu'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'MuTau_M_Tau'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'MuTau_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'ETau_M_SingleE'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'ETau_M_Tau'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'ETau_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'EMu_M_SingleE'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_M_SingleMu'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_M_MuonEG'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'TauTau_M_Tau'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'TauTau_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

h['NEvent'] = ROOT.TH1F ("NEvent","Number of Events;;N{event}",2,0,2)
h['JetPt'] = ROOT.TH1F ("JetPt",";p_{T};a.u.",1000,0,1000)

fcount = 0
for fname in fileList.readlines():
    #if fcount > 50: break
    fcount += 1
    #f=prefix+fname.replace('/eos/uscms', '').replace('\n', '')
    f=fname.replace('\n', '')
    print f
    events = Events(f)

    nevt = 0
    for event in events:
        nevt += 1
        
        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()

        event.getByLabel(labelGenMet, handleGenMet)
        genMet=handleGenMet.product()

        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()

        event.getByLabel(labelMet, handleMet)
        met=handleMet.product().front()

        event.getByLabel(labelJet, handleJet)
        jets=handleJet.product()

        event.getByLabel(labelMuon, handleMuon)
        muons=handleMuon.product()

        event.getByLabel(labelElectron, handleElectron)
        electrons=handleElectron.product()

        event.getByLabel(labelBoostedTau, handleBoostedTau)
        btaus = handleBoostedTau.product()

        event.getByLabel(labelTaus, handleTaus)
        taus = handleTaus.product()

        event.getByLabel(labelRho, handleRho)
        rho=handleRho.product()[0]

        h['NEvent'].Fill(0.5)
        h['NEvent'].Fill(1.5, genweight)

        event.getByLabel(labelHLT, handleHLT)
        triggerResults=handleHLT.product()
        names = event.object().triggerNames(triggerResults)

        if fcount == 1 and nevt == 1:
            for i in range(len(names)):
                print i, names.triggerName(i), triggerResults.accept(i)

        # gen information
        genMuons=[]
        genTauP=[]
        genTauM=[]
        genElectrons=[]
        genNutau=[]
        genNutauAnti=[]

        for particle in particles:
            if particle.pdgId()==15 and particle.mother().pdgId()==9999: 
                genTauM+=[particle]
            if particle.pdgId()==-15 and particle.mother().pdgId()==9999: 
                genTauP+=[particle]
            if abs(particle.pdgId())==11 and particle.isDirectHardProcessTauDecayProductFinalState(): 
                genElectrons+=[particle]
            if abs(particle.pdgId())==13 and particle.isDirectHardProcessTauDecayProductFinalState(): 
                genMuons+=[particle]
            if particle.pdgId()==16 and particle.isDirectHardProcessTauDecayProductFinalState():
                genNutau+=[particle]
            if particle.pdgId()==-16 and particle.isDirectHardProcessTauDecayProductFinalState():
                genNutauAnti+=[particle]

        event.getByLabel(labelGenJet, handleGenJet)
        genJets=[]
        for jet in handleGenJet.product():
            genJets+=[jet]
        genJets.sort(key=lambda x: x.pt(), reverse=True)

        selected_electrons=[]
        selected_isoelectrons=[]
        for electron in electrons:
            if electron.pt()<7: continue 
            if abs(electron.eta())>2.5: continue
            E_c = electron.superCluster().energy()
            if electron.isEB():
                eID = electron.full5x5_sigmaIetaIeta()<0.0112 \
                    and GsfEleEInverseMinusPInverse(electron)<0.193 \
                    and abs(dEtaInSeed(electron))<0.00377     \
                    and GsfEleMissingHitsCut(electron)<=1 \
                    and electron.passConversionVeto() \
                    and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.0884 \
                    and electron.hadronicOverEm()< (0.05 + 1.16/E_c + 0.0324*rho/E_c)
                if eID:
                    selected_electrons+=[electron]
                if eID and GsfEleEffAreaPFIsoCut(electron, rho) < 0.112+0.506/electron.pt():
                    selected_isoelectrons+=[electron]
            if electron.isEE():
                eID = electron.full5x5_sigmaIetaIeta() < 0.0425 \
                    and GsfEleEInverseMinusPInverse(electron) < 0.111 \
                    and abs(dEtaInSeed(electron)) < 0.00674 \
                    and GsfEleMissingHitsCut(electron)<=1 \
                    and electron.passConversionVeto() \
                    and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.169 \
                    and electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c)
                if eID:
                    selected_electrons+=[electron]
                if eID and GsfEleEffAreaPFIsoCut(electron, rho) < 0.108+0.963/electron.pt():
                    selected_isoelectrons+=[electron]
        selected_electrons.sort(key=lambda x: x.pt(), reverse=True)
        selected_isoelectrons.sort(key=lambda x: x.pt(), reverse=True)


        selected_muons=[]
        selected_isomuons=[]
        for muon in muons:
            if not muon.isLooseMuon(): continue
            if muon.pt()<3 or muon.eta()>2.4: continue
            selected_muons+=[muon]
            if muonIsoCut(muon)<0.25:
                selected_isomuons+=[muon]
        selected_muons.sort(key=lambda x: x.pt(), reverse=True)
        selected_isomuons.sort(key=lambda x: x.pt(), reverse=True)


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
            if MUF > 0.8: continue
            if CEMF > 0.9: continue
            if (NHF<0.90 and NEMF<0.90 and NumConst>1) and ((abs(jet.eta())<=2.4 and CHF>0 and CHM>0 and CEMF<0.99) or abs(jet.eta())>2.4) and abs(jet.eta())<=2.7:
                selected_jets+=[jet]
        selected_jets.sort(key=lambda x: x.pt(), reverse=True)

        if len(selected_jets)>0:
            h['JetPt'].Fill(selected_jets[0].pt(), genweight)

        selected_btaus=[]
        for tau in btaus:
            if tau.pt()<20 or abs(tau.eta())>2.3: continue
            if not tau.tauID("decayModeFinding"): continue
            if tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") < -0.5: continue
            if not tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"): continue
            selected_btaus+=[tau]
        selected_btaus.sort(key=lambda x: x.pt(), reverse=True)

        selected_taus=[]
        for tau in taus:
            if tau.pt()<20 or abs(tau.eta())>2.3: continue
            if not tau.tauID("decayModeFinding"): continue
            if tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") < -0.5: continue
            if not tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"): continue
            selected_taus+=[tau]
        selected_taus.sort(key=lambda x: x.pt(), reverse=True)


        # tau_e tau_mu
        if len(genElectrons) == 1 and len(genMuons) == 1 and met.pt()>100 and len(selected_isoelectrons) > 0 and len(selected_isomuons) > 0 and len(selected_jets) > 0:

            e = ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_isoelectrons[0].pt(), selected_isoelectrons[0].eta(), selected_isoelectrons[0].phi(), selected_isoelectrons[0].mass())

            mu = ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_isomuons[0].pt(), selected_isomuons[0].eta(), selected_isomuons[0].phi(), selected_isomuons[0].mass())

            if e.DeltaR(mu) < 0.4 and selected_isoelectrons[0].charge()*selected_isomuons[0].charge() < 0:
                if selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                    h['EMu_M_JetHT'].Fill((e+mu).M(), genweight)

                if (selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))):
                    h['EMu_M_SingleE'].Fill((e+mu).M(), genweight)

                if (selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                   or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                   or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))):
                    h['EMu_M_SingleMu'].Fill((e+mu).M(), genweight)

                if (selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                   or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                   or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                   or (mu.Pt()>8 and e.Pt()>23 and triggerResults.accept(names.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10"))) \
                   or (mu.Pt()>23 and e.Pt()>12 and triggerResults.accept(names.triggerIndex("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12"))):
                    h['EMu_M_MuonEG'].Fill((e+mu).M(), genweight)
                
        
        # tau_mu tau_h
        if len(genMuons) == 1 and len(genElectrons) == 0 and met.pt()>100 and len(selected_muons) > 0 and len(selected_taus) > 0 and len(selected_jets) > 0:
            genMu=ROOT.TLorentzVector()
            genMu.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

            if genMuons[0].pdgId() > 0:
                gentau = genTauP[0]
                gennu = genNutauAnti[0]
            else:
                gentau = genTauM[0]
                gennu = genNutau[0]
            genTau = ROOT.TLorentzVector()
            genTau.SetPtEtaPhiM(gentau.pt(), gentau.eta(), gentau.phi(), gentau.mass())
            genNu=ROOT.TLorentzVector()
            genNu.SetPtEtaPhiM(gennu.pt(), gennu.eta(), gennu.phi(), gennu.mass())

            mu = ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())

            tau = ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass())

            if mu.DeltaR(tau) < 0.4 and selected_muons[0].charge()*selected_taus[0].charge() < 0:
                if selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                    h['MuTau_M_JetHT'].Fill((mu+tau).M(), genweight)

                if (selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                   or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                   or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))):
                    h['MuTau_M_SingleMu'].Fill((mu+tau).M(), genweight)

                if (selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                   or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                   or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))) \
                   or (tau.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                   or (tau.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                    h['MuTau_M_Tau'].Fill((mu+tau).M(), genweight)

        # tau_e tau_h
        if len(genElectrons) == 1 and len(genMuons) == 0 and met.pt()>100 and len(selected_electrons) > 0 and len(selected_taus) > 0 and len(selected_jets) > 0:

            e = ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())

            tau = ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass())

            if e.DeltaR(tau) < 0.4 and selected_electrons[0].charge()*selected_taus[0].charge() < 0:
                if selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                    h['ETau_M_JetHT'].Fill((e+tau).M(), genweight)

                if (selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                   or (e.Pt()>24 and tau.Pt()>30 and triggerResults.accept(names.triggerIndex("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v9"))):
                    h['ETau_M_SingleE'].Fill((e+tau).M(), genweight)

                if (selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                   or (e.Pt()>24 and tau.Pt()>30 and triggerResults.accept(names.triggerIndex("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v9"))) \
                   or (tau.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                   or (tau.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                    h['ETau_M_Tau'].Fill((e+tau).M(), genweight)

        # tau_h tau_h
        if len(genMuons) == 0 and len(genElectrons) == 0 and len(genJets) > 0 and met.pt()>100 and len(selected_btaus) > 1 and len(selected_jets) > 0:
        
            tau1 = ROOT.TLorentzVector()
            tau1.SetPtEtaPhiM(selected_btaus[0].pt(), selected_btaus[0].eta(), selected_btaus[0].phi(), selected_btaus[0].mass())

            tau2 = ROOT.TLorentzVector()
            tau2.SetPtEtaPhiM(selected_btaus[1].pt(), selected_btaus[1].eta(), selected_btaus[1].phi(), selected_btaus[1].mass())

            if tau1.DeltaR(tau2) < 0.4 and selected_btaus[0].charge()*selected_btaus[1].charge() < 0:
                if selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                    h['TauTau_M_JetHT'].Fill((tau1+tau2).M(), genweight)

                if (selected_jets[0].pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (tau2.Pt()>40 and triggerResults.accept(names.triggerIndex("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v8"))) \
                   or (tau2.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v8"))) \
                   or (tau1.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                   or (tau1.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                    h['TauTau_M_Tau'].Fill((tau1+tau2).M(), genweight)

out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
    
        

    
