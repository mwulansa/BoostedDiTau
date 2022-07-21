import ROOT,sys, os
from DataFormats.FWLite import Events, Handle
import numpy as np
from looseElectron import *

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()
ROOT.TH1.SetDefaultSumw2()

fin = sys.argv[1]
whichSel = sys.argv[2]


def Selection1(l1s, l2s, js, isTauTau):
    j = ROOT.TLorentzVector()
    j.SetPtEtaPhiM(js[0].pt(), js[0].eta(), js[0].phi(), js[0].mass())
    
    if not isTauTau:
        l1 = ROOT.TLorentzVector()
        l1.SetPtEtaPhiM(l1s[0].pt(), l1s[0].eta(), l1s[0].phi(), l1s[0].mass())

        l2 = ROOT.TLorentzVector()
        l2.SetPtEtaPhiM(l2s[0].pt(), l2s[0].eta(), l2s[0].phi(), l2s[0].mass())

        if l1s[0].charge()*l2s[0].charge() < 0 and l1.DeltaR(l2) < 0.8 and l1.DeltaR(l2) > 0.05 and l1.DeltaR(j) > 0.8 and l2.DeltaR(j) > 0.8:
            return 0, 0
        else:
            return -9999, -9999
    else:
        l1 = ROOT.TLorentzVector()
        l1.SetPtEtaPhiM(l1s[0].pt(), l1s[0].eta(), l1s[0].phi(), l1s[0].mass())

        l2 = ROOT.TLorentzVector()
        l2.SetPtEtaPhiM(l2s[1].pt(), l2s[1].eta(), l2s[1].phi(), l2s[1].mass())

        if l1s[0].charge()*l2s[1].charge() < 0 and l1.DeltaR(l2) < 0.8 and l1.DeltaR(l2) > 0.05 and l1.DeltaR(j) > 0.8 and l2.DeltaR(j) > 0.8:
            return 0, 1
        else:
            return -9999, -9999

def Selection2(l1s, l2s, js, isTauTau):
    j = ROOT.TLorentzVector()
    j.SetPtEtaPhiM(js[0].pt(), js[0].eta(), js[0].phi(), js[0].mass())

    il1 = -9999
    il2 = -9999
    dR = 9999
    for l1 in l1s:
        for l2 in l2s:
            if isTauTau and l1s.index(l1) == l2s.index(l2): continue
            
            lep1 = ROOT.TLorentzVector()
            lep1.SetPtEtaPhiM(l1.pt(), l1.eta(), l1.phi(), l1.mass())

            lep2 = ROOT.TLorentzVector()
            lep2.SetPtEtaPhiM(l2.pt(), l2.eta(), l2.phi(), l2.mass())

            if l1.charge()*l2.charge() < 0 and lep1.DeltaR(lep2) < 0.8 and lep1.DeltaR(lep2) > 0.05 and lep1.DeltaR(j) > 0.8 and lep2.DeltaR(j) > 0.8:
                if lep1.DeltaR(lep2) < dR:
                    dR = lep1.DeltaR(lep2)
                    il1 = l1s.index(l1)
                    il2 = l2s.index(l2)

    return il1, il2




#prefix = "root://xrootd.unl.edu/"

fileList = open(fin, 'r')

if whichSel == '1':
    out=ROOT.TFile('h_trigStudy_Sel1_'+fin.split('/')[-1].replace('.txt','.root'),'recreate')
elif whichSel == '2':
    out=ROOT.TFile('h_trigStudy_Sel2_'+fin.split('/')[-1].replace('.txt','.root'),'recreate')
else:
    print "Missing argument!!!"
    sys.exit()

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
labelTaus = ('slimmedTausUnCleaned')

handleTausMuonCleaned = Handle ('vector<pat::Tau>')
labelTausMuonCleaned = ('slimmedTausMuonCleaned')

handleTausElectronCleaned = Handle ('vector<pat::Tau>')
labelTausElectronCleaned = ('slimmedTausElectronCleaned')

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
histname = 'MuTau_M_SingleMuonly'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'MuTau_Cleaned_M_SingleMu'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'MuTau_Cleaned_M_Tau'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'MuTau_Cleaned_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'MuTau_Cleaned_M_SingleMuonly'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'MuTau_Reco_Cleaned_M_SingleMu'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'MuTau_Reco_Cleaned_M_Tau'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'MuTau_Reco_Cleaned_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'MuTau_Reco_Cleaned_M_SingleMuonly'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'ETau_M_SingleE'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'ETau_M_Tau'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'ETau_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'ETau_Cleaned_M_SingleE'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'ETau_Cleaned_M_Tau'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'ETau_Cleaned_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'ETau_Reco_Cleaned_M_SingleE'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'ETau_Reco_Cleaned_M_Tau'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'ETau_Reco_Cleaned_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'EMu_M_SingleE'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_M_SingleMu'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_M_MuonEG'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_M_MuonEGonly'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_M_SingleMuonly'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_M_SingleEonly'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'EMu_Reco_M_SingleE'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_Reco_M_SingleMu'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_Reco_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_Reco_M_MuonEG'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_Reco_M_MuonEGonly'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_Reco_M_SingleMuonly'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'EMu_Reco_M_SingleEonly'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'TauTau_M_Tau'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'TauTau_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

histname = 'TauTau_Reco_M_Tau'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)
histname = 'TauTau_Reco_M_JetHT'
h[histname] = ROOT.TH1F (histname, ";M_{#tau#tau};a.u.", 100, 0, 100)

h['NEvent'] = ROOT.TH1F ("NEvent","Number of Events;;N{event}",2,0,2)
h['JetPt'] = ROOT.TH1F ("JetPt",";p_{T};a.u.",1000,0,1000)
h['genMtt'] = ROOT.TH1F ("genMtt",";m_{#tau#tau};a.u.",1000,0,100)

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

        event.getByLabel(labelTausMuonCleaned, handleTausMuonCleaned)
        mtaus = handleTausMuonCleaned.product()

        event.getByLabel(labelTausElectronCleaned, handleTausElectronCleaned)
        etaus = handleTausElectronCleaned.product()

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

        gentau1 = ROOT.TLorentzVector()
        gentau1.SetPtEtaPhiM(genTauM[0].pt(), genTauM[0].eta(), genTauM[0].phi(), genTauM[0].mass())

        gentau2 = ROOT.TLorentzVector()
        gentau2.SetPtEtaPhiM(genTauP[0].pt(), genTauP[0].eta(), genTauP[0].phi(), genTauP[0].mass())

        h['genMtt'].Fill((gentau1+gentau2).M(), genweight)
        

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

        selected_etaus=[]
        for tau in etaus:
            if tau.pt()<20 or abs(tau.eta())>2.3: continue
            if not tau.tauID("decayModeFinding"): continue
            if tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") < -0.5: continue
            if not tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"): continue
            selected_etaus+=[tau]
        selected_etaus.sort(key=lambda x: x.pt(), reverse=True)

        selected_mtaus=[]
        for tau in mtaus:
            if tau.pt()<20 or abs(tau.eta())>2.3: continue
            if not tau.tauID("decayModeFinding"): continue
            if tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") < -0.5: continue
            if not tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"): continue
            selected_mtaus+=[tau]
        selected_mtaus.sort(key=lambda x: x.pt(), reverse=True)


        # tau_e tau_mu
        if len(genElectrons) == 1 and len(genMuons) == 1 and met.pt()>100 and len(selected_isoelectrons) > 0 and len(selected_isomuons) > 0 and len(selected_jets) > 0:
            
            if whichSel == '1':
                il1, il2 = Selection1(selected_isoelectrons, selected_isomuons, selected_jets, False)
            if whichSel == '2':
                il1, il2 = Selection2(selected_isoelectrons, selected_isomuons, selected_jets, False)
                
            if il1 > -1 and il2 > -1:
                e = ROOT.TLorentzVector()
                e.SetPtEtaPhiM(selected_isoelectrons[il1].pt(), selected_isoelectrons[il1].eta(), selected_isoelectrons[il1].phi(), selected_isoelectrons[il1].mass())

                mu = ROOT.TLorentzVector()
                mu.SetPtEtaPhiM(selected_isomuons[il2].pt(), selected_isomuons[il2].eta(), selected_isomuons[il2].phi(), selected_isomuons[il2].mass())

                j = ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
            
                if j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                    h['EMu_M_JetHT'].Fill((e+mu).M(), genweight)

                if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))):
                    h['EMu_M_SingleE'].Fill((e+mu).M(), genweight)

                if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                   or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                   or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))):
                    h['EMu_M_SingleMu'].Fill((e+mu).M(), genweight)

                if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                   or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                   or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                   or (mu.Pt()>8 and e.Pt()>23 and triggerResults.accept(names.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10"))) \
                   or (mu.Pt()>23 and e.Pt()>12 and triggerResults.accept(names.triggerIndex("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12"))):
                    h['EMu_M_MuonEG'].Fill((e+mu).M(), genweight)

                if (mu.Pt()>8 and e.Pt()>23 and triggerResults.accept(names.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10"))) \
                   or (mu.Pt()>23 and e.Pt()>12 and triggerResults.accept(names.triggerIndex("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12"))):
                    h['EMu_M_MuonEGonly'].Fill((e+mu).M(), genweight)

                if (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                   or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))):
                    h['EMu_M_SingleMuonly'].Fill((e+mu).M(), genweight)

                if (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))):
                    h['EMu_M_SingleEonly'].Fill((e+mu).M(), genweight)

        if met.pt()>100 and len(selected_isoelectrons) > 0 and len(selected_isomuons) > 0 and len(selected_jets) > 0:

            if whichSel == '1':
                il1, il2 = Selection1(selected_isoelectrons, selected_isomuons, selected_jets, False)
            if whichSel == '2':
                il1, il2 = Selection2(selected_isoelectrons, selected_isomuons, selected_jets, False)
                
            if il1 > -1 and il2 > -1:
                e = ROOT.TLorentzVector()
                e.SetPtEtaPhiM(selected_isoelectrons[il1].pt(), selected_isoelectrons[il1].eta(), selected_isoelectrons[il1].phi(), selected_isoelectrons[il1].mass())

                mu = ROOT.TLorentzVector()
                mu.SetPtEtaPhiM(selected_isomuons[il2].pt(), selected_isomuons[il2].eta(), selected_isomuons[il2].phi(), selected_isomuons[il2].mass())

                j = ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())


                if j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                    h['EMu_Reco_M_JetHT'].Fill((e+mu).M(), genweight)

                if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))):
                    h['EMu_Reco_M_SingleE'].Fill((e+mu).M(), genweight)

                if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                   or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                   or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))):
                    h['EMu_Reco_M_SingleMu'].Fill((e+mu).M(), genweight)

                if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                   or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                   or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                   or (mu.Pt()>8 and e.Pt()>23 and triggerResults.accept(names.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10"))) \
                   or (mu.Pt()>23 and e.Pt()>12 and triggerResults.accept(names.triggerIndex("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12"))):
                    h['EMu_Reco_M_MuonEG'].Fill((e+mu).M(), genweight)

                if (mu.Pt()>8 and e.Pt()>23 and triggerResults.accept(names.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10"))) \
                   or (mu.Pt()>23 and e.Pt()>12 and triggerResults.accept(names.triggerIndex("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12"))):
                    h['EMu_Reco_M_MuonEGonly'].Fill((e+mu).M(), genweight)

                if (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                   or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))):
                    h['EMu_Reco_M_SingleMuonly'].Fill((e+mu).M(), genweight)

                if (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))):
                    h['EMu_Reco_M_SingleEonly'].Fill((e+mu).M(), genweight)
                
        
        # tau_mu tau_h
        if len(genMuons) == 1 and len(genElectrons) == 0 and met.pt()>100 and len(selected_jets) > 0:
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

            if len(selected_muons) > 0 and len(selected_taus) > 0:

                if whichSel == '1':
                    il1, il2 = Selection1(selected_muons, selected_taus, selected_jets, False)
                if whichSel == '2':
                    il1, il2 = Selection2(selected_muons, selected_taus, selected_jets, False)
                    
                if il1 > -1 and il2 > -1:
                    mu = ROOT.TLorentzVector()
                    mu.SetPtEtaPhiM(selected_muons[il1].pt(), selected_muons[il1].eta(), selected_muons[il1].phi(), selected_muons[il1].mass())
    
                    tau = ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(selected_taus[il2].pt(), selected_taus[il2].eta(), selected_taus[il2].phi(), selected_taus[il2].mass())
    
                    j = ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
    
                    if j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                        h['MuTau_M_JetHT'].Fill((mu+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                       or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                       or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))):
                        h['MuTau_M_SingleMu'].Fill((mu+tau).M(), genweight)
    
                    if (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                       or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                       or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))):
                        h['MuTau_M_SingleMuonly'].Fill((mu+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                       or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                       or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))) \
                       or (tau.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                       or (tau.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                        h['MuTau_M_Tau'].Fill((mu+tau).M(), genweight)

            if len(selected_muons) > 0 and len(selected_mtaus) > 0:

                if whichSel == '1':
                    il1, il2 = Selection1(selected_muons, selected_mtaus, selected_jets, False)
                if whichSel == '2':
                    il1, il2 = Selection2(selected_muons, selected_mtaus, selected_jets, False)
                    
                if il1 > -1 and il2 > -1:
                    mu = ROOT.TLorentzVector()
                    mu.SetPtEtaPhiM(selected_muons[il1].pt(), selected_muons[il1].eta(), selected_muons[il1].phi(), selected_muons[il1].mass())
    
                    tau = ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(selected_mtaus[il2].pt(), selected_mtaus[il2].eta(), selected_mtaus[il2].phi(), selected_mtaus[il2].mass())
    
                    j = ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
    
                    if j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                        h['MuTau_Cleaned_M_JetHT'].Fill((mu+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                       or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                       or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))):
                        h['MuTau_Cleaned_M_SingleMu'].Fill((mu+tau).M(), genweight)
    
                    if (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                       or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                       or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))):
                        h['MuTau_Cleaned_M_SingleMuonly'].Fill((mu+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                       or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                       or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))) \
                       or (tau.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                       or (tau.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                        h['MuTau_Cleaned_M_Tau'].Fill((mu+tau).M(), genweight)

        if met.pt()>100 and len(selected_jets) > 0:
            
            if len(selected_muons) > 0 and len(selected_mtaus) > 0:

                if whichSel == '1':
                    il1, il2 = Selection1(selected_muons, selected_mtaus, selected_jets, False)
                if whichSel == '2':
                    il1, il2 = Selection2(selected_muons, selected_mtaus, selected_jets, False)
                    
                if il1 > -1 and il2 > -1:
                    mu = ROOT.TLorentzVector()
                    mu.SetPtEtaPhiM(selected_muons[il1].pt(), selected_muons[il1].eta(), selected_muons[il1].phi(), selected_muons[il1].mass())
    
                    tau = ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(selected_mtaus[il2].pt(), selected_mtaus[il2].eta(), selected_mtaus[il2].phi(), selected_mtaus[il2].mass())
    
                    j = ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                    
                    if j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                        h['MuTau_Reco_Cleaned_M_JetHT'].Fill((mu+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                       or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                       or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))):
                        h['MuTau_Reco_Cleaned_M_SingleMu'].Fill((mu+tau).M(), genweight)
    
                    if (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                       or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                       or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))):
                        h['MuTau_Reco_Cleaned_M_SingleMuonly'].Fill((mu+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (mu.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_Mu50_v11"))) \
                       or (mu.Pt()>27 and triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13"))) \
                       or (mu.Pt()>24 and tau.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8"))) \
                       or (tau.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                       or (tau.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                        h['MuTau_Reco_Cleaned_M_Tau'].Fill((mu+tau).M(), genweight)

        # tau_e tau_h
        if len(genElectrons) == 1 and len(genMuons) == 0 and met.pt()>100 and len(selected_jets) > 0:
            
            if len(selected_electrons) > 0 and len(selected_taus) > 0:

                if whichSel == '1':
                    il1, il2 = Selection1(selected_electrons, selected_taus, selected_jets, False)
                if whichSel == '2':
                    il1, il2 = Selection2(selected_electrons, selected_taus, selected_jets, False)
                
                if il1 > -1 and il2 > -1:
                    e = ROOT.TLorentzVector()
                    e.SetPtEtaPhiM(selected_electrons[il1].pt(), selected_electrons[il1].eta(), selected_electrons[il1].phi(), selected_electrons[il1].mass())
    
                    tau = ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(selected_taus[il2].pt(), selected_taus[il2].eta(), selected_taus[il2].phi(), selected_taus[il2].mass())
    
                    j = ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
    
                    if j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                        h['ETau_M_JetHT'].Fill((e+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                       or (e.Pt()>24 and tau.Pt()>30 and triggerResults.accept(names.triggerIndex("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v9"))):
                        h['ETau_M_SingleE'].Fill((e+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                       or (e.Pt()>24 and tau.Pt()>30 and triggerResults.accept(names.triggerIndex("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v9"))) \
                       or (tau.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                       or (tau.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                        h['ETau_M_Tau'].Fill((e+tau).M(), genweight)

            if len(selected_electrons) > 0 and len(selected_etaus) > 0:

                if whichSel == '1':
                    il1, il2 = Selection1(selected_electrons, selected_etaus, selected_jets, False)
                if whichSel == '2':
                    il1, il2 = Selection2(selected_electrons, selected_etaus, selected_jets, False)
                    
                if il1 > -1 and il2 > -1:
                    e = ROOT.TLorentzVector()
                    e.SetPtEtaPhiM(selected_electrons[il1].pt(), selected_electrons[il1].eta(), selected_electrons[il1].phi(), selected_electrons[il1].mass())
    
                    tau = ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(selected_etaus[il2].pt(), selected_etaus[il2].eta(), selected_etaus[il2].phi(), selected_etaus[il2].mass())
    
                    j = ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                    
                    if j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                        h['ETau_Cleaned_M_JetHT'].Fill((e+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                       or (e.Pt()>24 and tau.Pt()>30 and triggerResults.accept(names.triggerIndex("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v9"))):
                        h['ETau_Cleaned_M_SingleE'].Fill((e+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                       or (e.Pt()>24 and tau.Pt()>30 and triggerResults.accept(names.triggerIndex("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v9"))) \
                       or (tau.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                       or (tau.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                        h['ETau_Cleaned_M_Tau'].Fill((e+tau).M(), genweight)

        if met.pt()>100 and len(selected_jets) > 0:

            if len(selected_electrons) > 0 and len(selected_etaus) > 0:

                if whichSel == '1':
                    il1, il2 = Selection1(selected_electrons, selected_etaus, selected_jets, False)
                if whichSel == '2':
                    il1, il2 = Selection2(selected_electrons, selected_etaus, selected_jets, False)
                    
                if il1 > -1 and il2 > -1:
                    e = ROOT.TLorentzVector()
                    e.SetPtEtaPhiM(selected_electrons[il1].pt(), selected_electrons[il1].eta(), selected_electrons[il1].phi(), selected_electrons[il1].mass())
    
                    tau = ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(selected_etaus[il2].pt(), selected_etaus[il2].eta(), selected_etaus[il2].phi(), selected_etaus[il2].mass())
    
                    j = ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                    
                    if j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                        h['ETau_Reco_Cleaned_M_JetHT'].Fill((e+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                       or (e.Pt()>24 and tau.Pt()>30 and triggerResults.accept(names.triggerIndex("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v9"))):
                        h['ETau_Reco_Cleaned_M_SingleE'].Fill((e+tau).M(), genweight)
    
                    if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                       or (e.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7"))) \
                       or (e.Pt()>24 and tau.Pt()>30 and triggerResults.accept(names.triggerIndex("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v9"))) \
                       or (tau.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                       or (tau.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                        h['ETau_Reco_Cleaned_M_Tau'].Fill((e+tau).M(), genweight)

        # tau_h tau_h
        if len(genMuons) == 0 and len(genElectrons) == 0 and len(genJets) > 0 and met.pt()>100 and len(selected_btaus) > 1 and len(selected_jets) > 0:
            if whichSel == '1':
                il1, il2 = Selection1(selected_btaus, selected_btaus, selected_jets, True)
            if whichSel == '2':
                il1, il2 = Selection2(selected_btaus, selected_btaus, selected_jets, True)
                
            if il1 > -1 and il2 > -1:
                tau1 = ROOT.TLorentzVector()
                tau1.SetPtEtaPhiM(selected_btaus[il1].pt(), selected_btaus[il1].eta(), selected_btaus[il1].phi(), selected_btaus[il1].mass())

                tau2 = ROOT.TLorentzVector()
                tau2.SetPtEtaPhiM(selected_btaus[il2].pt(), selected_btaus[il2].eta(), selected_btaus[il2].phi(), selected_btaus[il2].mass())

                j = ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

                if j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                    h['TauTau_M_JetHT'].Fill((tau1+tau2).M(), genweight)

                if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (tau2.Pt()>40 and triggerResults.accept(names.triggerIndex("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v8"))) \
                   or (tau2.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v8"))) \
                   or (tau1.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                   or (tau1.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                    h['TauTau_M_Tau'].Fill((tau1+tau2).M(), genweight)

        if len(genJets) > 0 and met.pt()>100 and len(selected_btaus) > 1 and len(selected_jets) > 0:
        
            if whichSel == '1':
                il1, il2 = Selection1(selected_btaus, selected_btaus, selected_jets, True)
            if whichSel == '2':
                il1, il2 = Selection2(selected_btaus, selected_btaus, selected_jets, True)
                
            if il1 > -1 and il2 > -1:
                tau1 = ROOT.TLorentzVector()
                tau1.SetPtEtaPhiM(selected_btaus[il1].pt(), selected_btaus[il1].eta(), selected_btaus[il1].phi(), selected_btaus[il1].mass())

                tau2 = ROOT.TLorentzVector()
                tau2.SetPtEtaPhiM(selected_btaus[il2].pt(), selected_btaus[il2].eta(), selected_btaus[il2].phi(), selected_btaus[il2].mass())

                j = ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                
                if j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")):
                    h['TauTau_Reco_M_JetHT'].Fill((tau1+tau2).M(), genweight)

                if (j.Pt() > 500 and triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17"))) \
                   or (tau2.Pt()>40 and triggerResults.accept(names.triggerIndex("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v8"))) \
                   or (tau2.Pt()>35 and triggerResults.accept(names.triggerIndex("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v8"))) \
                   or (tau1.Pt()>50 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"))) \
                   or (tau1.Pt()>180 and triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8"))):
                    h['TauTau_Reco_M_Tau'].Fill((tau1+tau2).M(), genweight)

out.cd()

for key in sorted(h.keys()):
    h[key].Write()

out.Close()
    
        

    
