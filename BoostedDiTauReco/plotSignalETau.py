import ROOT,sys, os
from DataFormats.FWLite import Events, Handle
import numpy as np

from looseElectron import *

jobs=np.linspace(100, 1, 100)
#jobs=np.linspace(50, 1, 50)
#jobs=[100]

if len(sys.argv)>1:
    mass=sys.argv[1]
else:
    mass='10'
    print "Use default signal mass: 10 GeV."

kin = "htjmin400"
#kin = "ptjmin100"


prefixDir = "/eos/uscms/store/user/zhangj/events/ALP/RunIISummer19UL17RECO/"
prefix = "root://cmseos.fnal.gov/"+prefixDir

#prefixDir = ""
#prefix = ""

out=ROOT.TFile("h_plotSignalETau_m"+mass+"_"+kin+"_v2.root",'recreate')

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('slimmedJets')

handleGenJet = Handle ('vector<reco::GenJet>')
labelGenJet = ('slimmedGenJets')

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

handleMuonCleanedTaus = Handle ('vector<pat::Tau>')
labelMuonCleanedTaus = ('slimmedTausMuonCleaned')

handleTausNewID = Handle ('vector<pat::Tau>')
labelTausNewID = ('slimmedTausNewID')

handleTausNewIDMuonCleaned = Handle ('vector<pat::Tau>')
labelTausNewIDMuonCleaned = ('slimmedTausNewIDMuonCleaned')

handleTausNewIDElectronCleaned = Handle ('vector<pat::Tau>')
labelTausNewIDElectronCleaned = ('slimmedTausNewIDElectronCleaned')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

h={}

# book histograms

histname = 'ETau_M'
h[histname] = ROOT.TH1F (histname, ";M_{e#tau};N_{events}", 100, 0, 100)
histname = 'ETau_M_deep'
h[histname] = ROOT.TH1F (histname, ";M_{e#tau};N_{events}", 100, 0, 100)

histname = 'ETau_M_eCleaned'
h[histname] = ROOT.TH1F (histname, ";M_{e#tau};N_{events}", 100, 0, 100)
histname = 'ETau_M_eCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";M_{e#tau};N_{events}", 100, 0, 100)

histname = 'ETau_M_baseline_eCleaned'
h[histname] = ROOT.TH1F (histname, ";M_{e#tau};N_{events}", 100, 0, 100)
histname = 'ETau_M_baseline_eCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";M_{e#tau};N_{events}", 100, 0, 100)

histname = 'ETau_M_eCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";M_{e#tau};N_{events}", 100, 0, 100)
histname = 'ETau_M_eCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";M_{e#tau};N_{events}", 100, 0, 100)

histname = 'ETau_M_baseline_eCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";M_{e#tau};N_{events}", 100, 0, 100)
histname = 'ETau_M_baseline_eCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";M_{e#tau};N_{events}", 100, 0, 100)

histname = 'ETau_tauPt_eCleaned'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)
histname = 'ETau_tauPt_eCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)

histname = 'ETau_tauPt_baseline_eCleaned'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)
histname = 'ETau_tauPt_baseline_eCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)

histname = 'ETau_tauPt_eCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)
histname = 'ETau_tauPt_eCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)

histname = 'ETau_tauPt_baseline_eCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)
histname = 'ETau_tauPt_baseline_eCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)


histname = 'ETau_ePt_baseline_eCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,e};N_{events}", 500, 0, 1000)
histname = 'ETau_ePt_baseline_eCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,e};N_{events}", 500, 0, 1000)

histname = 'ETau_ePt_eCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,e};N_{events}", 500, 0, 1000)
histname = 'ETau_ePt_eCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,e};N_{events}", 500, 0, 1000)

histname = 'ETau_dR_baseline_eCleaned'
h[histname] = ROOT.TH1F (histname, ";dR;N_{events}", 20, 0, 1)
histname = 'ETau_dR_baseline_eCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";dR;N_{events}", 20, 0, 1)
histname = 'ETau_dR_baseline_eCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";dR;N_{events}", 20, 0, 1)
histname = 'ETau_dR_baseline_eCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";dR;N_{events}", 20, 0, 1)

histname = 'ETau_IDraw_baseline_eCleaned'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)
histname = 'ETau_IDraw_baseline_eCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)

histname = 'ETau_IDraw_baseline_eCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)
histname = 'ETau_IDraw_baseline_eCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)

histname = 'ETau_IDraw_eCleaned'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)
histname = 'ETau_IDraw_eCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)

histname = 'ETau_IDraw_baseline'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)
histname = 'ETau_IDraw_baseline_deep'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)

histname = 'ETau_IDraw_baseline_genMatched'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)
histname = 'ETau_IDraw_baseline_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)

#h['MuTau_muPt_mCleaned_genMatched'] = ROOT.TH1F ("hMuTau_MuCleaned_MuPt_genMatched", "Muon Pt;P_{t,#mu};N_{events}", 50, 0, 500)
#h['MuTau_muPt_mCleaned'] = ROOT.TH1F ("hMuTau_MuCleaned_MuPt", "Muon Pt;P_{t,#mu};N_{events}", 50, 0, 500)


h['NEvent'] = ROOT.TH1F ("NEvent","Number of Events;;N{event}",1,0,2)
h['NMu'] = ROOT.TH1F ("NMu","Number of Events;;N{event}",1,0,2)
h['NE'] = ROOT.TH1F ("NE","Number of Events;;N{event}",1,0,2)
h['NTau'] = ROOT.TH1F ("NTau","Number of Events;;N{event}",1,0,2)
h['NmTau'] = ROOT.TH1F ("NmTau","Number of Events;;N{event}",1,0,2)
h['NeTau'] = ROOT.TH1F ("NeTau","Number of Events;;N{event}",1,0,2)


h['JetPt'] = ROOT.TH1F ("JetPt", "Jet Pt;P_{t,#tau};N_{events}", 100, 0, 1000)
h['GenJetPt'] = ROOT.TH1F ("GenJetPt", "Jet Pt;P_{t,#tau};N_{events}", 100, 0, 1000)

# loop over events
for job in jobs:
    filename = "ALP_m"+mass+"_w1_"+kin+"_RunIISummer19UL17RECO_MINIAODSIM_Cleaned_"+str(int(job))+".root"

    #filename = "ALP_m10_RunIISummer19UL17RECO_MINIAODSIM_Cleaned.root"
    
    print prefix+filename

    if not os.path.isfile(prefixDir+filename): continue
    
    events=Events(prefix+filename)

    ntot=0 
    nmatch1=0
    nmatch2=0
    
    for event in events:
        ntot+=1

        h['NEvent'].Fill(1)

        print "Processing event", ntot, "..."
    
        event.getByLabel(labelMuon, handleMuon)
        muons=handleMuon.product()

        event.getByLabel(labelElectron, handleElectron)
        electrons=handleElectron.product()

        event.getByLabel(labelRho, handleRho)
        rho=handleRho.product()[0]

        event.getByLabel(labelJet, handleJet)
        jets=handleJet.product()

        event.getByLabel(labelGenJet, handleGenJet)
        genjets=handleGenJet.product()

        event.getByLabel(labelTausNewID, handleTausNewID)
        taus=handleTausNewID.product()

        #event.getByLabel(labelTausNewIDMuonCleaned, handleTausNewIDMuonCleaned)
        #mtaus=handleTausNewIDMuonCleaned.product()

        event.getByLabel(labelTausNewIDElectronCleaned, handleTausNewIDElectronCleaned)
        etaus=handleTausNewIDElectronCleaned.product()

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

        selected_genjets = []
        for jet in genjets:
            if jet.pt() > 20 and abs(jet.eta()) < 2.5:
                selected_genjets+=[jet]
        selected_genjets.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_genjets) > 0: h['GenJetPt'].Fill(selected_genjets[0].pt())

        #Muon selection
        selected_muons=[]
        for muon in muons:
            if not muon.isLooseMuon(): continue
            #if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue 
            if muon.pt()<3 or muon.eta()>2.4: continue
            selected_muons+=[muon]
        selected_muons.sort(key=lambda x: x.pt(), reverse=True)

        if len(selected_muons) > 0: h['NMu'].Fill(1)

        selected_electrons=[]
        for electron in electrons:
            if electron.pt()<7: continue 
            if abs(electron.eta())>2.5: continue
            E_c = electron.superCluster().energy()
            if electron.isEB():
                if electron.full5x5_sigmaIetaIeta()<0.0112 \
                   and GsfEleEInverseMinusPInverse(electron)<0.193 \
                   and abs(dEtaInSeed(electron))<0.00377     \
                   and GsfEleMissingHitsCut(electron)<=1 \
                   and electron.passConversionVeto() \
                   and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.0884 \
                   and electron.hadronicOverEm()< (0.05 + 1.16/E_c + 0.0324*rho/E_c): 
                    selected_electrons+=[electron]

            if electron.isEE():
                if electron.full5x5_sigmaIetaIeta() < 0.0425 \
                   and GsfEleEInverseMinusPInverse(electron) < 0.111 \
                   and abs(dEtaInSeed(electron)) < 0.00674 \
                   and GsfEleMissingHitsCut(electron)<=1 \
                   and electron.passConversionVeto() \
                   and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.169 \
                   and electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c):
                    selected_electrons+=[electron]

        selected_electrons.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_electrons) > 0: h['NE'].Fill(1)
    

        #jet selection
        selected_jets=[]
        for jet in jets:
            if jet.pt() > 20 and abs(jet.eta()) < 2.5:
                NHF  = jet.neutralHadronEnergyFraction()
                NEMF = jet.neutralEmEnergyFraction()
                CHF  = jet.chargedHadronEnergyFraction()
                MUF  = jet.muonEnergyFraction()
                CEMF = jet.chargedEmEnergyFraction()
                NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity()
                NumNeutralParticles =jet.neutralMultiplicity()
                CHM      = jet.chargedMultiplicity()
                #if (NHF<0.99 and NEMF<0.99 and NumConst>1) and ((abs(jet.eta())<=2.4 and CHF>0 and CHM>0 and CEMF<0.99) or abs(jet.eta())>2.4) and abs(jet.eta())<=2.7:
                if CEMF<0.8 and CHM>0 and CHF>0 and NumConst>1 and NEMF<0.9 and MUF <0.8 and NHF < 0.9:
                    selected_jets+=[jet]                    
        selected_jets.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_jets) > 0: h['JetPt'].Fill(selected_jets[0].pt())

        #tau selection
        selected_taus=[]
        selected_etaus=[]
        selected_taus_deep=[]
        selected_etaus_deep=[]

        #print taus.size(), etaus.size()
        baseline_taus = []
        for tau in taus:
            if tau.pt() > 10 and abs(tau.eta()) < 2.3:
                baseline_taus += [tau]
                if tau.tauID('byMediumIsolationMVArun2017v2DBoldDMwLT2017'):
                    selected_taus+=[tau]
                if tau.tauID('byMediumDeepTau2017v2p1VSjet'): #and tau.tauID('byMediumDeepTau2017v2p1VSe'):
                    selected_taus_deep+=[tau]

        if len(selected_taus) > 0: h['NTau'].Fill(1)

        baseline_etaus = []
        for tau in etaus:
            if tau.pt() > 10 and abs(tau.eta()) < 2.3:
                baseline_etaus += [tau]
                if tau.tauID('byMediumIsolationMVArun2017v2DBoldDMwLT2017'):
                    selected_etaus+=[tau]   
            if tau.pt() > 10 and abs(tau.eta()) < 2.3:
                if tau.tauID('byMediumDeepTau2017v2p1ElectronCleanedVSjet'): #and tau.tauID('byMediumDeepTau2017v2p1ElectronCleanedVSe'):
                    selected_etaus_deep+=[tau]

        if len(selected_etaus) > 0: h['NeTau'].Fill(1)
                
        selected_taus.sort(key=lambda x: x.pt(), reverse=True)
        selected_etaus.sort(key=lambda x: x.pt(), reverse=True)
        selected_taus_deep.sort(key=lambda x: x.pt(), reverse=True)
        selected_etaus_deep.sort(key=lambda x: x.pt(), reverse=True)


        baseline_etaus.sort(key=lambda x: x.pt(), reverse=True)
        if len(baseline_etaus) > 0:
            h['ETau_IDraw_baseline_eCleaned'].Fill(baseline_etaus[0].tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017'))
            h['ETau_IDraw_baseline_eCleaned_deep'].Fill(baseline_etaus[0].tauID('byDeepTau2017v2p1ElectronCleanedVSjetraw'))

            if len(genElectrons)==1 and len(genMuons)==0: 
                genE=ROOT.TLorentzVector()
                genE.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                genTau=ROOT.TLorentzVector()
                genNt=ROOT.TLorentzVector()
                if genElectrons[0].pdgId()==11:
                    genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                    genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass()) 
                else:
                    genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                    genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())

                for t in baseline_etaus:
                    tau=ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(t.pt(), t.eta(), t.phi(), t.mass())
                    dR0 = 9999
                    rawId = -1
                    rawId_deep = -1
                    if tau.DeltaR(genTau-genNt) < dR0 and tau.DeltaR(genTau-genNt) < 0.4:
                        dR0 = tau.DeltaR(genTau-genNt)
                        rawId = t.tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017')
                        rawId_deep = t.tauID('byDeepTau2017v2p1ElectronCleanedVSjetraw')
                        h['ETau_IDraw_baseline_eCleaned_genMatched'].Fill(rawId)
                        h['ETau_IDraw_baseline_eCleaned_deep_genMatched'].Fill(rawId_deep)


        baseline_taus.sort(key=lambda x: x.pt(), reverse=True)
        if len(baseline_taus) > 0:
            h['ETau_IDraw_baseline'].Fill(baseline_taus[0].tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017'))
            h['ETau_IDraw_baseline_deep'].Fill(baseline_taus[0].tauID('byDeepTau2017v2p1VSjetraw'))

            if len(genElectrons)==1 and len(genMuons)==0: 
                genE=ROOT.TLorentzVector()
                genE.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                genTau=ROOT.TLorentzVector()
                genNt=ROOT.TLorentzVector()
                if genElectrons[0].pdgId()==11:
                    genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                    genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass()) 
                else:
                    genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                    genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())

                for t in baseline_taus:
                    tau=ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(t.pt(), t.eta(), t.phi(), t.mass())
                    dR0 = 9999
                    rawId = -1
                    rawId_deep = -1
                    if tau.DeltaR(genTau-genNt) < dR0 and tau.DeltaR(genTau-genNt) < 0.4:
                        dR0 = tau.DeltaR(genTau-genNt)
                        rawId = t.tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017')
                        rawId_deep = t.tauID('byDeepTau2017v2p1VSjetraw')
                        h['ETau_IDraw_baseline_genMatched'].Fill(rawId)
                        h['ETau_IDraw_baseline_deep_genMatched'].Fill(rawId_deep)
                        
        # tau_e tau_had
        if len(selected_electrons)>0 and len(selected_taus)>0 and len(selected_jets)>0: 
            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            if j.Pt() > 500 and e.DeltaR(tau)<0.4 and e.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                h['ETau_M'].Fill((e+tau).M())

        if len(selected_electrons)>0 and len(selected_taus_deep)>0 and len(selected_jets)>0: 
            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus_deep[0].pt(), selected_taus_deep[0].eta(), selected_taus_deep[0].phi(), selected_taus_deep[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            if j.Pt() > 500 and e.DeltaR(tau)<0.4 and e.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                h['ETau_M_deep'].Fill((e+tau).M())

        if len(selected_electrons)>0 and len(selected_etaus)>0 and len(selected_jets)>0: 
            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_etaus[0].pt(), selected_etaus[0].eta(), selected_etaus[0].phi(), selected_etaus[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            h['ETau_M_baseline_eCleaned'].Fill((e+tau).M())
            h['ETau_dR_baseline_eCleaned'].Fill(e.DeltaR(tau))
            h['ETau_tauPt_baseline_eCleaned'].Fill(tau.Pt())

            if len(genElectrons)==1 and len(genMuons)==0: 
                genE=ROOT.TLorentzVector()
                genE.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                genTau=ROOT.TLorentzVector()
                genNt=ROOT.TLorentzVector()
                if genElectrons[0].pdgId()==11:
                    genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                    genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass()) 
                else:
                    genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                    genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                        
                if e.DeltaR(genE)<0.4 and tau.DeltaR(genTau-genNt)<0.4:
                    h['ETau_M_baseline_eCleaned_genMatched'].Fill((e+tau).M())
                    h['ETau_dR_baseline_eCleaned_genMatched'].Fill(e.DeltaR(tau))

                else:
                    print(genE.Pt(), genE.Eta(), genE.Phi(), e.Pt(), e.Eta(), e.Phi())
                    print(genTau.Pt(), genTau.Eta(), genTau.Phi(), tau.Pt(), tau.Eta(), tau.Phi())

                if tau.DeltaR(genTau-genNt)<0.4:
                    h['ETau_tauPt_baseline_eCleaned_genMatched'].Fill(tau.Pt())

                if e.DeltaR(genE)<0.4:
                    h['ETau_ePt_baseline_eCleaned_genMatched'].Fill(e.Pt())
                

            
            if j.Pt() > 500 and e.DeltaR(tau)<0.4 and e.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                h['ETau_M_eCleaned'].Fill((e+tau).M())
                h['ETau_tauPt_eCleaned'].Fill(tau.Pt())
                h['ETau_IDraw_eCleaned'].Fill(selected_etaus[0].tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017'))

                if len(genElectrons)==1 and len(genMuons)==0: 
                    genE=ROOT.TLorentzVector()
                    genE.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                    genTau=ROOT.TLorentzVector()
                    genNt=ROOT.TLorentzVector()
                    if genElectrons[0].pdgId()==11:
                        genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                        genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass()) 
                    else:
                        genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                        genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                        
                    if e.DeltaR(genE)<0.4 and tau.DeltaR(genTau-genNt)<0.4:
                        h['ETau_M_eCleaned_genMatched'].Fill((e+tau).M())
                        h['ETau_tauPt_eCleaned_genMatched'].Fill(tau.Pt())

                    if tau.DeltaR(genTau-genNt)<0.4:
                        h['ETau_tauPt_eCleaned_genMatched'].Fill(tau.Pt())

                    if e.DeltaR(genE)<0.4:
                        h['ETau_ePt_eCleaned_genMatched'].Fill(e.Pt())
                
        if len(selected_electrons)>0 and len(selected_etaus_deep)>0 and len(selected_jets)>0: 
            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_etaus_deep[0].pt(), selected_etaus_deep[0].eta(), selected_etaus_deep[0].phi(), selected_etaus_deep[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            h['ETau_M_baseline_eCleaned_deep'].Fill((e+tau).M())
            h['ETau_dR_baseline_eCleaned_deep'].Fill(tau.DeltaR(e))
            h['ETau_tauPt_baseline_eCleaned_deep'].Fill(tau.Pt())
            
            if len(genElectrons)==1 and len(genMuons)==0: 
                genE=ROOT.TLorentzVector()
                genE.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                genTau=ROOT.TLorentzVector()
                genNt=ROOT.TLorentzVector()
                if genElectrons[0].pdgId()==11:
                    genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                    genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass()) 
                else:
                    genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                    genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                        
                if e.DeltaR(genE)<0.4 and tau.DeltaR(genTau-genNt)<0.4:
                    h['ETau_M_baseline_eCleaned_deep_genMatched'].Fill((e+tau).M())
                    h['ETau_dR_baseline_eCleaned_deep_genMatched'].Fill(tau.DeltaR(e))
                    
                if tau.DeltaR(genTau-genNt)<0.4:
                    h['ETau_tauPt_baseline_eCleaned_deep_genMatched'].Fill(tau.Pt())
                if e.DeltaR(genE)<0.4:
                    h['ETau_ePt_baseline_eCleaned_deep_genMatched'].Fill(tau.Pt())
                    
                
                    

            if j.Pt() > 500 and e.DeltaR(tau)<0.4 and e.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                h['ETau_M_eCleaned_deep'].Fill((e+tau).M())
                h['ETau_tauPt_eCleaned_deep'].Fill(tau.Pt())
                h['ETau_IDraw_eCleaned_deep'].Fill(selected_etaus_deep[0].tauID('byDeepTau2017v2p1ElectronCleanedVSjetraw'))

                if len(genElectrons)==1 and len(genMuons)==0: 
                    genE=ROOT.TLorentzVector()
                    genE.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                    genTau=ROOT.TLorentzVector()
                    genNt=ROOT.TLorentzVector()
                    if genElectrons[0].pdgId()==11:
                        genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                        genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass()) 
                    else:
                        genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                        genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                        
                    if e.DeltaR(genE)<0.4 and tau.DeltaR(genTau-genNt)<0.4:
                        h['ETau_M_eCleaned_deep_genMatched'].Fill((e+tau).M())
                        h['ETau_tauPt_eCleaned_deep_genMatched'].Fill(tau.Pt())

                    if tau.DeltaR(genTau-genNt)<0.4:
                        h['ETau_tauPt_eCleaned_deep_genMatched'].Fill(tau.Pt())
                    if e.DeltaR(genE)<0.4:
                        h['ETau_ePt_eCleaned_deep_genMatched'].Fill(tau.Pt())

out.cd()
for key in sorted(h.keys()):
    h[key].Write()
