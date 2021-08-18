import ROOT,sys, os
from DataFormats.FWLite import Events, Handle
import numpy as np

#from looseElectron import *

jobs=np.linspace(100, 1, 100)
#jobs=[100]

if len(sys.argv)>1:
    mass=sys.argv[1]
else:
    mass='10'
    print "Use default signal mass: 10 GeV."


prefixDir = "/eos/uscms/store/user/zhangj/events/ALP/RunIISummer19UL17RECO/"
prefix = "root://cmseos.fnal.gov/"+prefixDir

#prefixDir = ""
#prefix = ""

out=ROOT.TFile("h_plotSignalMuTau_m"+mass+"_v2.root",'recreate')

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

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

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

h={}

# book histograms

histname = 'MuTau_M'
h[histname] = ROOT.TH1F (histname, ";M_{#mu#tau};N_{events}", 100, 0, 100)
histname = 'MuTau_M_deep'
h[histname] = ROOT.TH1F (histname, ";M_{#mu#tau};N_{events}", 100, 0, 100)

histname = 'MuTau_M_mCleaned'
h[histname] = ROOT.TH1F (histname, ";M_{#mu#tau};N_{events}", 100, 0, 100)
histname = 'MuTau_M_mCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";M_{#mu#tau};N_{events}", 100, 0, 100)

histname = 'MuTau_M_baseline_mCleaned'
h[histname] = ROOT.TH1F (histname, ";M_{#mu#tau};N_{events}", 100, 0, 100)
histname = 'MuTau_M_baseline_mCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";M_{#mu#tau};N_{events}", 100, 0, 100)

histname = 'MuTau_M_mCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";M_{#mu#tau};N_{events}", 100, 0, 100)
histname = 'MuTau_M_mCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";M_{#mu#tau};N_{events}", 100, 0, 100)

histname = 'MuTau_M_baseline_mCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";M_{#mu#tau};N_{events}", 100, 0, 100)
histname = 'MuTau_M_baseline_mCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";M_{#mu#tau};N_{events}", 100, 0, 100)

histname = 'MuTau_tauPt_mCleaned'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)
histname = 'MuTau_tauPt_mCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)

histname = 'MuTau_tauPt_baseline_mCleaned'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)
histname = 'MuTau_tauPt_baseline_mCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)

histname = 'MuTau_tauPt_mCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)
histname = 'MuTau_tauPt_mCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)

histname = 'MuTau_tauPt_baseline_mCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)
histname = 'MuTau_tauPt_baseline_mCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";P_{t,#tau};N_{events}", 500, 0, 1000)

histname = 'MuTau_dR_baseline_mCleaned'
h[histname] = ROOT.TH1F (histname, ";dR;N_{events}", 20, 0, 1)
histname = 'MuTau_dR_baseline_mCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";dR;N_{events}", 20, 0, 1)
histname = 'MuTau_dR_baseline_mCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";dR;N_{events}", 20, 0, 1)
histname = 'MuTau_dR_baseline_mCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";dR;N_{events}", 20, 0, 1)

histname = 'MuTau_IDraw_baseline_mCleaned'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)
histname = 'MuTau_IDraw_baseline_mCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)

histname = 'MuTau_IDraw_baseline_mCleaned_genMatched'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)
histname = 'MuTau_IDraw_baseline_mCleaned_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)

histname = 'MuTau_IDraw_mCleaned'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)
histname = 'MuTau_IDraw_mCleaned_deep'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)

histname = 'MuTau_IDraw_baseline'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)
histname = 'MuTau_IDraw_baseline_deep'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)

histname = 'MuTau_IDraw_baseline_genMatched'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)
histname = 'MuTau_IDraw_baseline_deep_genMatched'
h[histname] = ROOT.TH1F (histname, ";id;N_{events}", 50, -1, 1)




#h['MuTau_muPt_mCleaned_genMatched'] = ROOT.TH1F ("hMuTau_MuCleaned_MuPt_genMatched", "Muon Pt;P_{t,#mu};N_{events}", 50, 0, 500)
#h['MuTau_muPt_mCleaned'] = ROOT.TH1F ("hMuTau_MuCleaned_MuPt", "Muon Pt;P_{t,#mu};N_{events}", 50, 0, 500)


h['NEvent'] = ROOT.TH1F ("NEvent","Number of Events;;N{event}",1,0,2)
h['NMu'] = ROOT.TH1F ("NMu","Number of Events;;N{event}",1,0,2)
h['NTau'] = ROOT.TH1F ("NTau","Number of Events;;N{event}",1,0,2)
h['NmTau'] = ROOT.TH1F ("NmTau","Number of Events;;N{event}",1,0,2)


h['JetPt'] = ROOT.TH1F ("JetPt", "Jet Pt;P_{t,#tau};N_{events}", 100, 0, 1000)
h['GenJetPt'] = ROOT.TH1F ("GenJetPt", "Jet Pt;P_{t,#tau};N_{events}", 100, 0, 1000)

# loop over events
for job in jobs:
    filename = "ALP_m"+mass+"_w1_htjmin400_RunIISummer19UL17RECO_MINIAODSIM_Cleaned_"+str(int(job))+".root"

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
    
        event.getByLabel(labelMuon, handleMuon)
        muons=handleMuon.product()

        event.getByLabel(labelJet, handleJet)
        jets=handleJet.product()

        event.getByLabel(labelGenJet, handleGenJet)
        genjets=handleGenJet.product()

        event.getByLabel(labelTausNewID, handleTausNewID)
        taus=handleTausNewID.product()

        #event.getByLabel(labelElectronCleanedTaus, handleElectronCleanedTaus)
        #etaus=handleElectronCleanedTaus.product()

        event.getByLabel(labelTausNewIDMuonCleaned, handleTausNewIDMuonCleaned)
        mtaus=handleTausNewIDMuonCleaned.product()

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

        #print muons.size(), len(selected_muons)

        if len(selected_muons) > 0: h['NMu'].Fill(1)
    

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
        selected_mtaus=[]
        selected_taus_deep=[]
        selected_mtaus_deep=[]

        #print taus.size(), mtaus.size()

        baseline_taus=[]
        for tau in taus:
            if tau.pt() > 10 and abs(tau.eta()) < 2.3:
                baseline_taus+=[tau]
                if tau.tauID('byMediumIsolationMVArun2017v2DBoldDMwLT2017'):
                    selected_taus+=[tau]
                if tau.tauID('byMediumDeepTau2017v2p1VSjet'):
                    selected_taus_deep+=[tau]

        if len(selected_taus) > 0: h['NTau'].Fill(1)

        
        baseline_mtaus=[]
        for tau in mtaus:
            if tau.pt() > 10 and abs(tau.eta()) < 2.3:
                #print tau.tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017')
                baseline_mtaus+=[tau]
                if tau.tauID('byMediumIsolationMVArun2017v2DBoldDMwLT2017'):
                    selected_mtaus+=[tau]   
                if tau.tauID('byMediumDeepTau2017v2p1MuonCleanedVSjet'):# and tau.tauID('byMediumDeepTau2017v2p1MuonCleanedVSmu'):
                    selected_mtaus_deep+=[tau]

        if len(selected_mtaus) > 0: h['NmTau'].Fill(1)
                
        selected_taus.sort(key=lambda x: x.pt(), reverse=True)
        selected_mtaus.sort(key=lambda x: x.pt(), reverse=True)
        selected_taus_deep.sort(key=lambda x: x.pt(), reverse=True)
        selected_mtaus_deep.sort(key=lambda x: x.pt(), reverse=True)

        baseline_mtaus.sort(key=lambda x: x.pt(), reverse=True)
        if len(baseline_mtaus) > 0:
            h['MuTau_IDraw_baseline_mCleaned'].Fill(baseline_mtaus[0].tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017'))
            h['MuTau_IDraw_baseline_mCleaned_deep'].Fill(baseline_mtaus[0].tauID('byDeepTau2017v2p1MuonCleanedVSjetraw'))

            if len(genElectrons)==0 and len(genMuons)==1: 
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

                for t in baseline_mtaus:
                    tau=ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(t.pt(), t.eta(), t.phi(), t.mass())
                    dR0 = 9999
                    rawId = -1
                    rawId_deep = -1
                    if tau.DeltaR(genTau-genNt) < dR0 and tau.DeltaR(genTau-genNt) < 0.4:
                        dR0 = tau.DeltaR(genTau-genNt)
                        rawId = t.tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017')
                        rawId_deep = t.tauID('byDeepTau2017v2p1MuonCleanedVSjetraw')
                        h['MuTau_IDraw_baseline_mCleaned_genMatched'].Fill(rawId)
                        h['MuTau_IDraw_baseline_mCleaned_deep_genMatched'].Fill(rawId_deep)
                

        baseline_taus.sort(key=lambda x: x.pt(), reverse=True)
        if len(baseline_taus) > 0:
            h['MuTau_IDraw_baseline'].Fill(baseline_taus[0].tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017'))
            h['MuTau_IDraw_baseline_deep'].Fill(baseline_taus[0].tauID('byDeepTau2017v2p1VSjetraw'))

            if len(genElectrons)==0 and len(genMuons)==1: 
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
                        h['MuTau_IDraw_baseline_genMatched'].Fill(rawId)
                        h['MuTau_IDraw_baseline_deep_genMatched'].Fill(rawId_deep)
        
        # tau_mu tau_had
        if len(selected_muons)>0 and len(selected_taus)>0 and len(selected_jets)>0: 
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            if j.Pt() > 500 and mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                h['MuTau_M'].Fill((mu+tau).M())

        if len(selected_muons)>0 and len(selected_taus_deep)>0 and len(selected_jets)>0: 
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus_deep[0].pt(), selected_taus_deep[0].eta(), selected_taus_deep[0].phi(), selected_taus_deep[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            if j.Pt() > 500 and mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                h['MuTau_M_deep'].Fill((mu+tau).M())

        if len(selected_muons)>0 and len(selected_mtaus)>0 and len(selected_jets)>0: 
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_mtaus[0].pt(), selected_mtaus[0].eta(), selected_mtaus[0].phi(), selected_mtaus[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            h['MuTau_M_baseline_mCleaned'].Fill((mu+tau).M())
            h['MuTau_dR_baseline_mCleaned'].Fill(mu.DeltaR(tau))
            h['MuTau_tauPt_baseline_mCleaned'].Fill(tau.Pt())

            if len(genElectrons)==0 and len(genMuons)==1: 
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
                    h['MuTau_M_baseline_mCleaned_genMatched'].Fill((mu+tau).M())
                    h['MuTau_dR_baseline_mCleaned_genMatched'].Fill(mu.DeltaR(tau))
                    h['MuTau_tauPt_baseline_mCleaned_genMatched'].Fill(tau.Pt())

            
            if j.Pt() > 500 and mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                h['MuTau_M_mCleaned'].Fill((mu+tau).M())
                h['MuTau_tauPt_mCleaned'].Fill(tau.Pt())
                h['MuTau_IDraw_mCleaned'].Fill(selected_mtaus[0].tauID('byIsolationMVArun2017v2DBoldDMwLTraw2017'))

                if len(genElectrons)==0 and len(genMuons)==1: 
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
                        h['MuTau_M_mCleaned_genMatched'].Fill((mu+tau).M())
                        h['MuTau_tauPt_mCleaned_genMatched'].Fill(tau.Pt())
                
        if len(selected_muons)>0 and len(selected_mtaus_deep)>0 and len(selected_jets)>0: 
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())

            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_mtaus_deep[0].pt(), selected_mtaus_deep[0].eta(), selected_mtaus_deep[0].phi(), selected_mtaus_deep[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            h['MuTau_M_baseline_mCleaned_deep'].Fill((mu+tau).M())
            h['MuTau_tauPt_baseline_mCleaned_deep'].Fill(tau.Pt())
            h['MuTau_dR_baseline_mCleaned_deep'].Fill(tau.DeltaR(mu))

            if len(genElectrons)==0 and len(genMuons)==1: 
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
                    h['MuTau_M_baseline_mCleaned_deep_genMatched'].Fill((mu+tau).M())
                    h['MuTau_tauPt_baseline_mCleaned_deep_genMatched'].Fill(tau.Pt())
                    h['MuTau_dR_baseline_mCleaned_deep_genMatched'].Fill(tau.DeltaR(mu))
                    

            if j.Pt() > 500 and mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8: 
                h['MuTau_M_mCleaned_deep'].Fill((mu+tau).M())
                h['MuTau_tauPt_mCleaned_deep'].Fill(tau.Pt())
                h['MuTau_IDraw_mCleaned_deep'].Fill(selected_mtaus_deep[0].tauID('byDeepTau2017v2p1MuonCleanedVSjetraw'))

                if len(genElectrons)==0 and len(genMuons)==1: 
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
                        h['MuTau_M_mCleaned_deep_genMatched'].Fill((mu+tau).M())
                        h['MuTau_tauPt_mCleaned_deep_genMatched'].Fill(tau.Pt())

out.cd()
for key in sorted(h.keys()):
    h[key].Write()
