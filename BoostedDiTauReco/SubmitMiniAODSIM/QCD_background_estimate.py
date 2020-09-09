import ROOT,sys,os
from DataFormats.FWLite import Events, Handle
from looseElectron import *
import numpy as np

sample = "QCD"

handleTausBoosted = Handle ('vector<pat::Tau>')
labelTausBoosted = ('slimmedTausBoosted', '', 'PAT')

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlineSlimmedPrimaryVertices')

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults','','HLT')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleBs = Handle ('reco::BeamSpot')
labelBs = ("offlineBeamSpot")

handleConv = Handle ('vector<reco::Conversion>')
labelConv = ('reducedEgamma', 'reducedConversions')

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('slimmedJets')

handlePatMETs = Handle("vector<pat::MET>")
labelPatMETs = ( 'slimmedMETs')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

prefix="root://cmseos.fnal.gov//eos/uscms/store/user/mwulansa/DIS/TCPAnalysis/Backgrounds/RunIIUL17/"

jobs=np.linspace(47, 0, 48)

outputFileName = "h_"+sample+"_backgroundEstimation.root"
print outputFileName

out=ROOT.TFile.Open(outputFileName,'recreate')

h = {}

h['hNEvent'] = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)

h['hTauTauBaseline_M'] = ROOT.TH1F ("hTauTau_Baseline_M", "#tau - #tau mass;M_{#tau#tau};", 200, 0, 100)
h['hTauTaudR_M'] = ROOT.TH1F ("hTauTau_dR_M", "#tau - #tau mass;M_{#tau#tau};", 200, 0, 100)
h['hTauTauTrig_M'] = ROOT.TH1F ("hTauTau_Trig_M", "#tau - #tau mass;M_{#tau#tau};", 200, 0, 100)
h['hTauTauMetcut_M'] = ROOT.TH1F ("hTauTau_Metcut_M", "#tau - #tau mass;M_{#tau#tau};", 200, 0, 100)
h['hTauTauInvMetcut_M'] = ROOT.TH1F ("hTauTau_InvMetcut_M", "#tau - #tau mass;M_{#tau#tau};", 200, 0, 100)

h['OS'] = ROOT.TH1F ("OS_M", "#tau - #tau mass;M_{#tau#tau};", 200, 0, 100)
h['SS'] = ROOT.TH1F ("SS_M", "#tau - #tau mass;M_{#tau#tau};", 200, 0, 100)
h['OS_invMET'] = ROOT.TH1F ("OS_invMET_M", "#tau - #tau mass;M_{#tau#tau};", 200, 0, 100)
h['SS_invMET'] = ROOT.TH1F ("SS_invMET_M", "#tau - #tau mass;M_{#tau#tau};", 200, 0, 100)

h['OS_lPt'] = ROOT.TH1F ("OS_lPt", ";p_{T};", 500, 0, 500)
h['OS_sPt'] = ROOT.TH1F ("OS_sPt", ";p_{T};", 500, 0, 500)

h['SS_lPt'] = ROOT.TH1F ("SS_lPt", ";p_{T};", 500, 0, 500)
h['SS_sPt'] = ROOT.TH1F ("SS_sPt", ";p_{T};", 500, 0, 500)

h['OS_invMET_lPt'] = ROOT.TH1F ("OS_invMET_lPt", ";p_{T};", 500, 0, 500)
h['OS_invMET_sPt'] = ROOT.TH1F ("OS_invMET_sPt", ";p_{T};", 500, 0, 500)

h['SS_invMET_lPt'] = ROOT.TH1F ("SS_invMET_lPt", ";p_{T};", 500, 0, 500)
h['SS_invMET_sPt'] = ROOT.TH1F ("SS_invMET_sPt", ";p_{T};", 500, 0, 500)

h['hMuPt'] = ROOT.TH1F ("hMu_Pt", ";p_{T};", 500, 0, 500)
h['hElePt'] = ROOT.TH1F ("hEle_Pt", ";p_{T};", 500, 0, 500)
h['hTauPt'] = ROOT.TH1F ("hTau_Pt", ";p_{T};", 500, 0, 500)
h['hJetPt'] = ROOT.TH1F ("hJet_Pt", ";p_{T};", 2000, 0, 2000)

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

        event.getByLabel(labelVertex, handleVertex)
        vertex=handleVertex.product()
        pv = vertex[0].position()

        event.getByLabel(labelMuon, handleMuon)
        muons=handleMuon.product()

        event.getByLabel(labelElectron, handleElectron)
        electrons=handleElectron.product()
        
        event.getByLabel(labelTausBoosted, handleTausBoosted)
        btaus = handleTausBoosted.product()

        event.getByLabel(labelBs, handleBs)
        bs=handleBs.product()

        event.getByLabel(labelConv, handleConv)
        convs=handleConv.product()

        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()
        
        event.getByLabel(labelJet, handleJet)
        jets=handleJet.product()

        event.getByLabel(labelRho, handleRho)
        rho=handleRho.product()[0]

        event.getByLabel(labelHLT, handleHLT)
        triggerResults=handleHLT.product()
        names = event.object().triggerNames(triggerResults)

        event.getByLabel(labelPatMETs, handlePatMETs)
        met=handlePatMETs.product().front()
        mets=[]
        mets+=[met]
        mets.sort(key=lambda x: x.pt(), reverse=True)
        
        h['hNEvent'].Fill(0.5, 1)
#        h['hNEvent'].Fill(1.5, genweight)


        selected_muons=[]
        for muon in muons:
            if not muon.isLooseMuon(): continue
            if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue 
            if muon.pt()<3 or muon.eta()>2.4: continue
            if muonIsoCut(muon)<0.25:
                selected_muons+=[muon]

        selected_muons.sort(key=lambda x: x.pt(), reverse=True) 
        if len(selected_muons)>0:
            h['hMuPt'].Fill(selected_muons[0].pt())

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
                and electron.hadronicOverEm()< (0.05 + 1.16/E_c + 0.0324*rho/E_c) \
                and abs(electron.gsfTrack().dz(pv))<0.1 \
                and abs(electron.gsfTrack().dxy(pv))<0.05:
                    selected_electrons+=[electron]

            if electron.isEE():
                if electron.full5x5_sigmaIetaIeta() < 0.0425 \
                and GsfEleEInverseMinusPInverse(electron) < 0.111 \
                and abs(dEtaInSeed(electron)) < 0.00674 \
                and GsfEleMissingHitsCut(electron)<=1 \
                and electron.passConversionVeto() \
                and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.169 \
                and electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c) \
                and abs(electron.gsfTrack().dz(pv)) < 0.2 \
                and abs(electron.gsfTrack().dxy(pv)) < 0.1:
                    selected_electrons+=[electron]

        selected_electrons.sort(key=lambda x: x.pt(), reverse=True) 
        if len(selected_electrons)>0:
            h['hElePt'].Fill(selected_electrons[0].pt())

        selected_btaus=[]
        
        for tau in btaus:
            if tau.pt() > 20 \
            and abs(tau.eta()) < 2.3 \
            and tau.leadChargedHadrCand().get().dz(vertex[0].position()) < 0.5 \
            and tau.leadChargedHadrCand().get().dxy(vertex[0].position()) < 0.2 \
            and tau.tauID("decayModeFinding") \
            and tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"):
                selected_btaus+=[tau]

        selected_btaus.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_btaus)>0:
            h['hTauPt'].Fill(selected_btaus[0].pt())

        selected_jets=[]
        selected_bjets=[]
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
                if jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.9535:
                    selected_bjets+=[jet]

        selected_jets.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_jets)>0:
            h['hJetPt'].Fill(selected_jets[0].pt())


        if len(selected_btaus)>=2 and len(selected_jets)>=1 and len(selected_muons)==0 and len(selected_electrons)==0 and len(selected_bjets)==0:

            btau1 = ROOT.TLorentzVector()
            btau1.SetPtEtaPhiM(selected_btaus[0].pt(), selected_btaus[0].eta(), selected_btaus[0].phi(), selected_btaus[0].mass())

            btau2 = ROOT.TLorentzVector()
            btau2.SetPtEtaPhiM(selected_btaus[1].pt(), selected_btaus[1].eta(), selected_btaus[1].phi(), selected_btaus[1].mass())
            
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            if (btau1+btau2).M() >= 10:

                h['hTauTauBaseline_M'].Fill((btau1+btau2).M())

                if btau1.DeltaR(btau2)<0.4 and btau1.DeltaR(jet)>0.8 and btau2.DeltaR(jet)>0.8:

                    h['hTauTaudR_M'].Fill((btau1+btau2).M())

                    m = ROOT.TLorentzVector()
                    m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

                    if jet.Pt()>500 \
                    and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) \
                    or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) \
                    or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5"))):

                        h['hTauTauTrig_M'].Fill((btau1+btau2).M())

                        if m.Pt()>100:

                            h['hTauTauMetcut_M'].Fill((btau1+btau2).M())

                            if selected_btaus[0].charge()*selected_btaus[1].charge()<0:
                                h['OS'].Fill((btau1+btau2).M())
                                h['OS_lPt'].Fill(btau1.Pt())
                                h['OS_sPt'].Fill(btau2.Pt())

                            if selected_btaus[0].charge()*selected_btaus[1].charge()>0:
                                h['SS'].Fill((btau1+btau2).M())
                                h['SS_lPt'].Fill(btau1.Pt())
                                h['SS_sPt'].Fill(btau2.Pt())

                        if m.Pt()<100:

                            h['hTauTauInvMetcut_M'].Fill((btau1+btau2).M())

                            if selected_btaus[0].charge()*selected_btaus[1].charge()<0:
                                h['OS_invMET'].Fill((btau1+btau2).M())
                                h['OS_invMET_lPt'].Fill(btau1.Pt())
                                h['OS_invMET_sPt'].Fill(btau2.Pt())
                            if selected_btaus[0].charge()*selected_btaus[1].charge()>0:
                                h['SS_invMET'].Fill((btau1+btau2).M())
                                h['SS_invMET_lPt'].Fill(btau1.Pt())
                                h['SS_invMET_sPt'].Fill(btau2.Pt())

out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
