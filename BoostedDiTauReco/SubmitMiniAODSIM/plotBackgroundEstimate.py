#Code for tau_mu tau_mu studies

import ROOT,sys,os
from DataFormats.FWLite import Events, Handle
from looseElectron import *
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

#inputFileListDir="./filelists/"+Sample+"/"
inputFileListName=sys.argv[1]

#inputFileList=inputFileListDir+inputFileListName
inputFileList=inputFileListName

if len(sys.argv)>2:
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./plots/"

outputFileName = outputFileDir+"h_BE_lowMET_charge_"+inputFileListName.split("/")[-1].replace(".txt",".root")

print outputFileName

pi = math.pi

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlineSlimmedPrimaryVertices')

handleBoostedTau = Handle ('vector<pat::Tau>')
labelBoostedTau = ('slimmedTausBoosted')

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults','','HLT')

handleBs = Handle ('reco::BeamSpot')
labelBs = ("offlineBeamSpot")

handleConv = Handle ('vector<reco::Conversion>')
labelConv = ('reducedEgamma', 'reducedConversions')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('slimmedJets')

handlePatMETs = Handle("vector<pat::MET>")
labelPatMETs = ( 'slimmedMETs')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

out=ROOT.TFile.Open(outputFileName,'recreate')

# <histograms>

pi = math.pi

h = {}

h['hNEvent'] = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)

h['hMuPt'] = ROOT.TH1F ("hMu_Pt", ";p_{T};", 500, 0, 500)
h['hElePt'] = ROOT.TH1F ("hEle_Pt", ";p_{T};", 500, 0, 500)
h['hTauPt_iso'] = ROOT.TH1F ("hTau_iso_Pt", ";p_{T};", 500, 0, 500)
h['hTauPt_aiso'] = ROOT.TH1F ("hTau_aiso_Pt", ";p_{T};", 500, 0, 500)
h['hJetPt'] = ROOT.TH1F ("hJet_Pt", ";p_{T};", 2000, 0, 2000)

h['hTauTauBaseline_oo_M'] = ROOT.TH1F ("hTauTau_Baseline_oo_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTaudR_oo_M'] = ROOT.TH1F ("hTauTau_dR_oo_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauTrig_oo_M'] = ROOT.TH1F ("hTauTau_Trig_oo_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauInvMetcut_oo_M'] = ROOT.TH1F ("hTauTau_Metcut_oo_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['hTauTauBaseline_oo_charge'] = ROOT.TH1F ("hTauTau_Baseline_oo_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)
h['hTauTauTrig_oo_charge'] = ROOT.TH1F ("hTauTau_Trig_oo_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)
h['hTauTaudR_oo_charge'] = ROOT.TH1F ("hTauTau_dR_oo_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)

h['OS_oo'] = ROOT.TH1F ("OS_oo_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['SS_oo'] = ROOT.TH1F ("SS_oo_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['OS_oo_lPt'] = ROOT.TH1F ("OS_oo_lPt", ";p_{T};", 500, 0, 500)
h['OS_oo_sPt'] = ROOT.TH1F ("OS_oo_sPt", ";p_{T};", 500, 0, 500)

h['SS_oo_lPt'] = ROOT.TH1F ("SS_oo_lPt", ";p_{T};", 500, 0, 500)
h['SS_oo_sPt'] = ROOT.TH1F ("SS_oo_sPt", ";p_{T};", 500, 0, 500)

h['hTauTauBaseline_ii_M'] = ROOT.TH1F ("hTauTau_Baseline_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTaudR_ii_M'] = ROOT.TH1F ("hTauTau_dR_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauTrig_ii_M'] = ROOT.TH1F ("hTauTau_Trig_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauInvMetcut_ii_M'] = ROOT.TH1F ("hTauTau_Metcut_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['hTauTauBaseline_ii_charge'] = ROOT.TH1F ("hTauTau_Baseline_ii_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)
h['hTauTauTrig_ii_charge'] = ROOT.TH1F ("hTauTau_Trig_ii_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)
h['hTauTaudR_ii_charge'] = ROOT.TH1F ("hTauTau_dR_ii_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)

h['OS_ii'] = ROOT.TH1F ("OS_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['SS_ii'] = ROOT.TH1F ("SS_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['OS_ii_lPt'] = ROOT.TH1F ("OS_ii_lPt", ";p_{T};", 500, 0, 500)
h['OS_ii_sPt'] = ROOT.TH1F ("OS_ii_sPt", ";p_{T};", 500, 0, 500)

h['SS_ii_lPt'] = ROOT.TH1F ("SS_ii_lPt", ";p_{T};", 500, 0, 500)
h['SS_ii_sPt'] = ROOT.TH1F ("SS_ii_sPt", ";p_{T};", 500, 0, 500)

h['hTauTauBaseline_io_M'] = ROOT.TH1F ("hTauTau_Baseline_io_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTaudR_io_M'] = ROOT.TH1F ("hTauTau_dR_io_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauTrig_io_M'] = ROOT.TH1F ("hTauTau_Trig_io_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauInvMetcut_io_M'] = ROOT.TH1F ("hTauTau_Metcut_io_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['hTauTauBaseline_io_charge'] = ROOT.TH1F ("hTauTau_Baseline_io_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)
h['hTauTauTrig_io_charge'] = ROOT.TH1F ("hTauTau_Trig_io_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)
h['hTauTaudR_io_charge'] = ROOT.TH1F ("hTauTau_dR_io_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)

h['OS_io'] = ROOT.TH1F ("OS_io_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['SS_io'] = ROOT.TH1F ("SS_io_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['OS_io_lPt'] = ROOT.TH1F ("OS_io_lPt", ";p_{T};", 500, 0, 500)
h['OS_io_sPt'] = ROOT.TH1F ("OS_io_sPt", ";p_{T};", 500, 0, 500)

h['SS_io_lPt'] = ROOT.TH1F ("SS_io_lPt", ";p_{T};", 500, 0, 500)
h['SS_io_sPt'] = ROOT.TH1F ("SS_io_sPt", ";p_{T};", 500, 0, 500)

h['hTauTauBaseline_oi_M'] = ROOT.TH1F ("hTauTau_Baseline_oi_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTaudR_oi_M'] = ROOT.TH1F ("hTauTau_dR_oi_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauTrig_oi_M'] = ROOT.TH1F ("hTauTau_Trig_oi_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauInvMetcut_oi_M'] = ROOT.TH1F ("hTauTau_Metcut_oi_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['hTauTauBaseline_oi_charge'] = ROOT.TH1F ("hTauTau_Baseline_oi_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)
h['hTauTauTrig_oi_charge'] = ROOT.TH1F ("hTauTau_Trig_oi_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)
h['hTauTaudR_oi_charge'] = ROOT.TH1F ("hTauTau_dR_oi_charge", "#tau - #tau charge; charge;", 3, -1.5, 1.5)

h['OS_oi'] = ROOT.TH1F ("OS_oi_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['SS_oi'] = ROOT.TH1F ("SS_oi_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['OS_oi_lPt'] = ROOT.TH1F ("OS_oi_lPt", ";p_{T};", 500, 0, 500)
h['OS_oi_sPt'] = ROOT.TH1F ("OS_oi_sPt", ";p_{T};", 500, 0, 500)

h['SS_oi_lPt'] = ROOT.TH1F ("SS_oi_lPt", ";p_{T};", 500, 0, 500)
h['SS_oi_sPt'] = ROOT.TH1F ("SS_oi_sPt", ";p_{T};", 500, 0, 500)



#</Histograms>

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    print inputFileName.replace("\n","")
    inputFileName=inputFileName.replace("\n","")
    f=ROOT.TFile.Open(inputFileName)

    if not f.IsZombie():
        events=Events(inputFileName)
    else:
        print "Can't Open File: "+inputFileName
        continue

    for event in events:
        
        event.getByLabel(labelVertex, handleVertex)
        vertex=handleVertex.product()
        
        event.getByLabel(labelMuon, handleMuon)
        muons=handleMuon.product()

        event.getByLabel(labelElectron, handleElectron)
        electrons=handleElectron.product()
        
        event.getByLabel(labelBoostedTau, handleBoostedTau)
        btaus = handleBoostedTau.product()

        event.getByLabel(labelBs, handleBs)
        bs=handleBs.product()

        event.getByLabel(labelConv, handleConv)
        convs=handleConv.product()

        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()
        
        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()

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
        h['hNEvent'].Fill(1.5, genweight)

        selected_muons=[]
        for muon in muons:
            if not muon.isLooseMuon(): continue
            #if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue 
            if muon.pt()<3 or muon.eta()>2.4: continue
            if muonIsoCut(muon)<0.25:
                selected_muons+=[muon]

        if len(selected_muons)>0:
            selected_muons.sort(key=lambda x: x.pt(), reverse=True)
            h['hMuPt'].Fill(selected_muons[0].pt(), genweight)

         #<electronSelection>
        selected_electrons=[]
        for electron in electrons:
            if electron.pt()<7: continue 
            if abs(electron.eta())>2.5: continue
            if electron.isEB(): 
                if electron.full5x5_sigmaIetaIeta()<0.011 \
                and electron.hadronicOverEm()<0.298 \
                and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.222 \
                and GsfEleEInverseMinusPInverse(electron)<0.241 \
                and abs(dEtaInSeed(electron))<0.00477 \
                and GsfEleMissingHitsCut(electron)<=1 \
                and electron.passConversionVeto() \
                and GsfEleEffAreaPFIsoCut(electron, rho)<0.0994:
                #and abs(electron.gsfTrack().dz(vertex[0].position()))<0.1 \
                #and abs(electron.gsfTrack().dxy(vertex[0].position()))<0.05:
                    selected_electrons+=[electron]
            if electron.isEE():
                if electron.full5x5_sigmaIetaIeta()<0.0314 \
                and electron.hadronicOverEm()<0.101 \
                and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.213 \
                and GsfEleEInverseMinusPInverse(electron)<0.14 \
                and abs(dEtaInSeed(electron))<0.00868 \
                and GsfEleMissingHitsCut(electron)<=1 \
                and electron.passConversionVeto() \
                and GsfEleEffAreaPFIsoCut(electron, rho)<0.107:
                #and abs(electron.gsfTrack().dz(vertex[0].position()))<0.2 \
                #and abs(electron.gsfTrack().dxy(vertex[0].position()))<0.1:
                    selected_electrons+=[electron]

        if len(selected_electrons)>0:
            selected_electrons.sort(key=lambda x: x.pt(), reverse=True) 
            h['hElePt'].Fill(selected_electrons[0].pt(), genweight)

        #<\electronSelection>

        selected_btaus_aiso=[]
        selected_btaus_iso=[]

                #if tau.leadChargedHadrCand().get().dz(vertex[0].position()) < 0.5
                #and tau.leadChargedHadrCand().get().dxy(vertex[0].position()) < 0.2
        
        for tau in btaus:
            if tau.pt() > 20 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") and tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"):
                if not tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT"): selected_btaus_aiso+=[tau]
                if tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT"): selected_btaus_iso+=[tau]

        if len(selected_btaus_aiso)>0:
            selected_btaus_aiso.sort(key=lambda x: x.pt(), reverse=True)
            h['hTauPt_aiso'].Fill(selected_btaus_aiso[0].pt(), genweight)

        if len(selected_btaus_iso)>0:
            selected_btaus_iso.sort(key=lambda x: x.pt(), reverse=True)
            h['hTauPt_iso'].Fill(selected_btaus_iso[0].pt(), genweight)

        leading_btaus=[]
        trailing_btaus=[]

        if len(selected_btaus_aiso)>0 and len(selected_btaus_iso)>0:
            if selected_btaus_iso[0].pt() > selected_btaus_aiso[0].pt(): leading_btaus.extend((selected_btaus_iso[0], selected_btaus_aiso[0]))
            if selected_btaus_iso[0].pt() < selected_btaus_aiso[0].pt(): trailing_btaus.extend((selected_btaus_aiso[0], selected_btaus_iso[0]))


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

        if len(selected_jets)>0:
            selected_jets.sort(key=lambda x: x.pt(), reverse=True)
            h['hJetPt'].Fill(selected_jets[0].pt(), genweight)


#Taus

        if len(selected_btaus_aiso)>=2 and len(selected_jets)>=1 and len(selected_muons)==0 and len(selected_electrons)==0 and len(selected_bjets)==0:

            btau1 = ROOT.TLorentzVector()
            btau1.SetPtEtaPhiM(selected_btaus_aiso[0].pt(), selected_btaus_aiso[0].eta(), selected_btaus_aiso[0].phi(), selected_btaus_aiso[0].mass())

            btau2 = ROOT.TLorentzVector()
            btau2.SetPtEtaPhiM(selected_btaus_aiso[1].pt(), selected_btaus_aiso[1].eta(), selected_btaus_aiso[1].phi(), selected_btaus_aiso[1].mass())
            
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            m = ROOT.TLorentzVector()
            m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

            if m.Pt()<100:

                h['hTauTauBaseline_oo_M'].Fill((btau1+btau2).M(), genweight)
                h['hTauTauBaseline_oo_charge'].Fill(selected_btaus_aiso[0].charge()*selected_btaus_aiso[1].charge(), genweight)
                
                if btau1.DeltaR(btau2)<0.4 and btau1.DeltaR(jet)>0.8 and btau2.DeltaR(jet)>0.8:

                    h['hTauTaudR_oo_M'].Fill((btau1+btau2).M(), genweight)
                    h['hTauTaudR_oo_charge'].Fill(selected_btaus_aiso[0].charge()*selected_btaus_aiso[1].charge(), genweight)

                    if jet.Pt()>500 \
                    and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) \
                    or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) \
                    or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5"))):

                        h['hTauTauTrig_oo_M'].Fill((btau1+btau2).M(), genweight)
                        h['hTauTauTrig_oo_charge'].Fill(selected_btaus_aiso[0].charge()*selected_btaus_aiso[1].charge(), genweight)

#                        h['hTauTauInvMetcut_oo_M'].Fill((btau1+btau2).M(), genweight)

                        if selected_btaus_aiso[0].charge()*selected_btaus_aiso[1].charge()<0:
                            h['OS_oo'].Fill((btau1+btau2).M(), genweight)
                            h['OS_oo_lPt'].Fill(btau1.Pt(), genweight)
                            h['OS_oo_sPt'].Fill(btau2.Pt(), genweight)
                        if selected_btaus_aiso[0].charge()*selected_btaus_aiso[1].charge()>0:
                            h['SS_oo'].Fill((btau1+btau2).M(), genweight)
                            h['SS_oo_lPt'].Fill(btau1.Pt(), genweight)
                            h['SS_oo_sPt'].Fill(btau2.Pt(), genweight)



        if len(selected_btaus_iso)>=2 and len(selected_jets)>=1 and len(selected_muons)==0 and len(selected_electrons)==0 and len(selected_bjets)==0:

            btau1 = ROOT.TLorentzVector()
            btau1.SetPtEtaPhiM(selected_btaus_iso[0].pt(), selected_btaus_iso[0].eta(), selected_btaus_iso[0].phi(), selected_btaus_iso[0].mass())

            btau2 = ROOT.TLorentzVector()
            btau2.SetPtEtaPhiM(selected_btaus_iso[1].pt(), selected_btaus_iso[1].eta(), selected_btaus_iso[1].phi(), selected_btaus_iso[1].mass())
            
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            m = ROOT.TLorentzVector()
            m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

            if m.Pt()<100:

                h['hTauTauBaseline_ii_M'].Fill((btau1+btau2).M(), genweight)
                h['hTauTauBaseline_ii_charge'].Fill(selected_btaus_iso[0].charge()*selected_btaus_iso[1].charge(), genweight)

                if btau1.DeltaR(btau2)<0.4 and btau1.DeltaR(jet)>0.8 and btau2.DeltaR(jet)>0.8:

                    h['hTauTaudR_ii_M'].Fill((btau1+btau2).M(), genweight)
                    h['hTauTaudR_ii_charge'].Fill(selected_btaus_iso[0].charge()*selected_btaus_iso[1].charge(), genweight)

                    if jet.Pt()>500 \
                    and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) \
                    or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) \
                    or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5"))):

                        h['hTauTauTrig_ii_M'].Fill((btau1+btau2).M(), genweight)
                        h['hTauTauTrig_ii_charge'].Fill(selected_btaus_iso[0].charge()*selected_btaus_iso[1].charge(), genweight)

#                        h['hTauTauInvMetcut_ii_M'].Fill((btau1+btau2).M(), genweight)

                        if selected_btaus_iso[0].charge()*selected_btaus_iso[1].charge()<0:
                            h['OS_ii'].Fill((btau1+btau2).M(), genweight)
                            h['OS_ii_lPt'].Fill(btau1.Pt(), genweight)
                            h['OS_ii_sPt'].Fill(btau2.Pt(), genweight)
                        if selected_btaus_iso[0].charge()*selected_btaus_iso[1].charge()>0:
                            h['SS_ii'].Fill((btau1+btau2).M(), genweight)
                            h['SS_ii_lPt'].Fill(btau1.Pt(), genweight) 
                            h['SS_ii_sPt'].Fill(btau2.Pt(), genweight)

        if len(leading_btaus)>=2 and len(selected_jets)>=1 and len(selected_muons)==0 and len(selected_electrons)==0 and len(selected_bjets)==0:

            btau1 = ROOT.TLorentzVector()
            btau1.SetPtEtaPhiM(leading_btaus[0].pt(), leading_btaus[0].eta(), leading_btaus[0].phi(), leading_btaus[0].mass())

            btau2 = ROOT.TLorentzVector()
            btau2.SetPtEtaPhiM(leading_btaus[1].pt(), leading_btaus[1].eta(), leading_btaus[1].phi(), leading_btaus[1].mass())
            
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            m = ROOT.TLorentzVector()
            m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

            if m.Pt()<100:
                
                h['hTauTauBaseline_io_M'].Fill((btau1+btau2).M(), genweight)
                h['hTauTauBaseline_io_charge'].Fill(leading_btaus[0].charge()*leading_btaus[1].charge(), genweight)

                if btau1.DeltaR(btau2)<0.4 and btau1.DeltaR(jet)>0.8 and btau2.DeltaR(jet)>0.8:

                    h['hTauTaudR_io_M'].Fill((btau1+btau2).M(), genweight)
                    h['hTauTaudR_io_charge'].Fill(leading_btaus[0].charge()*leading_btaus[1].charge(), genweight)
 
                    if jet.Pt()>500 \
                    and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) \
                    or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) \
                    or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5"))):

                        h['hTauTauTrig_io_M'].Fill((btau1+btau2).M(), genweight)
                        h['hTauTauTrig_io_charge'].Fill(leading_btaus[0].charge()*leading_btaus[1].charge(), genweight)

#                        h['hTauTauInvMetcut_io_M'].Fill((btau1+btau2).M(), genweight)

                        if leading_btaus[0].charge()*leading_btaus[1].charge()<0:
                            h['OS_io'].Fill((btau1+btau2).M(), genweight)
                            h['OS_io_lPt'].Fill(btau1.Pt(), genweight)
                            h['OS_io_sPt'].Fill(btau2.Pt(), genweight)
                        if leading_btaus[0].charge()*leading_btaus[1].charge()>0:
                            h['SS_io'].Fill((btau1+btau2).M(), genweight)
                            h['SS_io_lPt'].Fill(btau1.Pt(), genweight) 
                            h['SS_io_sPt'].Fill(btau2.Pt(), genweight)

        if len(trailing_btaus)>=2 and len(selected_jets)>=1 and len(selected_muons)==0 and len(selected_electrons)==0 and len(selected_bjets)==0:

            btau1 = ROOT.TLorentzVector()
            btau1.SetPtEtaPhiM(trailing_btaus[0].pt(), trailing_btaus[0].eta(), trailing_btaus[0].phi(), trailing_btaus[0].mass())

            btau2 = ROOT.TLorentzVector()
            btau2.SetPtEtaPhiM(trailing_btaus[1].pt(), trailing_btaus[1].eta(), trailing_btaus[1].phi(), trailing_btaus[1].mass())
            
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            m = ROOT.TLorentzVector()
            m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

            if m.Pt()<100:

                h['hTauTauBaseline_oi_M'].Fill((btau1+btau2).M(), genweight)
                h['hTauTauBaseline_oi_charge'].Fill(trailing_btaus[0].charge()*trailing_btaus[1].charge(), genweight)

                if btau1.DeltaR(btau2)<0.4 and btau1.DeltaR(jet)>0.8 and btau2.DeltaR(jet)>0.8:

                    h['hTauTaudR_oi_M'].Fill((btau1+btau2).M(), genweight)
                    h['hTauTaudR_oi_charge'].Fill(trailing_btaus[0].charge()*trailing_btaus[1].charge(), genweight)

                    if jet.Pt()>500 \
                    and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) \
                    or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) \
                    or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5"))):

                        h['hTauTauTrig_oi_M'].Fill((btau1+btau2).M(), genweight)
                        h['hTauTauTrig_oi_charge'].Fill(trailing_btaus[0].charge()*trailing_btaus[1].charge(), genweight)

#                        h['hTauTauInvMetcut_oi_M'].Fill((btau1+btau2).M(), genweight)

                        if trailing_btaus[0].charge()*trailing_btaus[1].charge()<0:
                            h['OS_oi'].Fill((btau1+btau2).M(), genweight)
                            h['OS_oi_lPt'].Fill(btau1.Pt(), genweight)
                            h['OS_oi_sPt'].Fill(btau2.Pt(), genweight)
                        if trailing_btaus[0].charge()*trailing_btaus[1].charge()>0:
                            h['SS_oi'].Fill((btau1+btau2).M(), genweight)
                            h['SS_oi_lPt'].Fill(btau1.Pt(), genweight) 
                            h['SS_oi_sPt'].Fill(btau2.Pt(), genweight)


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
