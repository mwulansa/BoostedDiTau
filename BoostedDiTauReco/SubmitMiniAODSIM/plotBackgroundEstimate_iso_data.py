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

outputFileName = outputFileDir+"h_BE_Iso_80_"+inputFileListName.split("/")[-1].replace(".txt",".root")

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

h['hTauTauBaseline_ii_M'] = ROOT.TH1F ("hTauTau_Baseline_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTaudR_ii_M'] = ROOT.TH1F ("hTauTau_dR_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauTrig_ii_M'] = ROOT.TH1F ("hTauTau_Trig_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauInvMetcut_ii_M'] = ROOT.TH1F ("hTauTau_InvMetcut_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauMetcut_ii_M'] = ROOT.TH1F ("hTauTau_Metcut_ii_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['OS_ii_lowMET'] = ROOT.TH1F ("OS_ii_lowMET_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['SS_ii_lowMET'] = ROOT.TH1F ("SS_ii_lowMET_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['OS_ii_lowMET_lPt'] = ROOT.TH1F ("OS_ii_lowMET_lPt", ";p_{T};", 500, 0, 500)
h['OS_ii_lowMET_sPt'] = ROOT.TH1F ("OS_ii_lowMET_sPt", ";p_{T};", 500, 0, 500)

h['SS_ii_lowMET_lPt'] = ROOT.TH1F ("SS_ii_lowMET_lPt", ";p_{T};", 500, 0, 500)
h['SS_ii_lowMET_sPt'] = ROOT.TH1F ("SS_ii_lowMET_sPt", ";p_{T};", 500, 0, 500)

h['OS_ii_hiMET'] = ROOT.TH1F ("OS_ii_hiMET_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['SS_ii_hiMET'] = ROOT.TH1F ("SS_ii_hiMET_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['OS_ii_hiMET_lPt'] = ROOT.TH1F ("OS_ii_hiMET_lPt", ";p_{T};", 500, 0, 500)
h['OS_ii_hiMET_sPt'] = ROOT.TH1F ("OS_ii_hiMET_sPt", ";p_{T};", 500, 0, 500)

h['SS_ii_hiMET_lPt'] = ROOT.TH1F ("SS_ii_hiMET_lPt", ";p_{T};", 500, 0, 500)
h['SS_ii_hiMET_sPt'] = ROOT.TH1F ("SS_ii_hiMET_sPt", ";p_{T};", 500, 0, 500)


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

        selected_muons=[]
        for muon in muons:
            if not muon.isLooseMuon(): continue
            #if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue 
            if muon.pt()<3 or muon.eta()>2.4: continue
            if muonIsoCut(muon)<0.25:
                selected_muons+=[muon]

        if len(selected_muons)>0:
            selected_muons.sort(key=lambda x: x.pt(), reverse=True)
            h['hMuPt'].Fill(selected_muons[0].pt())

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
            h['hElePt'].Fill(selected_electrons[0].pt())

        #<\electronSelection>

#        selected_btaus_aiso=[]
        selected_btaus_iso=[]

                #if tau.leadChargedHadrCand().get().dz(vertex[0].position()) < 0.5
                #and tau.leadChargedHadrCand().get().dxy(vertex[0].position()) < 0.2
        
        for tau in btaus:
            if tau.pt() > 20 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") :
#                if not tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"): selected_btaus_aiso+=[tau]
                if tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"): selected_btaus_iso+=[tau]

        # if len(selected_btaus_aiso)>0:
        #     selected_btaus_aiso.sort(key=lambda x: x.pt(), reverse=True)
        #     h['hTauPt_aiso'].Fill(selected_btaus_aiso[0].pt())

        # if len(selected_btaus_iso)>0:
        #     selected_btaus_iso.sort(key=lambda x: x.pt(), reverse=True)
        #     h['hTauPt_iso'].Fill(selected_btaus_iso[0].pt())

        # leading_btaus=[]
        # trailing_btaus=[]

        # if len(selected_btaus_aiso)>0 and len(selected_btaus_iso)>0:
        #     if selected_btaus_iso[0].pt() > selected_btaus_aiso[0].pt(): leading_btaus.extend((selected_btaus_iso[0], selected_btaus_aiso[0]))
        #     if selected_btaus_iso[0].pt() < selected_btaus_aiso[0].pt(): trailing_btaus.extend((selected_btaus_aiso[0], selected_btaus_iso[0]))

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
            h['hJetPt'].Fill(selected_jets[0].pt())

#Taus

        if len(selected_btaus_iso)>=2 and len(selected_jets)>=1 and len(selected_muons)==0 and len(selected_electrons)==0 and len(selected_bjets)==0:

            btau1 = ROOT.TLorentzVector()
            btau1.SetPtEtaPhiM(selected_btaus_iso[0].pt(), selected_btaus_iso[0].eta(), selected_btaus_iso[0].phi(), selected_btaus_iso[0].mass())

            btau2 = ROOT.TLorentzVector()
            btau2.SetPtEtaPhiM(selected_btaus_iso[1].pt(), selected_btaus_iso[1].eta(), selected_btaus_iso[1].phi(), selected_btaus_iso[1].mass())
            
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            h['hTauTauBaseline_ii_M'].Fill((btau1+btau2).M())

            if btau1.DeltaR(btau2)<0.4 and btau1.DeltaR(jet)>0.8 and btau2.DeltaR(jet)>0.8:

                h['hTauTaudR_ii_M'].Fill((btau1+btau2).M())

                m = ROOT.TLorentzVector()
                m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

                if jet.Pt()>500 \
                and (triggerResults.accept(5) \
                or triggerResults.accept(183) \
                or triggerResults.accept(212)):     
   
                    h['hTauTauTrig_ii_M'].Fill((btau1+btau2).M())

                    if m.Pt()<80:

                        h['hTauTauInvMetcut_ii_M'].Fill((btau1+btau2).M())

                        if selected_btaus_iso[0].charge()*selected_btaus_iso[1].charge()<0:
                            h['OS_ii_lowMET'].Fill((btau1+btau2).M())
                            h['OS_ii_lowMET_lPt'].Fill(btau1.Pt())
                            h['OS_ii_lowMET_sPt'].Fill(btau2.Pt())
                        if selected_btaus_iso[0].charge()*selected_btaus_iso[1].charge()>0:
                            h['SS_ii_lowMET'].Fill((btau1+btau2).M())
                            h['SS_ii_lowMET_lPt'].Fill(btau1.Pt())
                            h['SS_ii_lowMET_sPt'].Fill(btau2.Pt())

                    if m.Pt()>80:

                        h['hTauTauMetcut_ii_M'].Fill((btau1+btau2).M())

                        if selected_btaus_iso[0].charge()*selected_btaus_iso[1].charge()>0:
                            h['SS_ii_hiMET'].Fill((btau1+btau2).M())
                            h['SS_ii_hiMET_lPt'].Fill(btau1.Pt())
                            h['SS_ii_hiMET_sPt'].Fill(btau2.Pt())


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
