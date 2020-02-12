import ROOT,sys,os
from DataFormats.FWLite import Events, Handle
from looseElectron import *
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

inputFileListName=sys.argv[1]

inputFileList=inputFileListName

if len(sys.argv)>2:
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./plots/"
outputFileName = outputFileDir+"h_"+inputFileListName.split("/")[-1].replace(".txt",".root")
print outputFileName

out=ROOT.TFile.Open(outputFileName,'recreate')

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleBoostedTau = Handle ('vector<pat::Tau>')
labelBoostedTau = ('slimmedTausBoosted')

handleTaus = Handle ('vector<pat::Tau>')
labelTaus = ('slimmedTaus')

handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlineSlimmedPrimaryVertices')

handleBs = Handle ('reco::BeamSpot')
labelBs = ("offlineBeamSpot")

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults','','HLT')

handleConv = Handle ('vector<reco::Conversion>')
labelConv = ('reducedEgamma', 'reducedConversions')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('slimmedJets')

handleGenJet = Handle ('vector<reco::GenJet>')
labelGenJet = ('slimmedGenJets')

handlePatMETs = Handle("vector<pat::MET>")
labelPatMETs = ( 'slimmedMETs')

# Histograms
hNEvent = ROOT.TH1F ("hNEvent", "Number of Events;;N_{event}", 2, 0, 2)
hBTauBaseline_NbJets = ROOT.TH1F ("hBTau_Baseline_NbJets", "", 10, 0, 10)

hBTauBaseline_M = ROOT.TH1F ("hBTau_Baseline_M", ";M;N_{events}", 100, 0, 100)
hBTauBaseline_M_dR = ROOT.TH2F ("hBTau_Baseline_M_dR", ";M;#Delta R",100, 0, 100, 100, 0, 5)
hBTauBaseline_dRtj = ROOT.TH2F ("hBTau_Baseline_dRtj", ";#Delta R_{#tau1j}; #Delta R_{#tau2j}", 100, 0, 5, 100, 0, 5)
hBTauBaseline_pt1 = ROOT.TH1F ("hBTau_Baseline_pt1", ";p_{T};", 50, 0, 500)
hBTauBaseline_pt2 = ROOT.TH1F ("hBTau_Baseline_pt2", ";p_{T};", 50, 0, 500)
hBTauBaseline_jPt = ROOT.TH1F ("hBTau_Baseline_jPt", "jet pt;P_{t};", 2000, 0, 2000)
hBTauBaseline_METPt = ROOT.TH1F ("hBTau_Baseline_METPt", ";p_{T};", 500, 0, 500)
hBTauBaseline_M_METPt = ROOT.TH2F ("hBTau_Baseline_M_METPt", ";M;p_{T};", 100, 0, 100, 500, 0, 500)

hBTauTrig_M = ROOT.TH1F ("hBTau_Trig_M", ";M;N_{events}", 100, 0, 100)
hBTauTrig_M_dR = ROOT.TH2F ("hBTau_Trig_M_dR", ";M;#Delta R",100, 0, 100, 100, 0, 5)
hBTauTrig_dRtj = ROOT.TH2F ("hBTau_Trig_dRtj", ";#Delta R_{#tau1j}; #Delta R_{#tau2j}", 100, 0, 5, 100, 0, 5)
hBTauTrig_pt1 = ROOT.TH1F ("hBTau_Trig_pt1", ";p_{T};", 50, 0, 500)
hBTauTrig_pt2 = ROOT.TH1F ("hBTau_Trig_pt2", ";p_{T};", 50, 0, 500)
hBTauTrig_METPt = ROOT.TH1F ("hBTau_Trig_METPt", ";p_{T};", 500, 0, 500)
hBTauTrig_M_METPt = ROOT.TH2F ("hBTau_Trig_M_METPt", ";M;p_{T};", 100, 0, 100, 500, 0, 500)
hBTauTrig_jPt = ROOT.TH1F ("hBTau_Trig_jPt", "jet pt;P_{t};", 2000, 0, 2000)

hBTaudR_M = ROOT.TH1F ("hBTau_dR_M", ";M;N_{events}", 100, 0, 100)
hBTaudR_M_dR = ROOT.TH2F ("hBTau_dR_M_dR", ";M;#Delta R",100, 0, 100, 100, 0, 5)
hBTaudR_dRtj = ROOT.TH2F ("hBTau_dR_dRtj", ";#Delta R_{#tau1j}; #Delta R_{#tau2j}", 100, 0, 5, 100, 0, 5)
hBTaudR_pt1 = ROOT.TH1F ("hBTau_dR_pt1", ";p_{T};", 50, 0, 500)
hBTaudR_pt2 = ROOT.TH1F ("hBTau_dR_pt2", ";p_{T};", 50, 0, 500)
hBTaudR_METPt = ROOT.TH1F ("hBTau_dR_METPt", ";p_{T};", 500, 0, 500)
hBTaudR_M_METPt = ROOT.TH2F ("hBTau_dR_M_METPt", ";M;p_{T};", 100, 0, 100, 500, 0, 500)
hBTaudR_jPt = ROOT.TH1F ("hBTau_dR_jPt", "jet pt;P_{t};", 2000, 0, 2000)

hBTauMetcut_M = ROOT.TH1F ("hBTau_Metcut_M", ";M;N_{events}", 100, 0, 100)
hBTauMetcut_M_dR = ROOT.TH2F ("hBTau_Metcut_M_dR", ";M;#Delta R",100, 0, 100, 100, 0, 5)
hBTauMetcut_dRtj = ROOT.TH2F ("hBTau_Metcut_dRtj", ";#Delta R_{#tau1j}; #Delta R_{#tau2j}", 100, 0, 5, 100, 0, 5)
hBTauMetcut_pt1 = ROOT.TH1F ("hBTau_Metcut_pt1", ";p_{T};", 50, 0, 500)
hBTauMetcut_pt2 = ROOT.TH1F ("hBTau_Metcut_pt2", ";p_{T};", 50, 0, 500)
hBTauMetcut_METPt = ROOT.TH1F ("hBTau_Metcut_METPt", ";p_{T};", 500, 0, 500)
hBTauMetcut_M_METPt = ROOT.TH2F ("hBTau_Metcut_M_METPt", ";M;p_{T};", 100, 0, 100, 500, 0, 500)
hBTauMetcut_jPt = ROOT.TH1F ("hBTau_Metcut_jPt", "jet pt;P_{t};", 2000, 0, 2000)

# -----------

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

        event.getByLabel(labelMuon, handleMuon)
        muons=handleMuon.product()

        event.getByLabel(labelElectron, handleElectron)
        electrons=handleElectron.product()

        event.getByLabel(labelBoostedTau, handleBoostedTau)
        btaus = handleBoostedTau.product()

        event.getByLabel(labelTaus, handleTaus)
        taus = handleTaus.product()

        event.getByLabel(labelVertex, handleVertex)
        vertex=handleVertex.product()

        event.getByLabel(labelBs, handleBs)
        bs=handleBs.product()

        event.getByLabel(labelConv, handleConv)
        convs=handleConv.product()

        event.getByLabel(labelRho, handleRho)
        rho=handleRho.product()[0]

        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()

        event.getByLabel(labelHLT, handleHLT)
        triggerResults=handleHLT.product()
        names = event.object().triggerNames(triggerResults)

        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()

        event.getByLabel(labelJet, handleJet)
        jets=handleJet.product()

        event.getByLabel(labelGenJet, handleGenJet)

        event.getByLabel(labelPatMETs, handlePatMETs)
        met=handlePatMETs.product().front()
        mets=[]
        mets+=[met]
        mets.sort(key=lambda x: x.pt(), reverse=True)
        m=ROOT.TLorentzVector()
        m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

        hNEvent.Fill(0.5, 1)
        hNEvent.Fill(1.5, genweight)
     
        selected_btaus=[]

        for tau in btaus:
            if tau.pt()<20 or abs(tau.eta())>2.3 or tau.leadChargedHadrCand().get().dxy(vertex[0].position())>0.2 or tau.leadChargedHadrCand().get().dz(vertex[0].position())>0.5: continue
            if not tau.tauID("decayModeFinding"): continue
            if tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") < -0.5: continue
            if not tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"): continue
            selected_btaus+=[tau]

        selected_btaus.sort(key=lambda x: x.pt(), reverse=True)

        selected_muons=[]
        for muon in muons:
            if not muon.isLooseMuon(): continue
            if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue
            if muon.pt()<3 or muon.eta()>2.4: continue
            if muonIsoCut(muon)<0.25:
                selected_muons+=[muon]

        selected_muons.sort(key=lambda x: x.pt(), reverse=True)

        selected_electrons=[]
        for electron in electrons:
            if electron.pt()<7: continue
            if abs(electron.eta())>2.5: continue
            if electron.isEB():
                if not electron.full5x5_sigmaIetaIeta()<0.011 or not electron.hadronicOverEm()<0.298 or not abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.222 or not GsfEleEInverseMinusPInverse(electron)<0.241 or not abs(dEtaInSeed(electron))<0.00477 or not GsfEleMissingHitsCut(electron)<=1 or not electron.passConversionVeto() or not GsfEleEffAreaPFIsoCut(electron, rho)<0.0994 or not abs(electron.gsfTrack().dz(vertex[0].position()))<0.1 or not abs(electron.gsfTrack().dxy(vertex[0].position()))<0.05:
                    continue
            if electron.isEE():
                if not electron.full5x5_sigmaIetaIeta()<0.0314 or not electron.hadronicOverEm()<0.101 or not abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.213 or not GsfEleEInverseMinusPInverse(electron)<0.14 or not abs(dEtaInSeed(electron))<0.00868 or not GsfEleMissingHitsCut(electron)<=1 or not electron.passConversionVeto() or not GsfEleEffAreaPFIsoCut(electron, rho)<0.107 or not abs(electron.gsfTrack().dz(vertex[0].position()))<0.2 or not abs(electron.gsfTrack().dxy(vertex[0].position()))<0.1:
                    continue
            selected_electrons+=[electron]

        selected_electrons.sort(key=lambda x: x.pt(), reverse=True)

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
        selected_bjets.sort(key=lambda x: x.pt(), reverse=True)

        if len(selected_btaus)>=2 and selected_btaus[0].charge()*selected_btaus[1].charge()<0 and len(selected_jets)>=1 and len(selected_muons)==len(selected_electrons)==0:
            if len(selected_bjets)>=0:
                hBTauBaseline_NbJets.Fill(len(selected_bjets), genweight)

            if len(selected_bjets)==0:

                btau1 = ROOT.TLorentzVector()
                btau1.SetPtEtaPhiM(selected_btaus[0].pt(), selected_btaus[0].eta(), selected_btaus[0].phi(), selected_btaus[0].mass())

                btau2 = ROOT.TLorentzVector()
                btau2.SetPtEtaPhiM(selected_btaus[1].pt(), selected_btaus[1].eta(), selected_btaus[1].phi(), selected_btaus[1].mass())

                jet = ROOT.TLorentzVector()
                jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

                hBTauBaseline_M.Fill((btau1+btau2).M(), genweight)
                hBTauBaseline_M_dR.Fill((btau1+btau2).M(), btau1.DeltaR(btau2), genweight)
                hBTauBaseline_dRtj.Fill(btau1.DeltaR(jet), btau2.DeltaR(jet), genweight)
                hBTauBaseline_pt1.Fill(btau1.Pt(), genweight)
                hBTauBaseline_pt2.Fill(btau2.Pt(), genweight)
                hBTauBaseline_METPt.Fill(m.Pt(), genweight)
                hBTauBaseline_M_METPt.Fill((btau1+btau2).M(), m.Pt(), genweight)
                hBTauBaseline_jPt.Fill(jet.Pt())

                if (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5"))):
                    hBTauTrig_M.Fill((btau1+btau2).M(), genweight)
                    hBTauTrig_M_dR.Fill((btau1+btau2).M(), btau1.DeltaR(btau2), genweight)
                    hBTauTrig_dRtj.Fill(btau1.DeltaR(jet), btau2.DeltaR(jet), genweight)
                    hBTauTrig_pt1.Fill(btau1.Pt(), genweight)
                    hBTauTrig_pt2.Fill(btau2.Pt(), genweight)
                    hBTauTrig_METPt.Fill(m.Pt(), genweight)
#                    hBTauTrig_M_METPt.Fill((btau1+btau2).M(), m.Pt(), genweight)
                    hBTauTrig_jPt.Fill(jet.Pt())

                    if btau1.DeltaR(btau2)<0.4 and btau1.DeltaR(jet)>0.8 and btau2.DeltaR(jet)>0.8:
                        hBTaudR_M.Fill((btau1+btau2).M(), genweight)
#                        hBTaudR_M_dR.Fill((btau1+btau2).M(), btau1.DeltaR(btau2), genweight)
#                        hBTaudR_dRtj.Fill(btau1.DeltaR(jet), btau2.DeltaR(jet), genweight)
                        hBTaudR_pt1.Fill(btau1.Pt(), genweight)
                        hBTaudR_pt2.Fill(btau2.Pt(), genweight)
#                        hBTaudR_METPt.Fill(m.Pt(), genweight)
                        hBTaudR_M_METPt.Fill((btau1+btau2).M(), m.Pt(), genweight)
                        hBTaudR_jPt.Fill(jet.Pt())

                        if m.Pt()>100:
                            hBTauMetcut_M.Fill((btau1+btau2).M(), genweight)
                            hBTauMetcut_M_dR.Fill((btau1+btau2).M(), btau1.DeltaR(btau2), genweight)
#                            hBTauMetcut_dRtj.Fill(btau1.DeltaR(jet), btau2.DeltaR(jet), genweight)
                            hBTauMetcut_pt1.Fill(btau1.Pt(), genweight)
                            hBTauMetcut_pt2.Fill(btau2.Pt(), genweight)
#                            hBTauMetcut_METPt.Fill(m.Pt(), genweight)
#                            hBTauMetcut_M_METPt.Fill((btau1+btau2).M(), m.Pt(), genweight)
                            hBTauMetcut_jPt.Fill(jet.Pt())

out.cd()

hNEvent.Write()
hBTauBaseline_NbJets.Write()

hBTauBaseline_M.Write()
hBTauBaseline_M_dR.Write()
hBTauBaseline_dRtj.Write()
hBTauBaseline_pt1.Write()
hBTauBaseline_pt2.Write()
hBTauBaseline_METPt.Write()
hBTauBaseline_M_METPt.Write()
hBTauBaseline_jPt.Write()

hBTauTrig_M.Write()
hBTauTrig_M_dR.Write()
hBTauTrig_dRtj.Write()
hBTauTrig_pt1.Write()
hBTauTrig_pt2.Write()
hBTauTrig_METPt.Write()
hBTauTrig_M_METPt.Write()
hBTauTrig_jPt.Write()

hBTaudR_M.Write()
hBTaudR_M_dR.Write()
hBTaudR_dRtj.Write()
hBTaudR_pt1.Write()
hBTaudR_pt2.Write()
hBTaudR_METPt.Write()
hBTaudR_M_METPt.Write()
hBTaudR_jPt.Write()

hBTauMetcut_M.Write()
hBTauMetcut_M_dR.Write()
hBTauMetcut_dRtj.Write()
hBTauMetcut_pt1.Write()
hBTauMetcut_pt2.Write() 
hBTauMetcut_METPt.Write()
hBTauMetcut_M_METPt.Write()
hBTauMetcut_jPt.Write()

out.Close()
