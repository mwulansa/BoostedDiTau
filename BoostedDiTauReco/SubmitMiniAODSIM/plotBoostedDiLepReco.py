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
outputFileName = outputFileDir+"h_"+inputFileListName.split("/")[-1].replace(".txt",".root")
print outputFileName

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
labelHLT = ('TriggerResults','','HLT')

handleBs = Handle ('reco::BeamSpot')
labelBs = ("offlineBeamSpot")

handleConv = Handle ('vector<reco::Conversion>')
labelConv = ('reducedEgamma', 'reducedConversions')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handleGenJet = Handle ('vector<reco::GenJet>')
labelGenJet = ('slimmedGenJets')

handlePatMETs = Handle("vector<pat::MET>")
labelPatMETs = ( 'slimmedMETs')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

out=ROOT.TFile.Open(outputFileName,'recreate')

hNEvent = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)
hMET_pt = ROOT.TH1F ("hMET", "", 500, 0, 500)

hJet1Pt=ROOT.TH1F ("hJet1Pt", "jet pt;P_{t};", 2000, 0, 2000)
hJet2Pt=ROOT.TH1F ("hJet2Pt", "jet pt;P_{t};", 2000, 0, 2000)
hJet1PtTrig=ROOT.TH1F ("hJetTrigPt", "triggered jet pt;P_{t};", 2000, 0, 2000)
hJet1Pt_genMatched=ROOT.TH1F ("hJet1Pt_genMatched", "jet pt;P_{t};", 2000, 0, 2000)

# these 3 hists may not useful now
hMu1Pt=ROOT.TH1F ("hMu1Pt", "muon pt;P_{t};", 500, 0, 500)
hMu1PtIsoMuTrig= ROOT.TH1F ("hMu1PtIsoMuTrig", "muon pt;P_{t};", 500, 0, 500)
hMu1PtMuTrig= ROOT.TH1F ("hMu1PtMuTrig", "muon pt;P_{t};", 500, 0, 500)

# these 2 hists may not useful as well
hMuMuTrig_jPt = ROOT.TH1F ("hMuMuJetTrigPt", "triggered jet pt;P_{t};", 2000, 0, 2000)
hEMuTrig_jPt = ROOT.TH1F ("hEMuJetTrigPt", "triggered jet pt;P_{t};", 2000, 0, 2000)

hMuMuBaseline_M = ROOT.TH1F ("hMuMuBaseline_M", "#mu - #mu inv. mass;M_{#mu#mu};", 1000, 0, 200)
hMuMu_M = ROOT.TH1F ("hMuMu_M", "#mu - #mu inv. mass;M_{#mu#mu};", 1000, 0, 200)
hMuMuTrig_M = ROOT.TH1F ("hMuMuTrig_M", "#mu - #mu inv. mass;M_{#mu#mu};", 1000, 0, 200)
hMuMu_M_genMatched = ROOT.TH1F ("hMuMu_M_genMatched", "#mu - #mu inv. mass;M_{#mu#mu};", 1000, 0, 200)

hMuMuBaseline_lPt = ROOT.TH1F ("hMuMuBaseline_Mu1Pt", "leading muon pt;P_{t};", 500, 0, 500)
hMuMuTrig_lPt = ROOT.TH1F ("hMuMuTrig_Mu1Pt", "leading muon pt;P_{t};", 500, 0, 500)
hMuMu_lPt = ROOT.TH1F ("hMuMu_Mu1Pt", "leading muon pt;P_{t};", 500, 0, 500)
hMuMu_lPt_genMatched = ROOT.TH1F ("hMuMu_Mu1Pt_genMatched", "leading muon pt;P_{t};", 500, 0, 500)

hMuMuBaseline_lEta = ROOT.TH1F ("hMuMuBaseline_Mu1Eta", "leading muon eta;#eta;", 60, -3, 3)
hMuMuTrig_lEta = ROOT.TH1F ("hMuMuTrig_Mu1Eta", "leading muon eta;#eta;", 60, -3, 3)
hMuMu_lEta = ROOT.TH1F ("hMuMu_Mu1Eta", "leading muon eta;#eta;", 60, -3, 3)

hMuMuBaseline_tPt = ROOT.TH1F ("hMuMuBaseline_Mu2Pt", "sub-leading muon pt;P_{t};", 500, 0, 500)
hMuMu_tPt = ROOT.TH1F ("hMuMu_Mu2Pt", "sub-leading muon pt;P_{t};", 500, 0, 500)
hMuMuTrig_tPt = ROOT.TH1F ("hMuMuTrig_Mu2Pt", "sub-leading muon pt;P_{t};", 500, 0, 500)
hMuMu_tPt_genMatched = ROOT.TH1F ("hMuMu_Mu2Pt_genMatched", "sub-leading muon pt;P_{t};", 500, 0, 500)

hMuMuBaseline_ePt = ROOT.TH1F ("hMuMuBaseline_ePt", ";P_{t};", 500, 0, 500)
hMuMuBaseline_NJets = ROOT.TH1F ("hMuMuBaseline_NJets", ";N_{jets};", 10, 0, 10)
hMuMuBaseline_NbJets = ROOT.TH1F ("hMuMuBaseline_NbJets", ";N_{b-jets};", 10, 0, 10)

hMuMuBaseline_tEta = ROOT.TH1F ("hMuMuBaseline_Mu2Eta", "leading muon eta;#eta;", 60, -3, 3)
hMuMuTrig_tEta = ROOT.TH1F ("hMuMuTrig_Mu2Eta", "leading muon eta;#eta;", 60, -3, 3)
hMuMu_tEta = ROOT.TH1F ("hMuMu_Mu2Eta", "leading muon eta;#eta;", 60, -3, 3)
    
hEMuBaseline_M = ROOT.TH1F ("hEMuBaseline_M", "e - #mu inv. mass;M_{e#mu};", 1000, 0, 200)
hEMuTrig_M = ROOT.TH1F ("hEMuTrig_M", "e - #mu inv. mass;M_{e#mu};", 1000, 0, 200)
hEMu_M = ROOT.TH1F ("hEMu_M", "e - #mu inv. mass;M_{e#mu};", 1000, 0, 200)
hEMu_M_genMatched = ROOT.TH1F ("hEMu_M_genMatched", "e - #mu inv. mass;M_{e#mu};", 1000, 0, 200)
    
hEMuBaseline_muPt = ROOT.TH1F ("hEMuBaseline_MuPt", "muon pt;P_{t};", 500, 0, 500)
hEMuTrig_muPt = ROOT.TH1F ("hEMuTrig_MuPt", "muon pt;P_{t};", 500, 0, 500)
hEMu_muPt = ROOT.TH1F ("hEMu_MuPt", "muon pt;P_{t};", 500, 0, 500)
hEMu_muPt_genMatched = ROOT.TH1F ("hEMu_MuPt_genMatched", "muon pt;P_{t};", 500, 0, 500)

hEMuBaseline_tMuPt = ROOT.TH1F ("hEMuBaseline_tMuPt", ";P_{t};", 500, 0, 500)
hEMuBaseline_NJets = ROOT.TH1F ("hEMuBaseline_NJets", ";N_{jets};", 10, 0, 10)
hEMuBaseline_NbJets = ROOT.TH1F ("hEMuBaseline_NbJets", ";N_{b-jets};", 10, 0, 10)

hEMuBaseline_muEta = ROOT.TH1F ("hEMuBaseline_MuEta", "leading muon eta;#eta;", 60, -3, 3)
hEMuTrig_muEta = ROOT.TH1F ("hEMuTrig_MuEta", "leading muon eta;#eta;", 60, -3, 3)
hEMu_muEta = ROOT.TH1F ("hEMu_MuEta", "leading muon eta;#eta;", 60, -3, 3)
    
hEMuBaseline_ePt = ROOT.TH1F ("hEMuBaseline_EPt", "electron pt;P_{t};", 500, 0, 500)
hEMuTrig_ePt = ROOT.TH1F ("hEMuTrig_EPt", "electron pt;P_{t};", 500, 0, 500)
hEMu_ePt = ROOT.TH1F ("hEMu_EPt", "electron pt;P_{t};", 500, 0, 500)
hEMu_ePt_genMatched = ROOT.TH1F ("hEMu_EPt_genMatched", "electron pt;P_{t};", 500, 0, 500)

hEMuBaseline_eEta = ROOT.TH1F ("hEMuBaseline_EEta", "leading muon eta;#eta;", 60, -3, 3)
hEMuTrig_eEta = ROOT.TH1F ("hEMuTrig_EEta", "leading muon eta;#eta;", 60, -3, 3)
hEMu_eEta = ROOT.TH1F ("hEMu_EEta", "leading muon eta;#eta;", 60, -3, 3)

hMuMuBaseline_dR = ROOT.TH1F ("hMuMuBaseline_dR", ";#delta R;", 200, 0, 10)
hMuMuBaseline_dRjlMu = ROOT.TH1F ("hMuMuBaseline_dRjlMu", ";#delta R;", 200, 0, 10)
hMuMuBaseline_dRjtMu = ROOT.TH1F ("hMuMuBaseline_dRjtMu", ";#delta R;", 200, 0, 10)
hMuMuBaseline_M_dR = ROOT.TH2F ("hMuMuBaseline_M_dR", ";M_{#mu#mu};#delta R;", 1000, 0, 200, 200, 0, 10)
hMuMuBaseline_M_dRjlMu = ROOT.TH2F ("hMuMuBaseline_M_dRjlMu", ";M_{#mu#mu};#delta R;", 1000, 0, 200, 200, 0, 10)
hMuMuBaseline_M_dRjtMu = ROOT.TH2F ("hMuMuBaseline_M_dRjtMu", ";M_{#mu#mu};#delta R;", 1000, 0, 200, 200, 0, 10)
hMuMuBaseline_dRjlMu_dRjtMu = ROOT.TH2F ("hMuMuBaseline_dRjlMu_dRjtMu", ";#delta R;#delta R;", 200, 0, 10, 200, 0, 10)
hMuMuBaseline_lMuIso_dRjlMu = ROOT.TH2F ("hMuMuBaseline_lMuIso_dRjlMu", ";Isolation;#delta R;", 100, 0, 0.25, 200, 0, 10)
hMuMuBaseline_tMuIso_dRjtMu = ROOT.TH2F ("hMuMuBaseline_tMuIso_dRjtMu", ";Isolation;#delta R;", 100, 0, 0.25, 200, 0, 10)
hMuMuBaseline_lMuIso_lPt = ROOT.TH2F ("hMuMuBaseline_lMuIso_lPt", ";Isolation;P_{t};", 100, 0, 0.25, 500, 0, 500)
hMuMuBaseline_tMuIso_tPt = ROOT.TH2F ("hMuMuBaseline_tMuIso_tPt", ";Isolation;P_{t};", 100, 0, 0.25, 500, 0, 500)
hMuMuBaseline_jPt_dRjlMu = ROOT.TH2F ("hMuMuBaseline_jPt_dRjlMu", ";P_{t};#Delta R;", 2000, 0, 2000, 200, 0, 10)
hMuMuBaseline_jPt_dRjtMu = ROOT.TH2F ("hMuMuBaseline_jPt_dRjtMu", ";P_{t};#Delta R;", 2000, 0, 2000, 200, 0, 10)

hMuMuTrig_dR = ROOT.TH1F ("hMuMuTrig_dR", ";#delta R;", 200, 0, 10)
hMuMuTrig_dRjlMu = ROOT.TH1F ("hTrig_dRjlMu", ";#delta R;", 200, 0, 10)
hMuMuTrig_dRjtMu = ROOT.TH1F ("hTrig_dRjtMu", ";#delta R;", 200, 0, 10)
hMuMu_dR = ROOT.TH1F ("hMuMu_dR", ";#delta R;", 20, 0, 1)
hMuMu_dR_genMatched = ROOT.TH1F ("hMuMu_dR_genMatched", ";#delta R;", 20, 0, 1)

hEMuBaseline_dR = ROOT.TH1F ("hEMuBaseline_dR", ";#delta R;", 200, 0, 10)
hEMuBaseline_dRjE = ROOT.TH1F ("hEMuBaseline_dRjE", ";#delta R;", 200, 0, 10)
hEMuBaseline_dRjMu = ROOT.TH1F ("hEMuBaseline_dRjMu", ";#delta R;", 200, 0, 10)
hEMuBaseline_M_dR = ROOT.TH2F ("hEMuBaseline_M_dR", ";M_{e#mu};#delta R;", 1000, 0, 200, 200, 0, 10)
hEMuBaseline_M_dRjE = ROOT.TH2F ("hEMuBaseline_M_dRjE", ";M_{e#mu};#delta R;", 1000, 0, 200, 200, 0, 10)
hEMuBaseline_M_dRjMu = ROOT.TH2F ("hEMuBaseline_M_dRjMu", ";M_{e#mu};#delta R;", 1000, 0, 200, 200, 0, 10)
hEMuBaseline_dRjE_dRjMu = ROOT.TH2F ("hEMuBaseline_dRjE_dRjMu", ";#delta R;#delta R;", 200, 0, 10, 200, 0, 10)
hEMuBaseline_eIso_dRjE = ROOT.TH2F ("hEMuBaseline_eIso_dRjE", ";Isolation;#delta R;", 100, 0, 0.2, 200, 0, 10)
hEMuBaseline_muIso_dRjMu = ROOT.TH2F ("hEMuBaseline_muIso_dRjMu", ";Isolation;#delta R;", 100, 0, 0.25, 200, 0, 10)
hEMuBaseline_eIso_ePt = ROOT.TH2F ("hEMuBaseline_eIso_ePt", ";Isolation;P_{t};", 100, 0, 0.2, 500, 0, 500)
hEMuBaseline_muIso_muPt = ROOT.TH2F ("hEMuBaseline_muIso_muPt", ";Isolation;P_{t};", 100, 0, 0.25, 500, 0, 500)
hEMuBaseline_jPt_dRjMu = ROOT.TH2F ("hMuMuBaseline_jPt_dRjMu", ";P_{t};#Delta R;", 2000, 0, 2000, 200, 0, 10)
hEMuBaseline_jPt_dRjE = ROOT.TH2F ("hMuMuBaseline_jPt_dRjE", ";P_{t};#Delta R;", 2000, 0, 2000, 200, 0, 10)

hEMuTrig_dR = ROOT.TH1F ("hEMuTrig_dR", ";#delta R;", 200, 0, 10)
hEMuTrig_dRjE = ROOT.TH1F ("hEMuTrig_dRjE", ";#delta R;", 200, 0, 10)
hEMuTrig_dRjMu = ROOT.TH1F ("hEMuTrig_dRjMu", ";#delta R;", 200, 0, 10)
hEMu_dR = ROOT.TH1F ("hEMu_dR", ";#delta R;", 20, 0, 1)
hEMu_dR_genMatched = ROOT.TH1F ("hEMu_dR_genMatched", ";#delta R;", 20, 0, 1)

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

        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()

        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()

        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()
        
        event.getByLabel(labelPatMETs, handlePatMETs)
        met=handlePatMETs.product().front()
        hMET_pt.Fill(met.pt(), genweight)

        hNEvent.Fill(0.5, 1)
        hNEvent.Fill(1.5, genweight)
        
        # gen information
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

        event.getByLabel(labelGenJet, handleGenJet)
        #jets=handleGenJet.product()
    
        genJets=[]
        for jet in handleGenJet.product():
            genJets+=[jet]
            
        genJets.sort(key=lambda x: x.pt(), reverse=True) 

        # muon selection
        selected_muons=[]
        for muon in muons:
            if not muon.isLooseMuon(): continue
            if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue 
            if muon.pt()<3 or muon.eta()>2.4: continue
            if muonIsoCut(muon)<0.25:
                selected_muons+=[muon]

        selected_muons.sort(key=lambda x: x.pt(), reverse=True) 
        
        if len(selected_muons)>0: 
            hMu1Pt.Fill(selected_muons[0].pt(), genweight) 
            if triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")): 
                hMu1PtMuTrig.Fill(selected_muons[0].pt(), genweight)
            if triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")):
                hMu1PtIsoMuTrig.Fill(selected_muons[0].pt(), genweight)

        # electron selection
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

        # jet selection
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


        if len(selected_jets)>=1:
            #btags=selected_jets[0].tagInfoLabels()
            #print selected_jets[0].hasTagInfo("pfCombinedInclusiveSecondaryVertexV2BJetTags")
            #print selected_jets[0].hasTagInfo("pfCombinedSecondaryVertexV2BJetTags")
            #print selected_jets[0].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")
            #print selected_jets[0].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags")
            
            #for tag in btags:
            #    print tag
            #sys.exit()
            hJet1Pt.Fill(selected_jets[0].pt(), genweight) 
            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
            if len(genJets)>0:
                genJet1=ROOT.TLorentzVector(genJets[0].px(), genJets[0].py(), genJets[0].pz(), genJets[0].energy()) 
                if genJet1.DeltaR(j)<0.4: 
                    hJet1Pt_genMatched.Fill(selected_jets[0].pt(), genweight)
            if len(selected_jets)>=2:
                hJet2Pt.Fill(selected_jets[1].pt(), genweight)
            if triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")):
                hJet1PtTrig.Fill(selected_jets[0].pt(), genweight)

        # tau_mu tau_mu selection
        if len(selected_muons)>=2 and len(selected_jets)>=1 and selected_muons[0].charge()*selected_muons[1].charge()<0:
            hMuMuBaseline_NJets.Fill(len(selected_jets), genweight)
            hMuMuBaseline_NbJets.Fill(len(selected_bjets), genweight)
            if len(selected_electrons)>0:
                hMuMuBaseline_ePt.Fill(selected_electrons[0].pt())
                
            mu1=ROOT.TLorentzVector()
            mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass()) 
            iso1=muonIsoCut(selected_muons[0])
            
            mu2=ROOT.TLorentzVector()
            mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass())
            iso2=muonIsoCut(selected_muons[1])

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
            #if selected_muons[1].pt()<10:
            

            hMuMuBaseline_M.Fill((mu1+mu2).M(), genweight) 
            hMuMuBaseline_lPt.Fill(mu1.Pt(), genweight)
            hMuMuBaseline_tPt.Fill(mu2.Pt(), genweight)
            hMuMuBaseline_lEta.Fill(mu1.Eta(), genweight)
            hMuMuBaseline_tEta.Fill(mu2.Eta(), genweight)
            hMuMuBaseline_dR.Fill(mu1.DeltaR(mu2), genweight)
            hMuMuBaseline_dRjlMu.Fill(mu1.DeltaR(j), genweight)
            hMuMuBaseline_dRjtMu.Fill(mu2.DeltaR(j), genweight)
            hMuMuBaseline_M_dR.Fill((mu1+mu2).M(), mu1.DeltaR(mu2), genweight)
            hMuMuBaseline_M_dRjlMu.Fill((mu1+mu2).M(), mu1.DeltaR(j), genweight)
            hMuMuBaseline_M_dRjtMu.Fill((mu1+mu2).M(), mu2.DeltaR(j), genweight)
            hMuMuBaseline_dRjlMu_dRjtMu.Fill(mu1.DeltaR(j), mu2.DeltaR(j), genweight)
            hMuMuBaseline_lMuIso_dRjlMu.Fill(iso1, mu1.DeltaR(j), genweight)
            hMuMuBaseline_tMuIso_dRjtMu.Fill(iso2, mu2.DeltaR(j), genweight)
            hMuMuBaseline_lMuIso_lPt.Fill(iso1, mu1.Pt(), genweight)
            hMuMuBaseline_tMuIso_tPt.Fill(iso2, mu2.Pt(), genweight)
            hMuMuBaseline_jPt_dRjlMu.Fill(j.Pt(), mu1.DeltaR(j), genweight)
            hMuMuBaseline_jPt_dRjtMu.Fill(j.Pt(), mu2.DeltaR(j), genweight)

            if mu1.DeltaR(j)<0.02:
                NHF  = selected_jets[0].neutralHadronEnergyFraction() 
                NEMF = selected_jets[0].neutralEmEnergyFraction() 
                CHF  = selected_jets[0].chargedHadronEnergyFraction() 
                MUF  = selected_jets[0].muonEnergyFraction() 
                CEMF = selected_jets[0].chargedEmEnergyFraction()
                NumConst = selected_jets[0].chargedMultiplicity()+selected_jets[0].neutralMultiplicity()
                NumNeutralParticles =selected_jets[0].neutralMultiplicity()
                CHM      = selected_jets[0].chargedMultiplicity()
                #print "---"
                #print j.Pt(), j.Eta(), j.Phi()
                #print mu1.Pt(), mu1.Eta(), mu1.Phi()
                #print NHF, NEMF, CHF, MUF, CEMF, NumConst, NumNeutralParticles, CHM
                #print selected_muons[0].pfIsolationR04().sumChargedHadronPt, selected_muons[0].pfIsolationR04().sumPhotonEt, selected_muons[0].pfIsolationR04().sumNeutralHadronEt, selected_muons[0].pfIsolationR04().sumPUPt, selected_muons[0].pt()
                #print muonIsoCut(selected_muons[0]), selected_muons[0].isMediumMuon(), selected_muons[0].isTightMuon(vertex[0]), j.DeltaR(mu1)

            # trigger
            if (mu1.Pt()>26 and (triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")))) or (mu1.Pt()>50 and (triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")))) or (j.Pt()>500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))):
                hMuMuTrig_M.Fill((mu1+mu2).M(), genweight)
                hMuMuTrig_lPt.Fill(mu1.Pt(), genweight)
                hMuMuTrig_tPt.Fill(mu2.Pt(), genweight)
                hMuMuTrig_lEta.Fill(mu1.Eta(), genweight)
                hMuMuTrig_tEta.Fill(mu2.Eta(), genweight)
                hMuMuTrig_dR.Fill(mu1.DeltaR(mu2), genweight)
                hMuMuTrig_dRjlMu.Fill(mu1.DeltaR(j), genweight)
                hMuMuTrig_dRjtMu.Fill(mu2.DeltaR(j), genweight)
                 
                hMuMuTrig_jPt.Fill(j.Pt(), genweight)

                if mu1.DeltaR(mu2)<0.4 and mu1.DeltaR(j)>0.8 and mu2.DeltaR(j)>0.8:
                    #if mu1.DeltaR(mu2)<0.4 and mu1.DeltaR(j)>0.8 and mu2.DeltaR(j)>0.8 and (selected_muons[0].pfIsolationR04().sumChargedHadronPt+max(0,selected_muons[0].pfIsolationR04().sumPhotonEt+selected_muons[0].pfIsolationR04().sumNeutralHadronEt-0.5*selected_muons[0].pfIsolationR04().sumPUPt))/selected_muons[0].pt()<0.25: 
                    hMuMu_dR.Fill(mu1.DeltaR(mu2), genweight)
                    hMuMu_M.Fill((mu1+mu2).M(), genweight)
                    hMuMu_lPt.Fill(mu1.Pt(), genweight)
                    hMuMu_tPt.Fill(mu2.Pt(), genweight)
                    hMuMu_lEta.Fill(mu1.Eta(), genweight)
                    hMuMu_tEta.Fill(mu2.Eta(), genweight)

                    # gen match
                    evtMatched=False
                    if len(genMuons)==2 and len(genJets)>=1:
                        if not len(genElectrons)==0: print "Debug 1: Something is Wrong..."
                        genMu1=ROOT.TLorentzVector()
                        genMu1.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

                        genMu2=ROOT.TLorentzVector()
                        genMu2.SetPtEtaPhiM(genMuons[1].pt(), genMuons[1].eta(), genMuons[1].phi(), genMuons[1].mass())

                        if (mu1.DeltaR(genMu1)<0.3 or mu1.DeltaR(genMu2)<0.3) and (mu2.DeltaR(genMu1)<0.3 or mu2.DeltaR(genMu2)<0.3): 
                            hMuMu_M_genMatched.Fill((mu1+mu2).M(), genweight)
                            hMuMu_dR_genMatched.Fill(mu1.DeltaR(mu2), genweight)
                            hMuMu_lPt_genMatched.Fill(mu1.Pt(), genweight)
                            hMuMu_tPt_genMatched.Fill(mu2.Pt(), genweight)
                    

        # tau_mu tau_e selection
        if len(selected_muons)>0 and len(selected_electrons)>0 and len(selected_jets)>0 and selected_muons[0].charge()*selected_electrons[0].charge()<0:

            hEMuBaseline_NJets.Fill(len(selected_jets), genweight)
            hEMuBaseline_NbJets.Fill(len(selected_bjets), genweight)
            if len(selected_muons)>1:
                hEMuBaseline_tMuPt.Fill(selected_muons[1].pt())
            
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())
            isomu=muonIsoCut(selected_muons[0])

            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())
            isoe=GsfEleEffAreaPFIsoCut(selected_electrons[0], rho)

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
        
            hEMuBaseline_M.Fill((mu+e).M(), genweight) 
            hEMuBaseline_muPt.Fill(mu.Pt(), genweight)
            hEMuBaseline_ePt.Fill(e.Pt(), genweight)
            hEMuBaseline_muEta.Fill(mu.Eta(), genweight)
            hEMuBaseline_eEta.Fill(e.Eta(), genweight)
            hEMuBaseline_dR.Fill(mu.DeltaR(e), genweight)
            hEMuBaseline_dRjMu.Fill(mu.DeltaR(j), genweight)
            hEMuBaseline_dRjE.Fill(e.DeltaR(j), genweight)
            hEMuBaseline_M_dR.Fill((e+mu).M(), mu.DeltaR(e), genweight)
            hEMuBaseline_M_dRjE.Fill((e+mu).M(), e.DeltaR(j), genweight)
            hEMuBaseline_M_dRjMu.Fill((e+mu).M(), mu.DeltaR(j), genweight)
            hEMuBaseline_dRjE_dRjMu.Fill(e.DeltaR(j), mu.DeltaR(j), genweight)
            hEMuBaseline_eIso_dRjE.Fill(isoe, e.DeltaR(j), genweight)
            hEMuBaseline_muIso_dRjMu.Fill(isomu, mu.DeltaR(j), genweight)
            hEMuBaseline_eIso_ePt.Fill(isoe, e.Pt(), genweight)
            hEMuBaseline_muIso_muPt.Fill(isomu, mu.Pt(), genweight)
            hEMuBaseline_jPt_dRjE.Fill(j.Pt(), e.DeltaR(j), genweight)
            hEMuBaseline_jPt_dRjMu.Fill(j.Pt(), mu.DeltaR(j), genweight)

            # trigger
            if mu.Pt()>26 and (triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4"))) or mu.Pt()>50 and (triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5"))) or j.Pt()>500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5"))):
                hEMuTrig_M.Fill((mu+e).M(), genweight)
                hEMuTrig_muPt.Fill(mu.Pt(), genweight)
                hEMuTrig_ePt.Fill(e.Pt(), genweight)
                hEMuTrig_muEta.Fill(mu.Eta(), genweight)
                hEMuTrig_eEta.Fill(e.Eta(), genweight)
                hEMuTrig_dR.Fill(mu.DeltaR(e), genweight)
                hEMuTrig_dRjMu.Fill(mu.DeltaR(j), genweight)
                hEMuTrig_dRjE.Fill(e.DeltaR(j), genweight)
                
                hEMuTrig_jPt.Fill(j.Pt(), genweight)
            
                if mu.DeltaR(e)<0.4 and mu.DeltaR(j)>0.8 and e.DeltaR(j)>0.8: 
                    hEMu_M.Fill((mu+e).M(), genweight)
                    hEMu_muPt.Fill(mu.Pt(), genweight)
                    hEMu_ePt.Fill(e.Pt(), genweight)
                    hEMu_muEta.Fill(mu.Eta(), genweight)
                    hEMu_eEta.Fill(e.Eta(), genweight)
                    hEMu_dR.Fill(mu.DeltaR(e), genweight)

                    # gen match
                    if len(genMuons)==1 and len(genElectrons)==1: 
                        genMu=ROOT.TLorentzVector()
                        genMu.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

                        genEle=ROOT.TLorentzVector()
                        genEle.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())

                        if e.DeltaR(genEle)<0.3 and mu.DeltaR(genMu)<0.3: 
                            hEMu_M_genMatched.Fill((e+mu).M(), genweight)
                            hEMu_dR_genMatched.Fill(mu.DeltaR(e), genweight)
                            hEMu_muPt_genMatched.Fill(mu.Pt(), genweight)
                            hEMu_ePt_genMatched.Fill(e.Pt(), genweight)


out.cd()
hNEvent.Write()
hMET_pt.Write()

hJet1Pt.Write()
hJet1Pt_genMatched.Write()
hJet2Pt.Write()
hJet1PtTrig.Write()

hMuMuTrig_jPt.Write()
hEMuTrig_jPt.Write()

hMu1Pt.Write()
hMu1PtIsoMuTrig.Write()
hMu1PtMuTrig.Write()

hMuMuBaseline_ePt.Write()
hMuMuBaseline_NJets.Write()
hMuMuBaseline_NbJets.Write()

hMuMuBaseline_M.Write()
hMuMuBaseline_lPt.Write()
hMuMuBaseline_tPt.Write()
hMuMuBaseline_lEta.Write()
hMuMuBaseline_tEta.Write()
hMuMuBaseline_dR.Write()
hMuMuBaseline_dRjlMu.Write()
hMuMuBaseline_dRjtMu.Write()
hMuMuBaseline_M_dR.Write()
hMuMuBaseline_M_dRjlMu.Write()
hMuMuBaseline_M_dRjtMu.Write()
hMuMuBaseline_dRjlMu_dRjtMu.Write()
hMuMuBaseline_lMuIso_dRjlMu.Write()
hMuMuBaseline_tMuIso_dRjtMu.Write()
hMuMuBaseline_lMuIso_lPt.Write()
hMuMuBaseline_tMuIso_tPt.Write()
hMuMuBaseline_jPt_dRjlMu.Write()
hMuMuBaseline_jPt_dRjtMu.Write()

hMuMuTrig_M.Write()
hMuMuTrig_dR.Write()
hMuMuTrig_lPt.Write()
hMuMuTrig_tPt.Write()
hMuMuTrig_lEta.Write()
hMuMuTrig_tEta.Write()
hMuMuTrig_dRjlMu.Write()
hMuMuTrig_dRjtMu.Write()

hMuMu_M.Write()
hMuMu_lPt.Write()
hMuMu_tPt.Write()
hMuMu_lEta.Write()
hMuMu_tEta.Write()
hMuMu_dR.Write()

hMuMu_M_genMatched.Write()
hMuMu_lPt_genMatched.Write()
hMuMu_tPt_genMatched.Write()
hMuMu_dR_genMatched.Write()

hEMuBaseline_tMuPt.Write()
hEMuBaseline_NJets.Write()
hEMuBaseline_NbJets.Write()

hEMuBaseline_M.Write()
hEMuBaseline_muPt.Write()
hEMuBaseline_ePt.Write()
hEMuBaseline_muEta.Write()
hEMuBaseline_eEta.Write()
hEMuBaseline_dR.Write()
hEMuBaseline_dRjMu.Write()
hEMuBaseline_dRjE.Write()
hEMuBaseline_M_dR.Write()
hEMuBaseline_M_dRjE.Write()
hEMuBaseline_M_dRjMu.Write()
hEMuBaseline_dRjE_dRjMu.Write()
hEMuBaseline_eIso_dRjE.Write()
hEMuBaseline_muIso_dRjMu.Write()
hEMuBaseline_eIso_ePt.Write()
hEMuBaseline_muIso_muPt.Write()
hEMuBaseline_jPt_dRjE.Write()
hEMuBaseline_jPt_dRjMu.Write()

hEMuTrig_M.Write()
hEMuTrig_muPt.Write()
hEMuTrig_ePt.Write()
hEMuTrig_muEta.Write()
hEMuTrig_eEta.Write()
hEMuTrig_dR.Write()
hEMuTrig_dRjMu.Write()
hEMuTrig_dRjE.Write()

hEMu_M.Write()
hEMu_muPt.Write()
hEMu_ePt.Write()
hEMu_muEta.Write()
hEMu_eEta.Write()
hEMu_dR.Write()

hEMu_M_genMatched.Write()
hEMu_muPt_genMatched.Write()
hEMu_ePt_genMatched.Write()
hEMu_dR_genMatched.Write()

out.Close()    
