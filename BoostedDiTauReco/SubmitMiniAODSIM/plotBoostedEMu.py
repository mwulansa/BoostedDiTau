#Code for tau_mu tau_mu studies for 94X

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
outputFileName = outputFileDir+"h_EMu_OS_"+inputFileListName.split("/")[-1].replace(".txt",".root")
print outputFileName

pi = math.pi

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

# <histograms>

pi = math.pi

h = {}

h['hNEvent'] = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)
h['hEMuBaseline_NbJets'] = ROOT.TH1F ("hEMu_Baseline_NbJets","", 10, 0, 10)

#-----------------------EMu-----------------------                                                                                                                                                         
#-----Require 2 muons                                                                                                                                                                                          
h['hEMuBaseline_M'] = ROOT.TH1F ("hEMu_Baseline_M", "e - #mu mass;M_{e#mu};", 1000, 0, 200)
h['hEMuBaseline_muPt'] = ROOT.TH1F ("hEMu_Baseline_muPt", "#mu P_t; P_t;", 500, 0, 500)
h['hEMuBaseline_ePt'] = ROOT.TH1F ("hEMu_Baseline_ePt", "e P_t; P_t;", 500, 0, 500)
h['hEMuBaseline_jPt'] = ROOT.TH1F ("hEMu_Baseline_jPt", "jet pt;P_{t};", 2000, 0, 2000)
h['hEMuBaseline_dR'] = ROOT.TH1F ("hEMu_Baseline_dR", ";#Delta R;", 100, 0, 5)
h['hEMuBaseline_dRlj'] = ROOT.TH1F ("hEMu_Baseline_dRlj", ";#Delta R;", 100, 0, 5)
h['hEMuBaseline_METPt'] = ROOT.TH1F ("hEMu_Baseline_METPt", ";p_{T};", 500, 0, 500)
h['hEMuBaseline_dPhiMj'] = ROOT.TH1F ("hEMu_Baseline_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
h['hEMuBaseline_dPhiMmu'] = ROOT.TH1F ("hEMu_Baseline_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

#-----Passing Triggers                                                                                                                                                                       
h['hEMuTrig_M'] = ROOT.TH1F ("hEMu_Trig_M", "e - #mu mass;M_{e#mu};", 1000, 0, 200)
h['hEMuTrig_muPt'] = ROOT.TH1F ("hEMu_Trig_muPt", "#mu P_t;P_{t};", 500, 0, 500)
h['hEMuTrig_ePt'] = ROOT.TH1F ("hEMu_Trig_ePt", "e P_t;P_{t};", 500, 0, 500)
h['hEMuTrig_jPt'] = ROOT.TH1F ("hEMu_Trig_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
h['hEMuTrig_dR'] = ROOT.TH1F ("hEMu_Trig_dR", ";#Delta R;", 100, 0, 5)
h['hEMuTrig_dRlj'] = ROOT.TH1F ("hEMu_Trig_dRlj", ";#Delta R;", 100, 0, 5)
h['hEMuTrig_METPt'] = ROOT.TH1F ("hEMu_Trig_METPt", ";p_{T};", 500, 0, 500)
h['hEMuTrig_dPhiMj'] = ROOT.TH1F ("hEMu_Trig_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
h['hEMuTrig_dPhiMmu'] = ROOT.TH1F ("hEMu_Trig_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

#-----dR cuts                                                                                                                                                                                 
h['hEMudR_M'] = ROOT.TH1F ("hEMu_dR_M", "e - #mu mass;M_{e#mu};", 1000, 0, 200)
h['hEMudR_muPt'] = ROOT.TH1F ("hEMu_dR_muPt", "#mu P_t;P_{t};", 500, 0, 500)
h['hEMudR_ePt'] = ROOT.TH1F ("hEMu_dR_ePt", "e P_t;P_{t};", 500, 0, 500)
h['hEMudR_jPt'] = ROOT.TH1F ("hEMu_dR_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
h['hEMudR_dR'] = ROOT.TH1F ("hEMu_dR_dR", ";#Delta R;", 100, 0, 5)
h['hEMudR_dRlj'] = ROOT.TH1F ("hEMu_dR_dRlj", ";#Delta R;", 100, 0, 5)
h['hEMudR_METPt'] = ROOT.TH1F ("hEMu_dR_METPt", ";p_{T};", 500, 0, 500)
h['hEMudR_dPhiMj'] = ROOT.TH1F ("hEMu_dR_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
h['hEMudR_dPhiMmu'] = ROOT.TH1F ("hEMu_dR_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

#-----MET cuts                                                                                                                                                                                               
h['hEMuMetcut_M'] = ROOT.TH1F ("hEMu_Metcut_M", "e - #mu mass;M_{e#mu};", 1000, 0, 200)
h['hEMuMetcut_muPt'] = ROOT.TH1F ("hEMu_Metcut_muPt", "#mu P_t;P_{t};", 500, 0, 500)
h['hEMuMetcut_ePt'] = ROOT.TH1F ("hEMu_Metcut_ePt", "e P_t;P_{t};", 500, 0, 500)
h['hEMuMetcut_jPt'] = ROOT.TH1F ("hEMu_Metcut_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
h['hEMuMetcut_dR'] = ROOT.TH1F ("hEMu_Metcut_dR", ";#Delta R;", 100, 0, 5)
h['hEMuMetcut_dRlj'] = ROOT.TH1F ("hEMu_Metcut_dRlj", ";#Delta R;", 100, 0, 5)
h['hEMuMetcut_METPt'] = ROOT.TH1F ("hEMu_Metcut_METPt", ";p_{T};", 500, 0, 500)
h['hEMuMetcut_dPhiMj'] = ROOT.TH1F ("hEMu_Metcut_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
h['hEMuMetcut_dPhiMmu'] = ROOT.TH1F ("hEMu_Metcut_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

#-----dPhi cuts           
h['hEMudPhi_M'] = ROOT.TH1F ("hEMu_dPhi_M", "e - #mu mass;M_{e#mu};", 1000, 0, 200)
h['hEMudPhi_muPt'] = ROOT.TH1F ("hEMu_dPhi_muPt", "#mu P_t;P_{t};", 500, 0, 500)
h['hEMudPhi_ePt'] = ROOT.TH1F ("hEMu_dPhi_ePt", "e P_t;P_{t};", 500, 0, 500)
h['hEMudPhi_jPt'] = ROOT.TH1F ("hEMu_dPhi_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
h['hEMudPhi_dR'] = ROOT.TH1F ("hEMu_dPhi_dR", ";#Delta R;", 100, 0, 5)
h['hEMudPhi_dRlj'] = ROOT.TH1F ("hEMu_dPhi_dRlj", ";#Delta R;", 100, 0, 5)
h['hEMudPhi_METPt'] = ROOT.TH1F ("hEMu_dPhi_METPt", ";p_{T};", 500, 0, 500)
h['hEMudPhi_dPhiMj'] = ROOT.TH1F ("hEMu_dPhi_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
h['hEMudPhi_dPhiMmu'] = ROOT.TH1F ("hEMu_dPhi_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

#-----genMatching
h['hEMuGen_M'] = ROOT.TH1F ("hEMu_Gen_M", "e - #mu mass;M_{e#mu};", 1000, 0, 200)
h['hEMuGen_ePt'] = ROOT.TH1F ("hEMu_Gen_ePt", "e P_t;P_{t};", 500, 0, 500)
h['hEMuGen_muPt'] = ROOT.TH1F ("hEMu_Gen_muPt", "#mu P_t;P_{t};", 500, 0, 500)

#</histograms>

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
        mets=[]
        mets+=[met]
        mets.sort(key=lambda x: x.pt(), reverse=True)

        m=ROOT.TLorentzVector()
        m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())
        
        h['hNEvent'].Fill(0.5, 1)
        h['hNEvent'].Fill(1.5, genweight)

#<genInfo>
        
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

#</genInfo>

        #<muonSelection>
        selected_muons=[]
        for muon in muons:
            if not muon.isLooseMuon(): continue
            if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue 
            if muon.pt()<3 or muon.eta()>2.4: continue
            if muonIsoCut(muon)<0.25:
                selected_muons+=[muon]

        selected_muons.sort(key=lambda x: x.pt(), reverse=True) 
        #<\muonSelection>

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
                and GsfEleEffAreaPFIsoCut(electron, rho)<0.0994 \
                and abs(electron.gsfTrack().dz(vertex[0].position()))<0.1 \
                and abs(electron.gsfTrack().dxy(vertex[0].position()))<0.05:
                    selected_electrons+=[electron]
            if electron.isEE():
                if electron.full5x5_sigmaIetaIeta()<0.0314 \
                and electron.hadronicOverEm()<0.101 \
                and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.213 \
                and GsfEleEInverseMinusPInverse(electron)<0.14 \
                and abs(dEtaInSeed(electron))<0.00868 \
                and GsfEleMissingHitsCut(electron)<=1 \
                and electron.passConversionVeto() \
                and GsfEleEffAreaPFIsoCut(electron, rho)<0.107 \
                and abs(electron.gsfTrack().dz(vertex[0].position()))<0.2 \
                and abs(electron.gsfTrack().dxy(vertex[0].position()))<0.1:
                    selected_electrons+=[electron]

        selected_electrons.sort(key=lambda x: x.pt(), reverse=True) 
        #<\electronSelection>

        #<jetSelection>
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

#----------------------------- tau_e tau_mu selection ------------------------------

        if len(selected_muons)>=1 and len(selected_electrons)>=1 and len(selected_jets)>=1 and selected_muons[0].charge()*selected_electrons[0].charge()<0:
            h['hEMuBaseline_NbJets'].Fill(len(selected_bjets), genweight)

        if len(selected_muons)>=1 and len(selected_electrons)>=1 and len(selected_jets)>=1 and selected_muons[0].charge()*selected_electrons[0].charge()<0 and len(selected_bjets)==0:

            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass()) 
            
            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

            h['hEMuBaseline_M'].Fill((mu+e).M(), genweight)
            h['hEMuBaseline_muPt'].Fill(mu.Pt(), genweight)
            h['hEMuBaseline_ePt'].Fill(e.Pt(), genweight)
            h['hEMuBaseline_jPt'].Fill(j.Pt(), genweight)
            h['hEMuBaseline_dR'].Fill(mu.DeltaR(e), genweight)
            h['hEMuBaseline_dRlj'].Fill((mu+e).DeltaR(j), genweight)
            h['hEMuBaseline_METPt'].Fill(m.Pt(), genweight)
            h['hEMuBaseline_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            h['hEMuBaseline_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)
                
            if mu.DeltaR(e)<0.4 and mu.DeltaR(j)>0.8 and e.DeltaR(j)>0.8:

                h['hEMudR_M'].Fill((mu+e).M(), genweight)
                h['hEMudR_muPt'].Fill(mu.Pt(), genweight)
                h['hEMudR_ePt'].Fill(e.Pt(), genweight)
                h['hEMudR_jPt'].Fill(j.Pt(), genweight)
                h['hEMudR_dR'].Fill(mu.DeltaR(e), genweight)
                h['hEMudR_dRlj'].Fill((mu+e).DeltaR(j), genweight)
                h['hEMudR_METPt'].Fill(m.Pt(), genweight)
                h['hEMudR_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                h['hEMudR_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

                if m.Pt()>100:

                    h['hEMuMetcut_M'].Fill((mu+e).M(), genweight)
                    h['hEMuMetcut_muPt'].Fill(mu.Pt(), genweight)
                    h['hEMuMetcut_ePt'].Fill(e.Pt(), genweight)
                    h['hEMuMetcut_jPt'].Fill(j.Pt(), genweight)
                    h['hEMuMetcut_dR'].Fill(mu.DeltaR(e), genweight)
                    h['hEMuMetcut_dRlj'].Fill((mu+e).DeltaR(j), genweight)
                    h['hEMuMetcut_METPt'].Fill(m.Pt(), genweight)
                    h['hEMuMetcut_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                    h['hEMuMetcut_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

                    if abs(m.DeltaPhi(mu))<1 and abs(m.DeltaPhi(j))>2:
                        h['hEMudPhi_M'].Fill((mu+e).M(), genweight)
                        h['hEMudPhi_muPt'].Fill(mu.Pt(), genweight)
                        h['hEMudPhi_ePt'].Fill(e.Pt(), genweight)
                        h['hEMudPhi_jPt'].Fill(j.Pt(), genweight)
                        h['hEMudPhi_dR'].Fill(mu.DeltaR(e), genweight)
                        h['hEMudPhi_dRlj'].Fill((mu+e).DeltaR(j), genweight)
                        h['hEMudPhi_METPt'].Fill(m.Pt(), genweight)
                        h['hEMudPhi_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                        h['hEMudPhi_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

                        if len(genMuons)==1 and len(genJets)>=1 and len(genElectrons)==1:
                            genMu=ROOT.TLorentzVector()
                            genMu.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

                            genE=ROOT.TLorentzVector()
                            genE.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())

                            if mu.DeltaR(genMu)<0.3 and e.DeltaR(genE)<0.3:
                                h['hEMuGen_M'].Fill((mu+e).M(), genweight)
                                h['hEMuGen_ePt'].Fill(mu.Pt(), genweight)
                                h['hEMuGen_muPt'].Fill(e.Pt(), genweight)

                        if (mu.Pt()>26 and (triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")))) \
        or (mu.Pt()>50 and (triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")))) \
        or (j.Pt()>500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))):

                            h['hEMuTrig_M'].Fill((mu+e).M(), genweight)
                            h['hEMuTrig_muPt'].Fill(mu.Pt(), genweight)
                            h['hEMuTrig_ePt'].Fill(e.Pt(), genweight)
                            h['hEMuTrig_jPt'].Fill(j.Pt(), genweight)
                            h['hEMuTrig_dR'].Fill(mu.DeltaR(e), genweight)
                            h['hEMuTrig_dRlj'].Fill((mu+e).DeltaR(j), genweight)
                            h['hEMuTrig_METPt'].Fill(m.Pt(), genweight)
                            h['hEMuTrig_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                            h['hEMuTrig_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)


                        
out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
