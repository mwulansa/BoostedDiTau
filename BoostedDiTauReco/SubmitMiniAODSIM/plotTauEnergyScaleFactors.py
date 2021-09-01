#Code for generating up and down histograms for tau energy scale factor

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

outputFileName = outputFileDir+"h_TauScaleFactors_"+inputFileListName.split("/")[-1].replace(".txt",".root")

print outputFileName

pi = math.pi

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlineSlimmedPrimaryVertices')

handleBoostedTau = Handle ('vector<pat::Tau>')
#labelBoostedTau = ('slimmedTaus')
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
h['hNTau'] = ROOT.TH1F ("hNTau","Number of Events;;N_{events}", 5, 0, 5)
h['hNTauPair'] = ROOT.TH1F ("hNTauPair","Number of Events;;N_{events}", 5, 0, 5)

h['hMuPt'] = ROOT.TH1F ("hMu_Pt", ";p_{T};", 500, 0, 500)
h['hElePt'] = ROOT.TH1F ("hEle_Pt", ";p_{T};", 500, 0, 500)
h['hJetPt'] = ROOT.TH1F ("hJet_Pt", ";p_{T};", 2000, 0, 2000)
h['hTauPt'] = ROOT.TH1F ("hTau_Pt", ";p_{T};", 500, 0, 500)

h['hTauTauBaseline_M'] = ROOT.TH1F ("hTauTau_Baseline_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTaudR_M'] = ROOT.TH1F ("hTauTau_dR_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hTauTauMetcut_M'] = ROOT.TH1F ("hTauTau_Metcut_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['hTauTauTrig_M_DM0'] = ROOT.TH1F ('hTauTauTrig_M_DM0', '#tau - #tau Vis. Mass; M_{vis.};Events', 500, 0, 200)
h['hTauTauTrig_M_DM1'] = ROOT.TH1F ('hTauTauTrig_M_DM1', '#tau - #tau Vis. Mass; M_{vis.};Events', 500, 0, 200)
h['hTauTauTrig_M_DM10'] = ROOT.TH1F ('hTauTauTrig_M_DM10', '#tau - #tau Vis. Mass; M_{vis.};Events', 500, 0, 200)
h['hTauTauTrig_M_DM0_down'] = ROOT.TH1F ('hTauTauTrig_M_DM0_down', '#tau - #tau Vis.Mass down; M_{vis.};Events', 500, 0, 200)
h['hTauTauTrig_M_DM0_up'] = ROOT.TH1F ('hTauTauTrig_M_DM0_up', '#tau - #tau Vis.Mass up; M_{vis.};Events', 500, 0, 200)
h['hTauTauTrig_M_DM1_down'] = ROOT.TH1F ('hTauTauTrig_M_DM1_down', '#tau - #tau Vis.Mass down; M_{vis.};Events', 500, 0, 200)
h['hTauTauTrig_M_DM1_up'] = ROOT.TH1F ('hTauTauTrig_M_DM1_up', '#tau - #tau Vis.Mass up; M_{vis.};Events', 500, 0, 200)
h['hTauTauTrig_M_DM10_down'] = ROOT.TH1F ('hTauTauTrig_M_DM10_down', '#tau - #tau Vis.Mass down; M_{vis.};Events', 500, 0, 200)
h['hTauTauTrig_M_DM10_up'] = ROOT.TH1F ('hTauTauTrig_M_DM10_up', '#tau - #tau Vis.Mass up; M_{vis.};Events', 500, 0, 200)

# h['hTauTauTrig_lPt_DM0'] = ROOT.TH1F ('hTauTauTrig_lPt_DM0', '#tau - #tau p_{T}; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM1'] = ROOT.TH1F ('hTauTauTrig_lPt_DM1', '#tau - #tau p_{T}; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM2'] = ROOT.TH1F ('hTauTauTrig_lPt_DM2', '#tau - #tau p_{T}; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM10'] = ROOT.TH1F ('hTauTauTrig_lPt_DM10', '#tau - #tau p_{T}; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM0_down'] = ROOT.TH1F ('hTauTauTrig_lPt_DM0_down', '#tau - #tau p_{T} down; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM0_up'] = ROOT.TH1F ('hTauTauTrig_lPt_DM0_up', '#tau - #tau p_{T} up; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM1_down'] = ROOT.TH1F ('hTauTauTrig_lPt_DM1_down', '#tau - #tau p_{T} down; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM1_up'] = ROOT.TH1F ('hTauTauTrig_lPt_DM1_up', '#tau - #tau p_{T} up; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM2_down'] = ROOT.TH1F ('hTauTauTrig_lPt_DM2_down', '#tau - #tau p_{T} down; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM2_up'] = ROOT.TH1F ('hTauTauTrig_lPt_DM2_up', '#tau - #tau p_{T} up; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM10_down'] = ROOT.TH1F ('hTauTauTrig_lPt_DM10_down', '#tau - #tau p_{T} down; p_{T}(#tau);Events', 500, 0, 500)
# h['hTauTauTrig_lPt_DM10_up'] = ROOT.TH1F ('hTauTauTrig_lPt_DM10_up', '#tau - #tau p_{T} up; p_{T}(#tau);Events', 500, 0, 500)


for key in h.keys():
    h[key].Sumw2()

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
                    selected_electrons+=[electron]

        if len(selected_electrons)>0:
            selected_electrons.sort(key=lambda x: x.pt(), reverse=True) 
            h['hElePt'].Fill(selected_electrons[0].pt(), genweight)

        #<\electronSelection>

        selected_btaus_0 = [[],[],[]]
        selected_btaus_1 = [[],[],[]]
        selected_btaus_10 = [[],[],[]]

        btaus_0 = [0,0,0]
        btaus_1 = [0,0,0]
        btaus_10 = [0,0,0]

        DM0 = [0.984, 0.994, 1.004]
        DM1 = [0.986, 0.995, 1.004]
        DM10 = [0.989, 1, 1.011]
        

        for tau in btaus:
            if abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") and tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"):
                h['hNTau'].Fill(0.5)
                if tau.decayMode() == 0 :
                    h['hNTau'].Fill(1.5)
                    h['hNTau'].Fill(4.5)
                    for i in range(3):
                        btaus_0[i]=(tau.pt()*DM0[i])
                        if btaus_0[i] > 20 : 
                            selected_btaus_0[i]+=[tau]
                if tau.decayMode() == 1 :
                    h['hNTau'].Fill(2.5)
                    h['hNTau'].Fill(4.5)
                    for i in range(3):
                        btaus_1[i]=(tau.pt()*DM1[i])
                        if btaus_1[i] > 20 : 
                            selected_btaus_1[i]+=[tau]
                if tau.decayMode() == 10 :
                    h['hNTau'].Fill(3.5)
                    h['hNTau'].Fill(4.5)
                    for i in range(3):
                        btaus_10[i]=(tau.pt()*DM10[i])
                        if btaus_10[i] > 20 : 
                            selected_btaus_10[i]+=[tau]

        for j in range(3):
            if len(selected_btaus_0[j])>0 : selected_btaus_0[j].sort(key=lambda x: x.pt(), reverse=True)
            if len(selected_btaus_1[j])>0 : selected_btaus_1[j].sort(key=lambda x: x.pt(), reverse=True)
            if len(selected_btaus_10[j])>0 : selected_btaus_10[j].sort(key=lambda x: x.pt(), reverse=True)


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

#<Taus>

        hDM0 = ['hTauTauTrig_M_DM0_down', 'hTauTauTrig_M_DM0', 'hTauTauTrig_M_DM0_up']
        hDM1 = ['hTauTauTrig_M_DM1_down', 'hTauTauTrig_M_DM1', 'hTauTauTrig_M_DM1_up']
        hDM10 = ['hTauTauTrig_M_DM10_down', 'hTauTauTrig_M_DM10', 'hTauTauTrig_M_DM10_up']


        for k in range(3):
            if len(selected_btaus_0[k])>=2 and len(selected_jets)>=1 and len(selected_muons)==0 and len(selected_electrons)==0 and len(selected_bjets)==0 and selected_btaus_0[k][0].charge()*selected_btaus_0[k][1].charge()<0:

                btau1 = ROOT.TLorentzVector()
                btau1.SetPtEtaPhiM(selected_btaus_0[k][0].pt()*DM0[k], selected_btaus_0[k][0].eta(), selected_btaus_0[k][0].phi(), selected_btaus_0[k][0].mass()*DM0[k])

                btau2 = ROOT.TLorentzVector()
                btau2.SetPtEtaPhiM(selected_btaus_0[k][1].pt()*DM0[k], selected_btaus_0[k][1].eta(), selected_btaus_0[k][1].phi(), selected_btaus_0[k][1].mass()*DM0[k])

                jet = ROOT.TLorentzVector()
                jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

                if btau1.DeltaR(btau2)<0.4 and btau1.DeltaR(jet)>0.8 and btau2.DeltaR(jet)>0.8:
                    m = ROOT.TLorentzVector()
                    m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

                    if m.Pt()>100:
                        if jet.Pt()>500 \
                        and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) \
                        or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) \
                        or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5"))):

                           h[hDM0[k]].Fill((btau1+btau2).M(), genweight)
                            
            if len(selected_btaus_1[k])>=2 and len(selected_jets)>=1 and len(selected_muons)==0 and len(selected_electrons)==0 and len(selected_bjets)==0 and selected_btaus_1[k][0].charge()*selected_btaus_1[k][1].charge()<0:

                btau1 = ROOT.TLorentzVector()
                btau1.SetPtEtaPhiM(selected_btaus_1[k][0].pt()*DM1[k], selected_btaus_1[k][0].eta(), selected_btaus_1[k][0].phi(), selected_btaus_1[k][0].mass()*DM1[k])

                btau2 = ROOT.TLorentzVector()
                btau2.SetPtEtaPhiM(selected_btaus_1[k][1].pt()*DM1[k], selected_btaus_1[k][1].eta(), selected_btaus_1[k][1].phi(), selected_btaus_1[k][1].mass()*DM1[k])

                jet = ROOT.TLorentzVector()
                jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

#                h[hDM1[k]].Fill((btau1+btau2).M(), genweight)

                if btau1.DeltaR(btau2)<0.4 and btau1.DeltaR(jet)>0.8 and btau2.DeltaR(jet)>0.8:
                    m = ROOT.TLorentzVector()
                    m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

                    if m.Pt()>100:
                        if jet.Pt()>500 \
                        and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) \
                        or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) \
                        or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5"))):

                           h[hDM1[k]].Fill((btau1+btau2).M(), genweight)
                            

            if len(selected_btaus_10[k])>=2 and len(selected_jets)>=1 and len(selected_muons)==0 and len(selected_electrons)==0 and len(selected_bjets)==0 and selected_btaus_10[k][0].charge()*selected_btaus_10[k][1].charge()<0:

                btau1 = ROOT.TLorentzVector()
                btau1.SetPtEtaPhiM(selected_btaus_10[k][0].pt()*DM10[k], selected_btaus_10[k][0].eta(), selected_btaus_10[k][0].phi(), selected_btaus_10[k][0].mass()*DM10[k])

                btau2 = ROOT.TLorentzVector()
                btau2.SetPtEtaPhiM(selected_btaus_10[k][1].pt()*DM10[k], selected_btaus_10[k][1].eta(), selected_btaus_10[k][1].phi(), selected_btaus_10[k][1].mass()*DM10[k])

                jet = ROOT.TLorentzVector()
                jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

#                h[hDM10[k]].Fill((btau1+btau2).M(), genweight)

                if btau1.DeltaR(btau2)<0.4 and btau1.DeltaR(jet)>0.8 and btau2.DeltaR(jet)>0.8:
                    m = ROOT.TLorentzVector()
                    m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

                    if m.Pt()>100:
                        if jet.Pt()>500 \
                        and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) \
                        or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) \
                        or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5"))):

                           h[hDM10[k]].Fill((btau1+btau2).M(), genweight)
                            


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
