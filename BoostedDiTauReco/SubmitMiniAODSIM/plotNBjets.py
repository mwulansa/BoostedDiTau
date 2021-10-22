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
outputFileName = outputFileDir+"h_EMu_NBjets_"+inputFileListName.split("/")[-1].replace(".txt",".root")
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

h['hNbJets'] = ROOT.TH1F ("hNbJets", "", 10, 0, 10)

h['hJetPt'] = ROOT.TH1F ("hJetPt", "Jet Pt;;N_{events}", 2000, 0, 2000)

h['bDiscriminator'] = ROOT.TH1F ("hEMu_bDiscriminator_CSV","", 100, 0, 1)
h['bDiscriminator_flav'] = ROOT.TH1F ("hEMu_bDiscriminator_flav","", 100, 0, 1)

h['hEMuBaseline_NbJets_tight_flav'] = ROOT.TH1F ("hEMu_Baseline_NbJets_tight_flav","", 10, 0, 10)
h['hEMuBaseline_NbJets_medium_flav'] = ROOT.TH1F ("hEMu_Baseline_NbJets_medium_flav","", 10, 0, 10)
h['hEMuBaseline_NbJets_loose_flav'] = ROOT.TH1F ("hEMu_Baseline_NbJets_loose_flav","", 10, 0, 10)

h['hEMuBaseline_M_tight_flav'] = ROOT.TH1F ("hEMu_Baseline_M_tight_flav",'', 1000, 0, 200)
h['hEMuBaseline_M_medium_flav'] = ROOT.TH1F ("hEMu_Baseline_M_medium_flav",'', 1000, 0, 200)
h['hEMuBaseline_M_loose_flav'] = ROOT.TH1F ("hEMu_Baseline_M_loose_flav",'', 1000, 0, 200)

h['hEMuBaseline_NbJets_tight'] = ROOT.TH1F ("hEMu_Baseline_NbJets_tight","", 10, 0, 10)
h['hEMuBaseline_NbJets_medium'] = ROOT.TH1F ("hEMu_Baseline_NbJets_medium","", 10, 0, 10)
h['hEMuBaseline_NbJets_loose'] = ROOT.TH1F ("hEMu_Baseline_NbJets_loose","", 10, 0, 10)

h['hEMuBaseline_M_tight'] = ROOT.TH1F ("hEMu_Baseline_M_tight",'', 1000, 0, 200)
h['hEMuBaseline_M_medium'] = ROOT.TH1F ("hEMu_Baseline_M_medium",'', 1000, 0, 200)
h['hEMuBaseline_M_loose'] = ROOT.TH1F ("hEMu_Baseline_M_loose",'', 1000, 0, 200)


for key in h.keys():
    h[key].Sumw2()

#</histograms>

inputFileNames=open(inputFileList, 'r')

prefix = "root://cmseos.fnal.gov//eos/uscms/store/user/nbower/Events/TCP_m_50_w_1_htj_400toInf_slc6_amd64_gcc630_MINIAOD/"

for inputFileName in inputFileNames:
    print inputFileName.replace("\n","")

    inputFileName = prefix+inputFileName.replace("\n","")

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
        if len(handleRho.product())>0:
            rho=handleRho.product()[0]
        else:
            rho = 0

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
#            if abs(muon.innerTrack().dxy(vertex[0].position()))>0.2 or abs(muon.innerTrack().dz(vertex[0].position()))>0.5: continue 
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
            E_c = electron.superCluster().energy()
            if electron.isEB(): 
                if electron.full5x5_sigmaIetaIeta()<0.0112 \
                and electron.hadronicOverEm()<(0.05 + 1.16/E_c + 0.0324*rho/E_c) \
                and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.0884 \
                and GsfEleEInverseMinusPInverse(electron)<0.193 \
                and abs(dEtaInSeed(electron))<0.00377 \
                and GsfEleMissingHitsCut(electron)<=1 \
                and electron.passConversionVeto() \
                and GsfEleEffAreaPFIsoCut(electron, rho) < 0.112+0.506/electron.pt():
                    selected_electrons+=[electron]
            if electron.isEE():
                if electron.full5x5_sigmaIetaIeta()<0.0425 \
                and electron.hadronicOverEm()<(0.0441 + 2.54/E_c + 0.183*rho/E_c) \
                and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.169 \
                and GsfEleEInverseMinusPInverse(electron)<0.111 \
                and abs(dEtaInSeed(electron))<0.00674 \
                and GsfEleMissingHitsCut(electron)<=1 \
                and electron.passConversionVeto() \
                and GsfEleEffAreaPFIsoCut(electron, rho)<0.108+0.963/electron.pt():
                    selected_electrons+=[electron]

        selected_electrons.sort(key=lambda x: x.pt(), reverse=True) 
        #<\electronSelection>

        #<jetSelection>
        selected_jets=[]
        selected_bjets=[[],[],[]]
        selected_bjets_flav=[[],[],[]]
        for jet in jets:
            if jet.pt()<20 or abs(jet.eta())>=5: continue
            NHF  = jet.neutralHadronEnergyFraction() 
            NEMF = jet.neutralEmEnergyFraction() 
            CHF  = jet.chargedHadronEnergyFraction() 
            MUF  = jet.muonEnergyFraction() 
            CEMF = jet.chargedEmEnergyFraction()
            NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity()
            NumNeutralParticles =jet.neutralMultiplicity()
            NumChargedMult = jet.chargedMultiplicity()
            CHM      = jet.chargedMultiplicity()


            score = 0
            scoreflav = 0

            if abs(jet.eta())<=2.4:
                if NHF < 0.9 and NEMF < 0.9 and NumConst > 1.0 and MUF < 0.8 and CHM > 0.0 and NumChargedMult > 0.0 and CEMF < 0.8: 
                    selected_jets+=[jet]
                    if jet.bDiscriminator("pfDeepCSVJetTags:probb") >= 0 and jet.bDiscriminator("pfDeepCSVJetTags:probbb") >= 0:
                        h['bDiscriminator'].Fill(jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb"), genweight)
                        score = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb")
                        h['hNbJets'].Fill(1)
                        if score > 0.7738 : selected_bjets[0]+=[jet]
                        if score > 0.4506 : selected_bjets[1]+=[jet]
                        if score > 0.1355 : selected_bjets[2]+=[jet]

                    if jet.bDiscriminator("pfDeepFlavourJetTags:probb") >= 0 and jet.bDiscriminator("pfDeepFlavourJetTags:probbb") >= 0 and jet.bDiscriminator("pfDeepFlavourJetTags:problepb") >= 0:
                        scoreflav = jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb") + jet.bDiscriminator("pfDeepFlavourJetTags:problepb") 

                        h['bDiscriminator_flav'].Fill(scoreflav, genweight)
                        if scoreflav > 0.7476 : selected_bjets_flav[0]+=[jet]
                        if scoreflav > 0.3040 : selected_bjets_flav[1]+=[jet]
                        if scoreflav > 0.0532 : selected_bjets_flav[2]+=[jet]

            if abs(jet.eta())>2.4 and abs(jet.eta())<=2.7:
                if NHF < 0.9 and NEMF < 0.99:
                    selected_jets+=[jet]
                    if jet.bDiscriminator("pfDeepCSVJetTags:probb") >= 0 and jet.bDiscriminator("pfDeepCSVJetTags:probbb") >= 0: 
                        h['hNbJets'].Fill(1)
                        h['bDiscriminator'].Fill(jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb"), genweight)
                        score = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb")
                        if score > 0.7738 : selected_bjets[0]+=[jet]
                        if score > 0.4506 : selected_bjets[1]+=[jet]
                        if score > 0.1355 : selected_bjets[2]+=[jet]

                    if jet.bDiscriminator("pfDeepFlavourJetTags:probb") >= 0 and jet.bDiscriminator("pfDeepFlavourJetTags:probbb") >= 0 and jet.bDiscriminator("pfDeepFlavourJetTags:problepb") >= 0:
                        scoreflav = jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb") + jet.bDiscriminator("pfDeepFlavourJetTags:problepb")

                        h['bDiscriminator_flav'].Fill(scoreflav, genweight)
                        if scoreflav > 0.7476 : selected_bjets_flav[0]+=[jet]
                        if scoreflav > 0.3040 : selected_bjets_flav[1]+=[jet]
                        if scoreflav > 0.0532 : selected_bjets_flav[2]+=[jet]


            if abs(jet.eta())>2.7 and abs(jet.eta())<=3 :
                if NHF < 0.9 and NEMF > 0.0 and NEMF < 0.99 and NumNeutralParticles > 1.0:
                    selected_jets+=[jet]

            if abs(jet.eta())>3 and abs(jet.eta())<=5: 
                if NHF > 0.2 and NEMF < 0.9 and NumNeutralParticles > 10:
                    selected_jets+=[jet]

                
        selected_jets.sort(key=lambda x: x.pt(), reverse=True)       

#----------------------------- tau_e tau_mu selection ------------------------------

        if len(selected_muons)>=1 and len(selected_electrons)>=1 and len(selected_jets)>=1 and selected_muons[0].charge()*selected_electrons[0].charge()<0:
            h['hEMuBaseline_NbJets_tight'].Fill(len(selected_bjets[0]), genweight)
            h['hEMuBaseline_NbJets_medium'].Fill(len(selected_bjets[1]), genweight)
            h['hEMuBaseline_NbJets_loose'].Fill(len(selected_bjets[2]), genweight)

            h['hEMuBaseline_NbJets_tight_flav'].Fill(len(selected_bjets_flav[0]), genweight)
            h['hEMuBaseline_NbJets_medium_flav'].Fill(len(selected_bjets_flav[1]), genweight)
            h['hEMuBaseline_NbJets_loose_flav'].Fill(len(selected_bjets_flav[2]), genweight)

            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())

            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())


            if len(selected_bjets[0]) == 0 : h['hEMuBaseline_M_tight'].Fill((e+mu).M(), genweight)
            if len(selected_bjets[1]) == 0 : h['hEMuBaseline_M_medium'].Fill((e+mu).M(), genweight)
            if len(selected_bjets[2]) == 0 : h['hEMuBaseline_M_loose'].Fill((e+mu).M(), genweight)

            if len(selected_bjets_flav[0]) == 0 : h['hEMuBaseline_M_tight_flav'].Fill((e+mu).M(), genweight)
            if len(selected_bjets_flav[1]) == 0 : h['hEMuBaseline_M_medium_flav'].Fill((e+mu).M(), genweight)
            if len(selected_bjets_flav[2]) == 0 : h['hEMuBaseline_M_loose_flav'].Fill((e+mu).M(), genweight)

                        
out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
