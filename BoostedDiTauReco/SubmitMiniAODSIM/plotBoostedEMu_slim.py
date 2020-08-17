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

outputFileName = outputFileDir+"h_EMu_"+inputFileListName.split("/")[-1].replace(".txt",".root")
print outputFileName

pi = math.pi

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('selectedPATMuons','myMuons','PAT')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('selectedPATElectrons','myElectrons','PAT')

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('selectedPATJets')

handleTauElectronCleaned = Handle ('vector<pat::Tau>')
labelTauElectronCleaned = ('selectedPATTausElectronCleaned')

handleTau = Handle ('vector<pat::Tau>')
labelTau = ('selectedPATTaus')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handleHLT = Handle ('edm::TriggerResults')
#labelHLT = ('TriggerResults','','reMINIAOD')
labelHLT = ('TriggerResults','','PAT')
#labelHLT = ('TriggerResults','','HLT')

handleGenJet = Handle ('vector<reco::GenJet>')
labelGenJet = ('slimmedGenJets')

handlePatMETs = Handle("vector<pat::MET>")
labelPatMETs = ( 'slimmedMETs')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

out=ROOT.TFile.Open(outputFileName,'recreate')

# <histograms>

pi = math.pi

h = {}

h['hNEvent'] = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)
#h['hETauBaseline_NbJets'] = ROOT.TH1F ("hMuTau_Baseline_NbJets","", 10, 0, 10)
#h['hETauBaseline_muPt'] = ROOT.TH1F ("hMuTau_Baseline_muPt", "", 100, 0, 500)

#Objects
#-----------------------Taus------------------------

#h['hTauPt_raw'] = ROOT.TH1F ("hTauPt_raw", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_loose'] = ROOT.TH1F ("hTauPt_loose", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_medium'] = ROOT.TH1F ("hTauPt_medium", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_tight'] = ROOT.TH1F ("hTauPt_tight", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_vloose'] = ROOT.TH1F ("hTauPt_vloose", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_vtight'] = ROOT.TH1F ("hTauPt_vtight", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_vvloose'] = ROOT.TH1F ("hTauPt_vvloose", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_vvtight'] = ROOT.TH1F ("hTauPt_vvtight", "#tau P_t; P_t;", 100, 0, 500)

#-----------------------Muons------------------------

#h['hMuPt'] = ROOT.TH1F ("hMuPt", "#mu P_t; P_t;", 100, 0, 500)
h['hMuPt_iso'] = ROOT.TH1F ("hMuPt_iso", "#mu P_t; P_t;", 100, 0, 500)

#-----------------------Jets-------------------------

h['hJetPt'] = ROOT.TH1F ("hJetPt", "jet pt;P_{t};", 2000, 0, 2000)

#Event
#-----------------------ETau------------------------                                                                                                                                                         
#-----Require 1 Electron, 1 Tau                                                                                                                                                                                          
#h['hIsoMuTauBaseline_M'] = ROOT.TH1F ("hIsoMuTau_Baseline_M", "#tau - e mass;M_{#tau#tau};", 500, 0, 200)

h['hEMuBaseline_M'] = ROOT.TH1F ("hEMu_Baseline_M", "#mu - e mass;M_{e#mu};", 500, 0, 200)
# h['hEMuBaseline_ePt'] = ROOT.TH1F ("hEMu_Baseline_ePt", "e P_t; P_t;", 100, 0, 500)
# h['hEMuBaseline_muPt'] = ROOT.TH1F ("hEMu_Baseline_muPt", "#mu P_t; P_t;", 100, 0, 500)
# h['hEMuBaseline_jPt'] = ROOT.TH1F ("hEMu_Baseline_jPt", "jet pt;P_{t};", 2000, 0, 2000)
# h['hEMuBaseline_dR'] = ROOT.TH1F ("hEMu_Baseline_dR", ";#Delta R;", 100, 0, 5)
# h['hEMuBaseline_dRlj'] = ROOT.TH1F ("hEMu_Baseline_dRlj", ";#Delta R;", 100, 0, 5)
# h['hEMuBaseline_METPt'] = ROOT.TH1F ("hEMu_Baseline_METPt", ";p_{T};", 500, 0, 500)
# h['hEMuBaseline_dPhiMj'] = ROOT.TH1F ("hEMu_Baseline_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hEMuBaseline_dPhiME'] = ROOT.TH1F ("hEMu_Baseline_dPhiME", ";#delta#phi_{e};", 100, -pi, pi)
# h['hEMuBaseline_dPhiMtau'] = ROOT.TH1F ("hEMu_Baseline_dPhiMTau", ";#delta#phi_{#mu};", 100, -pi, pi)

h['hEMuBaseline_M_mod'] = ROOT.TH1F ("hEMu_Baseline_M_mod", "#mu - e mass;M_{e#mu};", 500, 0, 200)
# h['hEMuBaseline_ePt_mod'] = ROOT.TH1F ("hEMu_Baseline_ePt_mod", "e P_t; P_t;", 100, 0, 500)
# h['hEMuBaseline_muPt_mod'] = ROOT.TH1F ("hEMu_Baseline_muPt_mod", "#mu P_t; P_t;", 100, 0, 500)
# h['hEMuBaseline_jPt_mod'] = ROOT.TH1F ("hEMu_Baseline_jPt_mod", "jet pt;P_{t};", 2000, 0, 2000)
# h['hEMuBaseline_dR_mod'] = ROOT.TH1F ("hEMu_Baseline_dR_mod", ";#Delta R;", 100, 0, 5)
# h['hEMuBaseline_dRlj_mod'] = ROOT.TH1F ("hEMu_Baseline_dRlj_mod", ";#Delta R;", 100, 0, 5)
# h['hEMuBaseline_METPt_mod'] = ROOT.TH1F ("hEMu_Baseline_METPt_mod", ";p_{T};", 500, 0, 500)
# h['hEMuBaseline_dPhiMj_mod'] = ROOT.TH1F ("hEMu_Baseline_dPhiMj_mod", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hEMuBaseline_dPhiME_mod'] = ROOT.TH1F ("hEMu_Baseline_dPhiME_mod", ";#delta#phi_{e};", 100, -pi, pi)
# h['hEMuBaseline_dPhiMtau_mod'] = ROOT.TH1F ("hEMu_Baseline_dPhiMTau_mod", ";#delta#phi_{#mu};", 100, -pi, pi)

#-----Passing Triggers                                                                                                                                                                       
# h['hEMuTrig_M'] = ROOT.TH1F ("hEMu_Trig_M", "#mu - #mu mass;M_{#mu#mu};", 1000, 0, 200)
# h['hEMuTrig_ePt'] = ROOT.TH1F ("hEMu_Trig_ePt", "#mu_{1} P_t;P_{t};", 500, 0, 500)
# h['hEMuTrig_muPt'] = ROOT.TH1F ("hEMu_Trig_muPt", "#mu_{2} P_t;P_{t};", 500, 0, 500)
# h['hEMuTrig_jPt'] = ROOT.TH1F ("hEMu_Trig_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hEMuTrig_dR'] = ROOT.TH1F ("hEMu_Trig_dR", ";#Delta R;", 100, 0, 5)
# h['hEMuTrig_dRlj'] = ROOT.TH1F ("hEMu_Trig_dRlj", ";#Delta R;", 100, 0, 5)
# h['hEMuTrig_METPt'] = ROOT.TH1F ("hEMu_Trig_METPt", ";p_{T};", 500, 0, 500)
# h['hEMuTrig_dPhiMj'] = ROOT.TH1F ("hEMu_Trig_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hEMuTrig_dPhiMmu'] = ROOT.TH1F ("hEMu_Trig_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

# #-----dR cuts                                                                                                                                                                                 
h['hEMudR_M'] = ROOT.TH1F ("hEMu_dR_M", "#mu - #mu mass;M_{#mu#mu};", 500, 0, 200)
# h['hEMudR_ePt'] = ROOT.TH1F ("hEMu_dR_ePt", "#mu_{1} P_t;P_{t};", 500, 0, 500)
# h['hEMudR_muPt'] = ROOT.TH1F ("hEMu_dR_muPt", "#mu_{2} P_t;P_{t};", 500, 0, 500)
# h['hEMudR_jPt'] = ROOT.TH1F ("hEMu_dR_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hEMudR_dR'] = ROOT.TH1F ("hEMu_dR_dR", ";#Delta R;", 100, 0, 5)
# h['hEMudR_dRlj'] = ROOT.TH1F ("hEMu_dR_dRlj", ";#Delta R;", 100, 0, 5)
# h['hEMudR_METPt'] = ROOT.TH1F ("hEMu_dR_METPt", ";p_{T};", 500, 0, 500)
# h['hEMudR_dPhiMj'] = ROOT.TH1F ("hEMu_dR_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hEMudR_dPhiMmu'] = ROOT.TH1F ("hEMu_dR_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

h['hEMudR_M_mod'] = ROOT.TH1F ("hEMu_dR_M_mod", "#mu - #mu mass;M_{#mu#mu};", 500, 0, 200)

# #-----MET cuts                                                                                                                                                                                               
# h['hEMuMetcut_M'] = ROOT.TH1F ("hEMu_Metcut_M", "#mu - #mu mass;M_{#mu#mu};", 1000, 0, 200)
# h['hEMuMetcut_ePt'] = ROOT.TH1F ("hEMu_Metcut_ePt", "#mu_{1} P_t;P_{t};", 500, 0, 500)
# h['hEMuMetcut_muPt'] = ROOT.TH1F ("hEMu_Metcut_muPt", "#mu_{2} P_t;P_{t};", 500, 0, 500)
# h['hEMuMetcut_jPt'] = ROOT.TH1F ("hEMu_Metcut_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hEMuMetcut_dR'] = ROOT.TH1F ("hEMu_Metcut_dR", ";#Delta R;", 100, 0, 5)
# h['hEMuMetcut_dRlj'] = ROOT.TH1F ("hEMu_Metcut_dRlj", ";#Delta R;", 100, 0, 5)
# h['hEMuMetcut_METPt'] = ROOT.TH1F ("hEMu_Metcut_METPt", ";p_{T};", 500, 0, 500)
# h['hEMuMetcut_dPhiMj'] = ROOT.TH1F ("hEMu_Metcut_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hEMuMetcut_dPhiMmu'] = ROOT.TH1F ("hEMu_Metcut_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

# #-----dPhi cuts           
# h['hEMudPhi_M'] = ROOT.TH1F ("hEMu_dPhi_M", "#mu - #mu mass;M_{#mu#mu};", 1000, 0, 200)
# h['hEMudPhi_ePt'] = ROOT.TH1F ("hEMu_dPhi_ePt", "#mu_{1} P_t;P_{t};", 500, 0, 500)
# h['hEMudPhi_muPt'] = ROOT.TH1F ("hEMu_dPhi_muPt", "#mu_{2} P_t;P_{t};", 500, 0, 500)
# h['hEMudPhi_jPt'] = ROOT.TH1F ("hEMu_dPhi_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hEMudPhi_dR'] = ROOT.TH1F ("hEMu_dPhi_dR", ";#Delta R;", 100, 0, 5)
# h['hEMudPhi_dRlj'] = ROOT.TH1F ("hEMu_dPhi_dRlj", ";#Delta R;", 100, 0, 5)
# h['hEMudPhi_METPt'] = ROOT.TH1F ("hEMu_dPhi_METPt", ";p_{T};", 500, 0, 500)
# h['hEMudPhi_dPhiMj'] = ROOT.TH1F ("hEMu_dPhi_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hEMudPhi_dPhiMmu'] = ROOT.TH1F ("hEMu_dPhi_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

# #-----genMatching
h['hEMuGen_M'] = ROOT.TH1F ("hEMu_Gen_M", "#mu - #mu mass;M_{e#mu};", 500, 0, 200)
# h['hEMuGen_ePt'] = ROOT.TH1F ("hEMu_Gen_ePt", "#mu_{1} P_t;P_{t};", 500, 0, 500)
# h['hEMuGen_muPt'] = ROOT.TH1F ("hEMu_Gen_muPt", "#mu_{2} P_t;P_{t};", 500, 0, 500)

h['hEMuGen_M_mod'] = ROOT.TH1F ("hEMu_Gen_M_mod", "#mu - #mu mass;M_{e#mu};", 500, 0, 200)

#</histograms>

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    print inputFileName.replace("\n","")
    inputFileName = "root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCPAnalysis/RunIISummer17DR94Premix/"+inputFileName.replace("\n","")
#    inputFileName=inputFileName.replace("\n","")
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

        event.getByLabel(labelTauElectronCleaned, handleTauElectronCleaned)
        tausElectronCleaned=handleTauElectronCleaned.product()

        event.getByLabel(labelTau, handleTau)
        taus=handleTau.product()

        event.getByLabel(labelJet, handleJet)
        jets=handleJet.product()

        event.getByLabel(labelHLT, handleHLT)
        triggerResults=handleHLT.product()
#        print triggerResults
#        print triggerResults.size()
        names = event.object().triggerNames(triggerResults)
        
        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()

        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()

        event.getByLabel(labelRho, handleRho)
        if len(handleRho.product())>0:
            rho = handleRho.product()[0]
        else:
            rho = 0

        event.getByLabel(labelPatMETs, handlePatMETs)
        met=handlePatMETs.product().front()
        mets=[]
        mets+=[met]
        mets.sort(key=lambda x: x.pt(), reverse=True)
        
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
        selected_muons_iso=[]
        for muon in muons:
#            h['hMuPt'].Fill(muon.pt(),genweight)
            selected_muons+=[muon]
            if muonIsoCut(muon)<0.25:
#                h['hMuPt_iso'].Fill(muon.pt(),genweight)
                selected_muons_iso+=[muon]

        selected_muons.sort(key=lambda x: x.pt(), reverse=True) 
        selected_muons_iso.sort(key=lambda x: x.pt(), reverse=True) 

#<\muonSelection>

#<electronSelection>

        selected_electrons=[]
        selected_electrons_mod=[]
        for electron in electrons:
            E_c = electron.superCluster().energy()
            if electron.isEB():
                if electron.hadronicOverEm() < (0.05 + 1.16/E_c + 0.0324*rho/E_c) and GsfEleEffAreaPFIsoCut(electron, rho) < 0.112+0.506/electron.pt(): selected_electrons+=[electron]
                if electron.hadronicOverEm() < 0.215 and GsfEleEffAreaPFIsoCut(electron, rho) < 0.112+0.506/electron.pt(): selected_electrons_mod+=[electron]
            if electron.isEE():
                if electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c) and GsfEleEffAreaPFIsoCut(electron, rho) < 0.108+0.963/electron.pt(): selected_electrons+=[electron]
                if electron.hadronicOverEm() < 0.0984 and GsfEleEffAreaPFIsoCut(electron, rho) < 0.108+0.963/electron.pt(): selected_electrons_mod+=[electron]

        selected_electrons.sort(key=lambda x: x.pt(), reverse=True) 
        selected_electrons_mod.sort(key=lambda x: x.pt(), reverse=True) 

#<\electronSelection>

        # selected_eTaus=[]
        # for eTau in tausElectronCleaned:
        #     # if eTau.tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"): h['hTauPt_raw'].Fill(eTau.pt(),genweight)
        #     # if eTau.tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vvloose'].Fill(eTau.pt(),genweight)
        #     # if eTau.tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vloose'].Fill(eTau.pt(),genweight)
        #     # if eTau.tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_loose'].Fill(eTau.pt(),genweight)
        #     if eTau.tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"): 
        #         h['hTauPt_medium'].Fill(eTau.pt(),genweight)
        #         selected_eTaus+=[eTau]
        #     # if eTau.tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_tight'].Fill(eTau.pt(),genweight)
        #     # if eTau.tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vtight'].Fill(eTau.pt(),genweight)
        #     # if eTau.tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vvtight'].Fill(eTau.pt(),genweight)

        # selected_eTaus.sort(key=lambda x: x.pt(), reverse=True)

        # selected_taus=[]
        # for tau in taus:
        #     if tau.tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"): 
        #         selected_taus+=[tau]

        # selected_taus.sort(key=lambda x: x.pt(), reverse=True)

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
                h['hJetPt'].Fill(jet.pt(), genweight)
                selected_jets+=[jet] 
                if jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.9535:
                    selected_bjets+=[jet]
                
        selected_jets.sort(key=lambda x: x.pt(), reverse=True)       
        selected_bjets.sort(key=lambda x: x.pt(), reverse=True)


#----------------------------- tau_mu tau_e selection ------------------------------


        if len(selected_electrons)>=1 and len(selected_muons_iso)>=1 and len(selected_jets)>=1 and selected_electrons[0].charge()*selected_muons_iso[0].charge()<0:
#            h['hEMuBaseline_NbJets'].Fill(len(selected_bjets), genweight)

            if len(selected_bjets)==0:

#                if len(selected_muons)>=1: h['hEMuBaseline_ePt'].Fill(selected_muons[0].pt(), genweight)

                mu=ROOT.TLorentzVector()
                mu.SetPtEtaPhiM(selected_muons_iso[0].pt(), selected_muons_iso[0].eta(), selected_muons_iso[0].phi(), selected_muons_iso[0].mass())

                e=ROOT.TLorentzVector()
                e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass()) 

                m=ROOT.TLorentzVector()
                m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass()) 
                    
                j=ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                
                h['hEMuBaseline_M'].Fill((e+mu).M(), genweight)
                # h['hEMuBaseline_ePt'].Fill(e.Pt(), genweight)
                # h['hEMuBaseline_muPt'].Fill(mu.Pt(), genweight)
                # h['hEMuBaseline_jPt'].Fill(j.Pt(), genweight)
                # h['hEMuBaseline_dR'].Fill(e.DeltaR(mu), genweight)
                # h['hEMuBaseline_dRlj'].Fill((e+mu).DeltaR(j), genweight)
                # h['hEMuBaseline_METPt'].Fill(m.Pt(), genweight)
                # h['hEMuBaseline_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                # h['hEMuBaseline_dPhiME'].Fill(m.DeltaPhi(e), genweight)
                # h['hEMuBaseline_dPhiMtau'].Fill(m.DeltaPhi(mu), genweight)

                if e.DeltaR(mu)<0.4 and e.DeltaR(j)>0.8 and mu.DeltaR(j)>0.8:
                    h['hEMudR_M'].Fill((e+mu).M(), genweight)

                    if len(genMuons)==1 and len(genElectrons)==1: 
                        genMu=ROOT.TLorentzVector()
                        genMu.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

                        genEle=ROOT.TLorentzVector()
                        genEle.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())

                        if e.DeltaR(genEle)<0.3 and mu.DeltaR(genMu)<0.3:
                            h['hEMuGen_M'].Fill((e+mu).M(), genweight)


        if len(selected_electrons_mod)>=1 and len(selected_muons_iso)>=1 and len(selected_jets)>=1 and selected_electrons_mod[0].charge()*selected_muons_iso[0].charge()<0:
#            h['hEMuBaseline_NbJets'].Fill(len(selected_bjets), genweight)

            if len(selected_bjets)==0:

#                if len(selected_muons)>=1: h['hEMuBaseline_ePt'].Fill(selected_muons[0].pt(), genweight)

                mu=ROOT.TLorentzVector()
                mu.SetPtEtaPhiM(selected_muons_iso[0].pt(), selected_muons_iso[0].eta(), selected_muons_iso[0].phi(), selected_muons_iso[0].mass())

                e=ROOT.TLorentzVector()
                e.SetPtEtaPhiM(selected_electrons_mod[0].pt(), selected_electrons_mod[0].eta(), selected_electrons_mod[0].phi(), selected_electrons_mod[0].mass()) 

                m=ROOT.TLorentzVector()
                m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass()) 

                j=ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                
                h['hEMuBaseline_M_mod'].Fill((e+mu).M(), genweight)
                # h['hEMuBaseline_ePt_mod'].Fill(e.Pt(), genweight)
                # h['hEMuBaseline_muPt_mod'].Fill(mu.Pt(), genweight)
                # h['hEMuBaseline_jPt_mod'].Fill(j.Pt(), genweight)
                # h['hEMuBaseline_dR_mod'].Fill(e.DeltaR(mu), genweight)
                # h['hEMuBaseline_dRlj_mod'].Fill((e+mu).DeltaR(j), genweight)
                # h['hEMuBaseline_METPt_mod'].Fill(m.Pt(), genweight)
                # h['hEMuBaseline_dPhiMj_mod'].Fill(m.DeltaPhi(j), genweight)
                # h['hEMuBaseline_dPhiME_mod'].Fill(m.DeltaPhi(e), genweight)
                # h['hEMuBaseline_dPhiMtau_mod'].Fill(m.DeltaPhi(mu), genweight)

                if e.DeltaR(mu)<0.4 and e.DeltaR(j)>0.8 and mu.DeltaR(j)>0.8:
                    h['hEMudR_M_mod'].Fill((e+mu).M(), genweight)


                    if len(genMuons)==1 and len(genElectrons)==1: 
                        genMu=ROOT.TLorentzVector()
                        genMu.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

                        genEle=ROOT.TLorentzVector()
                        genEle.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())

                        if e.DeltaR(genEle)<0.3 and mu.DeltaR(genMu)<0.3: 
                            h['hEMuGen_M_mod'].Fill((e+mu).M(), genweight)


#            print triggerResults
#            print triggerResults.size()
#            print triggerResults.wasrun()
#            print names.triggerIndex("HLT_IsoTau4_v10")
#            print names.triggerIndex("HLT_Mu50_v11")


            # if (e.Pt()>26 and (triggerResults.accept(names.triggerIndex("HLT_IsoTkTau4_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoTau4_v4")))) \
            # or (e.Pt()>50 and (triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")))) \
            # or (j.Pt()>500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))):

            # if e.Pt()>26 or j.Pt()>500:# and triggerResults.accept(names.triggerIndex("HLT_IsoTau4_v10")):

            # if j.Pt() > 500 :

            #     h['hEMuTrig_M'].Fill((mu+tau).M(), genweight)
            #     h['hEMuTrig_ePt'].Fill(e.Pt(), genweight)
            #     h['hEMuTrig_muPt'].Fill(mu.Pt(), genweight)
            #     h['hEMuTrig_jPt'].Fill(j.Pt(), genweight)
            #     h['hEMuTrig_dR'].Fill(mu.DeltaR(tau), genweight)
            #     h['hEMuTrig_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #     h['hEMuTrig_METPt'].Fill(m.Pt(), genweight)
            #     h['hEMuTrig_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #     h['hEMuTrig_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)
                
            # if mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8:

            #     h['hEMudR_M'].Fill((mu+tau).M(), genweight)
            #     h['hEMudR_ePt'].Fill(e.Pt(), genweight)
            #     h['hEMudR_muPt'].Fill(mu.Pt(), genweight)
            #     h['hEMudR_jPt'].Fill(j.Pt(), genweight)
            #     h['hEMudR_dR'].Fill(mu.DeltaR(tau), genweight)
            #     h['hEMudR_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #     h['hEMudR_METPt'].Fill(m.Pt(), genweight)
            #     h['hEMudR_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #     h['hEMudR_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

            #     if m.Pt()>100:

            #         h['hEMuMetcut_M'].Fill((mu+tau).M(), genweight)
            #         h['hEMuMetcut_ePt'].Fill(e.Pt(), genweight)
            #         h['hEMuMetcut_muPt'].Fill(mu.Pt(), genweight)
            #         h['hEMuMetcut_jPt'].Fill(j.Pt(), genweight)
            #         h['hEMuMetcut_dR'].Fill(mu.DeltaR(tau), genweight)
            #         h['hEMuMetcut_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #         h['hEMuMetcut_METPt'].Fill(m.Pt(), genweight)
            #         h['hEMuMetcut_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #         h['hEMuMetcut_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

            #         if abs(m.DeltaPhi(mu))<1 and abs(m.DeltaPhi(j))>2:
            #             h['hEMudPhi_M'].Fill((mu+tau).M(), genweight)
            #             h['hEMudPhi_ePt'].Fill(e.Pt(), genweight)
            #             h['hEMudPhi_muPt'].Fill(mu.Pt(), genweight)
            #             h['hEMudPhi_jPt'].Fill(j.Pt(), genweight)
            #             h['hEMudPhi_dR'].Fill(mu.DeltaR(tau), genweight)
            #             h['hEMudPhi_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #             h['hEMudPhi_METPt'].Fill(m.Pt(), genweight)
            #             h['hEMudPhi_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #             h['hEMudPhi_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

            #             if len(genMuons)==2 and len(genJets)>=1 and len(genElectrons)==0:
            #                 genMu=ROOT.TLorentzVector()
            #                 genMu.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

            #                 genTau=ROOT.TLorentzVector()
            #                 genTau.SetPtEtaPhiM(genMuons[1].pt(), genMuons[1].eta(), genMuons[1].phi(), genMuons[1].mass())

            #                 if (mu.DeltaR(genMu)<0.3 or mu.DeltaR(genTau)<0.3) and (tau.DeltaR(genMu)<0.3 or tau.DeltaR(genTau)<0.3): 
            #                     h['hEMuGen_M'].Fill((mu+tau).M(), genweight)
            #                     h['hEMuGen_ePt'].Fill(e.Pt(), genweight)
            #                     h['hEMuGen_muPt'].Fill(mu.Pt(), genweight)

            #             if j.Pt() > 500 :
            #                 h['hEMuTrig_M'].Fill((mu+tau).M(), genweight)
            #                 h['hEMuTrig_ePt'].Fill(e.Pt(), genweight)
            #                 h['hEMuTrig_muPt'].Fill(mu.Pt(), genweight)
            #                 h['hEMuTrig_jPt'].Fill(j.Pt(), genweight)
            #                 h['hEMuTrig_dR'].Fill(mu.DeltaR(tau), genweight)
            #                 h['hEMuTrig_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #                 h['hEMuTrig_METPt'].Fill(m.Pt(), genweight)
            #                 h['hEMuTrig_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #                 h['hEMuTrig_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

                        
out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
