#Code for tau_mu tau_mu studies for 94X

import ROOT,sys,os
from DataFormats.FWLite import Events, Handle
from looseElectron import *
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

signal=True

#inputFileListDir="./filelists/"+Sample+"/"
inputFileListName=sys.argv[1]

#inputFileList=inputFileListDir+inputFileListName
inputFileList=inputFileListName

if len(sys.argv)>2:
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./plots/"

outputFileName = outputFileDir+"h_ETau_"+inputFileListName.split("/")[-1].replace(".txt",".root")
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

handleTauElectronCleanedFlat = Handle ('vector<pat::Tau>')
labelTauElectronCleanedFlat = ('selectedPATTausElectronCleanedFlat')

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

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

out=ROOT.TFile.Open(outputFileName,'recreate')

# <histograms>

pi = math.pi

h = {}

h['hNEvent'] = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)
#h['hETauBaseline_NbJets'] = ROOT.TH1F ("hMuTau_Baseline_NbJets","", 10, 0, 10)
#h['hETauBaseline_muPt'] = ROOT.TH1F ("hMuTau_Baseline_muPt", "", 100, 0, 500)

#Objects

#-----------------------Electrons------------------------

h['hEPt_EB'] = ROOT.TH1F ("hEPt_EB", "electron P_{t}; P_{t};", 100, 0, 500)
h['hEPt_EB_mod'] = ROOT.TH1F ("hEPt_EB_mod", "electron P_{t}; P_{t};", 100, 0, 500)
h['hEPt_EE'] = ROOT.TH1F ("hEPt_EE", "electron P_{t}; P_{t};", 100, 0, 500)
h['hEPt_EE_mod'] = ROOT.TH1F ("hEPt_EE_mod", "electron P_{t}; P_{t};", 100, 0, 500)

#-----------------------Taus------------------------

h['hTauPt'] = ROOT.TH1F ("hTauPt", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_mod'] = ROOT.TH1F ("hTauPt_mod", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_raw'] = ROOT.TH1F ("hTauPt_raw", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_loose'] = ROOT.TH1F ("hTauPt_loose", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_medium'] = ROOT.TH1F ("hTauPt_medium", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_medium_mod'] = ROOT.TH1F ("hTauPt_medium_mod", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_tight'] = ROOT.TH1F ("hTauPt_tight", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_vloose'] = ROOT.TH1F ("hTauPt_vloose", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_vtight'] = ROOT.TH1F ("hTauPt_vtight", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_vvloose'] = ROOT.TH1F ("hTauPt_vvloose", "#tau P_t; P_t;", 100, 0, 500)
#h['hTauPt_vvtight'] = ROOT.TH1F ("hTauPt_vvtight", "#tau P_t; P_t;", 100, 0, 500)

#-----------------------Muons------------------------

#h['hMuPt'] = ROOT.TH1F ("hMuPt", "#mu P_t; P_t;", 100, 0, 500)
#h['hMuPt_iso'] = ROOT.TH1F ("hMuPt_iso", "#mu P_t; P_t;", 100, 0, 500)

#-----------------------Jets-------------------------

h['hJetPt'] = ROOT.TH1F ("hJetPt", "jet pt;P_{t};", 2000, 0, 2000)

#Event
#-----------------------ETau------------------------                                                                                                                                                         
#-----Require 1 Electron, 1 Tau                                                                                                                                                                                          
#h['hIsoMuTauBaseline_M'] = ROOT.TH1F ("hIsoMuTau_Baseline_M", "#tau - e mass;M_{#tau#tau};", 500, 0, 200)

h['hETauBaseline_M'] = ROOT.TH1F ("hETau_Baseline_M", "#tau - e mass;M_{#tau#tau};", 500, 0, 200)
h['hETauBaseline_ePt'] = ROOT.TH1F ("hETau_Baseline_ePt", "e P_t; P_t;", 100, 0, 500)
h['hETauBaseline_tauPt'] = ROOT.TH1F ("hETau_Baseline_tauPt", "#tau P_t; P_t;", 100, 0, 500)
# h['hETauBaseline_jPt'] = ROOT.TH1F ("hETau_Baseline_jPt", "jet pt;P_{t};", 2000, 0, 2000)
# h['hETauBaseline_dR'] = ROOT.TH1F ("hETau_Baseline_dR", ";#Delta R;", 100, 0, 5)
# h['hETauBaseline_dRlj'] = ROOT.TH1F ("hETau_Baseline_dRlj", ";#Delta R;", 100, 0, 5)
# h['hETauBaseline_METPt'] = ROOT.TH1F ("hETau_Baseline_METPt", ";p_{T};", 500, 0, 500)
# h['hETauBaseline_dPhiMj'] = ROOT.TH1F ("hETau_Baseline_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hETauBaseline_dPhiME'] = ROOT.TH1F ("hETau_Baseline_dPhiME", ";#delta#phi_{e};", 100, -pi, pi)
# h['hETauBaseline_dPhiMtau'] = ROOT.TH1F ("hETau_Baseline_dPhiMTau", ";#delta#phi_{#tau};", 100, -pi, pi)

h['hETauBaseline_M_mod'] = ROOT.TH1F ("hETau_Baseline_M_mod", "#tau - e mass;M_{#tau#tau};", 500, 0, 200)
h['hETauBaseline_ePt_mod'] = ROOT.TH1F ("hETau_Baseline_ePt_mod", "e P_t; P_t;", 100, 0, 500)
h['hETauBaseline_tauPt_mod'] = ROOT.TH1F ("hETau_Baseline_tauPt_mod", "#tau P_t; P_t;", 100, 0, 500)
# h['hETauBaseline_jPt_mod'] = ROOT.TH1F ("hETau_Baseline_jPt_mod", "jet pt;P_{t};", 2000, 0, 2000)
# h['hETauBaseline_dR_mod'] = ROOT.TH1F ("hETau_Baseline_dR_mod", ";#Delta R;", 100, 0, 5)
# h['hETauBaseline_dRlj_mod'] = ROOT.TH1F ("hETau_Baseline_dRlj_mod", ";#Delta R;", 100, 0, 5)
# h['hETauBaseline_METPt_mod'] = ROOT.TH1F ("hETau_Baseline_METPt_mod", ";p_{T};", 500, 0, 500)
# h['hETauBaseline_dPhiMj_mod'] = ROOT.TH1F ("hETau_Baseline_dPhiMj_mod", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hETauBaseline_dPhiME_mod'] = ROOT.TH1F ("hETau_Baseline_dPhiME_mod", ";#delta#phi_{e};", 100, -pi, pi)
# h['hETauBaseline_dPhiMtau_mod'] = ROOT.TH1F ("hETau_Baseline_dPhiMTau_mod", ";#delta#phi_{#tau};", 100, -pi, pi)

h['hETauBaseline_M_nC'] = ROOT.TH1F ("hETau_Baseline_M_nC", "#tau - e mass;M_{#tau#tau};", 500, 0, 200)
h['hETauBaseline_M_nC_mod'] = ROOT.TH1F ("hETau_Baseline_M_nC_mod", "#tau - e mass;M_{#tau#tau};", 500, 0, 200)

#-----Passing Triggers                                                                                                                                                                       
# h['hETauTrig_M'] = ROOT.TH1F ("hETau_Trig_M", "#tau - #tau mass;M_{#tau#tau};", 1000, 0, 200)
# h['hETauTrig_muPt'] = ROOT.TH1F ("hETau_Trig_muPt", "#tau_{1} P_t;P_{t};", 500, 0, 500)
# h['hETauTrig_tauPt'] = ROOT.TH1F ("hETau_Trig_tauPt", "#tau_{2} P_t;P_{t};", 500, 0, 500)
# h['hETauTrig_jPt'] = ROOT.TH1F ("hETau_Trig_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hETauTrig_dR'] = ROOT.TH1F ("hETau_Trig_dR", ";#Delta R;", 100, 0, 5)
# h['hETauTrig_dRlj'] = ROOT.TH1F ("hETau_Trig_dRlj", ";#Delta R;", 100, 0, 5)
# h['hETauTrig_METPt'] = ROOT.TH1F ("hETau_Trig_METPt", ";p_{T};", 500, 0, 500)
# h['hETauTrig_dPhiMj'] = ROOT.TH1F ("hETau_Trig_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hETauTrig_dPhiMmu'] = ROOT.TH1F ("hETau_Trig_dPhiMmu", ";#delta#phi_{#tau};", 100, -pi, pi)

# #-----dR cuts                                                                                                                                                                                 
h['hETaudR_M'] = ROOT.TH1F ("hETau_dR_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hETaudR_ePt'] = ROOT.TH1F ("hETau_dR_ePt", "#tau_{1} P_t;P_{t};", 500, 0, 500)
# h['hETaudR_tauPt'] = ROOT.TH1F ("hETau_dR_tauPt", "#tau_{2} P_t;P_{t};", 500, 0, 500)
# h['hETaudR_jPt'] = ROOT.TH1F ("hETau_dR_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hETaudR_dR'] = ROOT.TH1F ("hETau_dR_dR", ";#Delta R;", 100, 0, 5)
# h['hETaudR_dRlj'] = ROOT.TH1F ("hETau_dR_dRlj", ";#Delta R;", 100, 0, 5)
# h['hETaudR_METPt'] = ROOT.TH1F ("hETau_dR_METPt", ";p_{T};", 500, 0, 500)
# h['hETaudR_dPhiMj'] = ROOT.TH1F ("hETau_dR_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hETaudR_dPhiMmu'] = ROOT.TH1F ("hETau_dR_dPhiMmu", ";#delta#phi_{#tau};", 100, -pi, pi)

h['hETaudR_M_mod'] = ROOT.TH1F ("hETau_dR_M_mod", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['hETaudR_M_nC'] = ROOT.TH1F ("hETau_dR_M_nC", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hETaudR_M_nC_mod'] = ROOT.TH1F ("hETau_dR_M_nC_mod", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

# #-----MET cuts                                                                                                                                                                                               
# h['hETauMetcut_M'] = ROOT.TH1F ("hETau_Metcut_M", "#tau - #tau mass;M_{#tau#tau};", 1000, 0, 200)
# h['hETauMetcut_muPt'] = ROOT.TH1F ("hETau_Metcut_muPt", "#tau_{1} P_t;P_{t};", 500, 0, 500)
# h['hETauMetcut_tauPt'] = ROOT.TH1F ("hETau_Metcut_tauPt", "#tau_{2} P_t;P_{t};", 500, 0, 500)
# h['hETauMetcut_jPt'] = ROOT.TH1F ("hETau_Metcut_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hETauMetcut_dR'] = ROOT.TH1F ("hETau_Metcut_dR", ";#Delta R;", 100, 0, 5)
# h['hETauMetcut_dRlj'] = ROOT.TH1F ("hETau_Metcut_dRlj", ";#Delta R;", 100, 0, 5)
# h['hETauMetcut_METPt'] = ROOT.TH1F ("hETau_Metcut_METPt", ";p_{T};", 500, 0, 500)
# h['hETauMetcut_dPhiMj'] = ROOT.TH1F ("hETau_Metcut_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hETauMetcut_dPhiMmu'] = ROOT.TH1F ("hETau_Metcut_dPhiMmu", ";#delta#phi_{#tau};", 100, -pi, pi)

# #-----dPhi cuts           
# h['hETaudPhi_M'] = ROOT.TH1F ("hETau_dPhi_M", "#tau - #tau mass;M_{#tau#tau};", 1000, 0, 200)
# h['hETaudPhi_muPt'] = ROOT.TH1F ("hETau_dPhi_muPt", "#tau_{1} P_t;P_{t};", 500, 0, 500)
# h['hETaudPhi_tauPt'] = ROOT.TH1F ("hETau_dPhi_tauPt", "#tau_{2} P_t;P_{t};", 500, 0, 500)
# h['hETaudPhi_jPt'] = ROOT.TH1F ("hETau_dPhi_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hETaudPhi_dR'] = ROOT.TH1F ("hETau_dPhi_dR", ";#Delta R;", 100, 0, 5)
# h['hETaudPhi_dRlj'] = ROOT.TH1F ("hETau_dPhi_dRlj", ";#Delta R;", 100, 0, 5)
# h['hETaudPhi_METPt'] = ROOT.TH1F ("hETau_dPhi_METPt", ";p_{T};", 500, 0, 500)
# h['hETaudPhi_dPhiMj'] = ROOT.TH1F ("hETau_dPhi_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hETaudPhi_dPhiMmu'] = ROOT.TH1F ("hETau_dPhi_dPhiMmu", ";#delta#phi_{#tau};", 100, -pi, pi)

# #-----genMatching
h['hETauGen_M'] = ROOT.TH1F ("hETau_Gen_M", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
h['hETauGen_M_nC'] = ROOT.TH1F ("hETau_Gen_M_nC", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)
# h['hETauGen_muPt'] = ROOT.TH1F ("hETau_Gen_muPt", "#tau_{1} P_t;P_{t};", 500, 0, 500)
# h['hETauGen_tauPt'] = ROOT.TH1F ("hETau_Gen_tauPt", "#tau_{2} P_t;P_{t};", 500, 0, 500)

h['hETauGen_M_mod'] = ROOT.TH1F ("hETau_Gen_M_mod", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

h['hETauGen_M_nC_mod'] = ROOT.TH1F ("hETau_Gen_M_nC_mod", "#tau - #tau mass;M_{#tau#tau};", 500, 0, 200)

#</histograms>

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    print inputFileName.replace("\n","")

    if signal==True:inputFileName = "root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCPAnalysis/RunIISummer17DR94Premix/"+inputFileName.replace("\n","")
    else:inputFileName = "root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCPAnalysis/Backgrounds/RunIIFall17DR94Premix/"+inputFileName.replace("\n","")

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

        event.getByLabel(labelTauElectronCleanedFlat, handleTauElectronCleanedFlat)
        tausElectronCleanedFlat=handleTauElectronCleanedFlat.product()

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


        if signal == True:

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
                if electron.hadronicOverEm() < (0.05 + 1.16/E_c + 0.0324*rho/E_c): 
                    h['hEPt_EB'].Fill(electron.pt(), genweight)
                    selected_electrons+=[electron]
                if electron.hadronicOverEm() < 0.215 :
                    h['hEPt_EB_mod'].Fill(electron.pt(), genweight) 
                    selected_electrons_mod+=[electron]
            if electron.isEE():
                if electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c): 
                    h['hEPt_EE'].Fill(electron.pt(), genweight)                    
                    selected_electrons+=[electron]
                if electron.hadronicOverEm() < 0.0984 : 
                    h['hEPt_EE_mod'].Fill(electron.pt(), genweight)
                    selected_electrons_mod+=[electron]


        selected_electrons.sort(key=lambda x: x.pt(), reverse=True) 
        selected_electrons_mod.sort(key=lambda x: x.pt(), reverse=True) 

#<\electronSelection>

        selected_eTaus=[]
        selected_eTaus_flat=[]

        for eTau in tausElectronCleanedFlat:
            h['hTauPt_mod'].Fill(eTau.pt(), genweight)
            if eTau.tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"):
                h['hTauPt_medium_mod'].Fill(eTau.pt(),genweight)
                selected_eTaus_flat+=[eTau]

        for eTau in tausElectronCleaned:
            h['hTauPt'].Fill(eTau.pt(), genweight)
            # if eTau.tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"): h['hTauPt_raw'].Fill(eTau.pt(),genweight)
            # if eTau.tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vvloose'].Fill(eTau.pt(),genweight)
            # if eTau.tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vloose'].Fill(eTau.pt(),genweight)
            # if eTau.tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_loose'].Fill(eTau.pt(),genweight)
            if eTau.tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"): 
                h['hTauPt_medium'].Fill(eTau.pt(),genweight)
                selected_eTaus+=[eTau]
            # if eTau.tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_tight'].Fill(eTau.pt(),genweight)
            # if eTau.tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vtight'].Fill(eTau.pt(),genweight)
            # if eTau.tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vvtight'].Fill(eTau.pt(),genweight)

        selected_eTaus.sort(key=lambda x: x.pt(), reverse=True)
        selected_eTaus_flat.sort(key=lambda x: x.pt(), reverse=True)

        selected_taus=[]
        for tau in taus:
            if tau.tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"): 
                selected_taus+=[tau]

        selected_taus.sort(key=lambda x: x.pt(), reverse=True)

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


#----------------------------- tau_h tau_e selection ------------------------------


        if len(selected_electrons)>=1 and len(selected_eTaus)>=1 and len(selected_jets)>=1 and selected_electrons[0].charge()*selected_eTaus[0].charge()<0:
#            h['hETauBaseline_NbJets'].Fill(len(selected_bjets), genweight)

            if len(selected_bjets)==0:

#                if len(selected_muons)>=1: h['hETauBaseline_muPt'].Fill(selected_muons[0].pt(), genweight)

                etau=ROOT.TLorentzVector()
                etau.SetPtEtaPhiM(selected_eTaus[0].pt(), selected_eTaus[0].eta(), selected_eTaus[0].phi(), selected_eTaus[0].mass())

                e=ROOT.TLorentzVector()
                e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass()) 

                m=ROOT.TLorentzVector()
                m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass()) 

                if len(selected_jets)==1:
                    j=ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                else:
                    j1=ROOT.TLorentzVector()
                    j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                    j2=ROOT.TLorentzVector()
                    j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())

                    if etau.DeltaR(j1) > 0.3:
                        j = j1
                    else:
                        j = j2

                h['hETauBaseline_M'].Fill((e+etau).M(), genweight)
                h['hETauBaseline_ePt'].Fill(e.Pt(), genweight)
                h['hETauBaseline_tauPt'].Fill(etau.Pt(), genweight)
                # h['hETauBaseline_jPt'].Fill(j.Pt(), genweight)
                # h['hETauBaseline_dR'].Fill(e.DeltaR(etau), genweight)
                # h['hETauBaseline_dRlj'].Fill((e+etau).DeltaR(j), genweight)
                # h['hETauBaseline_METPt'].Fill(m.Pt(), genweight)
                # h['hETauBaseline_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                # h['hETauBaseline_dPhiME'].Fill(m.DeltaPhi(e), genweight)
                # h['hETauBaseline_dPhiMtau'].Fill(m.DeltaPhi(etau), genweight)

                if e.DeltaR(etau)<0.4 and e.DeltaR(j)>0.8 and etau.DeltaR(j)>0.8:
                    h['hETaudR_M'].Fill((e+etau).M(), genweight)

                    if signal==True and len(genElectrons)==1 and len(genMuons)==0:
                    
                        genEle=ROOT.TLorentzVector()
                        genEle.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                        genTau=ROOT.TLorentzVector()
                        genNt=ROOT.TLorentzVector()

                        if genElectrons[0].pdgId()==11:
                            genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                            genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                        else:
                            genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                            genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())

                        if e.DeltaR(genEle)<0.2 and etau.DeltaR(genTau-genNt)<0.3:
                            h['hETauGen_M'].Fill((e+etau).M(), genweight)



        if len(selected_electrons_mod)>=1 and len(selected_eTaus_flat)>=1 and len(selected_jets)>=1 and selected_electrons_mod[0].charge()*selected_eTaus_flat[0].charge()<0:
#            h['hETauBaseline_NbJets'].Fill(len(selected_bjets), genweight)

            if len(selected_bjets)==0:

#                if len(selected_muons)>=1: h['hETauBaseline_muPt'].Fill(selected_muons[0].pt(), genweight)

                etau=ROOT.TLorentzVector()
                etau.SetPtEtaPhiM(selected_eTaus_flat[0].pt(), selected_eTaus_flat[0].eta(), selected_eTaus_flat[0].phi(), selected_eTaus_flat[0].mass())

                e=ROOT.TLorentzVector()
                e.SetPtEtaPhiM(selected_electrons_mod[0].pt(), selected_electrons_mod[0].eta(), selected_electrons_mod[0].phi(), selected_electrons_mod[0].mass()) 

                m=ROOT.TLorentzVector()
                m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass()) 

                if len(selected_jets)==1:
                    j=ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                else:
                    j1=ROOT.TLorentzVector()
                    j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                    j2=ROOT.TLorentzVector()
                    j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())

                    if etau.DeltaR(j1) > 0.3:
                        j = j1
                    else:
                        j = j2

                h['hETauBaseline_M_mod'].Fill((e+etau).M(), genweight)
                h['hETauBaseline_ePt_mod'].Fill(e.Pt(), genweight)
                h['hETauBaseline_tauPt_mod'].Fill(etau.Pt(), genweight)
                # h['hETauBaseline_jPt_mod'].Fill(j.Pt(), genweight)
                # h['hETauBaseline_dR_mod'].Fill(e.DeltaR(etau), genweight)
                # h['hETauBaseline_dRlj_mod'].Fill((e+etau).DeltaR(j), genweight)
                # h['hETauBaseline_METPt_mod'].Fill(m.Pt(), genweight)
                # h['hETauBaseline_dPhiMj_mod'].Fill(m.DeltaPhi(j), genweight)
                # h['hETauBaseline_dPhiME_mod'].Fill(m.DeltaPhi(e), genweight)
                # h['hETauBaseline_dPhiMtau_mod'].Fill(m.DeltaPhi(etau), genweight)

                if e.DeltaR(etau)<0.4 and e.DeltaR(j)>0.8 and etau.DeltaR(j)>0.8:
                    h['hETaudR_M_mod'].Fill((e+etau).M(), genweight)

                    if signal==True and len(genElectrons)==1 and len(genMuons)==0:
                    
                        genEle=ROOT.TLorentzVector()
                        genEle.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                        genTau=ROOT.TLorentzVector()
                        genNt=ROOT.TLorentzVector()

                        if genElectrons[0].pdgId()==11:
                            genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                            genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                        else:
                            genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                            genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())

                        if e.DeltaR(genEle)<0.2 and etau.DeltaR(genTau-genNt)<0.3:
                            h['hETauGen_M_mod'].Fill((e+etau).M(), genweight)



        if len(selected_electrons)>=1 and len(selected_eTaus_flat)>=1 and len(selected_jets)>=1 and selected_electrons[0].charge()*selected_eTaus_flat[0].charge()<0:
#            h['hETauBaseline_NbJets'].Fill(len(selected_bjets), genweight)

            if len(selected_bjets)==0:

#                if len(selected_muons)>=1: h['hETauBaseline_muPt'].Fill(selected_muons[0].pt(), genweight)

                tau=ROOT.TLorentzVector()
                tau.SetPtEtaPhiM(selected_eTaus_flat[0].pt(), selected_eTaus_flat[0].eta(), selected_eTaus_flat[0].phi(), selected_eTaus_flat[0].mass())

                e=ROOT.TLorentzVector()
                e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass()) 

                m=ROOT.TLorentzVector()
                m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass()) 

                if len(selected_jets)==1:
                    j=ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                else:
                    j1=ROOT.TLorentzVector()
                    j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                    j2=ROOT.TLorentzVector()
                    j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())

                    if etau.DeltaR(j1) > 0.3:
                        j = j1
                    else:
                        j = j2

                h['hETauBaseline_M_nC'].Fill((e+tau).M(), genweight)
                
                if e.DeltaR(tau)<0.4 and e.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8:
                    h['hETaudR_M_nC'].Fill((e+tau).M(), genweight)

                    if signal==True and len(genElectrons)==1 and len(genMuons)==1:
                    
                        genEle=ROOT.TLorentzVector()
                        genEle.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                        genTau=ROOT.TLorentzVector()
                        genNt=ROOT.TLorentzVector()

                        if genElectrons[0].pdgId()==11:
                            genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                            genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                        else:
                            genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                            genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())

                        if e.DeltaR(genEle)<0.2 and etau.DeltaR(genTau-genNt)<0.3:
                            h['hETauGen_M_nC'].Fill((e+etau).M(), genweight)


        if len(selected_electrons_mod)>=1 and len(selected_eTaus)>=1 and len(selected_jets)>=1 and selected_electrons_mod[0].charge()*selected_eTaus[0].charge()<0:
#            h['hETauBaseline_NbJets'].Fill(len(selected_bjets), genweight)

            if len(selected_bjets)==0:

#                if len(selected_muons)>=1: h['hETauBaseline_muPt'].Fill(selected_muons[0].pt(), genweight)

                tau=ROOT.TLorentzVector()
                tau.SetPtEtaPhiM(selected_eTaus[0].pt(), selected_eTaus[0].eta(), selected_eTaus[0].phi(), selected_eTaus[0].mass())

                e=ROOT.TLorentzVector()
                e.SetPtEtaPhiM(selected_electrons_mod[0].pt(), selected_electrons_mod[0].eta(), selected_electrons_mod[0].phi(), selected_electrons_mod[0].mass()) 

                m=ROOT.TLorentzVector()
                m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass()) 

                if len(selected_jets)==1:
                    j=ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                else:
                    j1=ROOT.TLorentzVector()
                    j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                    j2=ROOT.TLorentzVector()
                    j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())

                    if etau.DeltaR(j1) > 0.3:
                        j = j1
                    else:
                        j = j2

                h['hETauBaseline_M_nC_mod'].Fill((e+tau).M(), genweight)

                if e.DeltaR(tau)<0.4 and e.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8:
                    h['hETaudR_M_nC_mod'].Fill((e+tau).M(), genweight)

                    if signal==True and len(genElectrons)==1 and len(genMuons)==0:
                    
                        genEle=ROOT.TLorentzVector()
                        genEle.SetPtEtaPhiM(genElectrons[0].pt(), genElectrons[0].eta(), genElectrons[0].phi(), genElectrons[0].mass())
                        genTau=ROOT.TLorentzVector()
                        genNt=ROOT.TLorentzVector()

                        if genElectrons[0].pdgId()==11:
                            genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==-15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                            genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==-16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())
                        else:
                            genTau.SetPtEtaPhiM(genTaus[0].pt(), genTaus[0].eta(), genTaus[0].phi(), genTaus[0].mass()) if genTaus[0].pdgId()==15 else genTau.SetPtEtaPhiM(genTaus[1].pt(), genTaus[1].eta(), genTaus[1].phi(), genTaus[1].mass())
                            genNt.SetPtEtaPhiM(genNutaus[0].pt(), genNutaus[0].eta(), genNutaus[0].phi(), genNutaus[0].mass()) if genNutaus[0].pdgId()==16 else genNt.SetPtEtaPhiM(genNutaus[1].pt(), genNutaus[1].eta(), genNutaus[1].phi(), genNutaus[1].mass())

                        if e.DeltaR(genEle)<0.2 and etau.DeltaR(genTau-genNt)<0.3:
                            h['hETauGen_M_nC_mod'].Fill((e+etau).M(), genweight)


#            print triggerResults
#            print triggerResults.size()
#            print triggerResults.wasrun()
#            print names.triggerIndex("HLT_IsoTau4_v10")
#            print names.triggerIndex("HLT_Mu50_v11")


            # if (mu.Pt()>26 and (triggerResults.accept(names.triggerIndex("HLT_IsoTkTau4_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoTau4_v4")))) \
            # or (mu.Pt()>50 and (triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")))) \
            # or (j.Pt()>500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))):

            # if mu.Pt()>26 or j.Pt()>500:# and triggerResults.accept(names.triggerIndex("HLT_IsoTau4_v10")):

            # if j.Pt() > 500 :

            #     h['hETauTrig_M'].Fill((mu+tau).M(), genweight)
            #     h['hETauTrig_muPt'].Fill(mu.Pt(), genweight)
            #     h['hETauTrig_tauPt'].Fill(tau.Pt(), genweight)
            #     h['hETauTrig_jPt'].Fill(j.Pt(), genweight)
            #     h['hETauTrig_dR'].Fill(mu.DeltaR(tau), genweight)
            #     h['hETauTrig_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #     h['hETauTrig_METPt'].Fill(m.Pt(), genweight)
            #     h['hETauTrig_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #     h['hETauTrig_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)
                
            # if mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8:

            #     h['hETaudR_M'].Fill((mu+tau).M(), genweight)
            #     h['hETaudR_muPt'].Fill(mu.Pt(), genweight)
            #     h['hETaudR_tauPt'].Fill(tau.Pt(), genweight)
            #     h['hETaudR_jPt'].Fill(j.Pt(), genweight)
            #     h['hETaudR_dR'].Fill(mu.DeltaR(tau), genweight)
            #     h['hETaudR_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #     h['hETaudR_METPt'].Fill(m.Pt(), genweight)
            #     h['hETaudR_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #     h['hETaudR_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

            #     if m.Pt()>100:

            #         h['hETauMetcut_M'].Fill((mu+tau).M(), genweight)
            #         h['hETauMetcut_muPt'].Fill(mu.Pt(), genweight)
            #         h['hETauMetcut_tauPt'].Fill(tau.Pt(), genweight)
            #         h['hETauMetcut_jPt'].Fill(j.Pt(), genweight)
            #         h['hETauMetcut_dR'].Fill(mu.DeltaR(tau), genweight)
            #         h['hETauMetcut_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #         h['hETauMetcut_METPt'].Fill(m.Pt(), genweight)
            #         h['hETauMetcut_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #         h['hETauMetcut_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

            #         if abs(m.DeltaPhi(mu))<1 and abs(m.DeltaPhi(j))>2:
            #             h['hETaudPhi_M'].Fill((mu+tau).M(), genweight)
            #             h['hETaudPhi_muPt'].Fill(mu.Pt(), genweight)
            #             h['hETaudPhi_tauPt'].Fill(tau.Pt(), genweight)
            #             h['hETaudPhi_jPt'].Fill(j.Pt(), genweight)
            #             h['hETaudPhi_dR'].Fill(mu.DeltaR(tau), genweight)
            #             h['hETaudPhi_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #             h['hETaudPhi_METPt'].Fill(m.Pt(), genweight)
            #             h['hETaudPhi_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #             h['hETaudPhi_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

            #             if len(genMuons)==2 and len(genJets)>=1 and len(genElectrons)==0:
            #                 genMu=ROOT.TLorentzVector()
            #                 genMu.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

            #                 genTau=ROOT.TLorentzVector()
            #                 genTau.SetPtEtaPhiM(genMuons[1].pt(), genMuons[1].eta(), genMuons[1].phi(), genMuons[1].mass())

            #                 if (mu.DeltaR(genMu)<0.3 or mu.DeltaR(genTau)<0.3) and (tau.DeltaR(genMu)<0.3 or tau.DeltaR(genTau)<0.3): 
            #                     h['hETauGen_M'].Fill((mu+tau).M(), genweight)
            #                     h['hETauGen_muPt'].Fill(mu.Pt(), genweight)
            #                     h['hETauGen_tauPt'].Fill(tau.Pt(), genweight)

            #             if j.Pt() > 500 :
            #                 h['hETauTrig_M'].Fill((mu+tau).M(), genweight)
            #                 h['hETauTrig_muPt'].Fill(mu.Pt(), genweight)
            #                 h['hETauTrig_tauPt'].Fill(tau.Pt(), genweight)
            #                 h['hETauTrig_jPt'].Fill(j.Pt(), genweight)
            #                 h['hETauTrig_dR'].Fill(mu.DeltaR(tau), genweight)
            #                 h['hETauTrig_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #                 h['hETauTrig_METPt'].Fill(m.Pt(), genweight)
            #                 h['hETauTrig_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #                 h['hETauTrig_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

                        
out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
