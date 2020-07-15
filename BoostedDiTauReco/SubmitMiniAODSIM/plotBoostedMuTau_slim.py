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

outputFileName = outputFileDir+"h_MuTau_"+inputFileListName.split("/")[-1].replace(".txt",".root")
print outputFileName

pi = math.pi

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('selectedPATMuons','myMuons','PAT')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('selectedPATElectrons','myElectrons','PAT')

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('selectedPATJets')

handleTauMuonCleaned = Handle ('vector<pat::Tau>')
labelTauMuonCleaned = ('selectedPATTausMuonCleaned')

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

out=ROOT.TFile.Open(outputFileName,'recreate')

# <histograms>

pi = math.pi

h = {}

h['hNEvent'] = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)
h['hMuTauBaseline_NbJets'] = ROOT.TH1F ("hMuTau_Baseline_NbJets","", 10, 0, 10)
h['hMuTauBaseline_ePt'] = ROOT.TH1F ("hMuTau_Baseline_ePt", "", 100, 0, 500)

#Objects
#-----------------------Taus------------------------

h['hTauPt_raw'] = ROOT.TH1F ("hTauPt_raw", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_loose'] = ROOT.TH1F ("hTauPt_loose", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_medium'] = ROOT.TH1F ("hTauPt_medium", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_tight'] = ROOT.TH1F ("hTauPt_tight", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_vloose'] = ROOT.TH1F ("hTauPt_vloose", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_vtight'] = ROOT.TH1F ("hTauPt_vtight", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_vvloose'] = ROOT.TH1F ("hTauPt_vvloose", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_vvtight'] = ROOT.TH1F ("hTauPt_vvtight", "#tau P_t; P_t;", 100, 0, 500)

#-----------------------Muons------------------------

h['hMuPt'] = ROOT.TH1F ("hMuPt", "#mu P_t; P_t;", 100, 0, 500)
h['hMuPt_iso'] = ROOT.TH1F ("hMuPt_iso", "#mu P_t; P_t;", 100, 0, 500)

#-----------------------Jets-------------------------

h['hJetPt'] = ROOT.TH1F ("hJetPt", "jet pt;P_{t};", 2000, 0, 2000)

#Event
#-----------------------MuTau------------------------                                                                                                                                                         
#-----Require 2 muons                                                                                                                                                                                          
h['hIsoMuTauBaseline_M'] = ROOT.TH1F ("hIsoMuTau_Baseline_M", "#tau - #mu mass;M_{#tau#tau};", 500, 0, 200)

h['hMuTauBaseline_M'] = ROOT.TH1F ("hMuTau_Baseline_M", "#tau - #mu mass;M_{#tau#tau};", 500, 0, 200)
h['hMuTauBaseline_muPt'] = ROOT.TH1F ("hMuTau_Baseline_muPt", "#mu P_t; P_t;", 100, 0, 500)
h['hMuTauBaseline_tauPt'] = ROOT.TH1F ("hMuTau_Baseline_tauPt", "#tau P_t; P_t;", 100, 0, 500)
h['hMuTauBaseline_jPt'] = ROOT.TH1F ("hMuTau_Baseline_jPt", "jet pt;P_{t};", 2000, 0, 2000)
h['hMuTauBaseline_dR'] = ROOT.TH1F ("hMuTau_Baseline_dR", ";#Delta R;", 100, 0, 5)
h['hMuTauBaseline_dRlj'] = ROOT.TH1F ("hMuTau_Baseline_dRlj", ";#Delta R;", 100, 0, 5)
h['hMuTauBaseline_METPt'] = ROOT.TH1F ("hMuTau_Baseline_METPt", ";p_{T};", 500, 0, 500)
h['hMuTauBaseline_dPhiMj'] = ROOT.TH1F ("hMuTau_Baseline_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
h['hMuTauBaseline_dPhiMmu'] = ROOT.TH1F ("hMuTau_Baseline_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)
h['hMuTauBaseline_dPhiMtau'] = ROOT.TH1F ("hMuTau_Baseline_dPhiMTau", ";#delta#phi_{#tau};", 100, -pi, pi)

#-----Passing Triggers                                                                                                                                                                       
# h['hMuTauTrig_M'] = ROOT.TH1F ("hMuTau_Trig_M", "#tau - #tau mass;M_{#tau#tau};", 1000, 0, 200)
# h['hMuTauTrig_muPt'] = ROOT.TH1F ("hMuTau_Trig_muPt", "#tau_{1} P_t;P_{t};", 500, 0, 500)
# h['hMuTauTrig_tauPt'] = ROOT.TH1F ("hMuTau_Trig_tauPt", "#tau_{2} P_t;P_{t};", 500, 0, 500)
# h['hMuTauTrig_jPt'] = ROOT.TH1F ("hMuTau_Trig_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hMuTauTrig_dR'] = ROOT.TH1F ("hMuTau_Trig_dR", ";#Delta R;", 100, 0, 5)
# h['hMuTauTrig_dRlj'] = ROOT.TH1F ("hMuTau_Trig_dRlj", ";#Delta R;", 100, 0, 5)
# h['hMuTauTrig_METPt'] = ROOT.TH1F ("hMuTau_Trig_METPt", ";p_{T};", 500, 0, 500)
# h['hMuTauTrig_dPhiMj'] = ROOT.TH1F ("hMuTau_Trig_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hMuTauTrig_dPhiMmu'] = ROOT.TH1F ("hMuTau_Trig_dPhiMmu", ";#delta#phi_{#tau};", 100, -pi, pi)

# #-----dR cuts                                                                                                                                                                                 
# h['hMuTaudR_M'] = ROOT.TH1F ("hMuTau_dR_M", "#tau - #tau mass;M_{#tau#tau};", 1000, 0, 200)
# h['hMuTaudR_muPt'] = ROOT.TH1F ("hMuTau_dR_muPt", "#tau_{1} P_t;P_{t};", 500, 0, 500)
# h['hMuTaudR_tauPt'] = ROOT.TH1F ("hMuTau_dR_tauPt", "#tau_{2} P_t;P_{t};", 500, 0, 500)
# h['hMuTaudR_jPt'] = ROOT.TH1F ("hMuTau_dR_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hMuTaudR_dR'] = ROOT.TH1F ("hMuTau_dR_dR", ";#Delta R;", 100, 0, 5)
# h['hMuTaudR_dRlj'] = ROOT.TH1F ("hMuTau_dR_dRlj", ";#Delta R;", 100, 0, 5)
# h['hMuTaudR_METPt'] = ROOT.TH1F ("hMuTau_dR_METPt", ";p_{T};", 500, 0, 500)
# h['hMuTaudR_dPhiMj'] = ROOT.TH1F ("hMuTau_dR_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hMuTaudR_dPhiMmu'] = ROOT.TH1F ("hMuTau_dR_dPhiMmu", ";#delta#phi_{#tau};", 100, -pi, pi)

# #-----MET cuts                                                                                                                                                                                               
# h['hMuTauMetcut_M'] = ROOT.TH1F ("hMuTau_Metcut_M", "#tau - #tau mass;M_{#tau#tau};", 1000, 0, 200)
# h['hMuTauMetcut_muPt'] = ROOT.TH1F ("hMuTau_Metcut_muPt", "#tau_{1} P_t;P_{t};", 500, 0, 500)
# h['hMuTauMetcut_tauPt'] = ROOT.TH1F ("hMuTau_Metcut_tauPt", "#tau_{2} P_t;P_{t};", 500, 0, 500)
# h['hMuTauMetcut_jPt'] = ROOT.TH1F ("hMuTau_Metcut_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hMuTauMetcut_dR'] = ROOT.TH1F ("hMuTau_Metcut_dR", ";#Delta R;", 100, 0, 5)
# h['hMuTauMetcut_dRlj'] = ROOT.TH1F ("hMuTau_Metcut_dRlj", ";#Delta R;", 100, 0, 5)
# h['hMuTauMetcut_METPt'] = ROOT.TH1F ("hMuTau_Metcut_METPt", ";p_{T};", 500, 0, 500)
# h['hMuTauMetcut_dPhiMj'] = ROOT.TH1F ("hMuTau_Metcut_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hMuTauMetcut_dPhiMmu'] = ROOT.TH1F ("hMuTau_Metcut_dPhiMmu", ";#delta#phi_{#tau};", 100, -pi, pi)

# #-----dPhi cuts           
# h['hMuTaudPhi_M'] = ROOT.TH1F ("hMuTau_dPhi_M", "#tau - #tau mass;M_{#tau#tau};", 1000, 0, 200)
# h['hMuTaudPhi_muPt'] = ROOT.TH1F ("hMuTau_dPhi_muPt", "#tau_{1} P_t;P_{t};", 500, 0, 500)
# h['hMuTaudPhi_tauPt'] = ROOT.TH1F ("hMuTau_dPhi_tauPt", "#tau_{2} P_t;P_{t};", 500, 0, 500)
# h['hMuTaudPhi_jPt'] = ROOT.TH1F ("hMuTau_dPhi_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hMuTaudPhi_dR'] = ROOT.TH1F ("hMuTau_dPhi_dR", ";#Delta R;", 100, 0, 5)
# h['hMuTaudPhi_dRlj'] = ROOT.TH1F ("hMuTau_dPhi_dRlj", ";#Delta R;", 100, 0, 5)
# h['hMuTaudPhi_METPt'] = ROOT.TH1F ("hMuTau_dPhi_METPt", ";p_{T};", 500, 0, 500)
# h['hMuTaudPhi_dPhiMj'] = ROOT.TH1F ("hMuTau_dPhi_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hMuTaudPhi_dPhiMmu'] = ROOT.TH1F ("hMuTau_dPhi_dPhiMmu", ";#delta#phi_{#tau};", 100, -pi, pi)

# #-----genMatching
# h['hMuTauGen_M'] = ROOT.TH1F ("hMuTau_Gen_M", "#tau - #tau mass;M_{#tau#tau};", 1000, 0, 200)
# h['hMuTauGen_muPt'] = ROOT.TH1F ("hMuTau_Gen_muPt", "#tau_{1} P_t;P_{t};", 500, 0, 500)
# h['hMuTauGen_tauPt'] = ROOT.TH1F ("hMuTau_Gen_tauPt", "#tau_{2} P_t;P_{t};", 500, 0, 500)

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

        event.getByLabel(labelTauMuonCleaned, handleTauMuonCleaned)
        tausMuonCleaned=handleTauMuonCleaned.product()

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
        
        event.getByLabel(labelPatMETs, handlePatMETs)
        met=handlePatMETs.product().front()
        mets=[]
        mets+=[met]
        mets.sort(key=lambda x: x.pt(), reverse=True)
        
        h['hNEvent'].Fill(0.5, 1)
        h['hNEvent'].Fill(1.5, genweight)

#<genInfo>

        # genMuons=[]
        # genTaus=[]
        # genElectrons=[]
        # genNutaus=[]

        # for particle in particles:
        #     if abs(particle.pdgId())==15 and particle.mother().pdgId()==9999:
        #         genTaus+=[particle]
        #     if abs(particle.pdgId())==11 and particle.isDirectHardProcessTauDecayProductFinalState():
        #         genElectrons+=[particle]
        #     if abs(particle.pdgId())==13 and particle.isDirectHardProcessTauDecayProductFinalState():
        #         genMuons+=[particle]
        #     if abs(particle.pdgId())==16 and particle.isDirectHardProcessTauDecayProductFinalState():
        #         genNutaus+=[particle]
        
        # event.getByLabel(labelGenJet, handleGenJet)
        # #jets=handleGenJet.product()

        # genJets=[]
        # for jet in handleGenJet.product():
        #     genJets+=[jet]

        # genJets.sort(key=lambda x: x.pt(), reverse=True)


#</genInfo>
    
#<muonSelection>

        selected_muons=[]
        selected_muons_iso=[]
        for muon in muons:
            h['hMuPt'].Fill(muon.pt(),genweight)
            selected_muons+=[muon]
            if muonIsoCut(muon)<0.25:
                h['hMuPt_iso'].Fill(muon.pt(),genweight)
                selected_muons_iso+=[muon]

        selected_muons.sort(key=lambda x: x.pt(), reverse=True) 
        selected_muons_iso.sort(key=lambda x: x.pt(), reverse=True) 

#<\muonSelection>

#<electronSelection>

        selected_electrons=[]
        for electron in electrons:
            selected_electrons+=[electron]

        selected_electrons.sort(key=lambda x: x.pt(), reverse=True) 

#<\electronSelection>

        selected_mTaus=[]
        for mTau in tausMuonCleaned:
            if mTau.tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"): h['hTauPt_raw'].Fill(mTau.pt(),genweight)
            if mTau.tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vvloose'].Fill(mTau.pt(),genweight)
            if mTau.tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vloose'].Fill(mTau.pt(),genweight)
            if mTau.tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_loose'].Fill(mTau.pt(),genweight)
            if mTau.tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"): 
                h['hTauPt_medium'].Fill(mTau.pt(),genweight)
                selected_mTaus+=[mTau]
            if mTau.tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_tight'].Fill(mTau.pt(),genweight)
            if mTau.tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vtight'].Fill(mTau.pt(),genweight)
            if mTau.tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"): h['hTauPt_vvtight'].Fill(mTau.pt(),genweight)

        selected_mTaus.sort(key=lambda x: x.pt(), reverse=True)

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


#----------------------------- tau_h tau_mu selection ------------------------------


        if len(selected_muons)>=1 and len(selected_mTaus)>=1 and len(selected_jets)>=1 and selected_muons[0].charge()*selected_mTaus[0].charge()<0:
            h['hMuTauBaseline_NbJets'].Fill(len(selected_bjets), genweight)

            if len(selected_bjets)==0:

                if len(selected_electrons)>=1: h['hMuTauBaseline_ePt'].Fill(selected_electrons[0].pt(), genweight)

                mtau=ROOT.TLorentzVector()
                mtau.SetPtEtaPhiM(selected_mTaus[0].pt(), selected_mTaus[0].eta(), selected_mTaus[0].phi(), selected_mTaus[0].mass())

                mu=ROOT.TLorentzVector()
                mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass()) 

                m=ROOT.TLorentzVector()
                m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass()) 

                h['hMuTauBaseline_M'].Fill((mu+mtau).M(), genweight)

                if len(selected_jets)==1:
                    j=ROOT.TLorentzVector()
                    j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                else:
                    j1=ROOT.TLorentzVector()
                    j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                    j2=ROOT.TLorentzVector()
                    j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())

                    if mtau.DeltaR(j1) > 0.3:
                        j = j1
                    else:
                        j = j2

                if len(selected_muons_iso)>=1 and selected_muons_iso[0].charge()*selected_mTaus[0].charge()<0:
                    isomu=ROOT.TLorentzVector()
                    isomu.SetPtEtaPhiM(selected_muons_iso[0].pt(), selected_muons_iso[0].eta(), selected_muons_iso[0].phi(), selected_muons_iso[0].mass())
                    h['hIsoMuTauBaseline_M'].Fill((isomu+mtau).M(), genweight)

                h['hMuTauBaseline_M'].Fill((mu+mtau).M(), genweight)
                h['hMuTauBaseline_muPt'].Fill(mu.Pt(), genweight)
                h['hMuTauBaseline_tauPt'].Fill(mtau.Pt(), genweight)
                h['hMuTauBaseline_jPt'].Fill(j.Pt(), genweight)
                h['hMuTauBaseline_dR'].Fill(mu.DeltaR(mtau), genweight)
                h['hMuTauBaseline_dRlj'].Fill((mu+mtau).DeltaR(j), genweight)
                h['hMuTauBaseline_METPt'].Fill(m.Pt(), genweight)
                h['hMuTauBaseline_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                h['hMuTauBaseline_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)
                h['hMuTauBaseline_dPhiMtau'].Fill(m.DeltaPhi(mtau), genweight)

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

            #     h['hMuTauTrig_M'].Fill((mu+tau).M(), genweight)
            #     h['hMuTauTrig_muPt'].Fill(mu.Pt(), genweight)
            #     h['hMuTauTrig_tauPt'].Fill(tau.Pt(), genweight)
            #     h['hMuTauTrig_jPt'].Fill(j.Pt(), genweight)
            #     h['hMuTauTrig_dR'].Fill(mu.DeltaR(tau), genweight)
            #     h['hMuTauTrig_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #     h['hMuTauTrig_METPt'].Fill(m.Pt(), genweight)
            #     h['hMuTauTrig_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #     h['hMuTauTrig_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)
                
            # if mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8:

            #     h['hMuTaudR_M'].Fill((mu+tau).M(), genweight)
            #     h['hMuTaudR_muPt'].Fill(mu.Pt(), genweight)
            #     h['hMuTaudR_tauPt'].Fill(tau.Pt(), genweight)
            #     h['hMuTaudR_jPt'].Fill(j.Pt(), genweight)
            #     h['hMuTaudR_dR'].Fill(mu.DeltaR(tau), genweight)
            #     h['hMuTaudR_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #     h['hMuTaudR_METPt'].Fill(m.Pt(), genweight)
            #     h['hMuTaudR_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #     h['hMuTaudR_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

            #     if m.Pt()>100:

            #         h['hMuTauMetcut_M'].Fill((mu+tau).M(), genweight)
            #         h['hMuTauMetcut_muPt'].Fill(mu.Pt(), genweight)
            #         h['hMuTauMetcut_tauPt'].Fill(tau.Pt(), genweight)
            #         h['hMuTauMetcut_jPt'].Fill(j.Pt(), genweight)
            #         h['hMuTauMetcut_dR'].Fill(mu.DeltaR(tau), genweight)
            #         h['hMuTauMetcut_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #         h['hMuTauMetcut_METPt'].Fill(m.Pt(), genweight)
            #         h['hMuTauMetcut_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #         h['hMuTauMetcut_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

            #         if abs(m.DeltaPhi(mu))<1 and abs(m.DeltaPhi(j))>2:
            #             h['hMuTaudPhi_M'].Fill((mu+tau).M(), genweight)
            #             h['hMuTaudPhi_muPt'].Fill(mu.Pt(), genweight)
            #             h['hMuTaudPhi_tauPt'].Fill(tau.Pt(), genweight)
            #             h['hMuTaudPhi_jPt'].Fill(j.Pt(), genweight)
            #             h['hMuTaudPhi_dR'].Fill(mu.DeltaR(tau), genweight)
            #             h['hMuTaudPhi_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #             h['hMuTaudPhi_METPt'].Fill(m.Pt(), genweight)
            #             h['hMuTaudPhi_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #             h['hMuTaudPhi_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

            #             if len(genMuons)==2 and len(genJets)>=1 and len(genElectrons)==0:
            #                 genMu=ROOT.TLorentzVector()
            #                 genMu.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

            #                 genTau=ROOT.TLorentzVector()
            #                 genTau.SetPtEtaPhiM(genMuons[1].pt(), genMuons[1].eta(), genMuons[1].phi(), genMuons[1].mass())

            #                 if (mu.DeltaR(genMu)<0.3 or mu.DeltaR(genTau)<0.3) and (tau.DeltaR(genMu)<0.3 or tau.DeltaR(genTau)<0.3): 
            #                     h['hMuTauGen_M'].Fill((mu+tau).M(), genweight)
            #                     h['hMuTauGen_muPt'].Fill(mu.Pt(), genweight)
            #                     h['hMuTauGen_tauPt'].Fill(tau.Pt(), genweight)

            #             if j.Pt() > 500 :
            #                 h['hMuTauTrig_M'].Fill((mu+tau).M(), genweight)
            #                 h['hMuTauTrig_muPt'].Fill(mu.Pt(), genweight)
            #                 h['hMuTauTrig_tauPt'].Fill(tau.Pt(), genweight)
            #                 h['hMuTauTrig_jPt'].Fill(j.Pt(), genweight)
            #                 h['hMuTauTrig_dR'].Fill(mu.DeltaR(tau), genweight)
            #                 h['hMuTauTrig_dRlj'].Fill((mu+tau).DeltaR(j), genweight)
            #                 h['hMuTauTrig_METPt'].Fill(m.Pt(), genweight)
            #                 h['hMuTauTrig_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
            #                 h['hMuTauTrig_dPhiMmu'].Fill(m.DeltaPhi(mu), genweight)

                        
out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
