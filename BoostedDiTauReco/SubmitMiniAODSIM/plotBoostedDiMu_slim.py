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

outputFileName = outputFileDir+"h_MuMu_"+inputFileListName.split("/")[-1].replace(".txt",".root")
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
h['hMuMuBaseline_NbJets'] = ROOT.TH1F ("hMuMu_Baseline_NbJets","", 10, 0, 10)
h['hMuMuBaseline_ePt'] = ROOT.TH1F ("hMuMu_Baseline_ePt", "", 500, 0, 500)

#---------------------Tau IDs---------------------

h['hTauPt_raw'] = ROOT.TH1F ("hTauPt_raw", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_loose'] = ROOT.TH1F ("hTauPt_loose", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_medium'] = ROOT.TH1F ("hTauPt_medium", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_tight'] = ROOT.TH1F ("hTauPt_tight", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_vloose'] = ROOT.TH1F ("hTauPt_vloose", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_vtight'] = ROOT.TH1F ("hTauPt_vtight", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_vvloose'] = ROOT.TH1F ("hTauPt_vvloose", "#tau P_t; P_t;", 100, 0, 500)
h['hTauPt_vvtight'] = ROOT.TH1F ("hTauPt_vvtight", "#tau P_t; P_t;", 100, 0, 500)


#-----------------------MuMu-----------------------                                                                                                                                                         
#-----Require 2 muons                                                                                                                                                                                          
h['hMuMuBaseline_M'] = ROOT.TH1F ("hMuMu_Baseline_M", "#mu - #mu mass;M_{#mu#mu};", 500, 0, 200)
h['hMuTauBaseline_M'] = ROOT.TH1F ("hMuTau_Baseline_M", "#mu - #mu mass;M_{#mu#mu};", 500, 0, 200)
h['hMuMuBaseline_mu1Pt'] = ROOT.TH1F ("hMuMu_Baseline_mu1Pt", "#mu_{1} P_t; P_t;", 100, 0, 500)
h['hMuMuBaseline_mu2Pt'] = ROOT.TH1F ("hMuMu_Baseline_mu2Pt", "#mu_{2} P_t; P_t;", 100, 0, 500)
h['hMuMuBaseline_jPt'] = ROOT.TH1F ("hMuMu_Baseline_jPt", "jet pt;P_{t};", 2000, 0, 2000)
h['hMuMuBaseline_dR'] = ROOT.TH1F ("hMuMu_Baseline_dR", ";#Delta R;", 100, 0, 5)
h['hMuMuBaseline_dRlj'] = ROOT.TH1F ("hMuMu_Baseline_dRlj", ";#Delta R;", 100, 0, 5)
h['hMuMuBaseline_METPt'] = ROOT.TH1F ("hMuMu_Baseline_METPt", ";p_{T};", 500, 0, 500)
h['hMuMuBaseline_dPhiMj'] = ROOT.TH1F ("hMuMu_Baseline_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
h['hMuMuBaseline_dPhiMmu'] = ROOT.TH1F ("hMuMu_Baseline_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

#-----Passing Triggers                                                                                                                                                                       
# h['hMuMuTrig_M'] = ROOT.TH1F ("hMuMu_Trig_M", "#mu - #mu mass;M_{#mu#mu};", 1000, 0, 200)
# h['hMuMuTrig_mu1Pt'] = ROOT.TH1F ("hMuMu_Trig_mu1Pt", "#mu_{1} P_t;P_{t};", 500, 0, 500)
# h['hMuMuTrig_mu2Pt'] = ROOT.TH1F ("hMuMu_Trig_mu2Pt", "#mu_{2} P_t;P_{t};", 500, 0, 500)
# h['hMuMuTrig_jPt'] = ROOT.TH1F ("hMuMu_Trig_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hMuMuTrig_dR'] = ROOT.TH1F ("hMuMu_Trig_dR", ";#Delta R;", 100, 0, 5)
# h['hMuMuTrig_dRlj'] = ROOT.TH1F ("hMuMu_Trig_dRlj", ";#Delta R;", 100, 0, 5)
# h['hMuMuTrig_METPt'] = ROOT.TH1F ("hMuMu_Trig_METPt", ";p_{T};", 500, 0, 500)
# h['hMuMuTrig_dPhiMj'] = ROOT.TH1F ("hMuMu_Trig_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hMuMuTrig_dPhiMmu'] = ROOT.TH1F ("hMuMu_Trig_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

# #-----dR cuts                                                                                                                                                                                 
# h['hMuMudR_M'] = ROOT.TH1F ("hMuMu_dR_M", "#mu - #mu mass;M_{#mu#mu};", 1000, 0, 200)
# h['hMuMudR_mu1Pt'] = ROOT.TH1F ("hMuMu_dR_mu1Pt", "#mu_{1} P_t;P_{t};", 500, 0, 500)
# h['hMuMudR_mu2Pt'] = ROOT.TH1F ("hMuMu_dR_mu2Pt", "#mu_{2} P_t;P_{t};", 500, 0, 500)
# h['hMuMudR_jPt'] = ROOT.TH1F ("hMuMu_dR_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hMuMudR_dR'] = ROOT.TH1F ("hMuMu_dR_dR", ";#Delta R;", 100, 0, 5)
# h['hMuMudR_dRlj'] = ROOT.TH1F ("hMuMu_dR_dRlj", ";#Delta R;", 100, 0, 5)
# h['hMuMudR_METPt'] = ROOT.TH1F ("hMuMu_dR_METPt", ";p_{T};", 500, 0, 500)
# h['hMuMudR_dPhiMj'] = ROOT.TH1F ("hMuMu_dR_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hMuMudR_dPhiMmu'] = ROOT.TH1F ("hMuMu_dR_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

# #-----MET cuts                                                                                                                                                                                               
# h['hMuMuMetcut_M'] = ROOT.TH1F ("hMuMu_Metcut_M", "#mu - #mu mass;M_{#mu#mu};", 1000, 0, 200)
# h['hMuMuMetcut_mu1Pt'] = ROOT.TH1F ("hMuMu_Metcut_mu1Pt", "#mu_{1} P_t;P_{t};", 500, 0, 500)
# h['hMuMuMetcut_mu2Pt'] = ROOT.TH1F ("hMuMu_Metcut_mu2Pt", "#mu_{2} P_t;P_{t};", 500, 0, 500)
# h['hMuMuMetcut_jPt'] = ROOT.TH1F ("hMuMu_Metcut_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hMuMuMetcut_dR'] = ROOT.TH1F ("hMuMu_Metcut_dR", ";#Delta R;", 100, 0, 5)
# h['hMuMuMetcut_dRlj'] = ROOT.TH1F ("hMuMu_Metcut_dRlj", ";#Delta R;", 100, 0, 5)
# h['hMuMuMetcut_METPt'] = ROOT.TH1F ("hMuMu_Metcut_METPt", ";p_{T};", 500, 0, 500)
# h['hMuMuMetcut_dPhiMj'] = ROOT.TH1F ("hMuMu_Metcut_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hMuMuMetcut_dPhiMmu'] = ROOT.TH1F ("hMuMu_Metcut_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

# #-----dPhi cuts           
# h['hMuMudPhi_M'] = ROOT.TH1F ("hMuMu_dPhi_M", "#mu - #mu mass;M_{#mu#mu};", 1000, 0, 200)
# h['hMuMudPhi_mu1Pt'] = ROOT.TH1F ("hMuMu_dPhi_mu1Pt", "#mu_{1} P_t;P_{t};", 500, 0, 500)
# h['hMuMudPhi_mu2Pt'] = ROOT.TH1F ("hMuMu_dPhi_mu2Pt", "#mu_{2} P_t;P_{t};", 500, 0, 500)
# h['hMuMudPhi_jPt'] = ROOT.TH1F ("hMuMu_dPhi_jPt", "jet P_t;P_{t};", 2000, 0, 2000)
# h['hMuMudPhi_dR'] = ROOT.TH1F ("hMuMu_dPhi_dR", ";#Delta R;", 100, 0, 5)
# h['hMuMudPhi_dRlj'] = ROOT.TH1F ("hMuMu_dPhi_dRlj", ";#Delta R;", 100, 0, 5)
# h['hMuMudPhi_METPt'] = ROOT.TH1F ("hMuMu_dPhi_METPt", ";p_{T};", 500, 0, 500)
# h['hMuMudPhi_dPhiMj'] = ROOT.TH1F ("hMuMu_dPhi_dPhiMj", ";#delta#phi_{jet};", 100, -pi, pi)
# h['hMuMudPhi_dPhiMmu'] = ROOT.TH1F ("hMuMu_dPhi_dPhiMmu", ";#delta#phi_{#mu};", 100, -pi, pi)

# #-----genMatching
# h['hMuMuGen_M'] = ROOT.TH1F ("hMuMu_Gen_M", "#mu - #mu mass;M_{#mu#mu};", 1000, 0, 200)
# h['hMuMuGen_mu1Pt'] = ROOT.TH1F ("hMuMu_Gen_mu1Pt", "#mu_{1} P_t;P_{t};", 500, 0, 500)
# h['hMuMuGen_mu2Pt'] = ROOT.TH1F ("hMuMu_Gen_mu2Pt", "#mu_{2} P_t;P_{t};", 500, 0, 500)

#</histograms>

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    print inputFileName.replace("\n","")
#    inputFileName = "root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCPAnalysis/RunIISummer17DR94Premix/"+inputFileName.replace("\n","")
    inputFileName = "root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCPAnalysis/Backgrounds/RunIIFall17DR94Premix/"+inputFileName.replace("\n","")
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
        for muon in muons:
#            selected_muons+=[muon]
            if muonIsoCut(muon)<0.25:
                selected_muons+=[muon]

        selected_muons.sort(key=lambda x: x.pt(), reverse=True) 
        #<\muonSelection>

        #<electronSelection>
        selected_electrons=[]
        for electron in electrons:
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

#----------------------------- tau_mu tau_mu selection ------------------------------


        if len(selected_muons)>=2 and len(selected_jets)>=1 and selected_muons[0].charge()*selected_muons[1].charge()<0:
            h['hMuMuBaseline_NbJets'].Fill(len(selected_bjets), genweight)

            if len(selected_bjets)==0:

                m=ROOT.TLorentzVector()
                m.SetPtEtaPhiM(mets[0].pt(), mets[0].eta(), mets[0].phi(), mets[0].mass())

                mu1=ROOT.TLorentzVector()
                mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass()) 

                mu2=ROOT.TLorentzVector()
                mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass())

                j=ROOT.TLorentzVector()
                j.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())

                h['hMuMuBaseline_M'].Fill((mu1+mu2).M(), genweight)
                h['hMuMuBaseline_mu1Pt'].Fill(mu1.Pt(), genweight)
                h['hMuMuBaseline_mu2Pt'].Fill(mu2.Pt(), genweight)
                # h['hMuMuBaseline_jPt'].Fill(j.Pt(), genweight)
                # h['hMuMuBaseline_dR'].Fill(mu1.DeltaR(mu2), genweight)
                # h['hMuMuBaseline_dRlj'].Fill((mu1+mu2).DeltaR(j), genweight)
                # h['hMuMuBaseline_METPt'].Fill(m.Pt(), genweight)
                # h['hMuMuBaseline_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                # h['hMuMuBaseline_dPhiMmu'].Fill(m.DeltaPhi(mu1), genweight)

    #            print triggerResults
    #            print triggerResults.size()
    #            print triggerResults.wasrun()
    #            print names.triggerIndex("HLT_IsoMu24_v10")
    #            print names.triggerIndex("HLT_Mu50_v11")


                # if (mu1.Pt()>26 and (triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")) or triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")))) \
                # or (mu1.Pt()>50 and (triggerResults.accept(names.triggerIndex("HLT_TkMu50_v3")) or triggerResults.accept(names.triggerIndex("HLT_Mu50_v5")))) \
                # or (j.Pt()>500 and (triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) or triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) or triggerResults.accept(names.triggerIndex("HLT_CaloJet500_NoJetID_v5")))):

                # if mu1.Pt()>26 or j.Pt()>500:# and triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v10")):

                # if j.Pt() > 500 :

                #     h['hMuMuTrig_M'].Fill((mu1+mu2).M(), genweight)
                #     h['hMuMuTrig_mu1Pt'].Fill(mu1.Pt(), genweight)
                #     h['hMuMuTrig_mu2Pt'].Fill(mu2.Pt(), genweight)
                #     h['hMuMuTrig_jPt'].Fill(j.Pt(), genweight)
                #     h['hMuMuTrig_dR'].Fill(mu1.DeltaR(mu2), genweight)
                #     h['hMuMuTrig_dRlj'].Fill((mu1+mu2).DeltaR(j), genweight)
                #     h['hMuMuTrig_METPt'].Fill(m.Pt(), genweight)
                #     h['hMuMuTrig_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                #     h['hMuMuTrig_dPhiMmu'].Fill(m.DeltaPhi(mu1), genweight)

                # if mu1.DeltaR(mu2)<0.4 and mu1.DeltaR(j)>0.8 and mu2.DeltaR(j)>0.8:

                #     h['hMuMudR_M'].Fill((mu1+mu2).M(), genweight)
                #     h['hMuMudR_mu1Pt'].Fill(mu1.Pt(), genweight)
                #     h['hMuMudR_mu2Pt'].Fill(mu2.Pt(), genweight)
                #     h['hMuMudR_jPt'].Fill(j.Pt(), genweight)
                #     h['hMuMudR_dR'].Fill(mu1.DeltaR(mu2), genweight)
                #     h['hMuMudR_dRlj'].Fill((mu1+mu2).DeltaR(j), genweight)
                #     h['hMuMudR_METPt'].Fill(m.Pt(), genweight)
                #     h['hMuMudR_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                #     h['hMuMudR_dPhiMmu'].Fill(m.DeltaPhi(mu1), genweight)

                #     if m.Pt()>100:

                #         h['hMuMuMetcut_M'].Fill((mu1+mu2).M(), genweight)
                #         h['hMuMuMetcut_mu1Pt'].Fill(mu1.Pt(), genweight)
                #         h['hMuMuMetcut_mu2Pt'].Fill(mu2.Pt(), genweight)
                #         h['hMuMuMetcut_jPt'].Fill(j.Pt(), genweight)
                #         h['hMuMuMetcut_dR'].Fill(mu1.DeltaR(mu2), genweight)
                #         h['hMuMuMetcut_dRlj'].Fill((mu1+mu2).DeltaR(j), genweight)
                #         h['hMuMuMetcut_METPt'].Fill(m.Pt(), genweight)
                #         h['hMuMuMetcut_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                #         h['hMuMuMetcut_dPhiMmu'].Fill(m.DeltaPhi(mu1), genweight)

                #         if abs(m.DeltaPhi(mu1))<1 and abs(m.DeltaPhi(j))>2:
                #             h['hMuMudPhi_M'].Fill((mu1+mu2).M(), genweight)
                #             h['hMuMudPhi_mu1Pt'].Fill(mu1.Pt(), genweight)
                #             h['hMuMudPhi_mu2Pt'].Fill(mu2.Pt(), genweight)
                #             h['hMuMudPhi_jPt'].Fill(j.Pt(), genweight)
                #             h['hMuMudPhi_dR'].Fill(mu1.DeltaR(mu2), genweight)
                #             h['hMuMudPhi_dRlj'].Fill((mu1+mu2).DeltaR(j), genweight)
                #             h['hMuMudPhi_METPt'].Fill(m.Pt(), genweight)
                #             h['hMuMudPhi_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                #             h['hMuMudPhi_dPhiMmu'].Fill(m.DeltaPhi(mu1), genweight)

                #             if len(genMuons)==2 and len(genJets)>=1 and len(genElectrons)==0:
                #                 genMu1=ROOT.TLorentzVector()
                #                 genMu1.SetPtEtaPhiM(genMuons[0].pt(), genMuons[0].eta(), genMuons[0].phi(), genMuons[0].mass())

                #                 genMu2=ROOT.TLorentzVector()
                #                 genMu2.SetPtEtaPhiM(genMuons[1].pt(), genMuons[1].eta(), genMuons[1].phi(), genMuons[1].mass())

                #                 if (mu1.DeltaR(genMu1)<0.3 or mu1.DeltaR(genMu2)<0.3) and (mu2.DeltaR(genMu1)<0.3 or mu2.DeltaR(genMu2)<0.3): 
                #                     h['hMuMuGen_M'].Fill((mu1+mu2).M(), genweight)
                #                     h['hMuMuGen_mu1Pt'].Fill(mu1.Pt(), genweight)
                #                     h['hMuMuGen_mu2Pt'].Fill(mu2.Pt(), genweight)

                #             if j.Pt() > 500 :
                #                 h['hMuMuTrig_M'].Fill((mu1+mu2).M(), genweight)
                #                 h['hMuMuTrig_mu1Pt'].Fill(mu1.Pt(), genweight)
                #                 h['hMuMuTrig_mu2Pt'].Fill(mu2.Pt(), genweight)
                #                 h['hMuMuTrig_jPt'].Fill(j.Pt(), genweight)
                #                 h['hMuMuTrig_dR'].Fill(mu1.DeltaR(mu2), genweight)
                #                 h['hMuMuTrig_dRlj'].Fill((mu1+mu2).DeltaR(j), genweight)
                #                 h['hMuMuTrig_METPt'].Fill(m.Pt(), genweight)
                #                 h['hMuMuTrig_dPhiMj'].Fill(m.DeltaPhi(j), genweight)
                #                 h['hMuMuTrig_dPhiMmu'].Fill(m.DeltaPhi(mu1), genweight)

                        
out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
