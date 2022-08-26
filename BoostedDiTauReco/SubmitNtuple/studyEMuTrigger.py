import ROOT, sys, os
import numpy as np
import time

start_time = time.time()

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

outputTitle = "h_studyEMuTrigger"

if "-mc" in opts:
    isData == 0

ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/JetInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/MuonInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/ElectronInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/TauInfoDS.h"')

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

if len(sys.argv)>2 and not sys.argv[2].startswith("-"):
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./output/"

outputFileName = outputFileDir+outputTitle+"_"+inputFileListName.split("/")[-1].replace(".txt",".root")

out=ROOT.TFile.Open(outputFileName,'recreate')
print(outputFileName)

fchain = ROOT.TChain('tcpNtuples/analysisTree')
chain2 = ROOT.TChain('tcpTrigNtuples/triggerTree')
if isData == 0:
    chain3 = ROOT.TChain('lumiSummary/lumiTree')
    chain4 = ROOT.TChain('tcpGenNtuples/genTree')

pi = np.pi

h = {}

def define_event_histogram(region):

    h[region+"_Mass"] = ROOT.TH1F (region+"_Mass", region+"_Mass ; M_{vis.} (GeV) ; Events ", 100, 0, 100)
    h[region+"_Lepton1Pt"] = ROOT.TH1F (region+"_Lepton1Pt", region+"_Lepton1Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_Lepton2Pt"] = ROOT.TH1F (region+"_Lepton2Pt", region+"_Lepton2Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_JetPt"] = ROOT.TH1F (region+"_JetPt", region+"_JetPt ; JetP_{T} (GeV) ; Events ", 2000, 0, 2000)
    h[region+"_Mt"] = ROOT.TH1F (region+"_Mt", region+"_Mt ; M_{T} (GeV) ; Events ", 150, 0, 150)
    h[region+"_MetPt"] = ROOT.TH1F (region+"_MetPt", region+"_MetPt ; MET (GeV) ; Events ", 500, 0, 500)
    h[region+"_Nj"] = ROOT.TH1F (region+"_Nj", region+"_Nj ; N_{j} ; Events ", 10, 0, 10)
    h[region+"_dRl"] = ROOT.TH1F (region+"_dRl", region+"_dRl; dR(leptons); Events", 100, 0, 5)
    h[region+"_dRj"] = ROOT.TH1F (region+"_dRj", region+"t_dRj; dR(jet, ditau); Events", 100, 0, 5)


def define_general_histogram():

#-----Trigger bits-----

    h['isSingleJet'] = ROOT.TH1F("isSingleJet", "isSingleJet ; isSingleJet ; N", 4,-1.5,2.5)
    h['isHT'] = ROOT.TH1F("isHT", "isHT ; isHT ; N", 4,-1.5,2.5)

#-----Event counts-----

    h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)

    h['hMuMu_Events'] = ROOT.TH1F ("hMuMu_Events", "hMuMu_Events;;N", 6, 1, 7)
    h['hEMu_Events'] = ROOT.TH1F ("hEMu_Events", "hEMu_Events;;N", 6, 1, 7)

    h['hMuMu_Trigger'] = ROOT.TH1F ("hMuMu_Trigger", "hMuMu_Trigger;;N", 7, 1, 8)
    h['hEMu_Trigger'] = ROOT.TH1F ("hEMu_Trigger", "hEMu_Trigger;;N", 7, 1, 8)

    h['hMuMu_Inclusive_Trigger'] = ROOT.TH1F ("hMuMu_Inclusive_Trigger", "hMuMu_Inclusive_Trigger;;N", 5, 1, 6)
    h['hEMu_Inclusive_Trigger'] = ROOT.TH1F ("hEMu_Inclusive_Trigger", "hEMu_Inclusive_Trigger;;N", 5, 1, 6)

    h['hMuMu_SR_dRcut_highMET_dPhicut'] = ROOT.TH1F ("hMuMu_SR_dRcut_highMET_dPhicut", "hMuMu_SR_dRcut_highMET_dPhicut; ;N", 100, 0, 100)
    h['hEMu_SR_dRcut_highMET_dPhicut'] = ROOT.TH1F ("hEMu_SR_dRcut_highMET_dPhicut", "hEMu_SR_dRcut_highMET_dPhicut; ;N", 100, 0, 100)

#-----Objects-----

    h['hJetPt'] = ROOT.TH1F ("hJetPt", "Jet P_{T} ; P_{T} ; N", 1000, 0, 1000)
    h['hBtag'] = ROOT.TH1F ("hBtag", "DeepJet Score; score; N", 50, 0, 1)
    h['hBJetPt'] = ROOT.TH1F ("hBJetPt", "B-tagged Jet P_{T} ; P_{T} ; N", 1000, 0, 1000)
    h['hMuonPt'] = ROOT.TH1F ("hMuPt", "Muon P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hIsoMuonPt'] = ROOT.TH1F ("hIsoMuPt", "Isolated Muon P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hElectronPt'] = ROOT.TH1F ("hEPt", "Electron P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hIsoElectronPt'] = ROOT.TH1F ("hIsoEPt", "Isolated Electron P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hMetPt'] = ROOT.TH1F ("hMetPt", "MET P_{T} ; P_{T} ; N", 500, 0, 500)

    h['hTauUnCleanedPt'] = ROOT.TH1F ("hTauUnCleanedPt", "Uncleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hMuCleanedPt'] = ROOT.TH1F ("hMuCleanedPt", "Mu-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hMuCleanedPt_altered'] = ROOT.TH1F ("hMuCleanedPt_altered", "Mu-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hECleanedPt'] = ROOT.TH1F ("hECleanedPt", "E-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hTauBoostedPt'] = ROOT.TH1F ("hTauBoostedPt", "Boosted Tau P_{T}; P_{T}; a.u.", 500, 0, 500)



for key in h.keys():
    h[key].Sumw2()



def MuMu_Channel(mu):

   isMuMu = 0
   h['hMuMu_Events'].Fill(1, genweight)

   mu1 = ROOT.TLorentzVector()
   mu1.SetPtEtaPhiM(mu[0].pt, mu[0].eta, mu[0].phi, mu[0].mass)

   mu2 = ROOT.TLorentzVector()
   mu2.SetPtEtaPhiM(mu[1].pt, mu[1].eta, mu[1].phi, mu[1].mass)

   j = ROOT.TLorentzVector()
   j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

   m = ROOT.TLorentzVector()
   m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

   if j.Pt() > 100 and (mu1+mu2).M() > 1:   

       isJetHTEvent = 0
       isSingleMuonEvent = 0

       if met_pt < 100 and len(s_b) == 0:
           if j.Pt() > 500 and isSingleJet == 1 :
               h['hMuMu_Trigger'].Fill(1, genweight)
           if j.Pt() > 500 and isHT == 1 :
               h['hMuMu_Trigger'].Fill(2, genweight)
           if j.Pt() > 600 and isSingleJet == 1 :
               h['hMuMu_Trigger'].Fill(3, genweight)
           if j.Pt() > 600 and isHT == 1 :
               h['hMuMu_Trigger'].Fill(4, genweight)
           if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) :
               h['hMuMu_Trigger'].Fill(5, genweight)
           if ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) :
               h['hMuMu_Trigger'].Fill(6, genweight)

       if plotHTTrig == 1 :
           if ( j.Pt() > 500 and isHT == 1 ) : isJetHTEvent = 1

       if plotJetHT == 1 :
           if ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) : isJetHTEvent = 1

       if plotSingleMuon == 1 :
           if ( mu1.Pt() > 50 and isMu == 1 ) or ( mu1.Pt() > 27 and isIsoMu == 1 ) : isSingleMuonEvent = 1

       if met_pt < 100 and len(s_b) == 0:
           if isJetHTEvent == 1 and isSingleMuonEvent == 1:
               h['hMuMu_Inclusive_Trigger'].Fill(1, genweight)

           if isJetHTEvent == 1 and isSingleMuonEvent == 0:
               h['hMuMu_Inclusive_Trigger'].Fill(2, genweight)

           if isJetHTEvent == 0 and isSingleMuonEvent == 1:
               h['hMuMu_Inclusive_Trigger'].Fill(3, genweight)

           if isJetHTEvent == 1 or isSingleMuonEvent == 1:
               h['hMuMu_Inclusive_Trigger'].Fill(4, genweight)

       if ( isData == 0 and ( ( plotHTTrig == 1 and isJetHTEvent == 1 ) or ( plotJetHT == 1 and isJetHTEvent == 1 ) or ( plotSingleMuon == 1 and isSingleMuonEvent == 1 ) ) ) or \
          ( isData == 1 and ( isJetHTDataset == 1 and ( isJetHTEvent == 1 and isSingleMuonEvent == 0 ) ) ) or \
          ( isData == 1 and ( isSingleMuonDataset == 1 and isSingleMuonEvent == 1 ) ):

           # print("data? :", isData)
           # print("isJetHTDataset, isJetHTEvent, isSingleMuonEvent :", isJetHTDataset, isJetHTEvent, isSingleMuonEvent)
           # print("isSingleMuonDataset, isJetHTEvent, isSingleMuonEvent, ", isSingleMuonDataset, isJetHTEvent , isSingleMuonEvent)

           h['hMuMu_Events'].Fill(2, genweight)

           if len(s_b) == 0:
                if mu1.DeltaR(mu2) < 0.4 and mu1.DeltaR(j) > 0.8  and mu2.DeltaR(j) > 0.8:
                   if isData == 0 :
                       plot1Dhist('hMuMu_SR_dRcut', mu1, mu2, j, m)
                   if met_pt > 100:
                       if isData == 0 :
                           plot1Dhist('hMuMu_SR_dRcut_highMET', mu1, mu2, j, m)
                       if abs(m.DeltaPhi(mu1)) < 1 and abs(m.DeltaPhi(j)) > 2:
                           isMuMu = 1
                           h['ChOverlap'].Fill(7,1)
                           if isData == 0 :
                               h['hTauEvents'].Fill(2,genweight)
                               plot1Dhist('hMuMu_SR_dRcut_highMET_dPhicut', mu1, mu2, j, m)
#                           h['hMuMu_SR_dRcut_highMET_dPhicut'].Fill((mu1+mu2).M(), genweight)


   return isMuMu

def EMu_Channel(ele,mu_emu):

   isEMu = 0
   h['hEMu_Events'].Fill(1, genweight)

   mu = ROOT.TLorentzVector()
   mu.SetPtEtaPhiM(mu_emu[0].pt, mu_emu[0].eta, mu_emu[0].phi, mu_emu[0].mass)

   e = ROOT.TLorentzVector()
   e.SetPtEtaPhiM(ele[0].pt, ele[0].eta, ele[0].phi, ele[0].mass)

   j = ROOT.TLorentzVector()
   j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)

   m = ROOT.TLorentzVector()
   m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

   isJetHTEvent = 0
   isSingleMuonEvent = 0

   if met_pt < 100 and (e+mu).M() > 1:
       if j.Pt() > 500 and isSingleJet == 1 :
           h['hEMu_Trigger'].Fill(1, genweight)
       if j.Pt() > 500 and isHT == 1 :
           h['hEMu_Trigger'].Fill(2, genweight)
       if j.Pt() > 600 and isSingleJet == 1 :
           h['hEMu_Trigger'].Fill(3, genweight)
       if j.Pt() > 600 and isHT == 1 :
           h['hEMu_Trigger'].Fill(4, genweight)
       if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) :
           h['hEMu_Trigger'].Fill(5, genweight)
       if ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) :
           h['hEMu_Trigger'].Fill(6, genweight)

   if plotHTTrig == 1 and (e+mu).M() > 1 :
       if ( j.Pt() > 500 and isHT == 1 ) : isJetHTEvent = 1

   if plotJetHT == 1 and (e+mu).M() > 1:
       if ( j.Pt() > 600 and ( isSingleJet == 1 or isHT == 1 ) ) : isJetHTEvent = 1

   if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : isSingleMuonEvent = 1

   if met_pt < 100 and len(s_b) == 0 and j.Pt() > 100 and (e+mu).M() > 1:
       if isJetHTEvent == 1 and isSingleMuonEvent == 1:
           h['hEMu_Inclusive_Trigger'].Fill(1, genweight)

       if isJetHTEvent == 1 and isSingleMuonEvent == 0:
           h['hEMu_Inclusive_Trigger'].Fill(2, genweight)

       if isJetHTEvent == 0 and isSingleMuonEvent == 1:
           h['hEMu_Inclusive_Trigger'].Fill(3, genweight)

       if isJetHTEvent == 1 or isSingleMuonEvent == 1:
           h['hEMu_Inclusive_Trigger'].Fill(4, genweight)


   if isData == 0 and j.Pt() > 100 :

       if ( j.Pt() > 500 and ( isHT == 1 ) ) \
          or ( isMuonEG == 1 and ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) ) \
          or ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) \
          or ( e.Pt() > 35 and (isIsoEle == 1) ) :

           if (e+mu).M() > 1.0 :
               if mu.DeltaR(e) < 0.4 and mu.DeltaR(j)> 0.8  and e.DeltaR(j) > 0.8 :
                   plot1Dhist('hEMu_SR_dRcut_allTrig',e,mu,j,m)
                   if met_pt > 100.0 :
                       plot1Dhist('hEMu_SR_dRcut_highMET_allTrig',e,mu,j,m)
                       if abs(m.DeltaPhi(mu)) < 1.0 and abs( m.DeltaPhi(j) ) > 2.0 :
                           plot1Dhist('hEMu_SR_dRcut_highMET_dPhicut_allTrig',e,mu,j,m)
                           isEMu = 1

                    
   return isEMu


def Mt(lepton, met):

    cos = np.cos(met.DeltaPhi(lepton))
    Mt = np.sqrt(2*lepton.Pt()*met.Pt()*(1-cos))

    return Mt

def plot1Dhist(region,l1,l2,j,m):

    h[region+"_Lepton1Pt"].Fill(l1.Pt(), genweight)
    h[region+"_Mass"].Fill((l1+l2).M(), genweight)
    h[region+"_Lepton2Pt"].Fill(l2.Pt(), genweight)
    h[region+"_JetPt"].Fill(j.Pt(), genweight)
    h[region+"_Mt"].Fill(Mt(l2,m), genweight)
    h[region+"_MetPt"].Fill(m.Pt(), genweight)
    h[region+"_Nj"].Fill(len(s_j), genweight)



inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    inputFileName=inputFileName.replace("\n","")
    print(inputFileName.replace("\n",""))

    fchain.Add(inputFileName)
    chain2.Add(inputFileName)
    if isData == 0:
        chain3.Add(inputFileName)
        chain4.Add(inputFileName)

fchain.AddFriend(chain2)
if isData == 0:
    fchain.AddFriend(chain3)
    fchain.AddFriend(chain4)

jets = ROOT.JetInfoDS()
muons = ROOT.MuonInfoDS()
electrons = ROOT.ElectronInfoDS()
tausUnCleaned = ROOT.TauInfoDS()
tausECleaned = ROOT.TauInfoDS()
tausMCleaned = ROOT.TauInfoDS()
tausBoosted = ROOT.TauInfoDS()

fchain.SetBranchAddress("Jets", ROOT.AddressOf(jets))
fchain.SetBranchAddress("Muons", ROOT.AddressOf(muons))
fchain.SetBranchAddress("Electrons", ROOT.AddressOf(electrons))
fchain.SetBranchAddress("TausUnCleaned", ROOT.AddressOf(tausUnCleaned))
fchain.SetBranchAddress("TausECleaned", ROOT.AddressOf(tausECleaned))
fchain.SetBranchAddress("TausMCleaned", ROOT.AddressOf(tausMCleaned))
fchain.SetBranchAddress("TausBoosted", ROOT.AddressOf(tausBoosted))

if isData == 0:
    genParticle = ROOT.GenParticleInfoDS()
    fchain.SetBranchAddress("GenParticleInfo", ROOT.AddressOf(genParticle))

for iev in range(fchain.GetEntries()): # Be careful!!!                                                               

   fchain.GetEntry(iev)

   mets = fchain.GetBranch("Mets")
   met_pt = mets.GetLeaf('pt').GetValue()
   met_phi = mets.GetLeaf('phi').GetValue()

   h['hEvents'].Fill(0.5, 1)

   h['hMetPt'].Fill(met_pt)

   genweight = 1

   if isData == 0:
       weight = fchain.GetBranch("lumiInfo")
       genweight = weight.GetLeaf('weight').GetValue()

       s_genmu = []
       
       if genParticle.size()>0:
           for i in range(genParticle.size()):
               gen = genParticle.at(i)
               if abs(gen.pdgid) == 13 : s_genmu += [gen]

   h['hEvents'].Fill(1.5, genweight)

   isSingleJet = fchain.GetLeaf('isSingleJet').GetValue()
   isHT = fchain.GetLeaf('isHT').GetValue()
   isHTMHT = fchain.GetLeaf('isHTMHT').GetValue()
   isMu = fchain.GetLeaf('isMu').GetValue()
   isIsoMu = fchain.GetLeaf('isIsoMu').GetValue()
   isIsoMuTau = fchain.GetLeaf('isIsoMuTau').GetValue()
   isIsoEle = fchain.GetLeaf('isIsoEle').GetValue()
   isEleTau = fchain.GetLeaf('isEleTau').GetValue()
   isMuonEG = fchain.GetLeaf('isMuonEG').GetValue()

   h["isSingleJet"].Fill(isSingleJet, genweight)
   h["isHT"].Fill(isHT, genweight)

   s_j = []
   s_b = []
   s_e = []
   s_lowE = []
   s_mu = []
   s_isomu = []
   s_isoe = []

   unclean = []
   eclean = []
   mclean = []
   lowEclean = []

   boosted = []
   mclean_altered = []

   if jets.size()>0:
      for i in range(jets.size()):
         jet = jets.at(i)
         if jet.id == 2:
             h['hJetPt'].Fill(jet.pt, genweight)
             h['hBtag'].Fill(jet.deepjet)
             s_j+=[jet]
             if jet.deepjet > 0.7476:
                 h['hBJetPt'].Fill(jet.pt, genweight)
                 s_b+=[jet]

   if muons.size()>0:
      for i in range(muons.size()):
         muon = muons.at(i)
         if abs(muon.eta) < 2.4 : 
             if muon.id >= 1:
                 h['hMuonPt'].Fill(muon.pt, genweight) 
                 s_mu+=[muon]
                 if muon.iso < 0.25:
                     h['hIsoMuonPt'].Fill(muon.pt, genweight)
                     s_isomu+=[muon]

   if electrons.size()>0:
      for i in range(electrons.size()):
         electron = electrons.at(i)
         if abs(electron.eta) < 2.5 :
             if electron.id >= 1 :
                 h['hElectronPt'].Fill(electron.pt, genweight)
                 s_e+=[electron]
                 if electron.iso >= 1:
                     h['hIsoElectronPt'].Fill(electron.pt, genweight)
                     s_isoe+=[electron]

   s_j.sort(key=lambda x: x.pt, reverse=True)
   s_e.sort(key=lambda x: x.pt, reverse=True)
   s_isoe.sort(key=lambda x: x.pt, reverse=True)
   s_mu.sort(key=lambda x: x.pt, reverse=True)
   s_isomu.sort(key=lambda x: x.pt, reverse=True)

   isEMu = 0
   isMuMu = 0

   if len(s_isomu) > 1 and len(s_j) > 0 and s_isomu[0].charge*s_isomu[1].charge < 0 : 
       if MuMu_Channel(s_isomu) == 1: continue

   if len(s_isomu) > 0 and len(s_isoe) > 0 and len(s_j) > 0 and len(s_b) == 0 and s_isoe[0].charge*s_isomu[0].charge < 0 : 
       if EMu_Channel(s_isoe,s_isomu) == 1: continue



out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))
