import ROOT, sys, os
import numpy as np

ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/JetInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/MuonInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/ElectronInfoDS.h"')

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

if len(sys.argv)>2:
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./output/"

outputFileName = outputFileDir+"h_TriggerStudy_EMu_"+inputFileListName.split("/")[-1].replace(".txt",".root")

out=ROOT.TFile.Open(outputFileName,'recreate')

fchain = ROOT.TChain('tcpNtuples/analysisTree')

pi = np.pi

h = {}

h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)

h['hEMuJetHT'] = ROOT.TH1F ("hEMu_JetHT", "JetHT; M_{e#mu}; a.u.", 100, 0, 100)
h['hEMuSingleE'] = ROOT.TH1F ("hEMu_SingleE", "SingleE; M_{e#mu}; a.u.", 100, 0, 100)
h['hEMuSingleMu'] = ROOT.TH1F ("hEMu_SingleMu", "SingleMu; M_{e#mu}; a.u.", 100, 0, 100)
h['hEMuMuonEG'] = ROOT.TH1F ("hEMu_MuonEG", "MuonEG; M_{e#mu}; a.u.", 100, 0, 100)
h['hEMuMetHT'] = ROOT.TH1F ("hEMu_MetHT", "MetHT; M_{e#mu}; a.u.", 100, 0, 100)

h['hEMuJetHT_only'] = ROOT.TH1F ("hEMu_JetHT_only", "JetHT only; M_{e#mu}; a.u.", 100, 0, 100)
h['hEMuSingleE_only'] = ROOT.TH1F ("hEMu_SingleE_only", "SingleE only; M_{e#mu}; a.u.", 100, 0, 100)
h['hEMuSingleMu_only'] = ROOT.TH1F ("hEMu_SingleMu_only", "SingleMu only; M_{e#mu}; a.u.", 100, 0, 100)
h['hEMuMuonEG_only'] = ROOT.TH1F ("hEMu_MuonEG_only", "MuonEG only; M_{e#mu}; a.u.", 100, 0, 100)
h['hEMuMetHT_only'] = ROOT.TH1F ("hEMu_MetHT_only", "MetHT only; M_{e#mu}; a.u.", 100, 0, 100)

for key in h.keys():
    h[key].Sumw2()

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    inputFileName=inputFileName.replace("\n","")
    print(inputFileName.replace("\n",""))

    fchain.Add(inputFileName)
    fchain.AddFriend('tcpTrigNtuples/triggerTree', inputFileName)
    fchain.AddFriend('lumiSummary/lumiTree', inputFileName)

    jets = ROOT.JetInfoDS()
    muons = ROOT.MuonInfoDS()
    electrons = ROOT.ElectronInfoDS()

    fchain.SetBranchAddress("Jets", ROOT.AddressOf(jets))
    fchain.SetBranchAddress("Muons", ROOT.AddressOf(muons))
    fchain.SetBranchAddress("Electrons", ROOT.AddressOf(electrons))

    for iev in range(fchain.GetEntries()): # Be careful!!!                                                                                                   
       fchain.GetEntry(iev)

       mets = fchain.GetBranch("Mets")
       met_pt = mets.GetLeaf('pt').GetValue()
       met_phi = mets.GetLeaf('phi').GetValue()

       weight = fchain.GetBranch("lumiInfo")
       genweight = weight.GetLeaf('weight').GetValue()

       h['hEvents'].Fill(0.5, 1)
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

       s_j = []
       s_b = []
       s_e = []
       s_mu = []

       if jets.size()>0:
          for i in range(jets.size()):
             jet = jets.at(i)
             if jet.id == 2:
                 s_j+=[jet]
                 if jet.deepjet > 0.7476:
                     s_b+=[jet]

       if muons.size()>0:
          for i in range(muons.size()):
             muon = muons.at(i)
             if muon.id >= 1:
                 if muon.iso < 0.25:
                     s_mu+=[muon]

       if electrons.size()>0:
          for i in range(electrons.size()):
             electron = electrons.at(i)
             if electron.id >= 1 :
                 s_e+=[electron]

       s_j.sort(key=lambda x: x.pt, reverse=True)
       s_e.sort(key=lambda x: x.pt, reverse=True)
       s_mu.sort(key=lambda x: x.pt, reverse=True)


       if len(s_e) > 0 and len(s_j) > 0 and len(s_mu) > 0 and s_e[0].charge*s_mu[0].charge < 0 and len(s_b) == 0:

           mu = ROOT.TLorentzVector()
           mu.SetPtEtaPhiM(s_mu[0].pt, s_mu[0].eta, s_mu[0].phi, s_mu[0].mass)
           e = ROOT.TLorentzVector()
           e.SetPtEtaPhiM(s_e[0].pt, s_e[0].eta, s_e[0].phi, s_e[0].mass)
           j = ROOT.TLorentzVector()
           j.SetPtEtaPhiM(s_j[0].pt, s_j[0].eta, s_j[0].phi, s_j[0].mass)
           m = ROOT.TLorentzVector()
           m.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

           if mu.DeltaR(e) < 0.4 and mu.DeltaR(j) > 0.8 and e.DeltaR(j) > 0.8:
               if met_pt > 100:
                   if abs(m.DeltaPhi(mu))<1 and abs(m.DeltaPhi(j))>2:

                        if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                           or (mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ):
                            h['hEMuSingleMu'].Fill((e+mu).M(),  genweight)

                        if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                           or (mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) \
                           or (e.Pt() > 35 and (isIsoEle == 1 or isEleTau == 1 ) ):
                            h['hEMuSingleE'].Fill((e+mu).M(),  genweight)

                        if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                           or (mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) \
                           or (e.Pt() > 35 and (isIsoEle == 1 or isEleTau == 1 ) ) \
                           or ( isMuonEG == 1 and ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) ) \
                           or ( j.Pt() > 500 and isHTMHT == 1 ) :
                            h['hEMuMetHT'].Fill((e+mu).M(),  genweight)

                        if ( j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1 ) ) \
                           or ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) \
                           or ( e.Pt() > 35 and (isIsoEle == 1 or isEleTau == 1 ) ) \
                           or ( isMuonEG == 1 and ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) ) :
                            h['hEMuMuonEG'].Fill((e+mu).M(),  genweight)

                        if j.Pt() > 500 and ( isSingleJet == 1 or isHT == 1):
                            h['hEMuJetHT'].Fill((e+mu).M(),  genweight)

                        if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ):
                            h['hEMuSingleMu_only'].Fill((e+mu).M(),  genweight)

                        if e.Pt() > 35 and (isIsoEle == 1 or isEleTau == 1 ):
                            h['hEMuSingleE_only'].Fill((e+mu).M(),  genweight)

                        if isMuonEG == 1 and ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ):
                            h['hEMuMuonEG_only'].Fill((e+mu).M(),  genweight)

                        if j.Pt() > 500 and isHTMHT == 1 :
                            h['hEMuMetHT_only'].Fill((e+mu).M(),  genweight)


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
