import ROOT, sys, os
import numpy as np
import time

start_time = time.time()

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

outputTitle = "h_studyEMuDataMC"

if "-b" in opts:
    isData = 0

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

event_cut = {
    'jetPt': 100,
    'dRl': 0.4,
    'dRltau': 0.05,
    'dRlj': 0.8,
    'metcut': 100,
    'mtcut': 50,
    'dPhiml': 1,
    'dPhimj': 2,
    'mass' : 1
}

def define_event_histogram(region):

    h[region+"_Mass"] = ROOT.TH1F (region+"_Mass", region+"_Mass ; M_{vis.} (GeV) ; Events ", 100, 0, 100)
    h[region+"_Lepton1Pt"] = ROOT.TH1F (region+"_Lepton1Pt", region+"_Lepton1Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_Lepton2Pt"] = ROOT.TH1F (region+"_Lepton2Pt", region+"_Lepton2Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_JetPt"] = ROOT.TH1F (region+"_JetPt", region+"_JetPt ; JetP_{T} (GeV) ; Events ", 2000, 0, 2000)
    h[region+"_Mt"] = ROOT.TH1F (region+"_Mt", region+"_Mt ; M_{T} (GeV) ; Events ", 150, 0, 150)
    h[region+"_MetPt"] = ROOT.TH1F (region+"_MetPt", region+"_MetPt ; MET (GeV) ; Events ", 500, 0, 500)
    h[region+"_Nj"] = ROOT.TH1F (region+"_Nj", region+"_Nj ; N_{j} ; Events ", 10, 0, 10)
    h[region+"_dR"] = ROOT.TH2F (region+"_dR", region+"_dR ; dR(leptons) ; dR(jet, ditau)", 100, 0, 5, 100, 0, 5)
    h[region+"_dPhil"] = ROOT.TH2F (region+"_dPhil", region+"_dPhil ; dPhi(met,lepton) ; dPhi(met,lepton2)",  100, -pi, pi, 100, -pi, pi)
    h[region+"_dPhi"] = ROOT.TH2F (region+"_dPhi", region+"_dPhi ; dPhi(met,ditau) ; dPhi(met,jet)",  100, -pi, pi, 100, -pi, pi)


def define_general_histogram():

#-----Trigger bits-----

    h['isSingleJet'] = ROOT.TH1F("isSingleJet", "isSingleJet ; isSingleJet ; N", 4,-1.5,2.5)
    h['isHT'] = ROOT.TH1F("isHT", "isHT ; isHT ; N", 4,-1.5,2.5)

#-----Event counts-----

    h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)

    h['hMuMu_Events'] = ROOT.TH1F ("hMuMu_Events", "hMuMu_Events;;N", 6, 1, 7)
    h['hEMu_Events'] = ROOT.TH1F ("hEMu_Events", "hEMu_Events;;N", 6, 1, 7)

    h['hEMu_Trigger_Event'] = ROOT.TH1F ("hEMu_Trigger_Event", "hEMu_Trigger_Events ; ; N", 5, 0, 5)


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


hist_regions = [
    'hMuMu_SR_dRcut_highMET_dPhicut',
    'hMuMu_dRcut_highMET',
    'hMuMu_dRcut',
    'hMuMu_Baseline',
    'hEMu_Baseline',
    'hEMu_dRcut',
    'hEMu_JetHT+SingleMu+SingleE+MuonEG',
    'hEMu_JetHT+SingleMu',
    'hEMu_JetHT+SingleE',
    'hEMu_JetHT+MuonEG',
    'hEMu_JetHT+SingleMu+MuonEG',
    'hEMu_JetHT+SingleE+MuonEG',
    'hEMu_JetHT+MuonEG+SingleMu',
    'hEMu_JetHT+MuonEG+SingleE'
]

for r in hist_regions:
    define_event_histogram(r)

define_general_histogram()

for key in h.keys():
    h[key].Sumw2()



def get_TLorentzVector(l1, l2, js, met_pt, met_phi):

    vl1 = ROOT.TLorentzVector()
    vl1.SetPtEtaPhiM(l1.pt, l1.eta, l1.phi, l1.mass)

    vl2 = ROOT.TLorentzVector()
    vl2.SetPtEtaPhiM(l2.pt, l2.eta, l2.phi, l2.mass)

    vj = ROOT.TLorentzVector()
    vj.SetPtEtaPhiM(js.pt, js.eta, js.phi, js.mass)

    vm = ROOT.TLorentzVector()
    vm.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

    return vl1, vl2, vj, vm


def pass_deltaR(l1,l2,j, channel):

    if channel == "muTau" or channel == "eTau":
        if l1.DeltaR(l2) < event_cut["dRl"] and j.DeltaR(l1) > event_cut["dRlj"] and j.DeltaR(l2) > event_cut["dRlj"] and l1.DeltaR(l2) > event_cut["dRltau"]:
            return 1
        else:
            return -9999
    else:
        if l1.DeltaR(l2) < event_cut["dRl"] and j.DeltaR(l1) > event_cut["dRlj"] and j.DeltaR(l2) > event_cut["dRlj"]:
            return 1
        else:
            return -9999


def pass_baseline(l1, l2, j):

    if j.Pt() > event_cut['jetPt'] and (l1+l2).M() > event_cut["mass"]:
        return 1
    else:
        return -9999


def MuMu_Channel(mu, js, met_pt, met_phi):

    isMuMu = 0
    h['hMuMu_Events'].Fill(1, genweight)

    mu1, mu2, j, m = get_TLorentzVector(mu[0], mu[1], js[0], met_pt, met_phi)

    if pass_baseline(mu1, mu2, j):

        isJetHTEvent = 0
        isSingleMuonEvent = 0

        if ( j.Pt() > 500 and isHT == 1 ) : isJetHTEvent = 1

        if ( mu1.Pt() > 50 and isMu == 1 ) or ( mu1.Pt() > 27 and isIsoMu == 1 ) : isSingleMuonEvent = 1

        if isJetHTEvent or isSingleMuonEvent :
            plot_event_hist("hMuMu_Baseline", mu1, mu2, j, m)

            if pass_deltaR(mu1, mu2, j, "MuMu"):
                plot_event_hist('hMuMu_dRcut', mu1, mu2, j, m)

                if m.Pt() > 100:
                    plot_event_hist('hMuMu_dRcut_highMET', mu1, mu2, j, m)

                    if abs(m.DeltaPhi(mu1)) < 1 and abs(m.DeltaPhi(j)) > 2:
                        plot_event_hist('hMuMu_SR_dRcut_highMET_dPhicut', mu1, mu2, j, m)
                        isMuMu = 1

    return isMuMu


def EMu_Channel(ele,mu_emu, js, met_pt, met_phi):

    isEMu = 0
    h['hEMu_Events'].Fill(1, genweight)

    e, mu, j, m = get_TLorentzVector(ele[0], mu_emu[0], js[0], met_pt, met_phi)

    trigger = [0,0]

    if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : trigger[0] = 1
    
    if ( ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) and isMuonEG == 1 ) : trigger[1] = 1

    if trigger[0] == 1 or trigger[1] == 1:

        if pass_baseline(e, mu, j) == 1:
            if pass_deltaR(e, mu, j, "EMu") == 1:
                if m.Pt() > 100 :

                    h["hEMu_SR_Mtl1"].Fill(Mt(l1,m), genweight)
                    h["hEMu_SR_Mtl2"].Fill(Mt(l2,m), genweight)
                    h["hEMu_SR_Mtl"].Fill(Mt((l1+l2),m), genweight)

                    h['hEMu_SR_cosEm'].Fill(np.cos(m.DeltaPhi(e)), genweight)
                    h['hEMu_SR_cosMum'].Fill(np.cos(m.DeltaPhi(mu)), genweight)
                    
                if m.Pt() < 100 : #lowMET





def EMu_Channel_triggerStudy(ele,mu_emu, js, met_pt, met_phi):

    isEMu = 0
    h['hEMu_Events'].Fill(1, genweight)

    e, mu, j, m = get_TLorentzVector(ele[0], mu_emu[0], js[0], met_pt, met_phi)

    isJetHTEvent = 0
    isSingleMuonEvent = 0
    isSingleEEvent = 0
    isMuonEGEvent = 0

    trigger = [0, 0, 0, 0]

    if ( j.Pt() > 500 and isHT == 1 ) : trigger[0] = 1

    if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : trigger[1] = 1

    if ( e.Pt() > 35 and isIsoEle == 1 ) : trigger[2] = 1

    if ( ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) and isMuonEG == 1 ) : trigger[3] = 1

    if trigger[0] == 1 or trigger[1] == 1 or trigger[2] == 1 or trigger[3] == 1 :

        if pass_baseline(e, mu, j) == 1:
            plot_event_hist("hEMu_Baseline", e, mu, j ,m)

            if pass_deltaR(e, mu, j, "EMu") == 1:
                plot_event_hist("hEMu_dRcut", e, mu, j ,m)

                if m.Pt() > 100.0 :
                    plot_event_hist("hEMu_JetHT+SingleMu+SingleE+MuonEG", e, mu, j ,m)

                    h["hEMu_Trigger_Event"].Fill(4, genweight)

                    if trigger[0] == 1 and ( trigger[1] == trigger[2] == trigger[3] == 0 ):
                        h["hEMu_Trigger_Event"].Fill(0, genweight)
                    if trigger[1] == 1 and ( trigger[0]== trigger[2] == trigger[3] == 0 ):
                        h["hEMu_Trigger_Event"].Fill(1,genweight)
                    if trigger[2] == 1 and ( trigger[1]== trigger[0] == trigger[3] == 0 ):
                        h["hEMu_Trigger_Event"].Fill(2,genweight)
                    if trigger[3] == 1 and ( trigger[1]== trigger[2] == trigger[0] == 0 ):
                        h["hEMu_Trigger_Event"].Fill(3,genweight)

                    if trigger[0] == 1 or trigger[1] == 1:
                        plot_event_hist("hEMu_JetHT+SingleMu", e, mu, j, m)
                    if trigger[0] == 1 or trigger[1] == 1 or trigger[3] == 1:
                        plot_event_hist("hEMu_JetHT+SingleMu+MuonEG", e, mu, j, m)

                    if trigger[0] == 1 or trigger[2] == 1:
                        plot_event_hist("hEMu_JetHT+SingleE", e, mu, j, m)
                    if trigger[0] == 1 or trigger[2] == 1 or trigger[3] == 1:
                        plot_event_hist("hEMu_JetHT+SingleE+MuonEG", e, mu, j, m)

                    if trigger[0] == 1 or trigger[3] == 1:
                        plot_event_hist("hEMu_JetHT+MuonEG", e, mu, j, m)
                    if trigger[0] == 1 or trigger[3] == 1 or trigger[1] == 1:
                        plot_event_hist("hEMu_JetHT+MuonEG+SingleMu", e, mu, j, m)
                    if trigger[0] == 1 or trigger[3] == 1 or trigger[2] == 1: 
                        plot_event_hist("hEMu_JetHT+MuonEG+SingleE", e, mu, j, m)
                    
    return isEMu


def Mt(lepton, met):

    cos = np.cos(met.DeltaPhi(lepton))
    Mt = np.sqrt(2*lepton.Pt()*met.Pt()*(1-cos))

    return Mt

def plot_event_hist(region, l1, l2, j, m):
    
    h[region+"_Mass"].Fill((l1+l2).M(), genweight)
    h[region+"_Lepton1Pt"].Fill(l1.Pt(), genweight)
    h[region+"_Lepton2Pt"].Fill(l2.Pt(), genweight)
    h[region+"_JetPt"].Fill(j.Pt(), genweight)
    h[region+"_MetPt"].Fill(m.Pt(), genweight)
    h[region+"_Nj"].Fill(len(s_j), genweight)
    h[region+"_dR"].Fill(l1.DeltaR(l2), j.DeltaR(l1+l2), genweight)
    h[region+"_dPhil"].Fill(m.DeltaPhi(l1), m.DeltaPhi(l2), genweight)
    h[region+"_dPhi"].Fill(m.DeltaPhi(l1+l2), m.DeltaPhi(j), genweight)


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

   if len(s_isomu) > 1 and len(s_j) > 0 and s_isomu[0].charge*s_isomu[1].charge < 0 and len(s_b) == 0: 
       if MuMu_Channel(s_isomu, s_j, met_pt, met_phi) == 1: continue

   if len(s_isomu) > 0 and len(s_isoe) > 0 and len(s_j) > 0 and len(s_b) == 0 and s_isoe[0].charge*s_isomu[0].charge < 0 : 
       EMu_Channel_triggerStudy(s_isoe,s_isomu, s_j, met_pt, met_phi)

   if len(s_isomu) > 0 and len(s_isoe) > 0 and len(s_j) > 0 and len(s_b) == 0 and s_isoe[0].charge*s_isomu[0].charge < 0 :
       EMu_Channel(s_isoe,s_isomu, s_j, met_pt, met_phi)


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))
