import ROOT, sys, os
import numpy as np
import time

start_time = time.time()

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

outputTitle = "h_studyTTVarDependence"

isNonIso = False
triggerStudy = False
bEffStudy = False

isMuonEGDataset = 0
isSingleMuonDataset = 0

if "-b" in opts:
    isData = 0

if "-de" in opts :
    isData = 1
    isMuonEGDataset = 1
    isSingleMuonDataset = 0

if "-ds" in opts :
    isData = 1
    isMuonEGDataset = 0
    isSingleMuonDataset = 1

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
    'jetPt': 70,
    'dRl': 0.4,
    'dRltau': 0.05,
    'dRlj': 0.8,
    'metcut': 70.0,
    'mtcut': 50.0,
    'dPhiml': 1,
    'dPhimj': 2,
    'mass' : 5
}

met_range = [x for x in range(0,510,10)]
nb_range = [x for x in range(0,6)]
dR_range = [round(x*0.1, 1) for x in range(1, 16)]

h = {}

for metcut in met_range:
    h["hEMu_Nbjet_MetPt"+str(metcut)] = ROOT.TH1F ("hEMu_Nbjet_MetPt_"+str(metcut), "hEMu_Nbjet_MetPt_"+str(metcut)+"; Nbjet ; Events", 6, 0, 6)

for i in nb_range:
    h["hEMu_MetPt_Nbjet"+str(i)] = ROOT.TH1F ("hEMu_MetPt_Nbjet_"+str(i), "hEMu_MetPt_Nbjet_"+str(i)+"; MET (GeV) ; Events", 500, 0, 500)

for dr in dR_range:
    h["hEMu_Nbjet_MetPt_dR"+str(dr)] = ROOT.TH2F ("hEMu_Nbjet_MetPt_dR_"+str(dr), "hEMu_Nbjet_MetPt_dR_"+str(dr)+"; E_{T}^{Miss} ; N_{bjet}", 500, 0, 500, 6, 0, 6)
    h["hEMu_MetPt_bjet0_dR"+str(dr)] = ROOT.TH1F ("hEMu_MetPt_bjet0_dR_"+str(dr), "hEMu_MetPt_bjet0_dR_"+str(dr)+"; E_{T}^{Miss} ; Events", 500, 0, 500)
    h["hEMu_MetPt_bjet_dR"+str(dr)] = ROOT.TH1F ("hEMu_MetPt_bjet_dR_"+str(dr), "hEMu_MetPt_bjet_dR_"+str(dr)+"; E_{T}^{Miss} ; Events", 500, 0, 500)
    h["hEMu_Nbjet_highMET_dR"+str(dr)] = ROOT.TH1F ("hEMu_Nbjet_highMET_dR_"+str(dr), "hEMu_Nbjet_highMET_dR_"+str(dr)+"; N_{bjet} ; Events", 6, 0, 6)
    h["hEMu_Nbjet_lowMET_dR"+str(dr)] = ROOT.TH1F ("hEMu_Nbjet_lowMET_dR_"+str(dr), "hEMu_Nbjet_lowMET_dR_"+str(dr)+"; N_{bjet} ; Events", 6, 0, 6)

    h["hEMu_MetPt_dR"+str(dr)] = ROOT.TH1F ("hEMu_MetPt_dR_"+str(dr), "hEMu_MetPt_dR_"+str(dr)+"; E_{T}^{Miss} ; Event", 500, 0, 500)
    h["hEMu_Nbjet_dR"+str(dr)] = ROOT.TH1F ("hEMu_Nbjet_dR_"+str(dr), "hEMu_Nbjet_dR_"+str(dr)+"; N_{bjet} ; Event", 6, 0, 6)

def define_event_histogram(region):

    h[region+"_Count"] = ROOT.TH1F (region+"_Count", region+"_Count ; Events ; Events ", 1, 0, 1)

    h[region+"_Mass"] = ROOT.TH1F (region+"_Mass", region+"_Mass ; M_{vis.} (GeV) ; Events ", 100, 0, 100)
    h[region+"_Lepton1Pt"] = ROOT.TH1F (region+"_Lepton1Pt", region+"_Lepton1Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_Lepton2Pt"] = ROOT.TH1F (region+"_Lepton2Pt", region+"_Lepton2Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_JetPt"] = ROOT.TH1F (region+"_JetPt", region+"_JetPt ; JetP_{T} (GeV) ; Events ", 2000, 0, 2000)
    h[region+"_MetPt"] = ROOT.TH1F (region+"_MetPt", region+"_MetPt ; MET (GeV) ; Events ", 500, 0, 500)
    h[region+"_Nj"] = ROOT.TH1F (region+"_Nj", region+"_Nj ; N_{j} ; Events ", 10, 0, 10)
    h[region+"_dRl"] = ROOT.TH1F (region+"_dRl", region+"_dRl ; dR(leptons) ; Events", 100, 0, 5)
    h[region+"_dRj"] = ROOT.TH1F (region+"_dRj", region+"_dRj ; dR(jet, ditau) ; Events", 100, 0, 5)
    h[region+"_dPhil"] = ROOT.TH2F (region+"_dPhil", region+"_dPhil ; dPhi(met,lepton1) ; dPhi(met,lepton2)",  100, -pi, pi, 100, -pi, pi)
    h[region+"_dPhi"] = ROOT.TH2F (region+"_dPhi", region+"_dPhi ; dPhi(met,ditau) ; dPhi(met,jet)",  100, -pi, pi, 100, -pi, pi)


def define_general_histogram():

#-----Trigger bits-----

    h['isSingleJet'] = ROOT.TH1F("isSingleJet", "isSingleJet ; isSingleJet ; N", 4,-1.5,2.5)
    h['isHT'] = ROOT.TH1F("isHT", "isHT ; isHT ; N", 4,-1.5,2.5)

#-----Event counts-----

    h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)
    h['hEMu_Events'] = ROOT.TH1F ("hEMu_Events", "hEMu_Events;;N", 6, 1, 7)

#-----Objects----

    h['hJetPt'] = ROOT.TH1F ("hJetPt", "Jet P_{T} ; P_{T} ; N", 1000, 0, 1000)
    h['hBtag'] = ROOT.TH1F ("hBtag", "DeepJet Score; score; N", 50, 0, 1)
    h['hBJetPt'] = ROOT.TH1F ("hBJetPt", "B-tagged Jet P_{T} ; P_{T} ; N", 1000, 0, 1000)
    h['hMuonPt'] = ROOT.TH1F ("hMuPt", "Muon P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hIsoMuonPt'] = ROOT.TH1F ("hIsoMuPt", "Isolated Muon P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hNonIsoMuonPt'] = ROOT.TH1F ("hNonIsoMuPt", "Non-Isolated Muon P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hElectronPt'] = ROOT.TH1F ("hEPt", "Electron P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hIsoElectronPt'] = ROOT.TH1F ("hIsoEPt", "Isolated Electron P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hNonIsoElectronPt'] = ROOT.TH1F ("hNonIsoEPt", "Non-Isolated Electron P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hMetPt'] = ROOT.TH1F ("hMetPt", "MET P_{T} ; P_{T} ; N", 500, 0, 500)

    h['hTauUnCleanedPt'] = ROOT.TH1F ("hTauUnCleanedPt", "Uncleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hMuCleanedPt'] = ROOT.TH1F ("hMuCleanedPt", "Mu-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hMuCleanedPt_altered'] = ROOT.TH1F ("hMuCleanedPt_altered", "Mu-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hECleanedPt'] = ROOT.TH1F ("hECleanedPt", "E-cleaned Tau P_{T}; P_{T}; a.u.", 500, 0, 500)
    h['hTauBoostedPt'] = ROOT.TH1F ("hTauBoostedPt", "Boosted Tau P_{T}; P_{T}; a.u.", 500, 0, 500)


regions = [
    "hEMu_baseline",
    "hEMu_dRcut",
    "hEMu_dRjcut",
]

for r in regions:
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


def check_var_dependence(l1, l2, jets, met_pt, met_phi):

    e, mu, j, m = get_TLorentzVector(l1[0], l2[0], jets[0], met_pt, met_phi)

    trigger = [0,0]

    if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : trigger[0] = 1

    if ( ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) and isMuonEG == 1 ) : trigger[1] = 1
    
    if ( isData == 0 and ( trigger[0] == 1 or trigger[1] == 1 ) ) or \
       ( isData == 1 and ( isSingleMuonDataset == 1 and trigger[0] == 1 ) ) or \
       ( isData == 1 and ( isMuonEGDataset == 1 and ( trigger[0] == 0 and trigger[1] == 1 ) ) ) :

        if pass_baseline(e, mu, j) :            
            plot_event_hist("hEMu_baseline", e, mu, j, m)

            if j.DeltaR(e+mu) > 0.8 :

                plot_event_hist("hEMu_dRjcut", e, mu, j, m)

                for i in range(len(dR_range)-1) :
                    if dR_range[i] < 0.4 and isData == 1 : continue
                    if e.DeltaR(mu) > dR_range[i] and e.DeltaR(mu) < dR_range[i+1] :
                        h["hEMu_Nbjet_MetPt_dR"+str(dR_range[i])].Fill(m.Pt(), len(s_b), genweight)

                        h["hEMu_MetPt_dR"+str(dR_range[i])].Fill(m.Pt(),  genweight)
                        h["hEMu_Nbjet_dR"+str(dR_range[i])].Fill(len(s_b), genweight)

                        if len(s_b) == 0 :
                            h["hEMu_MetPt_bjet0_dR"+str(dR_range[i])].Fill(m.Pt(), genweight)

                        if len(s_b) > 0 :
                            h["hEMu_MetPt_bjet_dR"+str(dR_range[i])].Fill(m.Pt(), genweight)

                        if m.Pt() > 100 :
                            h["hEMu_Nbjet_highMET_dR"+str(dR_range[i])].Fill(len(s_b), genweight)

                        if m.Pt() < 100 :
                            h["hEMu_Nbjet_lowMET_dR"+str(dR_range[i])].Fill(len(s_b), genweight)


            if pass_deltaR(e, mu, j, "EMu") :

                plot_event_hist("hEMu_dRcut", e, mu, j, m)

                for i in range(len(met_range)-1):

                    if m.Pt() > met_range[i] and m.Pt() < met_range[i+1] :
                        h["hEMu_Nbjet_MetPt"+str(met_range[i])].Fill(len(s_b), genweight)
                    
                for i in nb_range:
                    
                    if len(s_b) == i :
                        h["hEMu_MetPt_Nbjet"+str(i)].Fill(m.Pt(), genweight)
        

def plot_event_hist(region, l1, l2, j, m):
    
    h[region+"_Count"].Fill(0, genweight)

    h[region+"_Mass"].Fill((l1+l2).M(), genweight)
    h[region+"_Lepton1Pt"].Fill(l1.Pt(), genweight)
    h[region+"_Lepton2Pt"].Fill(l2.Pt(), genweight)
    h[region+"_JetPt"].Fill(j.Pt(), genweight)
    h[region+"_MetPt"].Fill(m.Pt(), genweight)
    h[region+"_Nj"].Fill(len(s_j), genweight)
    h[region+"_dRl"].Fill(l1.DeltaR(l2), genweight)
    h[region+"_dRj"].Fill(j.DeltaR(l1+l2), genweight)
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

   if isData == 1: genweight = 1

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

   s_nonIsoMu = []
   s_nonIsoE = []

   unclean = []
   eclean = []
   mclean = []
   lowEclean = []

   boosted = []
   mclean_altered = []

   if jets.size()>0:
      for i in range(jets.size()):
         jet = jets.at(i)
         if abs(jet.eta) < 2.5 :
             if jet.id >= 2:
                 h['hJetPt'].Fill(jet.pt, genweight)
                 h['hBtag'].Fill(jet.deepjet, genweight)
                 s_j+=[jet]
                 if jet.deepjet > 0.7476:
                     h['hBJetPt'].Fill(jet.pt, genweight)
                     s_b+=[jet]

   if muons.size()>0:
       for i in range(muons.size()):
           muon = muons.at(i)
           if abs(muon.eta) < 2.4 and muon.pt > 8.0:
               if muon.id >= 1:
                   h['hMuonPt'].Fill(muon.pt, genweight) 
                   s_mu+=[muon]
                   if muon.iso < 0.25:
                       h['hIsoMuonPt'].Fill(muon.pt, genweight)
                       s_isomu+=[muon]
                   if muon.iso > 0.25:
                       h['hNonIsoMuonPt'].Fill(muon.pt, genweight)
                       s_nonIsoMu+=[muon]

   if electrons.size()>0:
      for i in range(electrons.size()):
         electron = electrons.at(i)
         if abs(electron.eta) < 2.5 and electron.pt > 12.0 :
             if electron.id >= 1 :
                 h['hElectronPt'].Fill(electron.pt, genweight)
                 s_e+=[electron]
                 if electron.iso >= 1:
                     h['hIsoElectronPt'].Fill(electron.pt, genweight)
                     s_isoe+=[electron]
                 if electron.iso == 0:
                     h['hNonIsoElectronPt'].Fill(electron.pt, genweight)
                     s_nonIsoE+=[electron]


   s_j.sort(key=lambda x: x.pt, reverse=True)
   s_e.sort(key=lambda x: x.pt, reverse=True)
   s_isoe.sort(key=lambda x: x.pt, reverse=True)
   s_mu.sort(key=lambda x: x.pt, reverse=True)
   s_isomu.sort(key=lambda x: x.pt, reverse=True)

   s_nonIsoMu.sort(key=lambda x: x.pt, reverse=True)
   s_nonIsoE.sort(key=lambda x: x.pt, reverse=True)

   isEMu = 0


   if len(s_isomu) > 0 and len(s_isoe) > 0 and len(s_j) > 0 and s_isomu[0].charge*s_isoe[0].charge < 0 :
       check_var_dependence(s_isoe, s_isomu, s_j, met_pt, met_phi)


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))
