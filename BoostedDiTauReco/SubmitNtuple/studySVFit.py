import ROOT, sys, os
import numpy as np
import time
import correctionlib
from array import array

start_time = time.time()

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

outputTitle = "h_studySVFit"

isData = 0

if "-b" in opts:
    isData = 0

if "-dj" in opts:
    isData = 1
    isJetHTSample = 1
    isSingleElectronSample = 0
    isHTMHTSample = 0

if "-de" in opts:
    isData = 1
    isJetHTSample = 0
    isSingleElectronSample = 1
    isHTMHTSample = 0

if "-dh" in opts:
    isData = 1
    isJetHTSample = 0
    isSingleElectronSample = 0
    isHTMHTSample = 1


ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/JetInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/MuonInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/ElectronInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/TauInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/TrigObjectInfoDS.h"')
ROOT.gInterpreter.ProcessLine('#include "../../../TauAnalysis/ClassicSVfit/test/testClassicSVfit.h"')

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
#    chain3 = ROOT.TChain('lumiSummary/lumiTree')
    chain4 = ROOT.TChain('tcpGenNtuples/genTree')
    chain5 = ROOT.TChain('tcpPrefiring/prefiringTree')
chain6 = ROOT.TChain('tcpMetfilter/metfilterTree')
chain7 = ROOT.TChain('testTrigObj/TriggerObjectTree')

pi = np.pi

h = {}

event_cut = {
    'jetpt': 100,
    'dRl': 0.4,
    'dRltau': 0.05,
    'dRlj': 0.8,
    'metcut': 100.0,
    'mtcut': 50.0,
    'dPhiml': 1,
    'dPhimj': 2,
    'mass' : 5
}


#define histograms here

def book_histogram():

    h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)
    h['hWeights'] = ROOT.TH1F ("hWeights", "Weights per events; weight; N", 100, 0, 2)
    h['hGenWeights'] = ROOT.TH1F ("hGenWeights", "Genweights per events; genweight; N", 100, 0, 2)
    h['hPuWeights'] = ROOT.TH1F ("hPuWeights", "PUweights per events; PUweight; N", 100, 0, 2)
    h['hPrWeights'] = ROOT.TH1F ("hPrWeights", "PreFiringWeights per events; PRweight; N", 100, 0, 2)

    # ---------- Objects ---------- #

    h['hJetPt'] = ROOT.TH1F ("hJetPt", "Jet P_{T} ; P_{T} ; N", 1500, 0, 1500)
    h['hDeepjet'] = ROOT.TH1F ("hDeepjet", "deepjet score ; score ; N", 100, 0, 1)
    h['hBJetPt'] = ROOT.TH1F ("hBJetPt", "BJet P_{T} ; P_{T} ; N", 1500, 0, 1500)

    h['hMuonPt'] = ROOT.TH1F ("hMuPt", "Muon P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hIsoMuonPt'] = ROOT.TH1F ("hIsoMuPt", "Isolated Muon P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hNonIsoMuonPt'] = ROOT.TH1F ("hNonIsoMuPt", "Non-Isolated Muon P_{T} ; P_{T} ; N", 500, 0, 500)

    h['hElectronPt'] = ROOT.TH1F ("hEPt", "Electron P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hIsoElectronPt'] = ROOT.TH1F ("hIsoEPt", "Isolated Electron P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hNonIsoElectronPt'] = ROOT.TH1F ("hNonIsoEPt", "Non-Isolated Electron P_{T} ; P_{T} ; N", 500, 0, 500)

    h['hTauECleanedPt'] = ROOT.TH1F ("hTauECleanedPt", "Electron-Cleaned Tau P_{T} ; P_{T} ; N", 500, 0, 500)
    h['hTauECleanedAlteredPt'] = ROOT.TH1F ("hTauECleanedAlteredPt", "Electron-Cleaned Tau Altered ID P_{T} ; P_{T} ; N", 500, 0, 500)

    h['hTauMuCleanedPt'] = ROOT.TH1F ("hTauMuCleanedPt", "Muon-Cleaned Tau P_{T} ; P_{T} ; N", 500, 0, 500)

    # ---------------

    h['MuTau_OS_dRcut_highMET_lowMt_SVFit'] = ROOT.TH1F ("MuTau_OS_dRcut_highMET_lowMt_SVFit", "MuTau Mass - SVFit ; Mass ; N", 150, 0, 150)

def book_event_histogram(region):

    h[region+"_Count"] = ROOT.TH1F (region+"_Count", region+"_Count ; Events ; Events ", 1, 0, 1)

    h[region+"_Mass"] = ROOT.TH1F (region+"_Mass", region+"_Mass ; M_{vis.} (GeV) ; Events ", 150, 0, 150)
    h[region+"_Lepton1Pt"] = ROOT.TH1F (region+"_Lepton1Pt", region+"_Lepton1Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_Lepton2Pt"] = ROOT.TH1F (region+"_Lepton2Pt", region+"_Lepton2Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_JetPt"] = ROOT.TH1F (region+"_JetPt", region+"_JetPt ; JetP_{T} (GeV) ; Events ", 2000, 0, 2000)
    h[region+"_MetPt"] = ROOT.TH1F (region+"_MetPt", region+"_MetPt ; MET (GeV) ; Events ", 500, 0, 500)
    h[region+"_Mt"] = ROOT.TH1F (region+"_Mt", region+"_Mt ; M_{T} (GeV) ; Events ", 150, 0, 150)
    h[region+"_Nj"] = ROOT.TH1F (region+"_Nj", region+"_Nj ; N_{j} ; Events ", 10, 0, 10)
    h[region+"_dRl"] = ROOT.TH1F (region+"_dRl", region+"_dRl ; dR(leptons) ; Events", 100, 0, 5)
    h[region+"_dRj"] = ROOT.TH1F (region+"_dRj", region+"_dRj ; dR(jet, ditau) ; Events", 100, 0, 5)
    h[region+"_dPhil"] = ROOT.TH2F (region+"_dPhil", region+"_dPhil ; dPhi(met,lepton1) ; dPhi(met,lepton2)",  100, -pi, pi, 100, -pi, pi)
    h[region+"_dPhi"] = ROOT.TH2F (region+"_dPhi", region+"_dPhi ; dPhi(met,ditau) ; dPhi(met,jet)",  100, -pi, pi, 100, -pi, pi)



def get_TLorentzVector(obj):
    
    v = ROOT.TLorentzVector()
    v.SetPtEtaPhiM(obj.pt, obj.eta, obj.phi, obj.mass)

    return v


def pass_deltaR(l1, l2, j, channel):

    if channel == "MuTau" or channel == "ETau":
        if l1.DeltaR(l2) <= event_cut["dRl"] and j.DeltaR(l1) >= event_cut["dRlj"] and j.DeltaR(l2) >= event_cut["dRlj"] and l1.DeltaR(l2) >= event_cut["dRltau"]:
            return 1
        else:
            return -9999
    if channel == "MuMu" or channel == "EMu" or channel == "EE":
        if l1.DeltaR(l2) <= event_cut["dRl"] and j.DeltaR(l1) >= event_cut["dRlj"] and j.DeltaR(l2) >= event_cut["dRlj"]:
            return 1
        else:
            return -9999

def Mt(lepton, met):

    cos = np.cos(met.DeltaPhi(lepton))
    Mt = np.sqrt(2*lepton.Pt()*met.Pt()*(1-cos))

    return Mt

def plot_variable(region, l1, l2, j, m, sf=1):

    h[region+"_Count"].Fill(0, weight*sf)

    h[region+"_Mass"].Fill((l1+l2).M(), weight*sf)
    h[region+"_Lepton1Pt"].Fill(l1.Pt(), weight*sf)
    h[region+"_Lepton2Pt"].Fill(l2.Pt(), weight*sf)
    h[region+"_JetPt"].Fill(j.Pt(), weight*sf)
    h[region+"_MetPt"].Fill(m.Pt(), weight*sf)
    h[region+"_Mt"].Fill(Mt(l1, m), weight*sf)
    h[region+"_Nj"].Fill(len(s_jet), weight*sf)
    h[region+"_dRl"].Fill(l1.DeltaR(l2), weight*sf)
    h[region+"_dRj"].Fill(j.DeltaR(l1+l2), weight*sf)
    h[region+"_dPhil"].Fill(m.DeltaPhi(l1), m.DeltaPhi(l2), weight*sf)
    h[region+"_dPhi"].Fill(m.DeltaPhi(l1+l2), m.DeltaPhi(j), weight*sf)


def mumu_channel():

    isMuMu = 0

    if s_isomuon[0].charge*s_isomuon[1].charge < 0 :

        isJetHTEvent = 0
        isSingleMuonEvent = 0

        mu1 = get_TLorentzVector(s_isomuon[0])
        mu2 = get_TLorentzVector(s_isomuon[1])
        jet = get_TLorentzVector(s_jet[0])

        if ( jet.Pt() > 510 and ( isHT == 1 or isSingleJet500 == 1 ) ) : isJetHTEvent = 1
        if ( mu1.Pt() > 52 and isMu == 1 ) or ( mu1.Pt() > 27 and isIsoMu == 1 ) : isSingleMuonEvent = 1

        if isData == 0 and ( isJetHTEvent == 1 ) \
           or ( isData == 1 and isJetHTSample == 1 and isJetHTEvent == 1 ) :

            if pass_deltaR(mu1, mu2, jet, 'MuMu') == 1 :

                if met.Pt() > event_cut['metcut'] :

                    isMuMu = 1

    return isMuMu


def mutau_channel():

    isMuTau = 0

    if s_muon[0].charge*s_tauMuclean[0].charge < 0:

        mu = get_TLorentzVector(s_muon[0])
        tau = get_TLorentzVector(s_tauMuclean[0])
        jet = get_TLorentzVector(s_jet[0])

        isMatchedMu = False
        for ito in tOisMu:
            if isMatchedMu == True: break
            trigObject = get_TLorentzVector(ito)
            if trigObject.DeltaR(mu) < 0.1 :
                isMatchedMu = True

        if isMatchedMu == True and mu.Pt() >= 50:
            plot_variable('MuTau_TriggerMatch_OS_0MuPtcut_isMu', mu, tau, jet , met)

            if pass_deltaR(mu, tau, jet, 'MuTau') == 1 :
                plot_variable('MuTau_TriggerMatch_OS_dRcut_0MuPtcut_isMu', mu, tau, jet , met)

                if met.Pt() > event_cut['metcut'] :
                    plot_variable('MuTau_TriggerMatch_OS_dRcut_highMET_0MuPtcut_isMu', mu, tau, jet , met)

                    if Mt(mu,met) < event_cut['mtcut'] : #Final event selection
                        plot_variable('MuTau_TriggerMatch_OS_dRcut_highMET_lowMt_0MuPtcut_isMu', mu, tau, jet , met)
                        isMuTau = 1

                        a=None
                        a = ROOT.svFit(met_x, met_y, mu.Pt(), mu.Eta(), mu.Phi(), mu.M(), tau.Pt(), tau.Eta(), tau.Phi(),  tau.M(), s_tauMuclean[0].decaymode, 1, met_covXX, met_covXY, met_covYY)

                        mSVFit = a.runSVFit()
                        print("Mass from SV fit = ", mSVFit)

                        h['MuTau_OS_dRcut_highMET_lowMt_SVFit'].Fill(mSVFit, weight)

    return isMuTau

def ee_channel():

    isEE = 0

    if s_isoelectron[0].charge*s_isoelectron[1].charge < 0 :
        
        e1 = get_TLorentzVector(s_isoelectron[0])
        e2 = get_TLorentzVector(s_isoelectron[1])
        jet = get_TLorentzVector(s_jet[0])

        isJetHTEvent = 0
        isSingleElectronEvent = 0

#        sf = get_sf(s_isoelectron[0], "ele") if isData == 0 else 1

        if ( jet.Pt() > 510 and ( isHT == 1  or isSingleJet500 == 1 ) ) : isJetHTEvent = 1
        if ( e1.Pt() > 37 and isIsoEle == 1 ) : isSingleElectronEvent = 1

        if ( isData == 0 and ( isJetHTEvent == 1 or isSingleElectronEvent == 1 ) ) \
           or ( isData == 1 and isJetHTSample == 1 and isJetHTEvent == 1 ) \
           or ( isData == 1 and ( isSingleElectronSample == 1 and isJetHTSample == 0 and isSingleElectronEvent == 1 ) ) :
             
            if pass_deltaR(e1, e2, jet, 'EE') == 1 :
                if met.Pt() > event_cut['metcut'] :

                    isEE = 1

    return isEE

def emu_channel():

    isEMu = 0

    if s_isomuon[0].charge*s_isoelectron[0].charge < 0:

        isJetHTEvent = 0
        isSingleMuonEvent = 0
        isMuonEGEvent = 0

        e = get_TLorentzVector(s_isoelectron[0])
        mu = get_TLorentzVector(s_isomuon[0])
        jet = get_TLorentzVector(s_jet[0])

        if ( jet.Pt() > 510 and ( isHT == 1 or isSingleJet500 == 1 ) ) : isJetHTEvent = 1
        
        if ( mu.Pt() > 52 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : isSingleMuonEvent = 1
        if ( ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) and isMuonEG == 1 ) : isMuonEGEvent = 1
       
        if ( isData == 0 and ( isJetHTEvent == 1 and ( isSingleMuonEvent == 1 or isMuonEGEvent == 1 ) ) ) \
           or ( isData == 1 and isJetHTSample == 1 and isJetHTEvent == 1 ):

            if pass_deltaR(e, mu, jet, 'EMu') == 1 :

                if met.Pt() > event_cut['metcut'] :
                    
                    isEMu = 1

    return isEMu



book_histogram()

regions = [
'MuTau_TriggerMatch_OS_0MuPtcut_isMu',
'MuTau_TriggerMatch_OS_dRcut_0MuPtcut_isMu',
'MuTau_TriggerMatch_OS_dRcut_highMET_0MuPtcut_isMu',
'MuTau_TriggerMatch_OS_dRcut_highMET_lowMt_0MuPtcut_isMu',
]

for r in regions:
    book_event_histogram(r)

for key in h.keys():
    h[key].Sumw2()

#-------- File loop --------#


inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    inputFileName=inputFileName.replace("\n","")
    print(inputFileName.replace("\n",""))

    fchain.Add(inputFileName)
    chain2.Add(inputFileName)
    if isData == 0:
#        chain3.Add(inputFileName)
        chain4.Add(inputFileName)
        chain5.Add(inputFileName)

    chain6.Add(inputFileName)
    chain7.Add(inputFileName)

#------- Adding friends to the main chain -------#


fchain.AddFriend(chain2)
if isData == 0:
#    fchain.AddFriend(chain3)
    fchain.AddFriend(chain4)
    fchain.AddFriend(chain5)
    
fchain.AddFriend(chain6)
fchain.AddFriend(chain7)

jets = ROOT.JetInfoDS()
muons = ROOT.MuonInfoDS()
electrons = ROOT.ElectronInfoDS()
tausUnCleaned = ROOT.TauInfoDS()
tausECleaned = ROOT.TauInfoDS()
tausMCleaned = ROOT.TauInfoDS()
tausBoosted = ROOT.TauInfoDS()
trigObj = ROOT.TrigObjectInfoDS()

fchain.SetBranchAddress("Jets", ROOT.AddressOf(jets))
fchain.SetBranchAddress("Muons", ROOT.AddressOf(muons))
fchain.SetBranchAddress("Electrons", ROOT.AddressOf(electrons))
fchain.SetBranchAddress("TausUnCleaned", ROOT.AddressOf(tausUnCleaned))
fchain.SetBranchAddress("TausECleaned", ROOT.AddressOf(tausECleaned))
fchain.SetBranchAddress("TausMCleaned", ROOT.AddressOf(tausMCleaned))
fchain.SetBranchAddress("TausBoosted", ROOT.AddressOf(tausBoosted))
fchain.SetBranchAddress("TriggerObjects", ROOT.AddressOf(trigObj))

if isData == 0:
    genParticle = ROOT.GenParticleInfoDS()
    fchain.SetBranchAddress("GenParticleInfo", ROOT.AddressOf(genParticle))

#----------- Event loop ----------#


for iev in range(fchain.GetEntries()): # Be careful!!!                                                               

    fchain.GetEntry(iev)

    mets = fchain.GetBranch("Mets")
    met_pt = mets.GetLeaf('pt').GetValue()
    met_phi = mets.GetLeaf('phi').GetValue()
    
    # ------- For SV fit

    met_covXX = mets.GetLeaf('covXX').GetValue()
    met_covXY = mets.GetLeaf('covXY').GetValue()
    met_covYY = mets.GetLeaf('covYY').GetValue()

    met_x = met_pt*np.cos(met_phi)
    met_y = met_pt*np.sin(met_phi)

    # ------- 

    met = ROOT.TLorentzVector()
    met.SetPtEtaPhiM(met_pt, 0, met_phi, 0)

    if isData == 0 :
        genweight = fchain.GetLeaf('genWeight').GetValue()
        puweight = fchain.GetLeaf('puWeight').GetValue()
        prweight = fchain.GetLeaf('prefiringWeight').GetValue()
    else :
        genweight = 1
        puweight = 1
        prweight = 1

    weight = genweight*puweight*prweight

    h['hEvents'].Fill(0.5, 1)
    h['hEvents'].Fill(1.5, weight)

    h['hWeights'].Fill(weight)
    h['hPuWeights'].Fill(puweight)
    h['hGenWeights'].Fill(genweight)
    h['hPrWeights'].Fill(prweight)

    #-------------- Gen particles -----------#

    if isData == 0 :
        gen_e = []
        if genParticle.size() > 0 :
            for i in range(genParticle.size()):
                igen = genParticle.at(i)
                if igen.isdirecthardprocesstaudecayproductfinalstate :
                    if abs(igen.pdgid) == 11 : gen_e+=[igen]

    #-------------- Trigger Objects ---------------#

    tOisEleJet = []
    tOisIsoMu = []
    tOisMu = []
    tOisSingleJet = []

    if trigObj.size() > 0:
        for i in range(trigObj.size()):
            iobj = trigObj.at(i)
            if iobj.isEleJet == 1 : 
                tOisEleJet+=[iobj]
                # h['hIsEleJetPt'].Fill(iobj.pt)
                # h['hIsEleJetEta'].Fill(iobj.eta)
            if iobj.isMu == 1 : tOisMu+=[iobj]
            if iobj.isIsoMu == 1 : tOisIsoMu +=[iobj]
            if iobj.isSingleJet == 1 : tOisSingleJet+=[iobj]
                
    tOisEleJet.sort(key=lambda x: x.pt, reverse=True)
    tOisIsoMu.sort(key=lambda x: x.pt, reverse=True)    
    tOisMu.sort(key=lambda x: x.pt, reverse=True)

    #-------------- Trigger definitions -----------#

    isSingleJet500 = fchain.GetLeaf('isSingleJet500').GetValue()
    isHT = fchain.GetLeaf('isHT').GetValue()
    isHTMHT = fchain.GetLeaf('isHTMHT').GetValue()
    isMu = fchain.GetLeaf('isMu').GetValue()
    isIsoMu = fchain.GetLeaf('isIsoMu').GetValue()
    isIsoMuTau = fchain.GetLeaf('isIsoMuTau').GetValue()
    isIsoEle = fchain.GetLeaf('isIsoEle').GetValue()
    isEleTau = fchain.GetLeaf('isEleTau').GetValue()
    isMuonEG = fchain.GetLeaf('isMuonEG').GetValue()
    isEle = fchain.GetLeaf('isEle').GetValue()
    isPhoton200 = fchain.GetLeaf('isPhoton200').GetValue()
    isEleJet = fchain.GetLeaf('isEleJet').GetValue()
    isSingleJet500 = fchain.GetLeaf('isSingleJet500').GetValue()
    isSingleJet450 = fchain.GetLeaf('isSingleJet450').GetValue()

    #------------ MET Filter --------------#

    pvFilter = fchain.GetLeaf('primaryvertexfilter').GetValue()
    haloFilter = fchain.GetLeaf('beamhalofilter').GetValue()
    hbheFilter = fchain.GetLeaf('hbhefilter').GetValue()
    hbheIsoFilter = fchain.GetLeaf('hbheisofilter').GetValue()
    eeBadSCFilter = fchain.GetLeaf('eebadscfilter').GetValue()
    ecalTPFilter = fchain.GetLeaf('ecaltpfilter').GetValue()
    badMuonFilter = fchain.GetLeaf('badpfmuonfilter').GetValue()
    ecalBadCalFilter = fchain.GetLeaf('ecalbadcalfilter').GetValue()

    #------------ Objects loop ------------#

    s_jet = []
    s_bjet = []

    iht = 0

    if jets.size() > 0:
        for i in range(jets.size()):
            ijet = jets.at(i)
            if abs(ijet.eta) < 2.5 :
                if ijet.id >= 2:
                    # h['hJetPt'].Fill(ijet.pt, weight)
                    # h['hDeepjet'].Fill(ijet.deepjet, weight)
                    s_jet+=[ijet]
                    iht = iht + ijet.pt
                    if ijet.deepjet >= 0.7476:
                        # h['hBJetPt'].Fill(ijet.pt, weight) 
                        s_bjet+=[ijet]

    s_muon = []
    s_isomuon = []
    s_nonisomuon = []

    if muons.size() > 0:
        for i in range(muons.size()):
            imuon = muons.at(i)
            if abs(imuon.eta) < 2.4 :
                if imuon.id >= 1: #loose Muons
                    # h['hMuonPt'].Fill(imuon.pt, weight) 
                    s_muon+=[imuon]
                    if imuon.iso <= 0.25:
                        # h['hIsoMuonPt'].Fill(imuon.pt, weight)
                        s_isomuon+=[imuon]
                    if imuon.iso > 0.25:
                        # h['hNonIsoMuonPt'].Fill(imuon.pt, weight)
                        s_nonisomuon+=[imuon]

    s_electron = []
    s_isoelectron = []
    s_nonisoelectron = []

    # if electrons.size() > 0:
    #     for i in range(electrons.size()):
    #         ielectron = electrons.at(i)
    #         if abs(ielectron.eta) < 2.5 :
    #             if ielectron.id >= 1 :
    #                 h['hElectronPt'].Fill(ielectron.pt, weight)
    #                 s_electron+=[ielectron]
    #                 if ielectron.iso >= 1:
    #                     h['hIsoElectronPt'].Fill(ielectron.pt, weight)
    #                     s_isoelectron+=[ielectron]
    #                 if ielectron.iso == 0:
    #                     h['hNonIsoElectronPt'].Fill(ielectron.pt, weight)
    #                     s_nonisoelectron+=[ielectron]

    s_tauEclean = []
    s_tauEcleanAltered = []

    # if tausECleaned.size()>0:
    #     for i in range(tausECleaned.size()):
    #         itau = tausECleaned.at(i)
    #         if abs(itau.eta) < 2.3 :
    #             if itau.mvaid >= 4:
    #                 h['hTauECleanedPt'].Fill(itau.pt, weight)
    #                 s_tauEclean+=[itau]
    #             if itau.mvaid < 4 and itau.mvaid >= 1 :
    #                 h['hTauECleanedAlteredPt'].Fill(itau.pt, weight)
    #                 s_tauEcleanAltered+=[itau]
                    

    s_tauMuclean = []

    if tausMCleaned.size()>0:
        for i in range(tausMCleaned.size()):
            itau = tausMCleaned.at(i)
            if abs(itau.eta) < 2.3 :
                if itau.mvaid >= 4 :
                    # h['hTauMuCleanedPt'].Fill(itau.pt, weight)
                    s_tauMuclean+=[itau]


    # ---------- Event Selections --------- #

    # if len(s_isomuon) >= 2 and len(s_jet) >= 1 : 
    #     if mumu_channel() == 1: continue

    # if len(s_isomuon) >= 1 and len(s_isoelectron) >= 1 and len(s_jet) >= 1 : 
    #     if emu_channel() == 1 : continue

    if len(s_muon) >= 1 and len(s_tauMuclean) >= 1 and len(s_jet) >= 1 : 
        mutau_channel()
        # mutau_trigger_study()
        # if mutau_channel() == 1 : continue


out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))

# result = np.array(regions)
# print(result)
# print(regions)
