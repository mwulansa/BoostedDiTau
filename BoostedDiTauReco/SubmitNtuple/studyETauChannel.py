import ROOT, sys, os
import numpy as np
import time
import correctionlib
from array import array

start_time = time.time()

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

sfDir = os.path.join('/cvmfs','cms.cern.ch','rsync','cms-nanoAOD','jsonpog-integration','POG')

btvjson = correctionlib.CorrectionSet.from_file(os.path.join(sfDir, 'BTV', '2017_UL', 'btagging.json.gz'))
elejson = correctionlib.CorrectionSet.from_file(os.path.join(sfDir, 'EGM', '2017_UL', 'electron.json.gz'))

outputTitle = "h_studyETauChannel"

isData = 0

if "-b" in opts:
    isData = 0
    isJetHTSample = 0
    isSingleElectronSample = 0
    isHTMHTSample = 0

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
    chain5 = ROOT.TChain('tcpPrefiring/prefiringTree')
chain6 = ROOT.TChain('tcpMetfilter/metfilterTree')

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

    h["ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet_JetPt_EPt"] = ROOT.TH2F ("ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet_JetPt_EPt", "JetPt vs. EPt ; JetPt (GeV) ; EPt (GeV)", 2000, 0,2000, 500, 0, 500)
    h["ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isSingleJet_JetPt_EPt"] = ROOT.TH2F ("ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isSingleJet_JetPt_EPt", "JetPt vs. EPt ; JetPt (GeV) ; EPt (GeV)", 2000, 0,2000, 500, 0, 500)
    h["ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet_isSingleJet_JetPt_EPt"] = ROOT.TH2F ("ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet_isSingleJet_JetPt_EPt", "JetPt vs. EPt ; JetPt (GeV) ; EPt (GeV)", 2000, 0,2000, 500, 0, 500)
    h["ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_Before_JetPt_EPt"] = ROOT.TH2F ("ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_Before_JetPt_EPt", "JetPt vs. EPt ; JetPt (GeV) ; EPt (GeV)", 2000, 0,2000, 500, 0, 500)

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


def book_trigger_histogram(region):

    h[region+"_Count"] = ROOT.TH1F (region+"_Count", region+"_Count ; Events ; Events ", 1, 0, 1)
    h[region+"_Lepton1Pt"] = ROOT.TH1F (region+"_Lepton1Pt", region+"_Lepton1Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_Lepton2Pt"] = ROOT.TH1F (region+"_Lepton2Pt", region+"_Lepton2Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_LeadingJetPt"] = ROOT.TH1F (region+"_LeadingJetPt", region+"_LeadingJetPt ; JetP_{T} (GeV) ; Events ", 2000, 0, 2000)
    h[region+"_HT"] = ROOT.TH1F (region+"_HT", region+"_HT ; HT (GeV) ; Events ", 4000, 0, 4000)
    h[region+"_Mass"] = ROOT.TH1F (region+"_Mass", region+"_Mass ; M_{vis.} (GeV) ; Events ", 150, 0, 150)
    h[region+"_dRl"] = ROOT.TH1F (region+"_dRl", region+"_dRl ; dR(leptons) ; Events", 100, 0, 5)
    h[region+"_dRj"] = ROOT.TH1F (region+"_dRj", region+"_dRj ; dR(jet, ditau) ; Events", 100, 0, 5)
    h[region+"_MetPt"] = ROOT.TH1F (region+"_MetPt", region+"_MetPt ; MET (GeV) ; Events ", 500, 0, 500)


def get_TLorentzVector(obj):
    
    v = ROOT.TLorentzVector()
    v.SetPtEtaPhiM(obj.pt, obj.eta, obj.phi, obj.mass)

    return v


def pass_deltaR(l1, l2, j, channel):

    if channel == "MuTau" or channel == "ETau":
        if l1.DeltaR(l2) < event_cut["dRl"] and j.DeltaR(l1) > event_cut["dRlj"] and j.DeltaR(l2) > event_cut["dRlj"] and l1.DeltaR(l2) > event_cut["dRltau"]:
            return 1
        else:
            return -9999
    if channel == "MuMu" or channel == "EMu" or channel == "EE":
        if l1.DeltaR(l2) < event_cut["dRl"] and j.DeltaR(l1) > event_cut["dRlj"] and j.DeltaR(l2) > event_cut["dRlj"]:
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



def plot_for_triggers(region, l1, l2, j, ht):

    global weight
    if weight<0 : weight = 1

    h[region+"_Count"].Fill(0, weight)

    h[region+"_Lepton1Pt"].Fill(l1.Pt(), weight)
    h[region+"_Lepton2Pt"].Fill(l2.Pt(), weight)
    h[region+"_LeadingJetPt"].Fill(j.Pt(), weight)
    h[region+"_HT"].Fill(ht, weight)
    h[region+"_Mass"].Fill((l1+l2).M(), weight)
    h[region+"_dRl"].Fill(l1.DeltaR(l2), weight)
    h[region+"_dRj"].Fill(j.DeltaR(l1+l2), weight)
    h[region+"_MetPt"].Fill(met.Pt(), weight)

def get_sf(objt, l):

    if l == "ele":

        try:
            valsf = elejson["UL-Electron-ID-SF"].evaluate("2017","sf","Loose",abs(objt.eta), objt.pt)
        except:
            valsf = 1

    elif l == "jet":

        if js.jetflavour != 0 :

            deepjet_col = 'deepJet_comb'

        if js.jetflavour == 0 :

            deepjet_col = 'deepJet_incl'

        try:
            valsf = btvjson[deepjet_col].evaluate("central", "tight", objt.jetflavour, abs(objt.eta), objt.pt)
        except:
            valsf = 1

    return valsf

def etau_trigger_study():

    e = get_TLorentzVector(s_electron[0])
    tau = get_TLorentzVector(s_tauEclean[0])
    jet = get_TLorentzVector(s_jet[0])

#    isEleJet = 0
#    isEle = 0
#    isPhoton200 = 0

    if isData == 0 :
        if len(gen_e) > 0 :
            genE = get_TLorentzVector(gen_e[0])
            if len(gen_e) > 1 :
                genE2 = get_TLorentzVector(gen_e[1])

    if s_electron[0].charge*s_tauEclean[0].charge < 0 :
        plot_for_triggers('ETau_TriggerStudy_OS_noGenMatch_Before', e, tau, jet, iht)

        if isHT == 1 and iht >= 1200 :
            plot_for_triggers('ETau_TriggerStudy_OS_noGenMatch_JetHT', e, tau, jet, iht)
        if isSingleJet500 == 1 and jet.Pt() >= 600 :
            plot_for_triggers('ETau_TriggerStudy_OS_noGenMatch_SingleJet', e, tau, jet, iht)
        if isSingleJet450 == 1 :
            plot_for_triggers('ETau_TriggerStudy_OS_noGenMatch_SingleJet450', e, tau, jet, iht)
        if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
            plot_for_triggers('ETau_TriggerStudy_OS_noGenMatch_SingleE', e, tau, jet, iht)
        if isEle == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
            plot_for_triggers('ETau_TriggerStudy_OS_noGenMatch_isEle', e, tau, jet, iht)
        if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
            plot_for_triggers('ETau_TriggerStudy_OS_noGenMatch_isEleJet', e, tau, jet, iht)
        if isHTMHT == 1 :
            plot_for_triggers('ETau_TriggerStudy_OS_noGenMatch_isHTMHT', e, tau, jet, iht)
            
        if isData == 0:

            if len(gen_e) == 1 and e.DeltaR(genE) < 0.1:

                plot_for_triggers('ETau_TriggerStudy_OS_yesGenMatch_Before', e, tau, jet, iht)

                if isHT == 1 and iht >= 1200 :
                    plot_for_triggers('ETau_TriggerStudy_OS_yesGenMatch_JetHT', e, tau, jet, iht)
                if isSingleJet500 == 1 and jet.Pt() >= 600 :
                    plot_for_triggers('ETau_TriggerStudy_OS_yesGenMatch_SingleJet', e, tau, jet, iht)
                if isSingleJet450 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_yesGenMatch_SingleJet450', e, tau, jet, iht)
                if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
                    plot_for_triggers('ETau_TriggerStudy_OS_yesGenMatch_SingleE', e, tau, jet, iht)
                if isEle == 1 or isPhoton200 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_yesGenMatch_isEle', e, tau, jet, iht)
                if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_yesGenMatch_isEleJet', e, tau, jet, iht)
                if isHTMHT == 1:
                    plot_for_triggers('ETau_TriggerStudy_OS_yesGenMatch_isHTMHT', e, tau, jet, iht)


        if e.DeltaR(tau) < event_cut["dRl"] and e.DeltaR(tau) > event_cut["dRltau"] :

            plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_noGenMatch_Before', e, tau, jet, iht)
            
            if isHT == 1 and iht >= 1200 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_noGenMatch_JetHT', e, tau, jet, iht)
            if isSingleJet500 == 1 and jet.Pt() >= 600 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_noGenMatch_SingleJet', e, tau, jet, iht)
            if isSingleJet450 == 1 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_noGenMatch_SingleJet450', e, tau, jet, iht)
            if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
                plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_noGenMatch_SingleE', e, tau, jet, iht)
            if isEle == 1 or isPhoton200 == 1 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_noGenMatch_isEle', e, tau, jet, iht)
            if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_noGenMatch_isEleJet', e, tau, jet, iht)
            if isHTMHT == 1:
                plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_noGenMatch_isHTMHT', e, tau, jet, iht)


            if isData == 0 and len(gen_e) == 1 and e.DeltaR(genE) < 0.1:

                plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_Before', e, tau, jet, iht)

                if isHT == 1 and iht >= 1200 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_JetHT', e, tau, jet, iht)
                if isSingleJet500 == 1 and jet.Pt() >= 600 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_SingleJet', e, tau, jet, iht)
                if isSingleJet450 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_SingleJet450', e, tau, jet, iht)
                if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
                    plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_SingleE', e, tau, jet, iht)
                if isEle == 1 or isPhoton200 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_isEle', e, tau, jet, iht)
                if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_isEleJet', e, tau, jet, iht)
                if isHTMHT == 1:
                    plot_for_triggers('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_isHTMHT', e, tau, jet, iht)


        if jet.DeltaR(tau) > event_cut["dRlj"] and jet.DeltaR(e) > event_cut["dRlj"] :

            plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_noGenMatch_Before', e, tau, jet, iht)

            if isHT == 1 and iht >= 1200 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_noGenMatch_JetHT', e, tau, jet, iht)
            if isSingleJet500 == 1 and jet.Pt() >= 600 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_noGenMatch_SingleJet', e, tau, jet, iht)
            if isSingleJet450 == 1 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_noGenMatch_SingleJet450', e, tau, jet, iht)
            if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
                plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_noGenMatch_SingleE', e, tau, jet, iht)
            if isEle == 1 or isPhoton200 == 1 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_noGenMatch_isEle', e, tau, jet, iht)
            if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_noGenMatch_isEleJet', e, tau, jet, iht)
            if isHTMHT == 1:
                plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_noGenMatch_isHTMHT', e, tau, jet, iht)


            if isData == 0 and len(gen_e) == 1 and e.DeltaR(genE) < 0.1:

                plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_Before', e, tau, jet, iht)

                if isHT == 1 and iht >= 1200 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_JetHT', e, tau, jet, iht)
                if isSingleJet500 == 1 and jet.Pt() >= 600 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_SingleJet', e, tau, jet, iht)
                if isSingleJet450 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_SingleJet450', e, tau, jet, iht)
                if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
                    plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_SingleE', e, tau, jet, iht)
                if isEle == 1 or isPhoton200 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_isEle', e, tau, jet, iht)
                if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_isEleJet', e, tau, jet, iht)
                if isHTMHT == 1:
                    plot_for_triggers('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_isHTMHT', e, tau, jet, iht)



        if pass_deltaR(e, tau, jet, 'ETau') == 1 :
            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_noGenMatch_Before', e, tau, jet, iht)

            if isHT == 1 and iht >= 1200 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_noGenMatch_JetHT', e, tau, jet, iht)
            if isSingleJet500 == 1 and jet.Pt() >= 600 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_noGenMatch_SingleJet', e, tau, jet, iht)
            if isSingleJet450 == 1 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_noGenMatch_SingleJet450', e, tau, jet, iht)
            if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_noGenMatch_SingleE', e, tau, jet, iht)
            if isEle == 1 or isPhoton200 == 1 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_noGenMatch_isEle', e, tau, jet, iht)
            if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_noGenMatch_isEleJet', e, tau, jet, iht)
            if isHTMHT == 1:
                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_noGenMatch_isHTMHT', e, tau, jet, iht)

            if isData == 0:

                if len(gen_e) == 1 and e.DeltaR(genE) < 0.1:

                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_yesGenMatch_Before', e, tau, jet, iht)

                    if isHT == 1 and iht >= 1200 :
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_yesGenMatch_JetHT', e, tau, jet, iht)
                    if isSingleJet500 == 1 and jet.Pt() >= 600 :
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_yesGenMatch_SingleJet', e, tau, jet, iht)
                    if isSingleJet450 == 1 :
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_yesGenMatch_SingleJet450', e, tau, jet, iht)
                    if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_yesGenMatch_SingleE', e, tau, jet, iht)
                    if isEle == 1 or isPhoton200 == 1 :
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_yesGenMatch_isEle', e, tau, jet, iht)
                    if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_yesGenMatch_isEleJet', e, tau, jet, iht)
                    if isHTMHT == 1:
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_yesGenMatch_isHTMHT', e, tau, jet, iht)

            if met.Pt() > 100 :
                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_Before', e, tau, jet, iht)

                if isHT == 1 and iht >= 1200 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_JetHT', e, tau, jet, iht)
                if isSingleJet500 == 1 and jet.Pt() >= 600 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_SingleJet', e, tau, jet, iht)
                if isSingleJet450 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_SingleJet450', e, tau, jet, iht)
                if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_SingleE', e, tau, jet, iht)
                if isEle == 1 or isPhoton200 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isEle', e, tau, jet, iht)
                if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isEleJet', e, tau, jet, iht)
                if isHTMHT == 1:
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isHTMHT', e, tau, jet, iht)

                if isIsoEle == 1 or ( ( isHT == 1 or isSingleJet500 == 1 ) and jet.Pt() > 600 ) :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_SingleE_pTcut', e, tau, jet, iht)
                if ( ( isEle == 1 or isPhoton200 == 1 ) and e.Pt() > 120 ) or ( isHT == 1 and iht > 1200 ) or ( isSingleJet500 == 1 and jet.Pt() > 600 ) :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isEle_pTcut', e, tau, jet, iht)
                if ( ( isEleJet == 1 or isPhoton200 == 1 ) and e.Pt() > 100 ) or isHT == 1 or isSingleJet500 == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isEleJet_EpTcut', e, tau, jet, iht)
                if ( ( isEleJet == 1  or  isPhoton200 == 1 ) and jet.Pt() > 200 ) or isSingleJet500 == 1 or isHT == 1 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isEleJet_JetpTcut', e, tau, jet, iht)

                if Mt(e, met) < 50 :
                    plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_Before', e, tau, jet, iht)
                    if isHT == 1 and iht >= 1200 :
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_JetHT', e, tau, jet, iht)
                    if isSingleJet500 == 1 and jet.Pt() >= 600 :
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_SingleJet', e, tau, jet, iht)
                    if isSingleJet450 == 1 :
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_SingleJet450', e, tau, jet, iht)
                    if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_SingleE', e, tau, jet, iht)
                    if isEle == 1 or isPhoton200 == 1 :
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_isEle', e, tau, jet, iht)
                    if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_isEleJet', e, tau, jet, iht)
                    if isHTMHT == 1:
                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_isHTMHT', e, tau, jet, iht)

                if isData == 0:

                    if len(gen_e) == 1 and e.DeltaR(genE) < 0.1:

                        plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_Before', e, tau, jet, iht)
                        if isHT == 1 and iht >= 1200 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_JetHT', e, tau, jet, iht)
                        if isSingleJet500 == 1 and jet.Pt() >= 600 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_SingleJet', e, tau, jet, iht)
                        if isSingleJet450 == 1 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_SingleJet450', e, tau, jet, iht)
                        if isIsoEle == 1 or isHT == 1 or isSingleJet500 == 1:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_SingleE', e, tau, jet, iht)
                        if isEle == 1 or isPhoton200 == 1 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEle', e, tau, jet, iht)
                        if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet', e, tau, jet, iht)
                        if isHTMHT == 1:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isHTMHT', e, tau, jet, iht)

                        if isIsoEle == 1 or ( ( isHT == 1 or isSingleJet500 == 1 ) and jet.Pt() > 600 ) :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_SingleE_pTcut', e, tau, jet, iht)
                        if ( ( isEle == 1 or isPhoton200 == 1 ) and e.Pt() > 120 ) or ( isHT == 1 and iht > 1200 ) or ( isSingleJet500 == 1 and jet.Pt() > 600 ) :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEle_pTcut', e, tau, jet, iht)

                        h["ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_Before_JetPt_EPt"].Fill(jet.Pt(), e.Pt(), weight)
                        if isEleJet == 1  or  isPhoton200 == 1 or isSingleJet500 == 1 or isHT == 1 :
                            h["ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet_isSingleJet_JetPt_EPt"].Fill(jet.Pt(), e.Pt(), weight)
                        if isEleJet == 1  or  isPhoton200 == 1 :
                            h["ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet_JetPt_EPt"].Fill(jet.Pt(), e.Pt(), weight)
                        if isSingleJet500 == 1 or isHT == 1 :
                            h["ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isSingleJet_JetPt_EPt"].Fill(jet.Pt(), e.Pt(), weight)

                        if e.Pt() > 100 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_EpTcut_Before', e, tau, jet, iht)
                        if jet.Pt() > 200:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_JetpTcut_Before', e, tau, jet, iht)
                        if ( isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 100:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet_EpTcut', e, tau, jet, iht)
                        if ( isEleJet == 1  or  isPhoton200 == 1 or isSingleJet500 == 1 or isHT == 1 ) and jet.Pt() > 200:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet_JetpTcut', e, tau, jet, iht)

                        if e.Pt() > 50 and jet.Pt() > 200 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_50EpTcut_Before', e, tau, jet, iht)
                        if ( isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 50 and jet.Pt() > 200:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_50EpTcut_isEleJet', e, tau, jet, iht)
                        if ( isIsoEle == 1 or isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 50 and jet.Pt() > 200:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_50EpTcut_isEleJet_isIsoEle',e, tau, jet, iht)

                        if e.Pt() > 60 and jet.Pt() > 180 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_180JetPtcut_60EpTcut_Before', e, tau, jet, iht)
                        if ( isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 60 and jet.Pt() > 180:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_180JetPtcut_60EpTcut_isEleJet', e, tau, jet, iht)
                        if ( isIsoEle == 1 or isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 60 and jet.Pt() > 180:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_180JetPtcut_60EpTcut_isEleJet_isIsoEle', e, tau, jet, iht)

                        if e.Pt() > 70 and jet.Pt() > 170 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_170JetPtcut_70EpTcut_Before', e, tau, jet, iht)
                        if ( isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 70 and jet.Pt() > 170:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_170JetPtcut_70EpTcut_isEleJet', e, tau, jet, iht)
                        if ( isIsoEle == 1 or isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 70 and jet.Pt() > 170:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_170JetPtcut_70EpTcut_isEleJet_isIsoEle', e, tau, jet, iht)

                        if e.Pt() > 80 and jet.Pt() > 165 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_165JetPtcut_80EpTcut_Before', e, tau, jet, iht)
                        if ( isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 80 and jet.Pt() > 165:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_165JetPtcut_80EpTcut_isEleJet', e, tau, jet, iht)
                        if ( isIsoEle == 1 or isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 80 and jet.Pt() > 165:
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_165JetPtcut_80EpTcut_isEleJet_isIsoEle', e, tau, jet, iht)

                        if e.Pt() > 90 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_90EpTcut_Before', e, tau, jet, iht)
                        if ( isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 90 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_90EpTcut_isEleJet', e, tau, jet, iht)
                        if ( isIsoEle == 1 or isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 90 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_90EpTcut_isEleJet_isIsoEle', e, tau, jet, iht)

                        if e.Pt() > 95 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_95EpTcut_Before', e, tau, jet, iht)
                        if ( isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 95 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_95EpTcut_isEleJet', e, tau, jet, iht)
                        if ( isIsoEle == 1 or isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 95 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_95EpTcut_isEleJet_isIsoEle', e, tau, jet, iht)

                        if e.Pt() > 100 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_100EpTcut_Before', e, tau, jet, iht)
                        if ( isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 100 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_100EpTcut_isEleJet', e, tau, jet, iht)
                        if ( isIsoEle == 1 or isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 100 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_100EpTcut_isEleJet_isIsoEle', e, tau, jet, iht)

                        if jet.Pt() > 200 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_0EpTcut_Before', e, tau, jet, iht)
                        if ( isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and jet.Pt() > 200 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_0EpTcut_isEleJet', e, tau, jet, iht)
                        if ( isIsoEle == 1 or isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and jet.Pt() > 200 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_0EpTcut_isEleJet_isIsoEle', e, tau, jet, iht)

                        if e.Pt() > 50 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_50EpTcut_Before', e, tau, jet, iht)
                        if ( isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 50 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_50EpTcut_isEleJet', e, tau, jet, iht)
                        if ( isIsoEle == 1 or isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 ) and e.Pt() > 50 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_50EpTcut_isEleJet_isIsoEle', e, tau, jet, iht)


                        if Mt(e, met) < 50 :
                            plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_Before', e, tau, jet, iht)
                            if isHT == 1 and iht >= 1200 :
                                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_JetHT', e, tau, jet, iht)
                            if isSingleJet500 == 1 and jet.Pt() >= 600 :
                                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_SingleJet', e, tau, jet, iht)
                            if isSingleJet450 == 1 :
                                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_SingleJet450', e, tau, jet, iht)
                            if isIsoEle == 1 :
                                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_SingleE', e, tau, jet, iht)
                            if isEle == 1 or isPhoton200 == 1 :
                                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_isEle', e, tau, jet, iht)
                            if isEleJet == 1 or isPhoton200 == 1 or isHT == 1 or isSingleJet500 == 1 :
                                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_isEleJet', e, tau, jet, iht)
                            if isHTMHT == 1:
                                plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_isHTMHT', e, tau, jet, iht)


                # if len(gen_e) > 1 and ( e.DeltaR(genE) < 0.1 or e.DeltaR(genE2) < 0.1 ) :

                #     plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_Before', e, tau, jet, iht)
                #     if isHT == 1 and iht >= 1200 :
                #         plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_JetHT', e, tau, jet, iht)
                #     if isSingleJet500 == 1 and jet.Pt() >= 600 :
                #         plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_SingleJet', e, tau, jet, iht)
                #     if isIsoEle == 1 :
                #         plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_SingleE', e, tau, jet, iht)

                #     if Mt(e, met) < 50 :
                #         plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_Before', e, tau, jet, iht)
                #         if isHT == 1 and iht >= 1200 :
                #             plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_JetHT', e, tau, jet, iht)
                #         if isSingleJet500 == 1 and jet.Pt() >= 600 :
                #             plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_SingleJet', e, tau, jet, iht)
                #         if isIsoEle == 1 :
                #             plot_for_triggers('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_SingleE', e, tau, jet, iht)



def etau_channel(eCleanedTau, tauid):

    if tauid == "nominal" : selectedTau = "ETau"
    if tauid == "altered" : selectedTau = "ETau_Altered"

    isJetHTEvent = 0
    isSingleElectronEvent = 0
    isHTMHTEvent = 0

    e = get_TLorentzVector(s_electron[0])
    tau = get_TLorentzVector(eCleanedTau[0])
    jet = get_TLorentzVector(s_jet[0])

    if ( jet.Pt() > 510 and ( isHT == 1 or isSingleJet500 == 1 ) ) : isJetHTEvent = 1
    if ( e.Pt() > 36 and isIsoEle == 1 ) : isSingleElectronEvent = 1
    if ( met.Pt() > 500 and isHTMHT == 1 ) : isHTMHTEvent = 1

    if (e+tau).M() > event_cut['mass'] and jet.Pt() > event_cut['jetpt']:

        if ( isData == 0 and ( isJetHTEvent == 1 or isSingleElectronEvent == 1 ) ) \
           or ( isData == 1 and ( isJetHTSample == 1 and isJetHTEvent == 1 ) ) \
           or ( isData == 1 and ( isSingleElectronSample == 1 and ( isJetHTEvent == 0 and isSingleElectronEvent == 1 ) ) ) :

            if s_electron[0].charge*eCleanedTau[0].charge < 0 : #OS
                plot_variable(selectedTau+'_OS', e, tau, jet, met)

                if e.DeltaR(tau) > 0.8 : #OS resolved
                    plot_variable(selectedTau+'_OS_resolved', e, tau, jet, met)

                if met.Pt() > event_cut['metcut'] : #OS highMET - nondRcut
                    plot_variable(selectedTau+'_OS_highMET', e, tau, jet, met)

                    if e.DeltaR(tau) > 0.8 : #OS highMET resolved
                        plot_variable(selectedTau+'_OS_highMET_resolved', e, tau, jet, met)

                if pass_deltaR(e, tau, jet, 'ETau') == 1 :
                    plot_variable(selectedTau+'_OS_dRcut', e, tau, jet, met)

                    if met.Pt() > event_cut['metcut'] : #OS highMET
                        plot_variable(selectedTau+'_OS_dRcut_highMET', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #OS highMET lowMt
                            plot_variable(selectedTau+'_OS_dRcut_highMET_lowMt', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #OS highMET highMt                                                               
                            plot_variable(selectedTau+'_OS_dRcut_highMET_highMt', e, tau, jet, met)

                    if met.Pt() < event_cut['metcut'] : #OS lowMET   
                        plot_variable(selectedTau+'_OS_dRcut_lowMET', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #OS lowMET lowMt                                                                                         
                            plot_variable(selectedTau+'_OS_dRcut_lowMET_lowMt', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #OS lowMET highMt          
                            plot_variable(selectedTau+'_OS_dRcut_lowMET_highMt', e, tau, jet, met)

            if s_electron[0].charge*eCleanedTau[0].charge > 0 : #SS
                plot_variable(selectedTau+'_SS', e, tau, jet, met)

                if met.Pt() > event_cut['metcut'] :
                    plot_variable(selectedTau+'_SS_highMET', e, tau, jet , met)

                if pass_deltaR(e, tau, jet, 'ETau') == 1 :
                    plot_variable(selectedTau+'_SS_dRcut', e, tau, jet, met)

                    if met.Pt() > event_cut['metcut'] : #SS highMET                                                                                   
                        plot_variable(selectedTau+'_SS_dRcut_highMET', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #SS highMET lowMt                                                     
                            plot_variable(selectedTau+'_SS_dRcut_highMET_lowMt', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #SS highMET highMt        
                            plot_variable(selectedTau+'_SS_dRcut_highMET_highMt', e, tau, jet, met)

                    if met.Pt() < event_cut['metcut'] : #SS lowMET                                                                          
                        plot_variable(selectedTau+'_SS_dRcut_lowMET', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #SS lowMET lowMt                                                           
                            plot_variable(selectedTau+'_SS_dRcut_lowMET_lowMt', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #SS lowMET highMt                        
                            plot_variable(selectedTau+'_SS_dRcut_lowMET_highMt', e, tau, jet, met)

        if ( isData == 0 and isJetHTEvent == 1 ) \
           or ( isData == 1 and ( isJetHTSample == 1 and isJetHTEvent == 1 ) ) : #JetHT trigger only

            if s_electron[0].charge*eCleanedTau[0].charge < 0 : #OS
                plot_variable(selectedTau+'_OS_JetHT', e, tau, jet, met)

                if e.DeltaR(tau) > 0.8 : #OS resolved
                    plot_variable(selectedTau+'_OS_resolved_JetHT', e, tau, jet, met)

                if met.Pt() > event_cut['metcut'] : #OS highMET - nondRcut
                    plot_variable(selectedTau+'_OS_highMET_JetHT', e, tau, jet, met)

                    if e.DeltaR(tau) > 0.8 : #OS highMET resolved
                        plot_variable(selectedTau+'_OS_highMET_resolved_JetHT', e, tau, jet, met)

                if pass_deltaR(e, tau, jet, 'ETau') == 1 :
                    plot_variable(selectedTau+'_OS_dRcut_JetHT', e, tau, jet, met)

                    if met.Pt() > event_cut['metcut'] : #OS highMET
                        plot_variable(selectedTau+'_OS_dRcut_highMET_JetHT', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #OS highMET lowMt
                            plot_variable(selectedTau+'_OS_dRcut_highMET_lowMt_JetHT', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #OS highMET highMt                                                               
                            plot_variable(selectedTau+'_OS_dRcut_highMET_highMt_JetHT', e, tau, jet, met)

                    if met.Pt() < event_cut['metcut'] : #OS lowMET   
                        plot_variable(selectedTau+'_OS_dRcut_lowMET_JetHT', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #OS lowMET lowMt                                                                                         
                            plot_variable(selectedTau+'_OS_dRcut_lowMET_lowMt_JetHT', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #OS lowMET highMt          
                            plot_variable(selectedTau+'_OS_dRcut_lowMET_highMt_JetHT', e, tau, jet, met)

            if s_electron[0].charge*eCleanedTau[0].charge > 0 : #SS
                plot_variable(selectedTau+'_SS_JetHT', e, tau, jet, met)

                if met.Pt() > event_cut['metcut'] :
                    plot_variable(selectedTau+'_SS_highMET_JetHT', e, tau, jet , met)

                if pass_deltaR(e, tau, jet, 'ETau') == 1 :
                    plot_variable(selectedTau+'_SS_dRcut_JetHT', e, tau, jet, met)

                    if met.Pt() > event_cut['metcut'] : #SS highMET                                                                                   
                        plot_variable(selectedTau+'_SS_dRcut_highMET_JetHT', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #SS highMET lowMt                                                     
                            plot_variable(selectedTau+'_SS_dRcut_highMET_lowMt_JetHT', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #SS highMET highMt        
                            plot_variable(selectedTau+'_SS_dRcut_highMET_highMt_JetHT', e, tau, jet, met)

                    if met.Pt() < event_cut['metcut'] : #SS lowMET                                                                          
                        plot_variable(selectedTau+'_SS_dRcut_lowMET_JetHT', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #SS lowMET lowMt                                                           
                            plot_variable(selectedTau+'_SS_dRcut_lowMET_lowMt_JetHT', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #SS lowMET highMt                        
                            plot_variable(selectedTau+'_SS_dRcut_lowMET_highMt_JetHT', e, tau, jet, met)


        if ( isData == 0 and ( isJetHTEvent == 0 and isSingleElectronEvent == 1 ) ) \
           or ( isData == 1 and ( isSingleElectronSample == 1 and ( isJetHTEvent == 0 and isSingleElectronEvent == 1 ) ) ) : #SingleElectron trigger only

            if s_electron[0].charge*eCleanedTau[0].charge < 0 : #OS
                plot_variable(selectedTau+'_OS_SingleElectron', e, tau, jet, met)

                if e.DeltaR(tau) > 0.8 : #OS resolved
                    plot_variable(selectedTau+'_OS_resolved_SingleElectron', e, tau, jet, met)

                if met.Pt() > event_cut['metcut'] : #OS highMET - nondRcut
                    plot_variable(selectedTau+'_OS_highMET_SingleElectron', e, tau, jet, met)

                    if e.DeltaR(tau) > 0.8 : #OS highMET resolved
                        plot_variable(selectedTau+'_OS_highMET_resolved_SingleElectron', e, tau, jet, met)

                if pass_deltaR(e, tau, jet, 'ETau') == 1 :
                    plot_variable(selectedTau+'_OS_dRcut_SingleElectron', e, tau, jet, met)

                    if met.Pt() > event_cut['metcut'] : #OS highMET
                        plot_variable(selectedTau+'_OS_dRcut_highMET_SingleElectron', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #OS highMET lowMt
                            plot_variable(selectedTau+'_OS_dRcut_highMET_lowMt_SingleElectron', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #OS highMET highMt                                                               
                            plot_variable(selectedTau+'_OS_dRcut_highMET_highMt_SingleElectron', e, tau, jet, met)

                    if met.Pt() < event_cut['metcut'] : #OS lowMET   
                        plot_variable(selectedTau+'_OS_dRcut_lowMET_SingleElectron', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #OS lowMET lowMt                                                                                         
                            plot_variable(selectedTau+'_OS_dRcut_lowMET_lowMt_SingleElectron', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #OS lowMET highMt          
                            plot_variable(selectedTau+'_OS_dRcut_lowMET_highMt_SingleElectron', e, tau, jet, met)

            if s_electron[0].charge*eCleanedTau[0].charge > 0 : #SS
                plot_variable(selectedTau+'_SS_SingleElectron', e, tau, jet, met)

                if met.Pt() > event_cut['metcut'] :
                    plot_variable(selectedTau+'_SS_highMET_SingleElectron', e, tau, jet , met)

                if pass_deltaR(e, tau, jet, 'ETau') == 1 :
                    plot_variable(selectedTau+'_SS_dRcut_SingleElectron', e, tau, jet, met)

                    if met.Pt() > event_cut['metcut'] : #SS highMET                                                                                   
                        plot_variable(selectedTau+'_SS_dRcut_highMET_SingleElectron', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #SS highMET lowMt                                                     
                            plot_variable(selectedTau+'_SS_dRcut_highMET_lowMt_SingleElectron', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #SS highMET highMt        
                            plot_variable(selectedTau+'_SS_dRcut_highMET_highMt_SingleElectron', e, tau, jet, met)

                    if met.Pt() < event_cut['metcut'] : #SS lowMET                                                                          
                        plot_variable(selectedTau+'_SS_dRcut_lowMET_SingleElectron', e, tau, jet, met)

                        if Mt(e, met) < event_cut['mtcut'] : #SS lowMET lowMt                                                           
                            plot_variable(selectedTau+'_SS_dRcut_lowMET_lowMt_SingleElectron', e, tau, jet, met)

                        if Mt(e, met) > event_cut['mtcut'] : #SS lowMET highMt                        
                            plot_variable(selectedTau+'_SS_dRcut_lowMET_highMt_SingleElectron', e, tau, jet, met)


def mutau_channel():

    isMuTau = 0

    if s_muon[0].charge*s_tauMuclean[0].charge < 0:

        mu = get_TLorentzVector(s_muon[0])
        tau = get_TLorentzVector(s_tauMuclean[0])
        jet = get_TLorentzVector(s_jet[0])

        isJetHTEvent = 0
        isSingleMuonEvent = 0

        if ( jet.Pt() > 510 and ( isHT == 1 or isSingleJet500 == 1 ) ) : isJetHTEvent = 1
        if ( mu.Pt() > 52 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : isSingleMuonEvent = 1

        if isData == 0 and ( isJetHTEvent == 1 ) \
           or ( isData == 1 and isJetHTSample == 1 and isJetHTEvent == 1 ) :

            if pass_deltaR(mu, tau, jet, 'MuTau') == 1 :

                if Mt(mu,met) < 50 :

                    if met.Pt() > event_cut['metcut'] :

                        plot_variable('MuTau', mu, tau, jet, met)

                        isMuTau = 1

    return isMuTau

def ee_channel():

    isEE = 0

    if s_isoelectron[0].charge*s_isoelectron[1].charge < 0 :
        
        e1 = get_TLorentzVector(s_isoelectron[0])
        e2 = get_TLorentzVector(s_isoelectron[1])
        jet = get_TLorentzVector(s_jet[0])

        isJetHTEvent = 0
        isSingleElectronEvent = 0

        sf = get_sf(s_isoelectron[0], "ele") if isData == 0 else 1

        if ( jet.Pt() > 510 and ( isHT == 1  or isSingleJet500 == 1 ) ) : isJetHTEvent = 1
        if ( e1.Pt() > 37 and isIsoEle == 1 ) : isSingleElectronEvent = 1

        if ( isData == 0 and ( isJetHTEvent == 1 or isSingleElectronEvent == 1 ) ) \
           or ( isData == 1 and isJetHTSample == 1 and isJetHTEvent == 1 ) \
           or ( isData == 1 and ( isSingleElectronSample == 1 and isJetHTSample == 0 and isSingleElectronEvent == 1 ) ) :
             
            if pass_deltaR(e1, e2, jet, 'EE') == 1 :
                if met.Pt() > event_cut['metcut'] :
                    plot_variable('EE', e1, e2, jet, met, sf)

                    isEE = 1

        if ( isData == 0 and ( isJetHTEvent == 1 or isSingleElectronEvent == 1 ) ) \
           or ( isData == 1 and ( isJetHTSample == 1 and isJetHTEvent == 1 ) ) \
           or ( isData == 1 and ( isSingleElectronSample == 1 and ( isJetHTEvent == 0 and isSingleElectronEvent == 1 ) ) ) :

            if e1.Pt() > 36 and e2.Pt() > 36 and jet.Pt() > event_cut['jetpt'] :

                if e1.DeltaR(e2) > 0.4 :
                    plot_variable('EE_resolved', e1, e2, jet, met, sf)

                    if e1.Pt() > 100 and e2.Pt() > 100 :
                        plot_variable('EE_resolved_highPtEle', e1, e2, jet, met, sf)

                if met.Pt() < event_cut['metcut'] :
                    plot_variable('EE_lowMET', e1, e2, jet, met, sf)

                    if e1.DeltaR(e2) > 0.4 :
                        plot_variable('EE_lowMET_resolved', e1, e2, jet, met, sf)

                if met.Pt() > event_cut['metcut'] :
                    plot_variable('EE_highMET', e1, e2, jet, met, sf)

                    if e1.DeltaR(e2) > 0.4 :
                        plot_variable('EE_highMET_resolved', e1, e2, jet, met, sf)

                        if e1.Pt() > 100 and e2.Pt() > 100 :
                            plot_variable('EE_highMET_resolved_highPtEle', e1, e2, jet, met, sf)

        if ( isData == 0 and ( isJetHTEvent == 1 ) ) or ( isData == 1 and ( isJetHTSample == 1 and isJetHTEvent == 1 ) ) :

            if e1.Pt() > 36 and e2.Pt() > 36 and jet.Pt() > event_cut['jetpt'] :

                if e1.DeltaR(e2) > 0.4 :
                    plot_variable('EE_resolved_JetHT', e1, e2, jet, met, sf)

                    if e1.Pt() > 100 and e2.Pt() > 100 :
                        plot_variable('EE_resolved_highPtEle_JetHT', e1, e2, jet, met, sf)

                if met.Pt() > event_cut['metcut'] :
                    plot_variable('EE_highMET_JetHT', e1, e2, jet, met, sf)

                    if e1.DeltaR(e2) > 0.4 :
                        plot_variable('EE_highMET_resolved_JetHT', e1, e2, jet, met, sf)

                        if e1.Pt() > 100 and e2.Pt() > 100 :
                            plot_variable('EE_highMET_resolved_highPtEle_JetHT', e1, e2, jet, met, sf)

        if ( isData == 0 and ( isSingleElectronEvent == 1 ) ) or ( isData == 1 and ( isSingleElectronSample == 1 and isSingleElectronEvent == 1 ) ) :

            if e1.Pt() > 36 and e2.Pt() > 36 and jet.Pt() > event_cut['jetpt'] :

                if e1.DeltaR(e2) > 0.4 :
                    plot_variable('EE_resolved_SingleElectron', e1, e2, jet, met, sf)

                    if e1.Pt() > 100 and e2.Pt() > 100 :
                        plot_variable('EE_resolved_highPtEle_SingleElectron', e1, e2, jet, met, sf)

                if met.Pt() > event_cut['metcut'] :
                    plot_variable('EE_highMET_SingleElectron', e1, e2, jet, met, sf)

                    if e1.DeltaR(e2) > 0.4 :
                        plot_variable('EE_highMET_resolved_SingleElectron', e1, e2, jet, met, sf)

                        if e1.Pt() > 100 and e2.Pt() > 100 :
                            plot_variable('EE_highMET_resolved_highPtEle_SingleElectron', e1, e2, jet, met, sf)

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

            plot_variable('EMu_Baseline', e, mu, jet, met)

            if pass_deltaR(e, mu, jet, 'EMu') == 1 :

                if met.Pt() > event_cut['metcut'] :
                    plot_variable('EMu', e, mu, jet, met)
                    
                    isEMu = 1

    return isEMu


regions = ['MuMu', 'MuTau', 'EMu', 'EE']
# etauR = ['ETau_OS','ETau_SS','ETau_Altered_OS', 'ETau_Altered_SS']
# etauM = ['highMt', 'lowMt']
# trig = ['SingleElectron', 'JetHT']

# for e in etauR :
#     regions.append(e)
#     regions.append(e+'_dRcut')
#     regions.append(e+'_dRcut_highMET')
#     regions.append(e+'_dRcut_lowMET')
#     for t in trig:
#         regions.append(e+'_'+t)
#         regions.append(e+'_dRcut_'+t)
#         regions.append(e+'_dRcut_highMET_'+t)
#         regions.append(e+'_dRcut_lowMET_'+t)
#     for m in etauM :
#         regions.append(e+'_dRcut_highMET_'+m)
#         regions.append(e+'_dRcut_lowMET_'+m)
#         for t in trig:
#             regions.append(e+'_dRcut_highMET_'+m+'_'+t)
#             regions.append(e+'_dRcut_lowMET_'+m+'_'+t)

# regions.append('ETau_OS_highMET')
# regions.append('ETau_OS_highMET_resolved')
# regions.append('ETau_OS_resolved')
# regions.append('ETau_SS_highMET')
# regions.append('ETau_Altered_OS_highMET')
# regions.append('ETau_Altered_OS_highMET_resolved')
# regions.append('ETau_Altered_OS_resolved')
# regions.append('ETau_Altered_SS_highMET')
# regions.append('ETau_OS_highMET_SingleElectron')
# regions.append('ETau_OS_highMET_resolved_SingleElectron')
# regions.append('ETau_OS_resolved_SingleElectron')
# regions.append('ETau_SS_highMET_SingleElectron')
# regions.append('ETau_Altered_OS_highMET_SingleElectron')
# regions.append('ETau_Altered_OS_highMET_resolved_SingleElectron')
# regions.append('ETau_Altered_OS_resolved_SingleElectron')
# regions.append('ETau_Altered_SS_highMET_SingleElectron')
# regions.append('ETau_OS_highMET_JetHT')
# regions.append('ETau_OS_highMET_resolved_JetHT')
# regions.append('ETau_OS_resolved_JetHT')
# regions.append('ETau_SS_highMET_JetHT')
# regions.append('ETau_Altered_OS_highMET_JetHT')
# regions.append('ETau_Altered_OS_highMET_resolved_JetHT')
# regions.append('ETau_Altered_OS_resolved_JetHT')
# regions.append('ETau_Altered_SS_highMET_JetHT')
regions.append('MuMu_Baseline')
regions.append('EMu_Baseline')
regions.append('EE_highMET')
regions.append('EE_highMET_resolved')
regions.append('EE_lowMET')
regions.append('EE_lowMET_resolved')
regions.append('EE_resolved')
regions.append('EE_highMET_resolved_highPtEle')
regions.append('EE_resolved_highPtEle')
regions.append('EE_highMET_JetHT')
regions.append('EE_highMET_resolved_JetHT')
regions.append('EE_resolved_JetHT')
regions.append('EE_highMET_resolved_highPtEle_JetHT')
regions.append('EE_resolved_highPtEle_JetHT')
regions.append('EE_highMET_SingleElectron')
regions.append('EE_highMET_resolved_SingleElectron')
regions.append('EE_resolved_SingleElectron')
regions.append('EE_highMET_resolved_highPtEle_SingleElectron')
regions.append('EE_resolved_highPtEle_SingleElectron')

book_trigger_histogram('ETau_TriggerStudy_OS_noGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_noGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_noGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_noGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_noGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_noGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_noGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_noGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_yesGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_yesGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_yesGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_yesGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_yesGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_yesGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_yesGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_yesGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_noGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_noGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_noGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_noGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_noGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_noGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_noGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_noGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRjcut_yesGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_noGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_noGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_noGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_noGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_noGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_noGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_noGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_noGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRlcut_yesGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_noGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_noGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_noGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_noGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_noGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_noGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_noGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_noGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_yesGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_yesGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_yesGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_yesGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_yesGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_yesGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_yesGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_yesGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_noGenMatch_isHTMHT')

book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_JetHT')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_SingleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_SingleJet450')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_SingleE')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_isEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_lowMt_yesGenMatch_isHTMHT')


book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_EpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_JetpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_SingleE_pTcut')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEle_pTcut')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet_EpTcut')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_isEleJet_JetpTcut')

book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_SingleE_pTcut')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isEle_pTcut')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isEleJet_EpTcut')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_noGenMatch_isEleJet_JetpTcut')

book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_50EpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_50EpTcut_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_50EpTcut_isEleJet_isIsoEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_180JetPtcut_60EpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_180JetPtcut_60EpTcut_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_180JetPtcut_60EpTcut_isEleJet_isIsoEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_170JetPtcut_70EpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_170JetPtcut_70EpTcut_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_170JetPtcut_70EpTcut_isEleJet_isIsoEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_165JetPtcut_80EpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_165JetPtcut_80EpTcut_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_165JetPtcut_80EpTcut_isEleJet_isIsoEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_90EpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_90EpTcut_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_90EpTcut_isEleJet_isIsoEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_95EpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_95EpTcut_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_95EpTcut_isEleJet_isIsoEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_0EpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_0EpTcut_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_200JetPtcut_0EpTcut_isEleJet_isIsoEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_100EpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_100EpTcut_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_100EpTcut_isEleJet_isIsoEle')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_50EpTcut_Before')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_50EpTcut_isEleJet')
book_trigger_histogram('ETau_TriggerStudy_OS_dRcut_highMET_yesGenMatch_0JetPtcut_50EpTcut_isEleJet_isIsoEle')
book_histogram()

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
        chain3.Add(inputFileName)
        chain4.Add(inputFileName)
        chain5.Add(inputFileName)

    chain6.Add(inputFileName)

#------- Adding friends to the main chain -------#


fchain.AddFriend(chain2)
if isData == 0:
    fchain.AddFriend(chain3)
    fchain.AddFriend(chain4)
    fchain.AddFriend(chain5)
    
fchain.AddFriend(chain6)

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

#----------- Event loop ----------#


for iev in range(fchain.GetEntries()): # Be careful!!!                                                               

    fchain.GetEntry(iev)

    mets = fchain.GetBranch("Mets")
    met_pt = mets.GetLeaf('pt').GetValue()
    met_phi = mets.GetLeaf('phi').GetValue()

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

    #-------------- Trigger definitions -----------#

    #isSingleJet500 = fchain.GetLeaf('isSingleJet500').GetValue()
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
                    h['hJetPt'].Fill(ijet.pt, weight)
                    h['hDeepjet'].Fill(ijet.deepjet, weight)
                    s_jet+=[ijet]
                    iht = iht + ijet.pt
                    if ijet.deepjet >= 0.7476:
                        h['hBJetPt'].Fill(ijet.pt, weight)
                        s_bjet+=[ijet]

    s_muon = []
    s_isomuon = []
    s_nonisomuon = []

    if muons.size() > 0:
        for i in range(muons.size()):
            imuon = muons.at(i)
            if abs(imuon.eta) < 2.4 :
                if imuon.id >= 1: #loose Muons
                    h['hMuonPt'].Fill(imuon.pt, weight) 
                    s_muon+=[imuon]
                    if imuon.iso <= 0.25:
                        h['hIsoMuonPt'].Fill(imuon.pt, weight)
                        s_isomuon+=[imuon]
                    if imuon.iso > 0.25:
                        h['hNonIsoMuonPt'].Fill(imuon.pt, weight)
                        s_nonisomuon+=[imuon]

    s_electron = []
    s_isoelectron = []
    s_nonisoelectron = []

    if electrons.size() > 0:
        for i in range(electrons.size()):
            ielectron = electrons.at(i)
            if abs(ielectron.eta) < 2.5 :
                if ielectron.id >= 1 :
                    h['hElectronPt'].Fill(ielectron.pt, weight)
                    s_electron+=[ielectron]
                    if ielectron.iso >= 1:
                        h['hIsoElectronPt'].Fill(ielectron.pt, weight)
                        s_isoelectron+=[ielectron]
                    if ielectron.iso == 0:
                        h['hNonIsoElectronPt'].Fill(ielectron.pt, weight)
                        s_nonisoelectron+=[ielectron]

    s_tauEclean = []
    s_tauEcleanAltered = []

    if tausECleaned.size()>0:
        for i in range(tausECleaned.size()):
            itau = tausECleaned.at(i)
            if abs(itau.eta) < 2.3 :
                if itau.mvaid >= 4:
                    h['hTauECleanedPt'].Fill(itau.pt, weight)
                    s_tauEclean+=[itau]
                if itau.mvaid < 4 and itau.mvaid >= 1 :
                    h['hTauECleanedAlteredPt'].Fill(itau.pt, weight)
                    s_tauEcleanAltered+=[itau]
                    

    s_tauMuclean = []

    if tausMCleaned.size()>0:
        for i in range(tausMCleaned.size()):
            itau = tausMCleaned.at(i)
            if abs(itau.eta) < 2.3 :
                if itau.mvaid >= 4 :
                    h['hTauMuCleanedPt'].Fill(itau.pt, weight)
                    s_tauMuclean+=[itau]


    # ---------- Event Selections --------- #

    ## b-jet veto not implemented

#    if pvFilter == haloFilter == haloFilter == hbheIsoFilter == ecalTPFilter == badMuonFilter == ecalBadCalFilter == 0 : continue

    if len(s_isomuon) >= 2 and len(s_jet) >= 1 : 
        if mumu_channel() == 1: continue

    if len(s_isomuon) >= 1 and len(s_isoelectron) >= 1 and len(s_jet) >= 1 : 
        if emu_channel() == 1 : continue

    if len(s_muon) >= 1 and len(s_tauMuclean) >= 1 and len(s_jet) >= 1 : 
        if mutau_channel() == 1 : continue

    if len(s_isoelectron) >=2 and len(s_jet) >= 1 :
        if ee_channel() == 1 : continue

    if len(s_electron) >= 1 and len(s_tauEclean) >= 1 and len(s_jet) >= 1 :
#        if isData == 0 : etau_trigger_study()
        etau_trigger_study()
#        etau_channel(s_tauEclean, "nominal")

    # if len(s_electron) >= 1 and len(s_tauEcleanAltered) >= 1 and len(s_jet) >= 1 and len(s_bjet) == 0 :
    #     etau_channel(s_tauEcleanAltered, "altered")

out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))

# result = np.array(regions)
# print(result)
# print(regions)
