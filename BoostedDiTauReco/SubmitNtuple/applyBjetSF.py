import ROOT, sys, os
import numpy as np
import time
import correctionlib

start_time = time.time()

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

outputTitle = "h_applyBjetSF"

sfDir = os.path.join('/cvmfs','cms.cern.ch','rsync','cms-nanoAOD','jsonpog-integration','POG','BTV','2017_UL')
btvjson = correctionlib.CorrectionSet.from_file(os.path.join(sfDir, 'btagging.json.gz'))

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

WP = 'T'

event_cut = {
    'jetPt': 70,
    'dRl': 0.4,
    'dRltau': 0.05,
    'dRlj': 0.8,
    'metcut': 100.0,
    'mtcut': 50.0,
    'dPhiml': 1,
    'dPhimj': 2,
    'mass' : 5
}

h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)

def define_eff_histogram(region):

    h[region+"_BFlavour_JetPtEta"] = ROOT.TH2F (region+"_BFlavour_JetPtEta", region+"_BFlavour_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)
    h[region+"_BFlavour_BTagged_JetPtEta"] = ROOT.TH2F (region+"_BFlavour_BTagged_JetPtEta", region+"_BFlavour_BTagged_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)

    h[region+"_CFlavour_JetPtEta"] = ROOT.TH2F (region+"_CFlavour_JetPtEta", region+"_CFlavour_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)
    h[region+"_CFlavour_BTagged_JetPtEta"] = ROOT.TH2F (region+"_CFlavour_BTagged_JetPtEta", region+"_CFlavour_BTagged_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)

    h[region+"_LFlavour_JetPtEta"] = ROOT.TH2F (region+"_LFlavour_JetPtEta", region+"_LFlavour_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)
    h[region+"_LFlavour_BTagged_JetPtEta"] = ROOT.TH2F (region+"_LFlavour_BTagged_JetPtEta", region+"_LFlavour_BTagged_JetPtEta ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)

    h[region+'_BFlavour_BTagging_Efficiency'] = ROOT.TH2F (region+"_BFlavour_BTagging_Efficiency", region+"_BFlavour_BTagging_Efficiency ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)
    h[region+'_CFlavour_BTagging_Efficiency'] = ROOT.TH2F (region+"_CFlavour_BTagging_Efficiency", region+"_CFlavour_BTagging_Efficiency ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)
    h[region+'_LFlavour_BTagging_Efficiency'] = ROOT.TH2F (region+"_LFlavour_BTagging_Efficiency", region+"_LFlavour_BTagging_Efficiency ; P_{T} (GeV) ; Eta", 2000, 0, 2000, 100, 0, 2.5)


def define_event_histogram(region):

    h[region+"_Count"] = ROOT.TH1F (region+"_Count", region+"_Count ; Events ; Events ", 1, 0, 1)

    h[region+"_Mass"] = ROOT.TH1F (region+"_Mass", region+"_Mass ; M_{vis.} (GeV) ; Events ", 100, 0, 100)
    h[region+"_Lepton1Pt"] = ROOT.TH1F (region+"_Lepton1Pt", region+"_Lepton1Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_Lepton2Pt"] = ROOT.TH1F (region+"_Lepton2Pt", region+"_Lepton2Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_JetPt"] = ROOT.TH1F (region+"_JetPt", region+"_JetPt ; JetP_{T} (GeV) ; Events ", 2000, 0, 2000)
    h[region+"_Jet2Pt"] = ROOT.TH1F (region+"_Jet2Pt", region+"_Jet2Pt ; Jet2P_{T} (GeV) ; Events ", 2000, 0, 2000)
    h[region+"_MetPt"] = ROOT.TH1F (region+"_MetPt", region+"_MetPt ; MET (GeV) ; Events ", 500, 0, 500)

    h[region+"_BJetPt"] = ROOT.TH1F (region+"_BJetPt", region+"_BJetPt ; BJetP_{T} (GeV) ; Events ", 2000, 0, 2000)
    h[region+"_BJet2Pt"] = ROOT.TH1F (region+"_BJet2Pt", region+"_BJet2Pt ; BJet2P_{T} (GeV) ; Events ", 2000, 0, 2000)

    h[region+"_Nj"] = ROOT.TH1F (region+"_Nj", region+"_Nj ; N_{j} ; Events ", 10, 0, 10)
    h[region+"_Nbj"] = ROOT.TH1F (region+"_Nbj", region+"_Nbj ; N_{bj} ; Events ", 10, 0, 10)

    h[region+"_dRl"] = ROOT.TH1F (region+"_dRl", region+"_dRl ; dR(leptons) ; Events", 100, 0, 5)
    h[region+"_dRj"] = ROOT.TH1F (region+"_dRj", region+"_dRj ; dR(jet, ditau) ; Events", 100, 0, 5)
    h[region+"_dRj2"] = ROOT.TH1F (region+"_dRj2", region+"_dRj2 ; dR(jet2, ditau) ; Events", 100, 0, 5)
    h[region+"_dRbj"] = ROOT.TH1F (region+"_dRbj", region+"_dRbj ; dR(bjet, ditau) ; Events", 100, 0, 5)
    h[region+"_dRbj2"] = ROOT.TH1F (region+"_dRbj2", region+"_dRbj2 ; dR(bjet2, ditau) ; Events", 100, 0, 5)

    h[region+'_deepjet'] = ROOT.TH1F (region+"_deepjet", region+"_deepjet ; deepjet score ; Events", 100, 0, 1)
    h[region+'_deepjet1'] = ROOT.TH1F (region+"_deepjet1", region+"_deepjet1 ; deepjet1 score ; Events", 100, 0, 1)
    h[region+'_deepjet2'] = ROOT.TH1F (region+"_deepjet2", region+"_deepjet2 ; deepjet2 score ; Events", 100, 0, 1)

plot_regions = [
    'lowMET_dRcut',
    'highMET_dRcut',
    'lowMET_invdRcut',
    'highMET_invdRcut',
    'Baseline',
    'lowMET',
    'highMET',
    'dRcut',
]

for r in plot_regions :

    define_eff_histogram(r)
    define_eff_histogram(r+"_0bjet")
    define_eff_histogram(r+"_bjet")
    define_event_histogram(r)
    define_event_histogram(r+"_noCorr")
    define_event_histogram(r+"_0bjet")
    define_event_histogram(r+"_0bjet_noCorr")
    define_event_histogram(r+"_bjet")
    define_event_histogram(r+"_bjet_noCorr")


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


def pasbjsaseline(l1, l2, j):

    if j.Pt() > event_cut['jetPt'] and (l1+l2).M() > event_cut["mass"]:
        return 1
    else:
        return -9999



def find_sf(js):

    if js.jetflavour != 0 :

        deepjet_col = 'deepJet_comb'

    if js.jetflavour == 0 :

        deepjet_col = 'deepJet_incl'

    jet_sf = btvjson[deepjet_col].evaluate("central", WP, js.jetflavour, abs(js.eta), js.pt)

    return jet_sf


def fill_efficiency(region):

    for i in range(len(js)) :

        if js[i].jetflavour == 5 :
            h[region+'_BFlavour_JetPtEta'].Fill(js[i].pt, abs(js[i].eta))
            if js[i].deepjet >= 0.7476 :
                h[region+'_BFlavour_BTagged_JetPtEta'].Fill(js[i].pt, abs(js[i].eta))

        if js[i].jetflavour == 4 :
            h[region+'_CFlavour_JetPtEta'].Fill(js[i].pt, abs(js[i].eta))
            if js[i].deepjet >= 0.7476:
                h[region+'_CFlavour_BTagged_JetPtEta'].Fill(js[i].pt, abs(js[i].eta))

        if js[i].jetflavour == 0 :
            h[region+'_LFlavour_JetPtEta'].Fill(js[i].pt, abs(js[i].eta))
            if js[i].deepjet >= 0.7476:
                h[region+'_LFlavour_BTagged_JetPtEta'].Fill(js[i].pt, abs(js[i].eta))


def measure_efficiency(region):

    h[region+'_BFlavour_BTagging_Efficiency'] = h[region+'_BFlavour_BTagged_JetPtEta'].Clone(region+'_BFlavour_BTagging_Efficiency')
    h[region+'_BFlavour_BTagging_Efficiency'].Divide(h[region+'_BFlavour_JetPtEta'])

    h[region+'_CFlavour_BTagging_Efficiency'] = h[region+'_CFlavour_BTagged_JetPtEta'].Clone(region+'_CFlavour_BTagging_Efficiency')
    h[region+'_CFlavour_BTagging_Efficiency'].Divide(h[region+'_CFlavour_JetPtEta'])

    h[region+'_LFlavour_BTagging_Efficiency'] = h[region+'_LFlavour_BTagged_JetPtEta'].Clone(region+'_LFlavour_BTagging_Efficiency')
    h[region+'_LFlavour_BTagging_Efficiency'].Divide(h[region+'_LFlavour_JetPtEta'])


def get_efficiency(region, jet):

    if jet.jetflavour == 5 :
        binPt = h[region+'_BFlavour_BTagging_Efficiency'].GetXaxis().FindBin(jet.pt)
        binEta = h[region+'_BFlavour_BTagging_Efficiency'].GetYaxis().FindBin(jet.eta)
        eff = h[region+'_BFlavour_BTagging_Efficiency'].GetBinContent(binPt, binEta)

    if jet.jetflavour == 4 :
        binPt = h[region+'_CFlavour_BTagging_Efficiency'].GetXaxis().FindBin(jet.pt)
        binEta = h[region+'_CFlavour_BTagging_Efficiency'].GetYaxis().FindBin(jet.eta)
        eff = h[region+'_CFlavour_BTagging_Efficiency'].GetBinContent(binPt, binEta)

    if jet.jetflavour == 0 :
        binPt = h[region+'_LFlavour_BTagging_Efficiency'].GetXaxis().FindBin(jet.pt)
        binEta = h[region+'_LFlavour_BTagging_Efficiency'].GetYaxis().FindBin(jet.eta)
        eff = h[region+'_LFlavour_BTagging_Efficiency'].GetBinContent(binPt, binEta)

    return eff


def EMu_regions_efficiency():

    e, mu, j, m = get_TLorentzVector(es[0], mus[0], js[0], met_pt, met_phi)

    trigger = [0,0]

    if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : trigger[0] = 1
    
    if ( ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) and isMuonEG == 1 ) : trigger[1] = 1

    if ( isData == 0 and ( trigger[0] == 1 or trigger[1] == 1 ) ) :

        if pasbjsaseline(e, mu, j) == 1 :

            region = 'Baseline'
            fill_efficiency(region)

            if m.Pt() > event_cut['metcut'] :

                region = 'highMET'
                fill_efficiency(region)

            if m.Pt() < event_cut['metcut'] :

                region = 'lowMET'
                fill_efficiency(region)

            if j.DeltaR(e+mu) > event_cut['dRlj'] :

                if e.DeltaR(mu) < 0.4 :

                    region = 'dRcut'
                    fill_efficiency(region)

                    if m.Pt() < event_cut['metcut'] :

                        region = 'lowMET_dRcut'
                        fill_efficiency(region)

                    if m.Pt() > event_cut['metcut'] :

                        region = 'highMET_dRcut'
                        fill_efficiency(region)

                if e.DeltaR(mu) > 0.4 and e.DeltaR(mu) < 0.6 :

                    if m.Pt() > event_cut['metcut'] :

                        region = 'highMET_invdRcut'
                        fill_efficiency(region)

                    if m.Pt() < event_cut['metcut'] :

                        region = 'lowMET_invdRcut'
                        fill_efficiency(region)



def get_weight(region, jets):

    pMC = 1.0
    pData = 1.0
    bweight = 1.0

    for js in jets:

        eff = get_efficiency(region, js)
        sf = find_sf(js)

        if js.deepjet >= 0.7476 :
            pMC = pMC*eff
            pData = pData*sf*eff

        else :
            pMC = pMC*(1.0-eff)
            pData = pData*sf*(1.0-eff)

    if pMC != 0.0 :
        bweight = pData/pMC

    return bweight


def EMu_events():

    e, mu, j, m = get_TLorentzVector(es[0], mus[0], js[0], met_pt, met_phi)

    trigger = [0,0]

    if ( mu.Pt() > 50 and isMu == 1 ) or ( mu.Pt() > 27 and isIsoMu == 1 ) : trigger[0] = 1

    if ( ( ( mu.Pt() > 8 and e.Pt() > 23 ) or ( mu.Pt() > 23 and e.Pt() > 12 ) ) and isMuonEG == 1 ) : trigger[1] = 1

    if ( isData == 0 and ( trigger[0] == 1 or trigger[1] == 1 ) ) :

        if pasbjsaseline(e, mu, j) == 1 :

            j2 = 0
            bj = 0
            bj2 = 0

            if len(js) > 1 :

                j2 = ROOT.TLorentzVector()
                j2.SetPtEtaPhiM(js[1].pt, js[1].eta, js[1].phi, js[1].mass)

            if len(bjs) > 0 :

                bj = ROOT.TLorentzVector()
                bj.SetPtEtaPhiM(bjs[0].pt, bjs[0].eta, bjs[0].phi, bjs[0].mass)

                if len(bjs) > 1 :

                    bj2 = ROOT.TLorentzVector()
                    bj2.SetPtEtaPhiM(bjs[1].pt, bjs[1].eta, bjs[1].phi, bjs[1].mass)

            region = 'Baseline'
            measure_efficiency(region)
            bweight = get_weight(region, js)

            plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
            plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)
        
            if m.Pt() < event_cut['metcut'] :

                region = 'lowMET'
                measure_efficiency(region)
                bweight = get_weight(region, js)

                plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                if len(bjs) == 0 :

                    region = 'lowMET_0bjet'
                    measure_efficiency(region)
                    bweight = get_weight(region, js)
                
                    plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                    plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                if len(bjs) > 0 :

                    region = 'lowMET_bjet'
                    measure_efficiency(region)
                    bweight = get_weight(region, js)
                
                    plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                    plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)


            if m.Pt() > event_cut['metcut'] :

                region = 'highMET'
                measure_efficiency(region)
                bweight = get_weight(region, js)

                plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                if len(bjs) == 0 :

                    region = 'highMET_0bjet'
                    measure_efficiency(region)
                    bweight = get_weight(region, js)
                
                    plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                    plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                if len(bjs) > 0 :

                    region = 'highMET_bjet'
                    measure_efficiency(region)
                    bweight = get_weight(region, js)

                    plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                    plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)


            if j.DeltaR(e+mu) > event_cut['dRlj'] :

                if e.DeltaR(mu) < 0.4 :

                    region = 'dRcut'
                    measure_efficiency(region)
                    bweight = get_weight(region, js)

                    plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                    plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                    if len(bjs) == 0 :

                        region = 'dRcut_0bjet'
                        measure_efficiency(region)
                        bweight = get_weight(region, js)
                
                        plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                        plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                    if len(bjs) > 0 :
                        
                        region = 'dRcut_bjet'
                        measure_efficiency(region)
                        bweight = get_weight(region, js)

                        plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                        plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)


                    if m.Pt() < event_cut['metcut'] :
                        
                        region = 'lowMET_dRcut'
                        measure_efficiency(region)
                        bweight = get_weight(region, js)

                        plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                        plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                        if len(bjs) == 0 :

                            region = 'lowMET_dRcut_0bjet'
                            measure_efficiency(region)
                            bweight = get_weight(region, js)

                            plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                            plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                        if len(bjs) > 0 :

                            region = 'lowMET_dRcut_bjet'
                            measure_efficiency(region)
                            bweight = get_weight(region, js)

                            plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                            plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)


                    if m.Pt() > event_cut['metcut'] :

                        region = 'highMET_dRcut'
                        measure_efficiency(region)
                        bweight = get_weight(region, js)

                        plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                        plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                        if len(bjs) == 0 :

                            region = 'highMET_dRcut_0bjet'
                            measure_efficiency(region)
                            bweight = get_weight(region, js)

                            plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                            plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                        if len(bjs) > 0 :

                            region = 'highMET_dRcut_bjet'
                            measure_efficiency(region)
                            bweight = get_weight(region, js)

                            plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                            plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)


                if e.DeltaR(mu) > 0.4 and e.DeltaR(mu) < 0.6 :

                    if m.Pt() > event_cut['metcut'] :

                        region = 'highMET_invdRcut'
                        measure_efficiency(region)
                        bweight = get_weight(region, js)

                        plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                        plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                        if len(bjs) == 0 :

                            region = 'highMET_invdRcut_0bjet'
                            measure_efficiency(region)
                            bweight = get_weight(region, js)

                            plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                            plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                        if len(bjs) > 0 :

                            region = 'highMET_invdRcut_bjet'
                            measure_efficiency(region)
                            bweight = get_weight(region, js)

                            plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                            plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)


                    if m.Pt() < event_cut['metcut'] :
                            
                        region = 'lowMET_invdRcut'
                        measure_efficiency(region)
                        bweight = get_weight(region, js)

                        plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                        plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)

                        if len(bjs) == 0 :

                            region = 'lowMET_invdRcut_0bjet'
                            measure_efficiency(region)
                            bweight = get_weight(region, js)

                            plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                            plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)


                        if len(bjs) > 0 :

                            region = 'lowMET_invdRcut_bjet'
                            measure_efficiency(region)
                            bweight = get_weight(region, js)

                            plot_event_hist(region+"_noCorr", e, mu, j, j2, m, bj, bj2, js, bjs, genweight)
                            plot_event_hist(region, e, mu, j, j2, m, bj, bj2, js, bjs, genweight*bweight)



def plot_event_hist(region, l1, l2, j, j2, m, bj, bj2, js, bjs, weight):
    
    h[region+"_Count"].Fill(0, weight)

    h[region+"_Mass"].Fill((l1+l2).M(), weight)
    h[region+"_Lepton1Pt"].Fill(l1.Pt(), weight)
    h[region+"_Lepton2Pt"].Fill(l2.Pt(), weight)
    h[region+"_JetPt"].Fill(j.Pt(), weight)
    h[region+"_MetPt"].Fill(m.Pt(), weight)

    if len(js) > 1 :
    
        h[region+"_Jet2Pt"].Fill(j2.Pt(), weight)
        h[region+"_dRj2"].Fill(j2.DeltaR(l1+l2), weight)

    if len(bjs) > 0 :

        h[region+"_BJetPt"].Fill(bj.Pt(), weight)
        h[region+"_dRbj"].Fill(bj.DeltaR(l1+l2), weight)

        if len(bjs) > 1 :

            h[region+"_BJet2Pt"].Fill(bj2.Pt(), weight)
            h[region+"_dRbj2"].Fill(bj2.DeltaR(l1+l2), weight)

    h[region+"_Nj"].Fill(len(js), weight)
    h[region+"_Nbj"].Fill(len(bjs), weight)

    h[region+"_dRl"].Fill(l1.DeltaR(l2), weight)
    h[region+"_dRj"].Fill(j.DeltaR(l1+l2), weight)

    h[region+"_deepjet1"].Fill(js[0].deepjet, weight)
    if len(js) > 1 :
        h[region+"_deepjet2"].Fill(js[1].deepjet, weight)

    for i in js:
        h[region+"_deepjet"].Fill(i.deepjet, weight)


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

fchain.SetBranchAddress("Jets", ROOT.AddressOf(jets))
fchain.SetBranchAddress("Muons", ROOT.AddressOf(muons))
fchain.SetBranchAddress("Electrons", ROOT.AddressOf(electrons))


##-------------Start Event loop----------------

for iev in range(fchain.GetEntries()): # Be careful!!!                                                               

    fchain.GetEntry(iev)

    mets = fchain.GetBranch("Mets")
    met_pt = mets.GetLeaf('pt').GetValue()
    met_phi = mets.GetLeaf('phi').GetValue()

    h['hEvents'].Fill(0.5, 1)

    if isData == 1: genweight = 1

    if isData == 0:
        weight = fchain.GetBranch("lumiInfo")
        gweight = weight.GetLeaf("weight").GetValue()
        puweight = fchain.GetLeaf("puWeight").GetValue()
        genweight = gweight*puweight

    h['hEvents'].Fill(1.5, genweight)

    isMu = fchain.GetLeaf('isMu').GetValue()
    isIsoMu = fchain.GetLeaf('isIsoMu').GetValue()
    isMuonEG = fchain.GetLeaf('isMuonEG').GetValue()

    js = []
    mus = []
    es = []
    bjs = []
    
    if jets.size()>0:
        for i in range(jets.size()):
            jet = jets.at(i)
            if abs(jet.eta) < 2.5 :
                if jet.id >= 2:
                    js+=[jet]
                    if jet.deepjet > 0.7476:
                        bjs+=[jet]

    if muons.size()>0:
        for i in range(muons.size()):
            muon = muons.at(i)
            if abs(muon.eta) < 2.4 and muon.pt > 8.0:
                if muon.id >= 1:
                    if muon.iso < 0.25:
                        mus+=[muon]

    if electrons.size()>0:
        for i in range(electrons.size()):
            electron = electrons.at(i)
            if abs(electron.eta) < 2.5 and electron.pt > 12.0 :
                if electron.id >= 1 :
                    if electron.iso >= 1:
                        es+=[electron]


    js.sort(key=lambda x: x.pt, reverse=True)
    bjs.sort(key=lambda x: x.pt, reverse=True)
    mus.sort(key=lambda x: x.pt, reverse=True)
    es.sort(key=lambda x: x.pt, reverse=True)

    if len(es) > 0 and len(mus) > 0 and es[0].charge*mus[0].charge < 0 and len(js) > 0 :
        EMu_regions_efficiency()
        EMu_events()

out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))
