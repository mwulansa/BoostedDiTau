import ROOT, sys, os
import numpy as np
import time

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/GenParticleInfoDS.h"')

start_time = time.time()

inputFileListName=sys.argv[1]
inputFileList=inputFileListName

if len(sys.argv)>2:
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./output/"

outputTitle = "h_GenDY"

outputFileName = outputFileDir+outputTitle+"_"+inputFileListName.split("/")[-1].replace(".txt",".root")

out=ROOT.TFile.Open(outputFileName,'recreate')
print(outputFileName)

chain4 = ROOT.TChain('tcpGenNtuples/genTree')

pi = np.pi

h = {}

def define_event_histogram(region):

    h[region+"_Count"] = ROOT.TH1F (region+"_Count", region+"_Count ; Events ; Events ", 1, 0, 1)

    h[region+"_Mass"] = ROOT.TH1F (region+"_Mass", region+"_Mass ; M_{vis.} (GeV) ; Events ", 100, 0, 100)
    h[region+"_Lepton1Pt"] = ROOT.TH1F (region+"_Lepton1Pt", region+"_Lepton1Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_Lepton2Pt"] = ROOT.TH1F (region+"_Lepton2Pt", region+"_Lepton2Pt ; P_{T} (GeV) ; Events ", 500, 0, 500)
    h[region+"_genJetPt"] = ROOT.TH1F (region+"_genJetPt", region+"_genJetPt ; P_{T} (GeV) ; Events ", 1500, 0, 1500)

def define_general_histogram():

    h['hEvents'] = ROOT.TH1F ("NEvents", "Number of Events; ;N", 2, 0, 2)
    h['genJet_pt'] = ROOT.TH1F ("genJet_pt", "genJet Pt; Pt (GeV); N", 1500, 0, 1500)
    h['nmu'] = ROOT.TH1F ("hnmu", "Muons Pt; Pt (GeV); N", 500, 0, 500)
    h['nmu_hardProcess'] = ROOT.TH1F ("hnmu_hardProcess", "Muons Pt; Pt (GeV); N", 500, 0, 500)
    h['nmu_taudecay'] = ROOT.TH1F ("hnmu_taudecay", "Muons Pt; Pt (GeV); N", 500, 0, 500)
    h['pmu'] = ROOT.TH1F ("hpmu", "Muons Pt; Pt (GeV); N", 500, 0, 500)
    h['pmu_hardProcess'] = ROOT.TH1F ("hpmu_hardProcess", "Muons Pt; Pt (GeV); N", 500, 0, 500)
    h['pmu_taudecay'] = ROOT.TH1F ("hpmu_taudecay", "Muons Pt; Pt (GeV); N", 500, 0, 500)
    h['ne'] = ROOT.TH1F ("hne", "Electrons Pt; Pt (GeV); N", 500, 0, 500)
    h['ne_hardProcess'] = ROOT.TH1F ("hne_hardProcess", "Electrons Pt; Pt (GeV); N", 500, 0, 500)
    h['ne_taudecay'] = ROOT.TH1F ("hne_taudecay", "Electrons Pt; Pt (GeV); N", 500, 0, 500)
    h['pe'] = ROOT.TH1F ("hpe", "Electrons Pt; Pt (GeV); N", 500, 0, 500)
    h['pe_hardProcess'] = ROOT.TH1F ("hpe_hardProcess", "Electrons Pt; Pt (GeV); N", 500, 0, 500)
    h['pe_taudecay'] = ROOT.TH1F ("hpe_taudecay", "Electrons Pt; Pt (GeV); N", 500, 0, 500)
    h['ntau'] = ROOT.TH1F ("hntau", "Taus Pt; Pt (GeV); N", 500, 0, 500)
    h['ntau_hardProcess'] = ROOT.TH1F ("hntau_hardProcess", "Taus Pt; Pt (GeV); N", 500, 0, 500)
    h['ptau'] = ROOT.TH1F ("hptau", "Tauss Pt; Pt (GeV); N", 500, 0, 500)
    h['ptau_hardProcess'] = ROOT.TH1F ("hptau_hardProcess", "Taus Pt; Pt (GeV); N", 500, 0, 500)


hist_regions = ['Muons',
                'Muons_taudecay',
                'Muons_taudecay_dRcut',
                'Muons_isHardProcess',
                'Muons_isHardProcess_hardJet100',
                'Muons_isHardProcess_hardJet200',
                'Muons_isHardProcess_hardJet300',
                'Muons_isHardProcess_hardJet400',
                'Muons_isHardProcess_hardJet500',
                'Muons_isHardProcess_dRcut',
                'Muons_lowMass',
                'Muons_lowMass_isHardProcess',
                'Electrons_taudecay',
                'Electrons_taudecay_dRcut',
                'Electrons_isHardProcess',
                'Electrons_isHardProcess_hardJet100',
                'Electrons_isHardProcess_hardJet200',
                'Electrons_isHardProcess_hardJet300',
                'Electrons_isHardProcess_hardJet400',
                'Electrons_isHardProcess_hardJet500',
                'Electrons_isHardProcess_dRcut',
                'Electrons_lowMass_isHardProcess',
                'Taus_isHardProcess',
                'Taus_isHardProcess_hardJet100',
                'Taus_isHardProcess_hardJet200',
                'Taus_isHardProcess_hardJet300',
                'Taus_isHardProcess_hardJet400',
                'Taus_isHardProcess_hardJet500',
                'Taus_isHardProcess_dRcut',
                'Taus_lowMass_isHardProcess',
                'Leptons_isHardProcess',
                'Leptons_isHardProcess_hardJet100',
                'Leptons_isHardProcess_hardJet200',
                'Leptons_isHardProcess_hardJet300',
                'Leptons_isHardProcess_hardJet400',
                'Leptons_isHardProcess_hardJet500',
                'Leptons_isHardProcess_dRcut',
                'Leptons_lowMass_isHardProcess',
]

for r in hist_regions:
    define_event_histogram(r)

define_general_histogram()


for key in h.keys():
    h[key].Sumw2()

def plot_hists(region,l1,l2,j):
    
    h[region+"_Count"].Fill(0, genweight)

    h[region+"_Mass"].Fill((l1+l2).M(), genweight)
    h[region+"_Lepton1Pt"].Fill(l1.Pt(), genweight)
    h[region+"_Lepton2Pt"].Fill(l2.Pt(), genweight)
    h[region+"_genJetPt"].Fill(j, genweight)


def get_TLorentzVector(l1, l2):

    vl1 = ROOT.TLorentzVector()
    vl1.SetPtEtaPhiM(l1.pt, l1.eta, l1.phi, l1.mass)

    vl2 = ROOT.TLorentzVector()
    vl2.SetPtEtaPhiM(l2.pt, l2.eta, l2.phi, l2.mass)

    return vl1, vl2


inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    inputFileName=inputFileName.replace("\n","")
    print(inputFileName.replace("\n",""))

    chain4.Add(inputFileName)


genParticle = ROOT.GenParticleInfoDS()
chain4.SetBranchAddress("GenParticleInfo", ROOT.AddressOf(genParticle))

for iev in range(chain4.GetEntries()): # Be careful!!!                                                                                                   
    chain4.GetEntry(iev)

#    print("NEvent", iev)

    genweight = chain4.GetLeaf('genWeight').GetValue()
    eventno = chain4.GetLeaf('event').GetValue()

    genJetInfo = chain4.GetBranch("genJetInfo")

    genJet_pt = genJetInfo.GetLeaf('pt').GetValue()
    genJet_eta = genJetInfo.GetLeaf('eta').GetValue()
    genJet_phi = genJetInfo.GetLeaf('phi').GetValue()
    genJet_mass = genJetInfo.GetLeaf('mass').GetValue()

    h['hEvents'].Fill(0.5, 1)
    h['hEvents'].Fill(1.5, genweight)

    pmu = []
    nmu = []

    pmu_hardProcess = []
    nmu_hardProcess = []
    pmu_taudecay = []
    nmu_taudecay = []

    pe = []
    ne = []

    pe_hardProcess = []
    ne_hardProcess = []
    pe_taudecay = []
    ne_taudecay = []

    ptau = []
    ntau = []

    ptau_hardProcess = []
    ntau_hardProcess = []

    h['genJet_pt'].Fill(genJet_pt)

    if genParticle.size()>0:
        for i in range(genParticle.size()):
            gen = genParticle.at(i)

            if gen.pdgid == 13 : 
                nmu += [gen]
                h['nmu'].Fill(gen.pt)
                if gen.ishardprocess:
                    nmu_hardProcess += [gen]
                    h['nmu_hardProcess'].Fill(gen.pt)
                if gen.isdirecthardprocesstaudecayproductfinalstate:
                    nmu_taudecay += [gen]
                    h['nmu_taudecay'].Fill(gen.pt)

            if gen.pdgid == -13 :
                pmu += [gen]
                h['pmu'].Fill(gen.pt)
                if gen.ishardprocess:
                    pmu_hardProcess += [gen]
                    h['pmu_hardProcess'].Fill(gen.pt)
                if gen.isdirecthardprocesstaudecayproductfinalstate:
                    pmu_taudecay += [gen]
                    h['pmu_taudecay'].Fill(gen.pt)

            if gen.pdgid == 11 :
                ne += [gen]
                h['ne'].Fill(gen.pt)
                if gen.ishardprocess:
                    ne_hardProcess += [gen]
                    h['ne_hardProcess'].Fill(gen.pt)
                if gen.isdirecthardprocesstaudecayproductfinalstate:
                    ne_taudecay += [gen]
                    h['ne_taudecay'].Fill(gen.pt)

            if gen.pdgid == -11 :
                pe += [gen]
                h['pe'].Fill(gen.pt)
                if gen.ishardprocess:
                    pe_hardProcess += [gen]
                    h['pe_hardProcess'].Fill(gen.pt)
                if gen.isdirecthardprocesstaudecayproductfinalstate:
                    pe_taudecay += [gen]
                    h['pe_taudecay'].Fill(gen.pt)

            if gen.pdgid == 15 :
                ntau += [gen]
                h['ntau'].Fill(gen.pt)
                if gen.ishardprocess:
                    ntau_hardProcess += [gen]
                    h['ntau_hardProcess'].Fill(gen.pt)

            if gen.pdgid == -15 :
                ptau += [gen]
                h['ptau'].Fill(gen.pt)
                if gen.ishardprocess:
                    ptau_hardProcess += [gen]
                    h['ptau_hardProcess'].Fill(gen.pt)


    pmu.sort(key=lambda x: x.pt, reverse=True)
    nmu.sort(key=lambda x: x.pt, reverse=True)
    pmu_hardProcess.sort(key=lambda x: x.pt, reverse=True)
    nmu_hardProcess.sort(key=lambda x: x.pt, reverse=True)

    pe.sort(key=lambda x: x.pt, reverse=True)
    ne.sort(key=lambda x: x.pt, reverse=True)
    pe_hardProcess.sort(key=lambda x: x.pt, reverse=True)
    ne_hardProcess.sort(key=lambda x: x.pt, reverse=True)

    ptau.sort(key=lambda x: x.pt, reverse=True)
    ntau.sort(key=lambda x: x.pt, reverse=True)
    ptau_hardProcess.sort(key=lambda x: x.pt, reverse=True)
    ntau_hardProcess.sort(key=lambda x: x.pt, reverse=True)


    if len(pmu_taudecay) == 1 and len(nmu_taudecay) == 1:
        
        mu1, mu2 = get_TLorentzVector(pmu_taudecay[0], nmu_taudecay[0])
        
        j = ROOT.TLorentzVector()
        j.SetPtEtaPhiM(genJet_pt, genJet_eta, genJet_phi, genJet_mass)

        plot_hists("Muons_taudecay", mu1, mu2, genJet_pt)

        if mu1.DeltaR(mu2) < 0.4 and j.DeltaR(mu1) > 0.8 and j.DeltaR(mu2) > 0.8 :
            
            plot_hists("Muons_taudecay_dRcut", mu1, mu2, genJet_pt)

    if len(pe_taudecay) == 1 and len(ne_taudecay) == 1:

        e1, e2 = get_TLorentzVector(pe_taudecay[0], ne_taudecay[0])

        j = ROOT.TLorentzVector()
        j.SetPtEtaPhiM(genJet_pt, genJet_eta, genJet_phi, genJet_mass)

        plot_hists("Electrons_taudecay", e1, e2, genJet_pt)

        if e1.DeltaR(e2) < 0.4 and j.DeltaR(e1) > 0.8 and j.DeltaR(e2) > 0.8 :

            plot_hists("Electrons_taudecay_dRcut", e1, e2, genJet_pt)


    if len(pmu) == 1 and len(nmu) == 1:

        mu1, mu2 = get_TLorentzVector(pmu[0], nmu[0])

        plot_hists("Muons", mu1, mu2, genJet_pt)

        if (mu1+mu2).M() < 4:

            plot_hists("Muons_lowMass", mu1, mu2, genJet_pt)

    if len(pmu_hardProcess) == 1 and len(nmu_hardProcess) == 1 :

        mu1, mu2 = get_TLorentzVector(pmu_hardProcess[0], nmu_hardProcess[0])

        j = ROOT.TLorentzVector()
        j.SetPtEtaPhiM(genJet_pt, genJet_eta, genJet_phi, genJet_mass)

        plot_hists("Muons_isHardProcess", mu1, mu2, genJet_pt)
        plot_hists("Leptons_isHardProcess", mu1, mu2, genJet_pt)

        if mu1.DeltaR(mu2) < 0.4 and j.DeltaR(mu1) > 0.8 and j.DeltaR(mu2) > 0.8 :

            plot_hists("Muons_isHardProcess_dRcut", mu1, mu2, genJet_pt)
            plot_hists("Leptons_isHardProcess_dRcut", mu1, mu2, genJet_pt)

        if genJet_pt > 100.0 :
            plot_hists("Muons_isHardProcess_hardJet100", mu1, mu2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet100", mu1, mu2, genJet_pt)

        if genJet_pt > 200.0 :
            plot_hists("Muons_isHardProcess_hardJet200", mu1, mu2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet200", mu1, mu2, genJet_pt)

        if genJet_pt > 300.0 :
            plot_hists("Muons_isHardProcess_hardJet300", mu1, mu2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet300", mu1, mu2, genJet_pt)

        if genJet_pt > 400.0 :
            plot_hists("Muons_isHardProcess_hardJet400", mu1, mu2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet400", mu1, mu2, genJet_pt)

        if genJet_pt > 500.0 : 
            plot_hists("Muons_isHardProcess_hardJet500", mu1, mu2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet500", mu1, mu2, genJet_pt)

        if (mu1+mu2).M() < 4:

            plot_hists("Muons_lowMass_isHardProcess", mu1, mu2, genJet_pt)
            plot_hists("Leptons_lowMass_isHardProcess", mu1, mu2, genJet_pt)

    if len(pe_hardProcess) == 1 and len(ne_hardProcess) == 1 :

        e1, e2 = get_TLorentzVector(pe_hardProcess[0], ne_hardProcess[0])

        j = ROOT.TLorentzVector()
        j.SetPtEtaPhiM(genJet_pt, genJet_eta, genJet_phi, genJet_mass)

        plot_hists("Electrons_isHardProcess", e1, e2, genJet_pt)
        plot_hists("Leptons_isHardProcess", e1, e2, genJet_pt)

        if e1.DeltaR(e2) < 0.4 and j.DeltaR(e1) > 0.8 and j.DeltaR(e2) > 0.8 :

            plot_hists("Electrons_isHardProcess_dRcut", e1, e2, genJet_pt)
            plot_hists("Leptons_isHardProcess_dRcut", e1, e2, genJet_pt)

        if genJet_pt > 100.0 :
            plot_hists("Electrons_isHardProcess_hardJet100", e1, e2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet100", e1, e2, genJet_pt)

        if genJet_pt > 200.0 :
            plot_hists("Electrons_isHardProcess_hardJet200", e1, e2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet200", e1, e2, genJet_pt)

        if genJet_pt > 300.0 :
            plot_hists("Electrons_isHardProcess_hardJet300", e1, e2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet300", e1, e2, genJet_pt)

        if genJet_pt > 400.0 :
            plot_hists("Electrons_isHardProcess_hardJet400", e1, e2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet400", e1, e2, genJet_pt)

        if genJet_pt > 500.0 :
            plot_hists("Electrons_isHardProcess_hardJet500", e1, e2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet500", e1, e2, genJet_pt)

        if (e1+e2).M() < 4:

            plot_hists("Electrons_lowMass_isHardProcess", e1, e2, genJet_pt)
            plot_hists("Leptons_lowMass_isHardProcess", e1, e2, genJet_pt)


    if len(ptau_hardProcess) == 1 and len(ntau_hardProcess) == 1 :

        tau1, tau2 = get_TLorentzVector(ptau_hardProcess[0], ntau_hardProcess[0])

        j = ROOT.TLorentzVector()
        j.SetPtEtaPhiM(genJet_pt, genJet_eta, genJet_phi, genJet_mass)

        plot_hists("Taus_isHardProcess", tau1, tau2, genJet_pt)
        plot_hists("Leptons_isHardProcess", tau1, tau2, genJet_pt)

        if tau1.DeltaR(tau2) < 0.4 and j.DeltaR(tau1) > 0.8 and j.DeltaR(tau2) > 0.8 :

            plot_hists("Taus_isHardProcess_dRcut", tau1, tau2, genJet_pt)
            plot_hists("Leptons_isHardProcess_dRcut", tau1, tau2, genJet_pt)

        if genJet_pt > 100.0 :
            plot_hists("Taus_isHardProcess_hardJet100", tau1, tau2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet100", tau1, tau2, genJet_pt)

        if genJet_pt > 200.0 :
            plot_hists("Taus_isHardProcess_hardJet200", tau1, tau2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet200", tau1, tau2, genJet_pt)

        if genJet_pt > 300.0 :
            plot_hists("Taus_isHardProcess_hardJet300", tau1, tau2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet300", tau1, tau2, genJet_pt)

        if genJet_pt > 400.0 :
            plot_hists("Taus_isHardProcess_hardJet400", tau1, tau2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet400", tau1, tau2, genJet_pt)

        if genJet_pt > 500.0 :
            plot_hists("Taus_isHardProcess_hardJet500", tau1, tau2, genJet_pt)
            plot_hists("Leptons_isHardProcess_hardJet500", tau1, tau2, genJet_pt)


        if (tau1+tau2).M() < 4:

            plot_hists("Taus_lowMass_isHardProcess", tau1, tau2, genJet_pt)
            plot_hists("Leptons_lowMass_isHardProcess", tau1, tau2, genJet_pt)

        

out.cd()

for key in h.keys():
    h[key].Write()

out.Close()

print("--- %s seconds ---" % (time.time() - start_time))
