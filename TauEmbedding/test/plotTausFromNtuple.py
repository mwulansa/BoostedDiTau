import ROOT

ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/TauInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../../MiniAODSkimmer/interface/GenParticleInfoDS.h"')

inputFile = 'TCPNtuple_Backgrounds_9.root'
out = ROOT.TFile('h_ntuple.root','recreate')

taus = ROOT.TauInfoDS()
genparticles = ROOT.GenParticleInfoDS()

fchain = ROOT.TChain('tcpNtuples/analysisTree')
fchain.Add(inputFile)

fchainGen = ROOT.TChain('tcpGenNtuples/genTree')
fchainGen.Add(inputFile)

fchain.AddFriend(fchainGen)

fchain.SetBranchAddress("TausUnCleaned", ROOT.AddressOf(taus))
fchain.SetBranchAddress("GenParticleInfo", ROOT.AddressOf(genparticles))

h1 = ROOT.TH1F ("genTauPt", "", 500, 0, 500)
h2 = ROOT.TH1F ("mtautau", "", 50, 0, 500)
h3 = ROOT.TH1F ("genTauDR", "", 70, 0, 7)
h4 = ROOT.TH1F ("tauPt", "", 500, 0, 500)

nthth = 0
for iev in range(fchain.GetEntries()): # Be careful!!!
    fchain.GetEntry(iev)
    if nthth == 1000: break
    gen_e = []
    gen_mu = []
    gen_tau = []
    for genparticle in genparticles:
        if abs(genparticle.pdgid) == 11: gen_e+=[genparticle]
        if abs(genparticle.pdgid) == 13: gen_mu+=[genparticle]
        if abs(genparticle.pdgid) == 15: gen_tau+=[genparticle]
        
    #if len(gen_e) == 0 and len(gen_mu) == 0 and len(gen_tau) > 0:
    if len(gen_tau) > 0:
        nthth += 1
        #print(iev, nthth, taus.size())
        selected_taus = []
        for tau in taus:
            if tau.deepid >=2: selected_taus+=[tau]

        if len(selected_taus) > 1:
            mu1=ROOT.TLorentzVector()
            mu1.SetPtEtaPhiM(selected_taus[0].pt, selected_taus[0].eta, selected_taus[0].phi, selected_taus[0].mass) 
    
            mu2=ROOT.TLorentzVector()
            mu2.SetPtEtaPhiM(selected_taus[1].pt, selected_taus[1].eta, selected_taus[1].phi, selected_taus[1].mass)

            h2.Fill((mu1+mu2).M())
            
            if mu1.Pt() > mu2.Pt(): h4.Fill(mu1.Pt())
            else: h4.Fill(mu2.Pt())

        if len(gen_tau) > 1:
            if gen_tau[0].pt > gen_tau[1].pt: h1.Fill(gen_tau[0].pt)
            else: h1.Fill(gen_tau[1].pt)

            t1=ROOT.TLorentzVector()
            t1.SetPtEtaPhiM(gen_tau[0].pt, gen_tau[0].eta, gen_tau[0].phi, gen_tau[0].mass) 
    
            t2=ROOT.TLorentzVector()
            t2.SetPtEtaPhiM(gen_tau[1].pt, gen_tau[1].eta, gen_tau[1].phi, gen_tau[1].mass)
            
            h3.Fill(t1.DeltaR(t2))

               
out.cd()
h1.Write()
h2.Write()
h3.Write()
h4.Write()
