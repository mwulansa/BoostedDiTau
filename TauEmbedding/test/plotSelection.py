import ROOT,sys,gc
from DataFormats.FWLite import Events, Handle

from looseElectron import *

#events=Events('condorCfg/embed_YMuMu_pth400_9.root')
#events=Events('outMuonSelection.root')
#events=Events('root://xrootd.unl.edu//store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0AB088EE-EA8A-E811-8636-0CC47A4C8E96.root')

out=ROOT.TFile('h_selection.root','recreate')


prefix = "root://cmseos.fnal.gov//eos/uscms/store/group/lpcsusyhiggs/Jingyu/events/UpsilonTauTau/UL2017/"
#prefix = "./condorCfg/"

h={}

hname = "mmumu1"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)

hname = "mmumu2"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)

hname = "mmumu3"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)

hname = "mmumu4"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)




handleGenParticles = Handle ('vector<reco::GenParticle>')
labelGenParticles = ('prunedGenParticles', '', 'PAT')

handleTaus = Handle ('vector<pat::Tau>')
#labelTaus = ('slimmedTausBoosted', '', 'PAT')
labelTaus = ("selectedPatTausBoostedReMINIAOD","","ReMINIAOD")


handleMuons = Handle ('vector<pat::Muon>')
labelMuons = ('slimmedMuons','','PAT')

handleElectrons = Handle ('vector<pat::Electron>')
labelElectrons = ('slimmedElectrons','','PAT')

handleCands = Handle ('vector<pat::PackedCandidate>')
labelCands = ('packedPFCandidates','','PAT')

handleMet = Handle ('vector<pat::MET>')
labelMet = ('slimmedMETs','','PAT')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator', '', 'GEN' )

handlePVInfo = Handle ('vector<reco::Vertex>')
labelPVInfo = ('offlineSlimmedPrimaryVertices', '', 'PAT')

handlePFJets = Handle ('vector<pat::Jet>')
labelPFJets = ('slimmedJets', '', 'PAT')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll', '', 'RECO')

handleCA8Jets = Handle ('vector<reco::PFJet>')
labelCA8Jets = ("ca8PFJetsCHSprunedForBoostedTausPATBoostedReMINIAOD","subJetsForSeedingBoostedTausPAT","ReMINIAOD")

handleZmumuCandidate = Handle ('vector<reco::CompositeCandidate>')
labelZmumuCandidate = ('ZmumuCandidates', '', 'SELECT')

nfail = 0
for i in range(1, 501):
#for i in range(1, 2):
    try:
        events = Events(prefix+"YMuMu_pth400_"+str(i)+".root")
        print(prefix+"YMuMu_pth400_"+str(i)+".root")
        nevt=0
        for event in events:
            nevt+=1
        
            event.getByLabel(labelMuons, handleMuons)
            muons=handleMuons.product()

            #event.getByLabel(labelZmumuCandidate, handleZmumuCandidate)
            #ZmumuCand = handleZmumuCandidate.product()

            genweight = 1
           
            selected_muons1=[]
            selected_muons2=[]
            selected_muons3=[]
            selected_muons4=[]

            selected_muons=[]
            for muon in muons:
                if not muon.isLooseMuon(): continue
                if muon.pt()<3 or muon.eta()>2.4: continue
                selected_muons1+=[muon]
                selected_muons2+=[muon]
                selected_muons3+=[muon]
                if (muon.pfIsolationR04().sumChargedHadronPt+max(0,muon.pfIsolationR04().sumPhotonEt+muon.pfIsolationR04().sumNeutralHadronEt-0.5*muon.pfIsolationR04().sumPUPt))/muon.pt()<0.4:
                    selected_muons4+=[muon]

            ## smallest dr
            if len(selected_muons1) > 1:
                drMin = 9999
                pt1 = -1
                pt2 = -1
                mass = -1
                for i, tau1 in enumerate(selected_muons1):
                    for j, tau2 in enumerate(selected_muons1):
                        if i > j or tau1.charge() * tau2.charge() > 0: continue
                        mu1=ROOT.TLorentzVector()
                        mu1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), tau1.mass()) 
            
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), tau2.mass())

                        if (mu1+mu2).M() < 3.6 or (mu1+mu2).M() > 20: continue
                        if mu1.DeltaR(mu2) < drMin:
                            mass = (mu1+mu2).M()
                            pt1 = mu1.Pt()
                            pt2 = mu2.Pt()
                            drMin = mu1.DeltaR(mu2)

                #if drMin > 0.8: continue
                h["mmumu1"].Fill(mass, genweight)

            ## smallest dr with minimum dr cut
            if len(selected_muons2) > 1:
                drMin = 9999
                pt1 = -1
                pt2 = -1
                mass = -1
                for i, tau1 in enumerate(selected_muons2):
                    for j, tau2 in enumerate(selected_muons2):
                        if i > j or tau1.charge() * tau2.charge() > 0: continue
                        mu1=ROOT.TLorentzVector()
                        mu1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), tau1.mass()) 
            
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), tau2.mass())

                        if (mu1+mu2).M() < 3.6 or (mu1+mu2).M() > 20: continue
                        if mu1.DeltaR(mu2) < drMin and mu1.DeltaR(mu2) > 0.05:
                            mass = (mu1+mu2).M()
                            pt1 = mu1.Pt()
                            pt2 = mu2.Pt()
                            drMin = mu1.DeltaR(mu2)

                #if drMin > 0.8: continue
                h["mmumu2"].Fill(mass, genweight)

            ## largest mass
            if len(selected_muons3) > 1:
                drMin = 9999
                pt1 = -1
                pt2 = -1
                mass = -1
                for i, tau1 in enumerate(selected_muons3):
                    for j, tau2 in enumerate(selected_muons3):
                        if i > j or tau1.charge() * tau2.charge() > 0: continue
                        mu1=ROOT.TLorentzVector()
                        mu1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), tau1.mass()) 
            
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), tau2.mass())


                        if (mu1+mu2).M() < 3.6 or (mu1+mu2).M() > 20: continue
                        if (mu1+mu2).M() > mass:
                            mass = (mu1+mu2).M()
                            pt1 = mu1.Pt()
                            pt2 = mu2.Pt()
                            drMin = mu1.DeltaR(mu2)

                if mass < 3.6: continue
                h["mmumu3"].Fill(mass, genweight)

            ## smallest dr isolated
            if len(selected_muons4) > 1:
                drMin = 9999
                pt1 = -1
                pt2 = -1
                mass = -1
                for i, tau1 in enumerate(selected_muons4):
                    for j, tau2 in enumerate(selected_muons4):
                        if i > j or tau1.charge() * tau2.charge() > 0: continue
                        mu1=ROOT.TLorentzVector()
                        mu1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), tau1.mass()) 
            
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), tau2.mass())


                        if (mu1+mu2).M() < 3.6 or (mu1+mu2).M() > 20: continue
                        if mu1.DeltaR(mu2) < drMin:
                            mass = (mu1+mu2).M()
                            pt1 = mu1.Pt()
                            pt2 = mu2.Pt()
                            drMin = mu1.DeltaR(mu2)

                if drMin > 0.8: continue
                h["mmumu4"].Fill(mass, genweight)
        
        del events
        gc.collect()
    except:
        nfail += 1
        print("Error...")

print("Nfail:", nfail)
out.cd()
for hist in sorted(h.keys()):
    h[hist].Write()
