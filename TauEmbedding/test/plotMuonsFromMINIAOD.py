import ROOT,sys,gc
from DataFormats.FWLite import Events, Handle

from looseElectron import *

#events=Events('condorCfg/embed_YMuMu_pth400_9.root')
#events=Events('outMuonSelection.root')
#events=Events('root://xrootd.unl.edu//store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0AB088EE-EA8A-E811-8636-0CC47A4C8E96.root')

out=ROOT.TFile('h_Ymumu.root','recreate')


prefix = "root://cmseos.fnal.gov//eos/uscms/store/group/lpcsusyhiggs/Jingyu/events/UpsilonTauTau/UL2017/"
#prefix = "./condorCfg/"

h={}
hname = "genTau1Pt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "genTau2Pt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)

hname = "genTaumu1Pt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "genTaumu2Pt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)

hname = "genTauDR"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = "genTauM"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)

hname = "pfJetPt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 1000)
hname = "ca8JetPt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 1000)

hname = "memu"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)
hname = "mmumu"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)

hname = "mmumu_genM"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 50)

hname = "mtautau"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)
hname = "mtautau_deep"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)

hname = 'dRtautau'
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = 'dRtautau_deep'
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = 'dRemu'
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = 'dRmumu'
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = "Weights"
h[hname] = ROOT.TH1F (hname, "", 20, 0, 2)
hname = "MET"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)

hname = "tau1Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau1Pt_deep"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2Pt_deep"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)

hname = "ePt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "muPt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)

hname = "mu1Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "mu2Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)


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

nfail = 0
#for i in range(1, 501):
for i in range(1, 2):
    try:
        events = Events(prefix+"YMuMu_pth400_"+str(i)+".root")
        print(prefix+"YMuMu_pth400_"+str(i)+".root")
        nevt=0
        for event in events:
            nevt+=1
            #print('-----', nevt)
        
            event.getByLabel(labelPVInfo, handlePVInfo)
            pv = handlePVInfo.product()[0]
            
            #event.getByLabel(labelTaus, handleTaus)
            #taus=handleTaus.product()
        
            event.getByLabel(labelMuons, handleMuons)
            muons=handleMuons.product()

            event.getByLabel(labelElectrons, handleElectrons)
            electrons=handleElectrons.product()
    
            event.getByLabel(labelCands, handleCands)
            cands=handleCands.product()
        
            event.getByLabel(labelMet, handleMet)
            met=handleMet.product()
        
            event.getByLabel(labelGenInfo, handleGenInfo)
            geninfo=handleGenInfo.product()
            genweight=geninfo.weight()
        
            event.getByLabel(labelGenParticles, handleGenParticles)
            genParticles=handleGenParticles.product()
        
            event.getByLabel(labelPFJets, handlePFJets)
            pfJets=handlePFJets.product()

            #event.getByLabel(labelCA8Jets, handleCA8Jets)
            #ca8Jets=handleCA8Jets.product()

            event.getByLabel(labelRho, handleRho)
            rho=handleRho.product()[0]
        
            if pfJets.size() > 0:
                h["pfJetPt"].Fill(pfJets[0].pt())

            #if ca8Jets.size() > 0:
            #    h["ca8JetPt"].Fill(ca8Jets[0].pt())
        
    
            h['MET'].Fill(met[0].pt())
            h['Weights'].Fill(genweight)
        
            selected_gentaus = []
            #print(genParticles.size())
            selected_gentaumus = []
            selected_ftaus=[]
            for particle in genParticles:
                #if particle.isHardProcess():
                #    print(particle.pdgId())
                if abs(particle.pdgId()) == 13 and (particle.mother().pdgId()==100553 or particle.mother().pdgId()==200553 or particle.mother().pdgId()==553):
                    #print("GenTau")
                    selected_gentaus += [particle]

                if abs(particle.pdgId()) == 13 and particle.status()==1 and abs(particle.mother().pdgId())==15:
                    selected_gentaumus += [particle]
                if particle.status()==1 and abs(particle.mother().pdgId())==15:
                    selected_ftaus +=[particle.pdgId()]

            #print(selected_ftaus)

            if len(selected_gentaumus) == 1:
                h["genTaumu1Pt"].Fill(selected_gentaumus[0].pt())

            if len(selected_gentaumus) == 2:
                if selected_gentaumus[0].pt() > selected_gentaumus[1].pt():
                    h["genTaumu1Pt"].Fill(selected_gentaumus[0].pt())
                    h["genTaumu2Pt"].Fill(selected_gentaumus[1].pt())
                else:
                    h["genTaumu1Pt"].Fill(selected_gentaumus[1].pt())
                    h["genTaumu2Pt"].Fill(selected_gentaumus[0].pt())

            genTauM = -1
            if len(selected_gentaus) > 1:
                if selected_gentaus[0].pt() > selected_gentaus[1].pt():
                    h["genTau1Pt"].Fill(selected_gentaus[0].pt(), genweight)
                    h["genTau2Pt"].Fill(selected_gentaus[1].pt(), genweight)
                else:
                    h["genTau2Pt"].Fill(selected_gentaus[0].pt(), genweight)
                    h["genTau1Pt"].Fill(selected_gentaus[1].pt(), genweight)
        
                t1=ROOT.TLorentzVector()
                t1.SetPtEtaPhiM(selected_gentaus[0].pt(), selected_gentaus[0].eta(), selected_gentaus[0].phi(), selected_gentaus[0].mass()) 
            
                t2=ROOT.TLorentzVector()
                t2.SetPtEtaPhiM(selected_gentaus[1].pt(), selected_gentaus[1].eta(), selected_gentaus[1].phi(), selected_gentaus[1].mass())
                    
                h["genTauDR"].Fill(t1.DeltaR(t2), genweight)
                h["genTauM"].Fill((t1+t2).M(), genweight)

                genTauM = (t1+t2).M()
        
            selected_muons=[]
            for muon in muons:
                if not muon.isLooseMuon(): continue
                if muon.pt()<3 or muon.eta()>2.4: continue
                #if (muon.pfIsolationR04().sumChargedHadronPt+max(0,muon.pfIsolationR04().sumPhotonEt+muon.pfIsolationR04().sumNeutralHadronEt-0.5*muon.pfIsolationR04().sumPUPt))/muon.pt()<0.4:
                selected_muons+=[muon]

            selected_electrons=[]
            for electron in electrons:
                if electron.pt()<7: continue 
                if abs(electron.eta())>2.5: continue
                E_c = electron.superCluster().energy()
                if electron.isEB():
                    if electron.full5x5_sigmaIetaIeta()<0.0112 \
                       and GsfEleEInverseMinusPInverse(electron)<0.193 \
                       and abs(dEtaInSeed(electron))<0.00377     \
                       and GsfEleMissingHitsCut(electron)<=1 \
                       and electron.passConversionVeto() \
                       and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.0884 \
                       and electron.hadronicOverEm()< (0.05 + 1.16/E_c + 0.0324*rho/E_c): 
                        selected_electrons+=[electron]

                if electron.isEE():
                    if electron.full5x5_sigmaIetaIeta() < 0.0425 \
                       and GsfEleEInverseMinusPInverse(electron) < 0.111 \
                       and abs(dEtaInSeed(electron)) < 0.00674 \
                       and GsfEleMissingHitsCut(electron)<=1 \
                       and electron.passConversionVeto() \
                       and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.169 \
                       and electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c):
                        selected_electrons+=[electron]

        

            if len(selected_muons) > 0 and len(selected_electrons) > 0:
                drMin = 9999
                pt1 = -1
                pt2 = -1
                mass = -1
                for i, tau1 in enumerate(selected_muons):
                    for j, tau2 in enumerate(selected_electrons):
                        if tau1.charge() * tau2.charge() > 0: continue
                        mu1=ROOT.TLorentzVector()
                        mu1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), tau1.mass()) 
            
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), tau2.mass())

                        if mu1.DeltaR(mu2) < drMin:
                            mass = (mu1+mu2).M()
                            pt1 = mu1.Pt()
                            pt2 = mu2.Pt()
                            drMin = mu1.DeltaR(mu2)

                if drMin > 0.8: continue
                h["memu"].Fill(mass, genweight)
                h["dRemu"].Fill(drMin, genweight)
                h["ePt"].Fill(pt2, genweight)
                h["muPt"].Fill(pt1, genweight)

            if len(selected_muons) > 1:
                drMin = 9999
                pt1 = -1
                pt2 = -1
                mass = -1
                for i, tau1 in enumerate(selected_muons):
                    for j, tau2 in enumerate(selected_muons):
                        if i > j or tau1.charge() * tau2.charge() > 0: continue
                        mu1=ROOT.TLorentzVector()
                        mu1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), tau1.mass()) 
            
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), tau2.mass())

                        if mu1.DeltaR(mu2) < drMin:
                            mass = (mu1+mu2).M()
                            pt1 = mu1.Pt()
                            pt2 = mu2.Pt()
                            drMin = mu1.DeltaR(mu2)

                if drMin > 0.8: continue
                h["mmumu"].Fill(mass, genweight)
                h["mmumu_genM"].Fill(genTauM, genweight)
                h["dRmumu"].Fill(drMin, genweight)
                if pt2 > pt1:
                    h["mu1Pt"].Fill(pt2, genweight)
                    h["mu2Pt"].Fill(pt1, genweight)
                else:
                    h["mu1Pt"].Fill(pt1, genweight)
                    h["mu2Pt"].Fill(pt2, genweight)
                
        
        del events
        gc.collect()
    except:
        nfail += 1
        print("Error...")

print("Nfail:", nfail)
out.cd()
for hist in sorted(h.keys()):
    h[hist].Write()
