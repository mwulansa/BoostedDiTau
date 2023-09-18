import ROOT,sys,gc, array
from DataFormats.FWLite import Events, Handle

#events=Events('condorCfg/embed_YMuMu_pth400_9.root')
#events=Events('outMuonSelection.root')
#events=Events('root://xrootd.unl.edu//store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0AB088EE-EA8A-E811-8636-0CC47A4C8E96.root')

out=ROOT.TFile('h_muon_embedded_boosted.root','recreate')


prefix = "root://cmseos.fnal.gov//eos/uscms/store/group/lpcsusyhiggs/Jingyu/embedding/YMuMu/UL2017/"
#prefix = "./condorCfg/"

h={}
hname = "genTauPt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "mmumu"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)
hname = "mu1Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "mu2Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "dRmumu"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = "mtautau"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)
hname = "mtautau_nocut"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)
hname = "mtautauEmbed"
h[hname ]= ROOT.TH1F (hname, "", 500, 0, 20)
hname = "mtautauEmbed_deep"
h[hname ]= ROOT.TH1F (hname, "", 500, 0, 20)

hname = "dRtautauEmbed"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = "dRtautauEmbed_deep"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = "Weights"
h[hname] = ROOT.TH1F (hname, "", 20, 0, 2)
hname = "MET"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "METEmbed"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "genTauDR"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = "genTauM"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)

hname = "tau1Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau1PtEmbed"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2PtEmbed"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)

hname = "tau1Pt_nocut"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2Pt_nocut"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)


hname = "tau1PtEmbed_deep"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2PtEmbed_deep"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)

hname = "pvxyEmbed"
h[hname] = ROOT.TH2F (hname, "", 80, -2, 2, 80, -2, 2)
hname = "pvzEmbed"
h[hname] = ROOT.TH1F (hname, "", 80, -20, 20)

hname = "pvxy"
h[hname] = ROOT.TH2F (hname, "", 80, -2, 2, 80, -2, 2)
hname = "pvz"
h[hname] = ROOT.TH1F (hname, "", 80, -20, 20)

hname = "pvxyGen"
h[hname] = ROOT.TH2F (hname, "", 80, -2, 2, 80, -2, 2)
hname = "pvzGen"
h[hname] = ROOT.TH1F (hname, "", 80, -20, 20)

hname = "nVtx"
h[hname] = ROOT.TH1F (hname, "", 5, 0, 5)

handleGenParticles = Handle ('vector<reco::GenParticle>')
labelGenParticles = ('prunedGenParticles', '', 'SIMembedding')

handleTaus = Handle ('vector<pat::Muon>')
#labelTaus = ('slimmedTausNewID','','SELECT')
labelTaus = ('slimmedMuons', '', 'SIMembedding')

handleMuons = Handle ('vector<pat::Muon>')
labelMuons = ('slimmedMuons','','PAT')

handleTausEmbed = Handle ('vector<pat::Muon>')
labelTausEmbed = ('embeddedPFCandidates','slimmedMuonsEmbedded','Embed')

handleCands = Handle ('vector<pat::PackedCandidate>')
labelCands = ('packedPFCandidates','','PAT')

handleCandsEmbed = Handle ('vector<pat::PackedCandidate>')
labelCandsEmbed = ('embeddedPFCandidates','packedPFCandidatesEmbedded','Embed')

handleMet = Handle ('vector<pat::MET>')
labelMet = ('slimmedMETs','','PAT')

handleMetEmbed = Handle ('vector<pat::MET>')
labelMetEmbed = ('slimmedMETsTEST','','Embed')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator', '', 'SIMembedding' )

handlePVInfo = Handle ('vector<reco::Vertex>')
labelPVInfo = ('offlineSlimmedPrimaryVertices', '', 'PAT')

handlePVInfoEmbed = Handle ('vector<reco::Vertex>')
labelPVInfoEmbed = ('embeddedPFCandidates','offlineSlimmedPrimaryVerticesEmbedded','Embed')

handlePFJets = Handle ('vector<reco::PFJet>')
labelPFJets = ('ak4PFJets', '', 'SIMembedding')

handlePFJetsEmbed = Handle ('vector<reco::PFJet>')
labelPFJetsEmbed = ('ak4PFJets', '', 'Embed')

handlePFTaus = Handle ('vector<reco::PFTau>')
#labelPFTaus = ('combinatoricRecoTaus', '', 'SIMembedding')
labelPFTaus = ('hpsPFTauProducerSansRefs', '', 'SIMembedding')

handlePFTausEmbed = Handle ('vector<reco::PFTau>')
#labelPFTausEmbed = ('combinatoricRecoTausEmbed', '', 'Embed')
labelPFTausEmbed = ('hpsPFTauProducerSansRefsEmbed', '', 'Embed')

handlePosition = Handle ('ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >')
labelPosition = ('externalLHEProducer','vertexPosition', 'LHEembedding')

handleZmumuCandidate = Handle ('vector<reco::CompositeCandidate>')
labelZmumuCandidate = ('ZmumuCandidates', '', 'SELECT')

handleCA8Jets = Handle ('vector<reco::PFJet>')
labelCA8Jets = ("ca8PFJetsCHSprunedForBoostedTausPATBoostedEmbed","subJetsForSeedingBoostedTausPAT","Embed")

nfail = 0
for i in range(1, 501):
#for i in range(1, 2):
    try:
        events = Events(prefix+"embed_muon_v1_YMuMu_pth400_"+str(i)+".root")
        print(prefix+"embed_muon_v1_YMuMu_pth400_"+str(i)+".root")
        nevt=0
        for event in events:
            nevt+=1
        
            event.getByLabel(labelPVInfo, handlePVInfo)
            pv = handlePVInfo.product()[0]
            h["pvxy"].Fill(pv.position().x(), pv.position().y())
            h["pvz"].Fill(pv.position().z())
        
            event.getByLabel(labelPVInfoEmbed, handlePVInfoEmbed)
            pvEmbed = handlePVInfoEmbed.product()[0]
            h["pvxyEmbed"].Fill(pvEmbed.position().x(), pvEmbed.position().y())
            h["pvzEmbed"].Fill(pvEmbed.position().z())

            #h["nVtx"]
        
            event.getByLabel(labelPosition, handlePosition)
            position = handlePosition.product()
            #print(position.Coordinates().x())
            h["pvxyGen"].Fill(position.Coordinates().x(), position.Coordinates().y())
            h["pvzGen"].Fill(position.Coordinates().z())
            
            event.getByLabel(labelTaus, handleTaus)
            taus=handleTaus.product()
        
            event.getByLabel(labelMuons, handleMuons)
            muons=handleMuons.product()
            
            event.getByLabel(labelTausEmbed, handleTausEmbed)
            tausEmbed=handleTausEmbed.product()
         
            event.getByLabel(labelCands, handleCands)
            cands=handleCands.product()
         
            event.getByLabel(labelCandsEmbed, handleCandsEmbed)
            candsEmbed=handleCandsEmbed.product()
            
            event.getByLabel(labelMet, handleMet)
            met=handleMet.product()
         
            event.getByLabel(labelMetEmbed, handleMetEmbed)
            metEmbed=handleMetEmbed.product()
        
            event.getByLabel(labelGenInfo, handleGenInfo)
            geninfo=handleGenInfo.product()
            genweight=geninfo.weight()
        
            event.getByLabel(labelGenParticles, handleGenParticles)
            genParticles=handleGenParticles.product()
        
            event.getByLabel(labelPFJets, handlePFJets)
            pfJets=handlePFJets.product()

            event.getByLabel(labelCA8Jets, handleCA8Jets)
            ca8Jets=handleCA8Jets.product()

            event.getByLabel(labelZmumuCandidate, handleZmumuCandidate)
            ZmumuCand = handleZmumuCandidate.product()
            #print(ZmumuCand)
        
            h['MET'].Fill(met[0].pt())
            h['METEmbed'].Fill(metEmbed[0].pt(), genweight)
            h['Weights'].Fill(genweight)
        
            selected_gentaus = []
            for particle in genParticles:
                if abs(particle.pdgId()) == 13 and particle.isHardProcess():
                    selected_gentaus += [particle]
            if len(selected_gentaus) > 1:
                #if selected_gentaus[0].pt() > selected_gentaus[1].pt(): h["genTauPt"].Fill(selected_gentaus[0].pt())
                #else: h["genTauPt"].Fill(selected_gentaus[1].pt())
                h["genTauPt"].Fill(selected_gentaus[0].pt(), genweight)
                h["genTauPt"].Fill(selected_gentaus[1].pt(), genweight)
        
                t1=ROOT.TLorentzVector()
                t1.SetPtEtaPhiM(selected_gentaus[0].pt(), selected_gentaus[0].eta(), selected_gentaus[0].phi(), selected_gentaus[0].mass()) 
            
                t2=ROOT.TLorentzVector()
                t2.SetPtEtaPhiM(selected_gentaus[1].pt(), selected_gentaus[1].eta(), selected_gentaus[1].phi(), selected_gentaus[1].mass())
                    
                h["genTauDR"].Fill(t1.DeltaR(t2), genweight)
                h["genTauM"].Fill((t1+t2).M(), genweight)
        
                #h2.Fill((t1+t2).M(), genweight)

            selected_muons=[]
            dR = 9999
            ipair = -1
            for i, cand in enumerate(ZmumuCand):
                p1 = cand.daughter(0)
                p2 = cand.daughter(1)
                m1=ROOT.TLorentzVector()
                m1.SetPtEtaPhiM(p1.pt(), p1.eta(), p1.phi(), p1.mass()) 
            
                m2=ROOT.TLorentzVector()
                m2.SetPtEtaPhiM(p2.pt(), p2.eta(), p2.phi(), p2.mass())

                if m1.DeltaR(m2)<dR:
                    dR = m1.DeltaR(m2)
                    ipair = i

            selected_muons+=[ZmumuCand[ipair].daughter(0)]
            selected_muons+=[ZmumuCand[ipair].daughter(1)]
        
            if len(selected_muons) > 1:
                mu1=ROOT.TLorentzVector()
                mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass()) 
            
                mu2=ROOT.TLorentzVector()
                mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass())

                h["mmumu"].Fill((mu1+mu2).M())
                h["dRmumu"].Fill(mu1.DeltaR(mu2))
                if selected_muons[0].pt() > selected_muons[1].pt():
                    h["mu1Pt"].Fill(mu1.Pt(), genweight)
                    h["mu2Pt"].Fill(mu2.Pt(), genweight)
                else:
                    h["mu1Pt"].Fill(mu2.Pt(), genweight)
                    h["mu2Pt"].Fill(mu1.Pt(), genweight)

            selected_taus=[]
            for tau in taus:
                if tau.pt()>3 and abs(tau.eta())<2.3 and tau.isLooseMuon():
                    selected_taus+=[tau]
        
            selected_taus_nocut=[]
            for tau in taus:
                if tau.pt()>3 and abs(tau.eta())<2.3 and tau.isMediumMuon():
                    selected_taus_nocut+=[tau]
        
            selected_taus_embed=[]
            for tau in tausEmbed:
                if tau.pt()>3 and abs(tau.eta())<2.3 and tau.isLooseMuon():
                    selected_taus_embed+=[tau]
        
            selected_taus_embed_deep=[]
            for tau in tausEmbed:
                if tau.pt()>3 and abs(tau.eta())<2.3 and tau.isLooseMuon():
                    selected_taus_embed_deep+=[tau]
        
            if len(selected_taus) > 1:
                
                mu1=ROOT.TLorentzVector()
                mu1.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass())
            
                mu2=ROOT.TLorentzVector()
                mu2.SetPtEtaPhiM(selected_taus[1].pt(), selected_taus[1].eta(), selected_taus[1].phi(), selected_taus[1].mass())

                drMin = 9999
                pt1 = -1
                pt2 = -1
                mass = -1
                for i, tau1 in enumerate(selected_taus):
                    for j, tau2 in enumerate(selected_taus):
                        if i > j or tau1.charge() * tau2.charge() > 0: continue
                        mu1=ROOT.TLorentzVector()
                        mu1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), tau1.mass()) 
            
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), tau2.mass())

                        if (mu1+mu2).M()<=3.56: continue
                        if mu1.DeltaR(mu2) < drMin:
                            mass = (mu1+mu2).M()
                            pt1 = mu1.Pt()
                            pt2 = mu2.Pt()
                            drMin = mu1.DeltaR(mu2)
                
                h["mtautau"].Fill(mass, genweight)
                if pt1 > pt2:
                    h["tau1Pt"].Fill(pt1, genweight)
                    h["tau2Pt"].Fill(pt2, genweight)
                else:
                    h["tau1Pt"].Fill(pt2, genweight)
                    h["tau2Pt"].Fill(pt1, genweight)
        
            if len(selected_taus_nocut) > 1:
                
                mu1=ROOT.TLorentzVector()
                mu1.SetPtEtaPhiM(selected_taus_nocut[0].pt(), selected_taus_nocut[0].eta(), selected_taus_nocut[0].phi(), selected_taus_nocut[0].mass())
            
                mu2=ROOT.TLorentzVector()
                mu2.SetPtEtaPhiM(selected_taus_nocut[1].pt(), selected_taus_nocut[1].eta(), selected_taus_nocut[1].phi(), selected_taus_nocut[1].mass())

                drMin = 9999
                pt1 = -1
                pt2 = -1
                mass = -1
                for i, tau1 in enumerate(selected_taus):
                    for j, tau2 in enumerate(selected_taus):
                        if i > j or tau1.charge() * tau2.charge() > 0: continue
                        mu1=ROOT.TLorentzVector()
                        mu1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), tau1.mass()) 
            
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), tau2.mass())

                        if (mu1+mu2).M()<=3.56: continue
                        if mu1.DeltaR(mu2) < drMin:
                            mass = (mu1+mu2).M()
                            pt1 = mu1.Pt()
                            pt2 = mu2.Pt()
                            drMin = mu1.DeltaR(mu2)
                
                h["mtautau_nocut"].Fill(mass, genweight)
                if pt1 > pt2:
                    h["tau1Pt_nocut"].Fill(pt1, genweight)
                    h["tau2Pt_nocut"].Fill(pt2, genweight)
                else:
                    h["tau1Pt_nocut"].Fill(pt2, genweight)
                    h["tau2Pt_nocut"].Fill(pt1, genweight)
        
            if len(selected_taus_embed) > 1:
                drMin = 9999
                pt1 = -1
                pt2 = -1
                mass = -1
                for i, tau1 in enumerate(selected_taus_embed):
                    for j, tau2 in enumerate(selected_taus_embed):
                        if i > j or tau1.charge() * tau2.charge() > 0: continue
                        mu1=ROOT.TLorentzVector()
                        mu1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), tau1.mass()) 
            
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), tau2.mass())

                        if (mu1+mu2).M()<=3.56: continue
                        if mu1.DeltaR(mu2) < drMin:
                            mass = (mu1+mu2).M()
                            pt1 = mu1.Pt()
                            pt2 = mu2.Pt()
                            drMin = mu1.DeltaR(mu2)
    
                #if drMin > 0.8: continue
                h["mtautauEmbed"].Fill(mass, genweight)
                if pt1 > pt2:
                    h["tau1PtEmbed"].Fill(pt1, genweight)
                    h["tau2PtEmbed"].Fill(pt2, genweight)
                else:
                    h["tau1PtEmbed"].Fill(pt2, genweight)
                    h["tau2PtEmbed"].Fill(pt1, genweight)
                h["dRtautauEmbed"].Fill(drMin, genweight)
        
            if len(selected_taus_embed_deep) > 1:
                drMin = 9999
                pt1 = -1
                pt2 = -1
                mass = -1
                for i, tau1 in enumerate(selected_taus_embed_deep):
                    for j, tau2 in enumerate(selected_taus_embed_deep):
                        if i > j or tau1.charge() * tau2.charge() > 0: continue
                        mu1=ROOT.TLorentzVector()
                        mu1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), tau1.mass()) 
            
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), tau2.mass())

                        if (mu1+mu2).M()<=3.56: continue
                        if mu1.DeltaR(mu2) < drMin:
                            mass = (mu1+mu2).M()
                            pt1 = mu1.Pt()
                            pt2 = mu2.Pt()
                            drMin = mu1.DeltaR(mu2)
    
                #if drMin > 0.8: continue
                h["mtautauEmbed_deep"].Fill(mass, genweight)
                if pt1 > pt2:
                    h["tau1PtEmbed_deep"].Fill(pt1, genweight)
                    h["tau2PtEmbed_deep"].Fill(pt2, genweight)
                else:
                    h["tau1PtEmbed_deep"].Fill(pt2, genweight)
                    h["tau2PtEmbed_deep"].Fill(pt1, genweight)
                h["dRtautauEmbed_deep"].Fill(drMin, genweight)
        del events
        gc.collect()
        print('------')
    except:
        nfail += 1
        print("Error...")

print("Nfail:", nfail)
out.cd()
for hist in sorted(h.keys()):
    h[hist].Write()
