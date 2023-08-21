import ROOT,sys,gc
from DataFormats.FWLite import Events, Handle

#events=Events('condorCfg/embed_YMuMu_pth400_9.root')
#events=Events('outMuonSelection.root')
#events=Events('root://xrootd.unl.edu//store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0AB088EE-EA8A-E811-8636-0CC47A4C8E96.root')

out=ROOT.TFile('h_lhe.root','recreate')


#prefix = "root://cmseos.fnal.gov//store/user/zhangj/embedding/YMuMu/UL2017/"
#prefix = "./condorCfg/"
prefix = ""

h={}
hname = "genTauPt"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "mmumu"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)
hname = "muPt"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "dRmumu"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)


hname = "mu1Pt"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "mu2Pt"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "mu1Iso"
h[hname] = ROOT.TH1F (hname, "", 70, 0, 7)
hname = "mu2Iso"
h[hname] = ROOT.TH1F (hname, "", 70, 0, 7)

hname = "mu2Charge"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "mu2Neutral"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "mu2Photon"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)

hname = "mtautau"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "mtautauEmbed"
h[hname ]= ROOT.TH1F (hname, "", 500, 0, 500)

hname = "mtautau_deep"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "mtautauEmbed_deep"
h[hname ]= ROOT.TH1F (hname, "", 500, 0, 500)

hname = "Weights"
h[hname] = ROOT.TH1F (hname, "", 20, 0, 2)
hname = "MET"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "METEmbed"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "genTauDR"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = "genTauM"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)


handleGenParticles = Handle ('vector<reco::GenParticle>')
labelGenParticles = ('prunedGenParticles', '', 'SIMembedding')

handleTaus = Handle ('vector<pat::Tau>')
#labelTaus = ('slimmedTausNewID','','SELECT')
labelTaus = ('slimmedTausBoosted', '', 'SIMembedding')

handleMuons = Handle ('vector<pat::Muon>')
labelMuons = ('slimmedMuons','','PAT')

handleJets = Handle ('vector<pat::Jet>')
labelJets = ("slimmedJets","","PAT")

handleTausEmbed = Handle ('vector<pat::Tau>')
labelTausEmbed = ('selectedPatTausBoostedEmbed','','Embed')

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

nfail = 0
for i in range(1, 2):
    #events = Events(prefix+"embed_YMuMu_pth400_"+str(i)+".root")
    #print(events)
    #try:
        events = Events(prefix+"simulated_muon_YMuMu_pth400_"+str(i)+".root")
        print(prefix+"simulated_muon_YMuMu_pth400_"+str(i)+".root")
        nevt=0
        for event in events:
            nevt+=1
            #print('-----', nevt)
        
            genweight=1.
        
            event.getByLabel(labelGenParticles, handleGenParticles)
            genParticles=handleGenParticles.product()

            event.getByLabel(labelMuons, handleMuons)
            muons=handleMuons.product()

            event.getByLabel(labelJets, handleJets)
            jets=handleJets.product()

            event.getByLabel(labelZmumuCandidate, handleZmumuCandidate)
            ZmumuCand = handleZmumuCandidate.product()
            #print(ZmumuCand)
        
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

            #print(ZmumuCand.size())
            
            #selected_muons+=ZmumuCand.daughter(0)
            #selected_muons+=ZmumuCand.daughter(1)
            

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

            #print(ZmumuCand[ipair].numberOfDaughters())
            #print(type(ZmumuCand[ipair].daughter(0)))

            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(jets[0].pt(), jets[0].eta(), jets[0].phi(), jets[0].mass()) 
            
            #if dR < 0.4:
            selected_muons+=[ZmumuCand[ipair].daughter(0)]
            selected_muons+=[ZmumuCand[ipair].daughter(1)]
            muon1 = ZmumuCand[ipair].daughter(0)
            muon2 = ZmumuCand[ipair].daughter(1)
        
            if len(selected_muons) > 1:
                mu1=ROOT.TLorentzVector()
                mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass()) 
            
                mu2=ROOT.TLorentzVector()
                mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass())

                h["mmumu"].Fill((mu1+mu2).M())
                h["dRmumu"].Fill(mu1.DeltaR(mu2))
                h["muPt"].Fill(mu1.Pt(), genweight)
                h["muPt"].Fill(mu2.Pt(), genweight)

                if mu1.DeltaR(jet) > 1 and mu2.DeltaR(jet) > 1:

                    for muon in muons:
                        if not muon.isLooseMuon(): continue
                        if muon.pt()<3 or muon.eta()>2.4: continue
                        if muon.pt()-muon1.pt() < 0.01 and muon.charge() == muon1.charge():
                            iso =  (muon.pfIsolationR04().sumChargedHadronPt+max(0,muon.pfIsolationR04().sumPhotonEt+muon.pfIsolationR04().sumNeutralHadronEt-0.5*muon.pfIsolationR04().sumPUPt))/muon.pt()
                                                   
                            h['mu1Iso'].Fill(iso)
                            h['mu1Pt'].Fill(muon1.pt())
                        if muon.pt()-muon2.pt() < 0.01 and muon.charge() == muon2.charge():
                            iso =  (muon.pfIsolationR04().sumChargedHadronPt+max(0,muon.pfIsolationR04().sumPhotonEt+muon.pfIsolationR04().sumNeutralHadronEt-0.5*muon.pfIsolationR04().sumPUPt))/muon.pt()
                           
                            h['mu2Iso'].Fill(iso)
                            h['mu2Pt'].Fill(muon2.pt())
                            h['mu2Charge'].Fill(muon.pfIsolationR04().sumChargedHadronPt)
                            h['mu2Neutral'].Fill(muon.pfIsolationR04().sumPhotonEt)
                            h['mu2Photon'].Fill(muon.pfIsolationR04().sumNeutralHadronEt)
                
        
        

        del events
        gc.collect()
    #except:
    #    nfail += 1
    #    print("File Not Found...")

print("Nfail:", nfail)
out.cd()
for hist in sorted(h.keys()):
    h[hist].Write()
