import ROOT,sys
from DataFormats.FWLite import Events, Handle

events=Events('embedded.root')
#events=Events('outMuonSelection.root')
#events=Events('root://xrootd.unl.edu//store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0AB088EE-EA8A-E811-8636-0CC47A4C8E96.root')

out=ROOT.TFile('h_embedded.root','recreate')


h={}
hname = "genTauPt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "mmumu"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "mtautau"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "mtautauEmbed"
h[hname ]= ROOT.TH1F (hname, "", 50, 0, 500)

hname = "mtautau_deep"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "mtautauEmbed_deep"
h[hname ]= ROOT.TH1F (hname, "", 50, 0, 500)

hname = "Weights"
h[hname] = ROOT.TH1F (hname, "", 20, 0, 2)
hname = "MET"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "METEmbed"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "genTauDR"
h[hname] = ROOT.TH1F (hname, "", 70, 0, 7)
hname = "genTauM"
h[hname] = ROOT.TH1F (hname, "", 150, 0, 150)

hname = "tauPt"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 100)
hname = "tauPtEmbed"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 100)

hname = "tauPt_deep"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 100)
hname = "tauPtEmbed_deep"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 100)

hname = 'pfJetPt'
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = 'pfJetPtEmbed'
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)

hname = 'pfTauPt'
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = 'pfTauPtEmbed'
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles', '', 'SIMembedding')

handleTaus = Handle ('vector<pat::Tau>')
#labelTaus = ('slimmedTausNewID','','SELECT')
labelTaus = ('slimmedTaus', '', 'SIMembedding')

handleMuons = Handle ('vector<pat::Muon>')
labelMuons = ('slimmedMuons','','DQM')

handleTausEmbed = Handle ('vector<pat::Tau>')
labelTausEmbed = ('selectedPatTausEmbed','','Embed')

handleGenParticles = Handle ('vector<reco::GenParticle>')
labelGenParticles = ('prunedGenParticles','','SIMembedding')

handleCands = Handle ('vector<pat::PackedCandidate>')
labelCands = ('packedPFCandidates','','DQM')

handleCandsEmbed = Handle ('vector<pat::PackedCandidate>')
labelCandsEmbed = ('embeddedPFCandidates','packedPFCandidatesEmbedded','Embed')

handleMet = Handle ('vector<pat::MET>')
labelMet = ('slimmedMETs','','DQM')

handleMetEmbed = Handle ('vector<pat::MET>')
labelMetEmbed = ('slimmedMETsTEST','','Embed')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator', '', 'SIMembedding' )

handlePVInfo = Handle ('vector<reco::Vertex>')
labelPVInfo = ('offlineSlimmedPrimaryVertices', '', 'DQM')

handlePVInfoEmbed = Handle ('vector<reco::Vertex>')
labelPVInfoEmbed = ('embeddedPFCandidates','offlineSlimmedPrimaryVerticesEmbedded','Embed')

handlePFJets = Handle ('vector<reco::PFJet>')
labelPFJets = ('ak4PFJets', '', 'SIMembedding')

handlePFJetsEmbed = Handle ('vector<reco::PFJet>')
labelPFJetsEmbed = ('ak4PFJetsPATEmbed', '', 'Embed')

handlePFTaus = Handle ('vector<reco::PFTau>')
#labelPFTaus = ('combinatoricRecoTaus', '', 'SIMembedding')
labelPFTaus = ('hpsPFTauProducerSansRefs', '', 'SIMembedding')

handlePFTausEmbed = Handle ('vector<reco::PFTau>')
#labelPFTausEmbed = ('combinatoricRecoTausEmbed', '', 'Embed')
labelPFTausEmbed = ('hpsPFTauProducerSansRefsEmbed', '', 'Embed')

handlePosition = Handle ('ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >')
labelPosition = ('externalLHEProducer','vertexPosition', 'LHE')

nevt=0
for event in events:
    nevt+=1
    print('-----', nevt)

    event.getByLabel(labelPVInfo, handlePVInfo)
    pv = handlePVInfo.product()[0]

    event.getByLabel(labelPVInfoEmbed, handlePVInfoEmbed)
    pvEmbed = handlePVInfoEmbed.product()[0]

    event.getByLabel(labelPosition, handlePosition)
    position = handlePosition.product()

    print(position.x(), position.y(), position.z())
    print(handlePVInfo.product().size(), pv.x(), pv.y(), pv.z())
    print(handlePVInfoEmbed.product().size(), pvEmbed.x(), pvEmbed.y(), pvEmbed.z())

    #for track in pv.tracks():
    #    print(track)
    
    event.getByLabel(labelTaus, handleTaus)
    taus=handleTaus.product()

    event.getByLabel(labelMuons, handleMuons)
    muons=handleMuons.product()
    
    event.getByLabel(labelTausEmbed, handleTausEmbed)
    tausEmbed=handleTausEmbed.product()
  
    event.getByLabel(labelGenParticles, handleGenParticles)
    genParticle=handleGenParticles.product()
 
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

    if pfJets.size() > 0:
        h["pfJetPt"].Fill(pfJets[0].pt())

    event.getByLabel(labelPFJetsEmbed, handlePFJetsEmbed)
    pfJetsEmbed=handlePFJetsEmbed.product()

    if pfJetsEmbed.size() > 0:
        h["pfJetPtEmbed"].Fill(pfJetsEmbed[0].pt())

    event.getByLabel(labelPFTaus, handlePFTaus)
    pfTaus=handlePFTaus.product()

    if pfTaus.size() > 0:
        h["pfTauPt"].Fill(pfTaus[0].pt())

    event.getByLabel(labelPFTausEmbed, handlePFTausEmbed)
    pfTausEmbed=handlePFTausEmbed.product()

    if pfTausEmbed.size() > 0:
        h["pfTauPtEmbed"].Fill(pfTausEmbed[0].pt())
    
    #print(taus.size(), tausEmbed.size(), tausMC.size())

    h['MET'].Fill(met[0].pt())
    h['METEmbed'].Fill(metEmbed[0].pt(), genweight)
    h['Weights'].Fill(genweight)

    selected_gentaus = []
    for particle in genParticles:
        if abs(particle.pdgId()) == 15 and particle.isHardProcess():
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
    for muon in muons:
        if not muon.isLooseMuon(): continue
        if muon.pt()<3 or muon.eta()>2.4: continue
        if (muon.pfIsolationR04().sumChargedHadronPt+max(0,muon.pfIsolationR04().sumPhotonEt+muon.pfIsolationR04().sumNeutralHadronEt-0.5*muon.pfIsolationR04().sumPUPt))/muon.pt()<0.4:
            selected_muons+=[muon]

    selected_taus=[]
    for tau in taus:
        #print(tau.pt())
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding"):
        if tau.pt()>18 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") > 0.5 and tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"):
            selected_taus+=[tau]

    selected_taus_embed=[]
    for tau in tausEmbed:
        #print("embed:",tau.pt())
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding"):
        if tau.pt()>18 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") and tau.tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") > 0.5 and tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"):
            selected_taus_embed+=[tau]

    selected_taus_deep=[]
    for tau in taus:
        #print(tau.pt())
        if tau.pt()>18 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") > 0.5 and tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"):
            selected_taus_deep+=[tau]

    selected_taus_embed_deep=[]
    for tau in tausEmbed:
        #print("embed:",tau.pt())
        if tau.pt()>18 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") and tau.tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") > 0.5 and tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"):
            selected_taus_embed_deep+=[tau]

    if len(selected_muons) > 1:
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass())

        h["mmumu"].Fill((mu1+mu2).M())

    if len(selected_taus) > 1:
        if selected_taus[0].charge() * selected_taus[1].charge() > 0: continue
        
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_taus[1].pt(), selected_taus[1].eta(), selected_taus[1].phi(), selected_taus[1].mass())

        h["mtautau"].Fill((mu1+mu2).M(), genweight)

        #if mu1.Pt() > mu2.Pt(): h["tauPt"].Fill(mu1.Pt(), genweight)
        #else: h["tauPt"].Fill(mu2.Pt(), genweight)
        h["tauPt"].Fill(mu1.Pt(), genweight)
        h["tauPt"].Fill(mu2.Pt(), genweight)

    if len(selected_taus_embed) > 1:
        if selected_taus_embed[0].charge() * selected_taus_embed[1].charge() > 0: continue
        
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_taus_embed[0].pt(), selected_taus_embed[0].eta(), selected_taus_embed[0].phi(), selected_taus_embed[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_taus_embed[1].pt(), selected_taus_embed[1].eta(), selected_taus_embed[1].phi(), selected_taus_embed[1].mass())

        h["mtautauEmbed"].Fill((mu1+mu2).M(), genweight)

        #if mu1.Pt() > mu2.Pt(): h["tauPtEmbed"].Fill(mu1.Pt(), genweight)
        #else: h["tauPtEmbed"].Fill(mu2.Pt(), genweight)
        h["tauPtEmbed"].Fill(mu1.Pt(), genweight)
        h["tauPtEmbed"].Fill(mu2.Pt(), genweight)

    if len(selected_taus_deep) > 1:
        if selected_taus_deep[0].charge() * selected_taus_deep[1].charge() > 0: continue
        
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_taus_deep[0].pt(), selected_taus_deep[0].eta(), selected_taus_deep[0].phi(), selected_taus_deep[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_taus_deep[1].pt(), selected_taus_deep[1].eta(), selected_taus_deep[1].phi(), selected_taus_deep[1].mass())

        h["mtautau_deep"].Fill((mu1+mu2).M(), genweight)
        h["tauPt_deep"].Fill(mu1.Pt(), genweight)
        h["tauPt_deep"].Fill(mu2.Pt(), genweight)

    if len(selected_taus_embed_deep) > 1:
        if selected_taus_embed_deep[0].charge() * selected_taus_embed_deep[1].charge() > 0: continue
        
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_taus_embed_deep[0].pt(), selected_taus_embed_deep[0].eta(), selected_taus_embed_deep[0].phi(), selected_taus_embed_deep[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_taus_embed_deep[1].pt(), selected_taus_embed_deep[1].eta(), selected_taus_embed_deep[1].phi(), selected_taus_embed_deep[1].mass())

        h["mtautauEmbed_deep"].Fill((mu1+mu2).M(), genweight)
        h["tauPtEmbed_deep"].Fill(mu1.Pt(), genweight)
        h["tauPtEmbed_deep"].Fill(mu2.Pt(), genweight)


out.cd()
for hist in sorted(h.keys()):
    h[hist].Write()
