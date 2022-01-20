import ROOT,sys
from DataFormats.FWLite import Events, Handle

events=Events('embedded.root')
#events=Events('outMuonSelection.root')
#events=Events('root://xrootd.unl.edu//store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0AB088EE-EA8A-E811-8636-0CC47A4C8E96.root')

out=ROOT.TFile('h_embedded.root','recreate')

h1 = ROOT.TH1F ("mmumu", "", 50, 0, 500)
h2 = ROOT.TH1F ("mtautau", "", 50, 0, 500)
h3 = ROOT.TH1F ("mtautauEmbed", "", 50, 0, 500)



handleTaus = Handle ('vector<pat::Tau>')
labelTaus = ('slimmedTausNewID','','SELECT')
#labelTaus = ('slimmedTaus')

handleMuons = Handle ('vector<pat::Muon>')
labelMuons = ('slimmedMuons','','DQM')

handleTausEmbed = Handle ('vector<pat::Tau>')
labelTausEmbed = ('selectedPatTausEmbed','','Embed')

handleTausMC = Handle ('vector<pat::Tau>')
labelTausMC = ('slimmedTaus','','SIMembedding')

handleGenParticles = Handle ('vector<reco::GenParticle>')
labelGenParticles = ('prunedGenParticles','','SIMembedding')

handleCands = Handle ('vector<pat::PackedCandidate>')
labelCands = ('packedPFCandidates','','DQM')

handleCandsEmbed = Handle ('vector<pat::PackedCandidate>')
labelCandsEmbed = ('embeddedPFCandidates','packedPFCandidatesEmbedded','Embed')

for event in events:
    #print('-----')
    
    event.getByLabel(labelTaus, handleTaus)
    taus=handleTaus.product()

    event.getByLabel(labelMuons, handleMuons)
    muons=handleMuons.product()
    
    event.getByLabel(labelTausEmbed, handleTausEmbed)
    tausEmbed=handleTausEmbed.product()
 
    event.getByLabel(labelTausMC, handleTausMC)
    tausMC=handleTausMC.product()
 
    event.getByLabel(labelGenParticles, handleGenParticles)
    genParticle=handleGenParticles.product()
 
    event.getByLabel(labelCands, handleCands)
    cands=handleCands.product()
 
    event.getByLabel(labelCandsEmbed, handleCandsEmbed)
    candsEmbed=handleCandsEmbed.product()
 
    #print(taus.size(), tausEmbed.size(), tausMC.size())

    selected_muons=[]
    for muon in muons:
        if not muon.isLooseMuon(): continue
        if muon.pt()<3 or muon.eta()>2.4: continue
        if (muon.pfIsolationR04().sumChargedHadronPt+max(0,muon.pfIsolationR04().sumPhotonEt+muon.pfIsolationR04().sumNeutralHadronEt-0.5*muon.pfIsolationR04().sumPUPt))/muon.pt()<0.4:
            selected_muons+=[muon]

    selected_taus=[]
    for tau in taus:
        #print(tau.pt())
        if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet") and tau.tauID("byVLooseDeepTau2017v2p1VSmu"):
            selected_taus+=[tau]

    selected_taus_embed=[]
    for tau in tausEmbed:
        #print("embed:",tau.pt())
        if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
            selected_taus_embed+=[tau]

    if len(selected_muons) > 1:
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass())

        h1.Fill((mu1+mu2).M())

    if len(selected_taus) > 1:
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_taus[1].pt(), selected_taus[1].eta(), selected_taus[1].phi(), selected_taus[1].mass())

        h2.Fill((mu1+mu2).M())

    if len(selected_taus_embed) > 1:
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_taus_embed[0].pt(), selected_taus_embed[0].eta(), selected_taus_embed[0].phi(), selected_taus_embed[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_taus_embed[1].pt(), selected_taus_embed[1].eta(), selected_taus_embed[1].phi(), selected_taus_embed[1].mass())

        h3.Fill((mu1+mu2).M())

        

    #print("GenParticle")
    #for p in genParticle:
        #if abs(p.pdgId()) == 15 and p.isHardProcess():
            #print(p.pt())
        
    #print(len(selected_taus), len(selected_taus_embed))


out.cd()
h1.Write()
h2.Write()
h3.Write()
