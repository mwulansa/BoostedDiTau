import ROOT,sys
from DataFormats.FWLite import Events, Handle

events=Events('embedded.root')
#events=Events('outMuonSelection.root')
#events=Events('root://xrootd.unl.edu//store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0AB088EE-EA8A-E811-8636-0CC47A4C8E96.root')

out=ROOT.TFile('h_embedded.root','recreate')

h0 = ROOT.TH1F ("genTauPt", "", 500, 0, 500)
h1 = ROOT.TH1F ("mmumu", "", 50, 0, 500)
h2 = ROOT.TH1F ("mtautau", "", 50, 0, 500)
h3 = ROOT.TH1F ("mtautauEmbed", "", 50, 0, 500)
h4 = ROOT.TH1F ("Weights", "", 20, 0, 2)
h6 = ROOT.TH1F ("MET", "", 50, 0, 500)
h7 = ROOT.TH1F ("METEmbed", "", 50, 0, 500)
h8 = ROOT.TH1F ("genTauDR", "", 70, 0, 7)


handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles', '', 'SIMembedding')

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

handleMet = Handle ('vector<pat::MET>')
labelMet = ('slimmedMETs','','DQM')

handleMetEmbed = Handle ('vector<pat::MET>')
labelMetEmbed = ('slimmedMETsTEST','','Embed')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator', '', 'SIMembedding' )

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

    event.getByLabel(labelMet, handleMet)
    met=handleMet.product()
 
    event.getByLabel(labelMetEmbed, handleMetEmbed)
    metEmbed=handleMetEmbed.product()

    event.getByLabel(labelGenInfo, handleGenInfo)
    geninfo=handleGenInfo.product()
    genweight=geninfo.weight()

    event.getByLabel(labelGenParticles, handleGenParticles)
    genParticles=handleGenParticles.product()
    
    #print(taus.size(), tausEmbed.size(), tausMC.size())

    h6.Fill(met[0].pt())
    h7.Fill(metEmbed[0].pt(), genweight)
    h4.Fill(genweight)

    selected_gentaus = []
    for particle in genParticles:
        if abs(particle.pdgId()) == 15 and particle.isHardProcess():
            selected_gentaus += [particle]
    if len(selected_gentaus) > 1:
        if selected_gentaus[0].pt() > selected_gentaus[1].pt(): h0.Fill(selected_gentaus[0].pt())
        else: h0.Fill(selected_gentaus[1].pt())

        t1=ROOT.TLorentzVector()
        t1.SetPtEtaPhiM(selected_gentaus[0].pt(), selected_gentaus[0].eta(), selected_gentaus[0].phi(), selected_gentaus[0].mass()) 
    
        t2=ROOT.TLorentzVector()
        t2.SetPtEtaPhiM(selected_gentaus[1].pt(), selected_gentaus[1].eta(), selected_gentaus[1].phi(), selected_gentaus[1].mass())
            
        h8.Fill(t1.DeltaR(t2))

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
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
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

        h3.Fill((mu1+mu2).M(), genweight)

        

    #print("GenParticle")
    #for p in genParticle:
        #if abs(p.pdgId()) == 15 and p.isHardProcess():
            #print(p.pt())
        
    #print(len(selected_taus), len(selected_taus_embed))


out.cd()
h0.Write()
h1.Write()
h2.Write()
h3.Write()
h4.Write()
h6.Write()
h7.Write()
h8.Write()
