import ROOT,sys,gc
from DataFormats.FWLite import Events, Handle

from looseElectron import *

#events=Events('condorCfg/embed_YMuMu_pth400_9.root')
#events=Events('outMuonSelection.root')
#events=Events('root://xrootd.unl.edu//store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0AB088EE-EA8A-E811-8636-0CC47A4C8E96.root')

i = sys.argv[1]

out=ROOT.TFile('./tmp/h_embedded_boosted_'+str(i)+'.root','recreate')

#prefix = "root://cmseos.fnal.gov//eos/uscms/store/group/lpcsusyhiggs/Jingyu/embedding/YMuMu/UL2017/"
prefix = "./condorCfg/"


def calcDRWeight(dR):
    return 1.04-0.34*dR

h={}
hname = "genTau1Pt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "genTau2Pt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "genTauDR"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = "genTauM"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)

hname = "genTaumu1Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "genTaumu2Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "genTaumuDR"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = "genTaumuM"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)

hname = "mmumu"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)
hname = "mu1Pt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "mu2Pt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 500)
hname = "dRmumu"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = "mmumuEmbed"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)

hname = "mmumuEmbed_cor"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)

hname = "mmumuEmbed_genM"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 20)

hname = "mu1PtEmbed"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "mu2PtEmbed"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)

hname = "mu1PtEmbed_cor"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "mu2PtEmbed_cor"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)

hname = "mmumuEmbed_tight"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)
hname = "mu1PtEmbed_tight"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "mu2PtEmbed_tight"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "dRmumuEmbed_tight"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = "mmumuEmbed_med"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)
hname = "mu1PtEmbed_med"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "mu2PtEmbed_med"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "dRmumuEmbed_med"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = "dRmumuEmbed"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = "dRmumuEmbed_cor"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)


hname = "ca8JetPt"
h[hname] = ROOT.TH1F (hname, "", 500, 0, 1000)

hname = "mtautau"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)
hname = "mtautauEmbed"
h[hname ]= ROOT.TH1F (hname, "", 50, 0, 50)

hname = "dRtautauEmbed"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = "dRtautauEmbed_deep"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)
hname = "dRtautauEmbed_deep_cor"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = "mtautau_deep"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)
hname = "mtautauEmbed_deep"
h[hname ]= ROOT.TH1F (hname, "", 50, 0, 50)
hname = "mtautauEmbed_deep_cor"
h[hname ]= ROOT.TH1F (hname, "", 50, 0, 50)

hname = 'memuEmbed'
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)
hname = 'memuEmbed_cor'
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)

hname = 'ePtEmbed'
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = 'muPtEmbed'
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = 'dRemuEmbed'
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = 'ePtEmbed_cor'
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = 'muPtEmbed_cor'
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = 'dRemuEmbed_cor'
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = "Weights"
h[hname] = ROOT.TH1F (hname, "", 20, 0, 2)
hname = "MET"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)
hname = "METEmbed"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 500)

hname = "tau1Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2Pt"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau1PtEmbed"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2PtEmbed"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)

hname = "tau1Pt_deep"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2Pt_deep"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau1PtEmbed_deep"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2PtEmbed_deep"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)

hname = "tau1PtEmbed_deep_cor"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "tau2PtEmbed_deep_cor"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)

hname = "mmumuSim"
h[hname] = ROOT.TH1F (hname, "", 50, 0, 50)
hname = "mu1PtSim"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "mu2PtSim"
h[hname] = ROOT.TH1F (hname, "", 100, 0, 500)
hname = "dRmumuSim"
h[hname] = ROOT.TH1F (hname, "", 80, 0, 2)

hname = "pvxyEmbed"
h[hname] = ROOT.TH2F (hname, "", 80, -2, 2, 80, -2, 2)
hname = "pvzEmbed"
h[hname] = ROOT.TH1F (hname, "", 80, -20, 20)

hname = "pvxySim"
h[hname] = ROOT.TH2F (hname, "", 80, -2, 2, 80, -2, 2)
hname = "pvzSim"
h[hname] = ROOT.TH1F (hname, "", 80, -20, 20)

hname = "pvxySimBS"
h[hname] = ROOT.TH2F (hname, "", 80, -2, 2, 80, -2, 2)
hname = "pvzSimBS"
h[hname] = ROOT.TH1F (hname, "", 80, -20, 20)

hname = "pvxyReco"
h[hname] = ROOT.TH2F (hname, "", 80, -2, 2, 80, -2, 2)
hname = "pvzReco"
h[hname] = ROOT.TH1F (hname, "", 80, -20, 20)

hname = "pvxy"
h[hname] = ROOT.TH2F (hname, "", 80, -2, 2, 80, -2, 2)
hname = "pvz"
h[hname] = ROOT.TH1F (hname, "", 80, -20, 20)

hname = "bsxy"
h[hname] = ROOT.TH2F (hname, "", 80, -2, 2, 80, -2, 2)
hname = "bsz"
h[hname] = ROOT.TH1F (hname, "", 80, -20, 20)

hname = "pvxyGen"
h[hname] = ROOT.TH2F (hname, "", 80, -2, 2, 80, -2, 2)
hname = "pvzGen"
h[hname] = ROOT.TH1F (hname, "", 80, -20, 20)

hname = "nVtxSim"
h[hname] = ROOT.TH1F (hname, "", 5, 0, 5)

handleGenParticles = Handle ('vector<reco::GenParticle>')
labelGenParticles = ('prunedGenParticles', '', 'SIMembedding')

handleTaus = Handle ('vector<pat::Tau>')
#labelTaus = ('slimmedTausNewID','','SELECT')
labelTaus = ('slimmedTausBoosted', '', 'SIMembedding')

handleMuons = Handle ('vector<pat::Muon>')
labelMuons = ('slimmedMuons','','PAT')

handleMuonsSim = Handle ('vector<pat::Muon>')
labelMuonsSim = ('slimmedMuons', '', 'SIMembedding')

handleMuonsEmbed = Handle ('vector<pat::Muon>')
labelMuonsEmbed = ("embeddedPFCandidates","slimmedMuonsEmbedded","Embed")

handleElectronsEmbed = Handle ('vector<pat::Electron>')
labelElectronsEmbed = ("embeddedPFCandidates","slimmedElectronsEmbedded","Embed")

handleTausEmbed = Handle ('vector<pat::Tau>')
labelTausEmbed = ('selectedPatTausBoostedEmbed','','Embed')
#labelTausEmbed = ("selectedPatTausNoNewIDsBoostedEmbed","","Embed")

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

handlePVInfoSimBS = Handle ('vector<reco::Vertex>')
labelPVInfoSimBS = ('offlineSlimmedPrimaryVerticesWithBS', '', 'SIMembedding')

handlePVInfoSim = Handle ('vector<reco::Vertex>')
labelPVInfoSim = ('offlineSlimmedPrimaryVertices', '', 'SIMembedding')

handlePVInfoReco = Handle ('vector<reco::Vertex>')
labelPVInfoReco = ('offlinePrimaryVertices', '', 'SIMembedding')

handlePosition = Handle ('ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >')
labelPosition = ('externalLHEProducer','vertexPosition', 'LHEembedding')

handlePVInfoEmbed = Handle ('vector<reco::Vertex>')
labelPVInfoEmbed = ('embeddedPFCandidates','offlineSlimmedPrimaryVerticesEmbedded','Embed')

handlePVInfoGen = Handle ('edm::HepMCProduct')
labelPVInfoGen = ('VtxSmeared', '', 'SIMembedding')

handleBS = Handle ('reco::BeamSpot')
labelBS = ('offlineBeamSpot', '', 'SIMembedding')

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

handleZmumuCandidate = Handle ('vector<reco::CompositeCandidate>')
labelZmumuCandidate = ('ZmumuCandidates', '', 'SELECT')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll', '', 'RECO')

handleCA8Jets = Handle ('vector<reco::PFJet>')
labelCA8Jets = ("ca8PFJetsCHSprunedForBoostedTausPATBoostedEmbed","subJetsForSeedingBoostedTausPAT","Embed")

nfail = 0
#for i in range(1, 501):
#for i in range(149, 150):
#    try:
#events = Events(prefix+"embed_v1_YMuMu_pth400_"+str(i)+".root")
#print(prefix+"embed_v1_YMuMu_pth400_"+str(i)+".root")
events = Events(prefix+"embed_v1_YMuMu_pth400_"+str(i)+".root")
print(prefix+"embed_v1_YMuMu_pth400_"+str(i)+".root")
nevt=0
for event in events:
    nevt+=1
    #print('-----', nevt)

    event.getByLabel(labelPVInfo, handlePVInfo)
    pv = handlePVInfo.product()[0]

    h["pvxy"].Fill(pv.position().x(), pv.position().y())
    h["pvz"].Fill(pv.position().z())

    event.getByLabel(labelPVInfoEmbed, handlePVInfoEmbed)
    pvEmbed = handlePVInfoEmbed.product()[0]

    h["pvxyEmbed"].Fill(pvEmbed.position().x(), pvEmbed.position().y())
    h["pvzEmbed"].Fill(pvEmbed.position().z())

    event.getByLabel(labelPVInfoSim, handlePVInfoSim)
    pvSim = handlePVInfoSim.product()[0]

    h["pvxySim"].Fill(pvSim.position().x(), pvSim.position().y())
    h["pvzSim"].Fill(pvSim.position().z())
    h["nVtxSim"].Fill(handlePVInfoSim.product().size())

    event.getByLabel(labelPVInfoSimBS, handlePVInfoSimBS)
    pvSimBS = handlePVInfoSimBS.product()[0]

    h["pvxySimBS"].Fill(pvSimBS.position().x(), pvSimBS.position().y())
    h["pvzSimBS"].Fill(pvSimBS.position().z())

    event.getByLabel(labelPVInfoReco, handlePVInfoReco)
    pvReco = handlePVInfoReco.product()[0]

    h["pvxyReco"].Fill(pvReco.position().x(), pvReco.position().y())
    h["pvzReco"].Fill(pvReco.position().z())
    
    event.getByLabel(labelPVInfoGen, handlePVInfoGen)
    pvGen = handlePVInfoGen.product()

    #print(pvGen.GetEvent().vertices_begin().position().z())
    h["pvxyGen"].Fill(pvGen.GetEvent().vertices_begin().position().x()/10, pvGen.GetEvent().vertices_begin().position().y()/10)
    h["pvzGen"].Fill(pvGen.GetEvent().vertices_begin().position().z()/10)
    
    event.getByLabel(labelTaus, handleTaus)
    taus=handleTaus.product()

    event.getByLabel(labelMuons, handleMuons)
    muons=handleMuons.product()

    event.getByLabel(labelMuonsSim, handleMuonsSim)
    muonsSim=handleMuonsSim.product()

    event.getByLabel(labelMuonsEmbed, handleMuonsEmbed)
    muonsEmbed=handleMuonsEmbed.product()

    event.getByLabel(labelElectronsEmbed, handleElectronsEmbed)
    electronsEmbed=handleElectronsEmbed.product()
    
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

    event.getByLabel(labelRho, handleRho)
    rho=handleRho.product()[0]

    if ca8Jets.size() > 0:
        h["ca8JetPt"].Fill(ca8Jets[0].pt())

    event.getByLabel(labelZmumuCandidate, handleZmumuCandidate)
    ZmumuCand = handleZmumuCandidate.product()
    #print(ZmumuCand)

    event.getByLabel(labelPosition, handlePosition)
    position = handlePosition.product()
    #print(position.Coordinates().x())
    #h["pvxyGen"].Fill(position.Coordinates().x(), position.Coordinates().y())
    #h["pvzGen"].Fill(position.Coordinates().z())

    event.getByLabel(labelBS, handleBS)
    bs = handleBS.product()
    h["bsxy"].Fill(bs.position().x(), bs.position().y())
    h["bsz"].Fill(bs.position().z())
    

    h['MET'].Fill(met[0].pt())
    h['METEmbed'].Fill(metEmbed[0].pt(), genweight)
    h['Weights'].Fill(genweight)

    selected_gentaus = []
    decay_products = []
    selected_gentaumus = []
    for particle in genParticles:
        if abs(particle.pdgId()) == 15 and particle.isHardProcess():
            selected_gentaus += [particle]
        if particle.isDirectHardProcessTauDecayProductFinalState():
            decay_products +=[particle.pdgId()]
            if abs(particle.pdgId()) == 13:
                selected_gentaumus+=[particle]

    if len(selected_gentaumus) == 1:
        h["genTaumu1Pt"].Fill(selected_gentaumus[0].pt())

    if len(selected_gentaumus) == 2:
        if selected_gentaumus[0].pt() > selected_gentaumus[1].pt():
            h["genTaumu1Pt"].Fill(selected_gentaumus[0].pt())
            h["genTaumu2Pt"].Fill(selected_gentaumus[1].pt())
        else:
            h["genTaumu1Pt"].Fill(selected_gentaumus[1].pt())
            h["genTaumu2Pt"].Fill(selected_gentaumus[0].pt())

        t1=ROOT.TLorentzVector()
        t1.SetPtEtaPhiM(selected_gentaumus[0].pt(), selected_gentaumus[0].eta(), selected_gentaumus[0].phi(), selected_gentaumus[0].mass()) 
    
        t2=ROOT.TLorentzVector()
        t2.SetPtEtaPhiM(selected_gentaumus[1].pt(), selected_gentaumus[1].eta(), selected_gentaumus[1].phi(), selected_gentaumus[1].mass())
        
        h["genTaumuM"].Fill((t1+t2).M())
        h["genTaumuDR"].Fill(t1.DeltaR(t2))
    
    genTauM = -1
    genTau1Pt = -1
    genTau2Pt = -1
    genDR = -1
    if len(selected_gentaus) > 1:
        #if selected_gentaus[0].pt() > selected_gentaus[1].pt(): h["genTauPt"].Fill(selected_gentaus[0].pt())
        #else: h["genTauPt"].Fill(selected_gentaus[1].pt())
        if selected_gentaus[0].pt() > selected_gentaus[1].pt():
            h["genTau1Pt"].Fill(selected_gentaus[0].pt(), genweight)
            h["genTau2Pt"].Fill(selected_gentaus[1].pt(), genweight)
        else:
            h["genTau1Pt"].Fill(selected_gentaus[1].pt(), genweight)
            h["genTau2Pt"].Fill(selected_gentaus[0].pt(), genweight)

        t1=ROOT.TLorentzVector()
        t1.SetPtEtaPhiM(selected_gentaus[0].pt(), selected_gentaus[0].eta(), selected_gentaus[0].phi(), selected_gentaus[0].mass()) 
    
        t2=ROOT.TLorentzVector()
        t2.SetPtEtaPhiM(selected_gentaus[1].pt(), selected_gentaus[1].eta(), selected_gentaus[1].phi(), selected_gentaus[1].mass())
            
        h["genTauDR"].Fill(t1.DeltaR(t2), genweight)
        h["genTauM"].Fill((t1+t2).M(), genweight)
        genTauM = (t1+t2).M()
        genTau1Pt = t1.Pt()
        genTau2Pt = t2.Pt()
        genDR = t1.DeltaR(t2)
        #h2.Fill((t1+t2).M(), genweight)

    if genDR < 0.01: continue
    drWeight = calcDRWeight(genDR)

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
    selected_muons+=[ZmumuCand[ipair].daughter(0)]
    selected_muons+=[ZmumuCand[ipair].daughter(1)]

    selected_electrons_embed=[]
    for electron in electronsEmbed:
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
                selected_electrons_embed+=[electron]

        if electron.isEE():
            if electron.full5x5_sigmaIetaIeta() < 0.0425 \
               and GsfEleEInverseMinusPInverse(electron) < 0.111 \
               and abs(dEtaInSeed(electron)) < 0.00674 \
               and GsfEleMissingHitsCut(electron)<=1 \
               and electron.passConversionVeto() \
               and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.169 \
               and electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c):
                selected_electrons_embed+=[electron]

    selected_muons_embed=[]
    for muon in muonsEmbed:
        if not muon.isLooseMuon(): continue
        if muon.pt()<3 or muon.eta()>2.4: continue
        selected_muons_embed+=[muon]

    selected_muons_sim=[]
    for muon in muonsSim:
        if not muon.isLooseMuon(): continue
        if muon.pt()<3 or muon.eta()>2.4: continue
        selected_muons_sim+=[muon]

    #print("N muons:", len(selected_muons_embed))

    selected_muons_embed_tight=[]
    for muon in muonsEmbed:
        if not muon.isTightMuon(pvEmbed): continue
        if muon.pt()<3 or muon.eta()>2.4: continue
        selected_muons_embed_tight+=[muon]

    selected_muons_embed_med=[]
    for muon in muonsEmbed:
        if not muon.isMediumMuon(): continue
        if muon.pt()<3 or muon.eta()>2.4: continue
        selected_muons_embed_med+=[muon]

    #print("N muons:", len(selected_muons_embed), len(selected_muons_sim), len(selected_gentaumus))
    
    selected_taus=[]
    for tau in taus:
        #print(tau.pt())
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding"):
        #if tau.pt()>10 and abs(tau.eta())<2.3:
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") > 0.5 and tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"):
            selected_taus+=[tau]

    selected_taus_deep=[]
    for tau in taus:
        #print(tau.pt())
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding"):
        if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") > 0.5 and tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
            selected_taus_deep+=[tau]

    selected_taus_embed=[]
    for tau in tausEmbed:
        #print("embed:",tau.pt())
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding"):
        #if tau.pt()>10 and abs(tau.eta())<2.3:
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") and tau.tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") > 0.5 and tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"):
            selected_taus_embed+=[tau]

    selected_taus_embed_deep=[]
    for tau in tausEmbed:
        #print("embed:",tau.pt())
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFindingNewDMs") and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("byVLooseDeepTau2017v2p1VSjet"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding"):
        if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") and tau.tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"):
        #if tau.pt()>10 and abs(tau.eta())<2.3 and tau.tauID("decayModeFinding") > 0.5 and tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"):
            selected_taus_embed_deep+=[tau]

    if len(selected_muons) > 1:
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass())

        h["mmumu"].Fill((mu1+mu2).M())
        h["dRmumu"].Fill(mu1.DeltaR(mu2))
        if mu1.Pt()>mu2.Pt():
            h["mu1Pt"].Fill(mu1.Pt(), genweight)
            h["mu2Pt"].Fill(mu2.Pt(), genweight)
        else:
            h["mu1Pt"].Fill(mu2.Pt(), genweight)
            h["mu2Pt"].Fill(mu1.Pt(), genweight)

    if len(selected_taus) > 1:
        #print("Simulated:", len(selected_taus))
        if selected_taus[0].charge() * selected_taus[1].charge() > 0: continue
        
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_taus[1].pt(), selected_taus[1].eta(), selected_taus[1].phi(), selected_taus[1].mass())

        if mu1.DeltaR(mu2) > 0.8: continue
        h["mtautau"].Fill((mu1+mu2).M(), genweight)
        if mu1.Pt()>mu2.Pt():
            h["tau1Pt"].Fill(mu1.Pt(), genweight)
            h["tau2Pt"].Fill(mu2.Pt(), genweight)
        else:
            h["tau1Pt"].Fill(mu2.Pt(), genweight)
            h["tau2Pt"].Fill(mu1.Pt(), genweight)

    if len(selected_taus_deep) > 1:
        if selected_taus_deep[0].charge() * selected_taus_deep[1].charge() > 0: continue
        
        mu1=ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selected_taus_deep[0].pt(), selected_taus_deep[0].eta(), selected_taus_deep[0].phi(), selected_taus_deep[0].mass()) 
    
        mu2=ROOT.TLorentzVector()
        mu2.SetPtEtaPhiM(selected_taus_deep[1].pt(), selected_taus_deep[1].eta(), selected_taus_deep[1].phi(), selected_taus_deep[1].mass())

        h["mtautau_deep"].Fill((mu1+mu2).M(), genweight)
        if mu1.Pt()>mu2.Pt():
            h["tau1Pt_deep"].Fill(mu1.Pt(), genweight)
            h["tau2Pt_deep"].Fill(mu2.Pt(), genweight)
        else:
            h["tau1Pt_deep"].Fill(mu2.Pt(), genweight)
            h["tau2Pt_deep"].Fill(mu1.Pt(), genweight)

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

                if mu1.DeltaR(mu2) < drMin:
                    mass = (mu1+mu2).M()
                    pt1 = mu1.Pt()
                    pt2 = mu2.Pt()
                    drMin = mu1.DeltaR(mu2)

        if drMin > 0.8: continue
        h["mtautauEmbed"].Fill(mass, genweight)
        h["dRtautauEmbed"].Fill(drMin, genweight)
        if pt1>pt2:
            h["tau1PtEmbed"].Fill(pt1, genweight)
            h["tau2PtEmbed"].Fill(pt2, genweight)
        else:
            h["tau1PtEmbed"].Fill(pt2, genweight)
            h["tau2PtEmbed"].Fill(pt1, genweight)
        

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

                if mu1.DeltaR(mu2) < drMin:
                    mass = (mu1+mu2).M()
                    pt1 = mu1.Pt()
                    pt2 = mu2.Pt()
                    drMin = mu1.DeltaR(mu2)

        if drMin > 0.8: continue
        h["mtautauEmbed_deep"].Fill(mass, genweight)
        h["dRtautauEmbed_deep"].Fill(drMin, genweight)
        #if genTauM < 10.4 and genTauM > 0:
        if genTau1Pt > 50 and genTau2Pt > 20:
            h["mtautauEmbed_deep_cor"].Fill(mass, genweight*drWeight)
            h["dRtautauEmbed_deep_cor"].Fill(drMin, genweight*drWeight)
        if pt1>pt2:
            h["tau1PtEmbed_deep"].Fill(pt1, genweight)
            h["tau2PtEmbed_deep"].Fill(pt2, genweight)
            if genTau1Pt > 50 and genTau2Pt > 20:
                h["tau1PtEmbed_deep_cor"].Fill(pt1, genweight*drWeight)
                h["tau2PtEmbed_deep_cor"].Fill(pt2, genweight*drWeight)
        else:
            h["tau1PtEmbed_deep"].Fill(pt2, genweight)
            h["tau2PtEmbed_deep"].Fill(pt1, genweight)
            if genTau1Pt > 50 and genTau2Pt > 20:
                h["tau1PtEmbed_deep_cor"].Fill(pt1, genweight*drWeight)
                h["tau2PtEmbed_deep_cor"].Fill(pt2, genweight*drWeight)
            
        #h["mtautauEmbed_deep_cor"].Fill(mass, genweight)
        #h["tau1PtEmbed_deep_cor"].Fill(pt1, genweight)
        #h["tau2PtEmbed_deep_cor"].Fill(pt2, genweight)

    if len(selected_muons_embed) > 0 and len(selected_electrons_embed) > 0:
        drMin = 9999
        pt1 = -1
        pt2 = -1
        mass = -1
        for i, tau1 in enumerate(selected_muons_embed):
            for j, tau2 in enumerate(selected_electrons_embed):
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
        h["memuEmbed"].Fill(mass, genweight)
        h["dRemuEmbed"].Fill(drMin, genweight)
        h["ePtEmbed"].Fill(pt2, genweight)
        h["muPtEmbed"].Fill(pt1, genweight)
        if genTau1Pt > 50 and genTau2Pt > 20:
            h["memuEmbed_cor"].Fill(mass, genweight*drWeight)
            h["dRemuEmbed_cor"].Fill(drMin, genweight*drWeight)
            h["ePtEmbed_cor"].Fill(pt2, genweight*drWeight)
            h["muPtEmbed_cor"].Fill(pt1, genweight*drWeight)

    if len(selected_muons_embed) > 1:
        drMin = 9999
        pt1 = -1
        pt2 = -1
        mass = -1
        for i, tau1 in enumerate(selected_muons_embed):
            for j, tau2 in enumerate(selected_muons_embed):
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
        h["mmumuEmbed"].Fill(mass, genweight)
        h["mmumuEmbed_genM"].Fill(genTauM, genweight)
        h["dRmumuEmbed"].Fill(drMin, genweight)
        #if genTauM < 10.4 and genTauM > 0:
        if genTau1Pt > 50 and genTau2Pt > 20:
            h["mmumuEmbed_cor"].Fill(mass, genweight*drWeight)
            h["dRmumuEmbed_cor"].Fill(drMin, genweight*drWeight)
        if pt2 > pt1:
            h["mu1PtEmbed"].Fill(pt2, genweight)
            h["mu2PtEmbed"].Fill(pt1, genweight)
            if genTau1Pt > 50 and genTau2Pt > 20:
                h["mu1PtEmbed_cor"].Fill(pt2, genweight*drWeight)
                h["mu2PtEmbed_cor"].Fill(pt1, genweight*drWeight)
        else:
            h["mu1PtEmbed"].Fill(pt1, genweight)
            h["mu2PtEmbed"].Fill(pt2, genweight)
            if genTau1Pt > 50 and genTau2Pt > 20:
                h["mu1PtEmbed_cor"].Fill(pt2, genweight*drWeight)
                h["mu2PtEmbed_cor"].Fill(pt1, genweight*drWeight)

    if len(selected_muons_embed_tight) > 1:
        drMin = 9999
        pt1 = -1
        pt2 = -1
        mass = -1
        for i, tau1 in enumerate(selected_muons_embed_tight):
            for j, tau2 in enumerate(selected_muons_embed_tight):
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
        if genTau1Pt > 50 and genTau2Pt > 20:
            h["mmumuEmbed_tight"].Fill(mass, genweight)
            h["dRmumuEmbed_tight"].Fill(drMin, genweight)
            if pt2 > pt1:
                h["mu1PtEmbed_tight"].Fill(pt2, genweight)
                h["mu2PtEmbed_tight"].Fill(pt1, genweight)
            else:
                h["mu1PtEmbed_tight"].Fill(pt1, genweight)
                h["mu2PtEmbed_tight"].Fill(pt2, genweight)

    if len(selected_muons_embed_med) > 1:
        drMin = 9999
        pt1 = -1
        pt2 = -1
        mass = -1
        for i, tau1 in enumerate(selected_muons_embed_med):
            for j, tau2 in enumerate(selected_muons_embed_med):
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
        if genTau1Pt > 50 and genTau2Pt > 20:
            h["mmumuEmbed_med"].Fill(mass, genweight)
            h["dRmumuEmbed_med"].Fill(drMin, genweight)
            if pt2 > pt1:
                h["mu1PtEmbed_med"].Fill(pt2, genweight)
                h["mu2PtEmbed_med"].Fill(pt1, genweight)
            else:
                h["mu1PtEmbed_med"].Fill(pt1, genweight)
                h["mu2PtEmbed_med"].Fill(pt2, genweight)

    if len(selected_muons_sim) > 1:
        drMin = 9999
        pt1 = -1
        pt2 = -1
        mass = -1
        for i, tau1 in enumerate(selected_muons_sim):
            for j, tau2 in enumerate(selected_muons_sim):
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
        h["mmumuSim"].Fill(mass, genweight)
        h["dRmumuSim"].Fill(drMin, genweight)
        if pt2 > pt1:
            h["mu1PtSim"].Fill(pt2, genweight)
            h["mu2PtSim"].Fill(pt1, genweight)
        else:
            h["mu1PtSim"].Fill(pt1, genweight)
            h["mu2PtSim"].Fill(pt2, genweight) 
        
#        del events
#        gc.collect()
#    except:
#        nfail += 1
#        print("Error...")

print("Nfail:", nfail)
out.cd()
for hist in sorted(h.keys()):
    h[hist].Write()
