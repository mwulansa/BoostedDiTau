import ROOT, sys, math, gc, os
from DataFormats.FWLite import Events, Handle
import numpy as np

from looseElectron import *

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('slimmedJets')

handleTaus = Handle ('vector<pat::Tau>')
labelTaus = ('slimmedTaus')

handleElectronCleanedTaus = Handle ('vector<pat::Tau>')
labelElectronCleanedTaus = ('slimmedTausElectronCleaned', '', 'PAT')

handleMuonCleanedTaus = Handle ('vector<pat::Tau>')
labelMuonCleanedTaus = ('slimmedTausMuonCleaned', '', 'PAT')

handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlineSlimmedPrimaryVertices')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('genParticles')

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

# book histograms
h={}
h['NEvent'] = ROOT.TH1F ("NEvent", "", 2, 0, 2)
h['NElectrons'] = ROOT.TH1F ("NElectrons", "", 4, 0, 4)
h['NIsoElectrons'] = ROOT.TH1F ("NIsoElectrons", "", 4, 0, 4)
h['ePt'] = ROOT.TH1F("ePt", "", 500, 0, 500)
h['isoEPt'] = ROOT.TH1F("isoEPt", "", 500, 0, 500)

h['gen_ePt'] = ROOT.TH1F ('gen_ePt', '', 500, 0, 500)

h['gen_matched_ePt'] = ROOT.TH1F ('gen_matched_ePt', '', 500, 0, 500)


prefix="root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer17DR94Premix/"

#mass = "50_120"
mass = "10"

out=ROOT.TFile("h_plot_electronIds_m_"+mass+"_94X.root",'recreate')

#searchString = "/ZToEE*"+mass+"/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"

#os.system('dasgoclient --query "'+searchString+'"')
#query = 'dasgoclient --query "file dataset='+searchString+'"'
#print query
#files=os.popen(query).read().split()
#print files
#sys.exit()

#files = ["/store/mc/RunIIFall17MiniAODv2/ZToEE_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/12E59949-4C44-E811-A054-0CC47A4C8E56.root"]


jobs=np.linspace(100, 1, 100)
files = []
for job in jobs:
    filename = "ALP_m"+mass+"_w1_htjmin400_RunIISummer17DR94Premix_MINIAODSIM_Cleaned_"+str(int(job))+".root"
    files+=[filename]

ntot=0
for fil in files:

    print prefix+fil
    events=Events(prefix+fil)

    for event in events:
    
        ntot+=1
        #if not ntot==879: continue
        #print "Begin processing Event", ntot
            
        event.getByLabel(labelVertex, handleVertex)
        vertex=handleVertex.product()
        pv = vertex[0].position()
    
        event.getByLabel(labelElectron, handleElectron)
        electrons=handleElectron.product()
    
        event.getByLabel(labelRho, handleRho)
        if len(handleRho.product())>0:
            rho = handleRho.product()[0]
        else:
            rho = 0
    
        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()
    
        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()

        h['NEvent'].Fill(0.5, 1)
        h['NEvent'].Fill(1.5, genweight)

        #gen electrons
        gen_electrons = []
        for particle in particles:
            if abs(particle.pdgId()) == 11 and particle.isDirectHardProcessTauDecayProductFinalState():
                gen_electrons+=[particle]
    
        #if not len(gen_electrons) == 2:
        #    print "ERROR!", len(gen_electrons)
    
        selected_gen_electrons = []
        for electron in gen_electrons:
            if electron.pt()>7 and abs(electron.eta())<2.5:
                selected_gen_electrons+=[electron]

        selected_gen_electrons.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_gen_electrons)>0:
            h['gen_ePt'].Fill(selected_gen_electrons[0].pt(), genweight)
    
        #Electron Selection
        selected_electrons=[]
        if electrons.size()>0:
            for electron in electrons:
                h['NElectrons'].Fill(0.5, genweight)
                if electron.pt()>7 and abs(electron.eta())<2.5:
                    h['NElectrons'].Fill(1.5, genweight)
                    E_c = electron.superCluster().energy()
                    if electron.isEB():
                        if electron.full5x5_sigmaIetaIeta()<0.0112 \
                        and GsfEleEInverseMinusPInverse(electron)<0.193 \
                        and abs(dEtaInSeed(electron))<0.00377     \
                        and GsfEleMissingHitsCut(electron)<=1 \
                        and electron.passConversionVeto() \
                        and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.0884 \
                        and electron.hadronicOverEm()< (0.05 + 1.16/E_c + 0.0324*rho/E_c):
                            h['NElectrons'].Fill(2.5, genweight)
                            if abs(electron.gsfTrack().dz(pv))<0.1 \
                            and abs(electron.gsfTrack().dxy(pv))<0.05:
                                h['NElectrons'].Fill(3.5, genweight)
                                selected_electrons+=[electron]
                                
                    if electron.isEE():
                        if electron.full5x5_sigmaIetaIeta() < 0.0425 \
                        and GsfEleEInverseMinusPInverse(electron) < 0.111 \
                        and abs(dEtaInSeed(electron)) < 0.00674 \
                        and GsfEleMissingHitsCut(electron)<=1 \
                        and electron.passConversionVeto() \
                        and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.169 \
                        and electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c):
                            h['NElectrons'].Fill(2.5, genweight)
                            if abs(electron.gsfTrack().dz(pv)) < 0.2 \
                            and abs(electron.gsfTrack().dxy(pv)) < 0.1:
                                h['NElectrons'].Fill(3.5, genweight)
                                selected_electrons+=[electron]
        selected_electrons.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_electrons) > 0:
            h['ePt'].Fill(selected_electrons[0].pt(), genweight)
    
        #Isolated electron Selection
        selected_electrons_iso=[]
        if electrons.size()>0:
            for electron in electrons:
                h['NIsoElectrons'].Fill(0.5, genweight)
                if electron.pt()>7 and abs(electron.eta())<2.5:
                    h['NIsoElectrons'].Fill(1.5, genweight)
                    E_c = electron.superCluster().energy()
                    if electron.isEB():
                        if electron.full5x5_sigmaIetaIeta()<0.0112 \
                        and GsfEleEInverseMinusPInverse(electron)<0.193 \
                        and abs(electron.deltaPhiSuperClusterTrackAtVtx())<0.0884 \
                        and abs(dEtaInSeed(electron))<0.00377     \
                        and GsfEleMissingHitsCut(electron)<=1 \
                        and GsfEleEffAreaPFIsoCut(electron, rho) < 0.112+0.506/electron.pt() \
                        and electron.passConversionVeto() \
                        and electron.hadronicOverEm()< (0.05 + 1.16/E_c + 0.0324*rho/E_c) :
                            h['NIsoElectrons'].Fill(2.5, genweight)
                            if abs(electron.gsfTrack().dz(pv))<0.1 \
                            and abs(electron.gsfTrack().dxy(pv))<0.05:
                                h['NIsoElectrons'].Fill(3.5, genweight)
                                selected_electrons_iso+=[electron]
                                
                    if electron.isEE():
                        if electron.full5x5_sigmaIetaIeta() < 0.0425 \
                        and GsfEleEInverseMinusPInverse(electron) < 0.111 \
                        and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.169 \
                        and abs(dEtaInSeed(electron)) < 0.00674 \
                        and GsfEleMissingHitsCut(electron)<=1 \
                        and GsfEleEffAreaPFIsoCut(electron, rho) < 0.108+0.963/electron.pt() \
                        and electron.passConversionVeto() \
                        and electron.hadronicOverEm() < (0.0441 + 2.54/E_c + 0.183*rho/E_c) :
                            h['NIsoElectrons'].Fill(2.5, genweight)
                            if abs(electron.gsfTrack().dz(pv)) < 0.2 \
                            and abs(electron.gsfTrack().dxy(pv)) < 0.1:
                                h['NIsoElectrons'].Fill(3.5, genweight)
                                selected_electrons_iso+=[electron]
        selected_electrons.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_electrons_iso) > 0:
            h['isoEPt'].Fill(selected_electrons_iso[0].pt(), genweight)

        #gen match
        if len(selected_gen_electrons)>0 and len(selected_electrons)>0:
            leading_gen_e = selected_gen_electrons[0]
            leading_e = selected_electrons[0]
            if leading_gen_e.charge() == leading_e.charge():
                genE = ROOT.TLorentzVector()
                genE.SetPtEtaPhiM(leading_gen_e.pt(), leading_gen_e.eta(), leading_gen_e.phi(), leading_gen_e.mass())

                e = ROOT.TLorentzVector()
                e.SetPtEtaPhiM(leading_e.pt(), leading_e.eta(), leading_e.phi(), leading_e.mass())
                if genE.DeltaR(e) < 0.1:
                    h['gen_matched_ePt'].Fill(e.Pt(), genweight)
            

out.cd()
for key in h.keys():
    h[key].Write()
out.Close()
