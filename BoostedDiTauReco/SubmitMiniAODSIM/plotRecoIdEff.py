##TODO Next

#Implement impact parameter
 

import ROOT,sys,os
from DataFormats.FWLite import Events, Handle
from looseElectron import *
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

isCheckVar = False
isTCP = True

inputFileListName=sys.argv[1]

inputFileList=inputFileListName

if len(sys.argv)>2:
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./plots/"
outputFileName = outputFileDir+"h_"+inputFileListName.split("/")[-1].replace(".txt","_eff.root")
print outputFileName

out=ROOT.TFile.Open(outputFileName,'recreate')

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

# Histograms

h={}

h['NEvent'] = ROOT.TH1F ("NEvent", "", 2, 0, 2)
h['NElectron'] = ROOT.TH1F ("NElectron", "", 6, 0, 6)
h['NElectron_matched'] = ROOT.TH1F ("NElectron_matched", "", 4, 0, 4)
h['NElectronID'] = ROOT.TH1F ("NElectronID", "", 4, 0, 4)
h['NElectronID_matched'] = ROOT.TH1F ("NElectronID_matched", "", 4, 0, 4)
h['NElectronIso'] = ROOT.TH1F ("NElectronIso", "", 4, 0, 4)

h['ePt'] = ROOT.TH1F("ePt", "", 500, 0, 500)
h['ePt_veto'] = ROOT.TH1F("ePt_veto", "", 500, 0, 500)
h['ePt_loose'] = ROOT.TH1F("ePt_loose", "", 500, 0, 500)
h['ePt_medium'] = ROOT.TH1F("ePt_medium", "", 500, 0, 500)
h['ePt_tight'] = ROOT.TH1F("ePt_tight", "", 500, 0, 500)
h['matched_ePt'] = ROOT.TH1F("matched_ePt", "", 500, 0, 500)
h['matched_ePt_veto'] = ROOT.TH1F("matched_ePt_veto", "", 500, 0, 500)
h['matched_ePt_loose'] = ROOT.TH1F("matched_ePt_loose", "", 500, 0, 500)
h['matched_ePt_medium'] = ROOT.TH1F("matched_ePt_medium", "", 500, 0, 500)
h['matched_ePt_tight'] = ROOT.TH1F("matched_ePt_tight", "", 500, 0, 500)

h['gen_ePt'] = ROOT.TH1F("gen_ePt", "", 500, 0, 500)
h['matched_gen_ePt'] = ROOT.TH1F("matched_gen_ePt", "", 500, 0, 500)
h['matched_gen_ePt_veto'] = ROOT.TH1F("matched_gen_ePt_veto", "", 500, 0, 500)
h['matched_gen_ePt_loose'] = ROOT.TH1F("matched_gen_ePt_loose", "", 500, 0, 500)
h['matched_gen_ePt_medium'] = ROOT.TH1F("matched_gen_ePt_medium", "", 500, 0, 500)
h['matched_gen_ePt_tight'] = ROOT.TH1F("matched_gen_ePt_tight", "", 500, 0, 500)

h['EE_sigmaIetaIeta'] = ROOT.TH1F ("EE_sigmaIetaIeta", "", 40, 0, 20E-3)
h['EE_dEtaSeed'] = ROOT.TH1F ("EE_dEtaSeed", "", 30, -15E-3, 15E-3)
h['EE_dPhiIn'] = ROOT.TH1F ("EE_dPhiIn", "", 20, -0.1, 0.1)
h['EE_hOverE'] = ROOT.TH1F ("EE_hOverE", "", 40, 0, 0.2)
h['EE_relIso'] = ROOT.TH1F ("EE_relIso", "", 50, 0, 0.5)
h['EE_eInvpInv'] = ROOT.TH1F ("EE_eInvpInv", "", 40, 0, 0.2)
h['EE_missHits'] = ROOT.TH1F ("EE_missHits", "", 5, -0.5, 4.5)               
h['EE_dz'] = ROOT.TH1F ("EE_dz", "", 40,-0.2, 0.2)
h['EE_dxdy'] = ROOT.TH1F ("EE_dxdy", "", 40,-0.2, 0.2)

h['EB_sigmaIetaIeta'] = ROOT.TH1F ("EB_sigmaIetaIeta", "", 40, 0, 20E-3)
h['EB_dEtaSeed'] = ROOT.TH1F ("EB_dEtaSeed", "", 30, -15E-3, 15E-3)
h['EB_dPhiIn'] = ROOT.TH1F ("EB_dPhiIn", "", 20, -0.1, 0.1)
h['EB_hOverE'] = ROOT.TH1F ("EB_hOverE", "", 40, 0, 0.2)
h['EB_relIso'] = ROOT.TH1F ("EB_relIso", "", 50, 0, 0.5)
h['EB_eInvpInv'] = ROOT.TH1F ("EB_eInvpInv", "", 40, 0, 0.2)
h['EB_missHits'] = ROOT.TH1F ("EB_missHits", "", 5, -0.5, 4.5)               
h['EB_dz'] = ROOT.TH1F ("EB_dz", "", 40,-0.2, 0.2)
h['EB_dxdy'] = ROOT.TH1F ("EB_dxdy", "", 40,-0.2, 0.2)

workingPoint={
    'barrel':{
        'veto':{'sigmaIEIE':0.0115,'dEtaSeed':0.00749,'dPhi':0.228,'invEinvP':0.299,'miss':2},
        'loose':{'sigmaIEIE':0.011,'dEtaSeed':0.00477,'dPhi':0.222,'invEinvP':0.241,'miss':1},
        'medium':{'sigmaIEIE':0.00998,'dEtaSeed':0.00311,'dPhi':0.103,'invEinvP':0.134,'miss':1},
        'tight':{'sigmaIEIE':0.00998,'dEtaSeed':0.00308,'dPhi':0.0816,'invEinvP':0.0129,'miss':1}
    },
    'endcap':{
        'veto':{'sigmaIEIE':0.037,'dEtaSeed':0.00895,'dPhi':0.213,'invEinvP':0.15,'miss':3},
        'loose':{'sigmaIEIE':0.0314,'dEtaSeed':0.00868,'dPhi':0.213,'invEinvP':0.14,'miss':1},
        'medium':{'sigmaIEIE':0.0298,'dEtaSeed':0.00609,'dPhi':0.045,'invEinvP':0.13,'miss':1},
        'tight':{'sigmaIEIE':0.0292,'dEtaSeed':0.00605,'dPhi':0.0394,'invEinvP':0.0129,'miss':1}
    }
}

# -----------

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    if isTCP == True:
#        print "root://cmseos.fnal.gov//store/user/zhangj/events/ALP/RunIISummer17DR94Premix/"+inputFileName.replace("\n","")
#        inputFileName="root://cmseos.fnal.gov//store/user/zhangj/events/ALP/RunIISummer17DR94Premix/"+inputFileName.replace("\n","")
        print inputFileName.replace("\n","")
        inputFileName=inputFileName.replace("\n","")
    else:
        print "root://cmxrootd.fnal.gov/"+inputFileName.replace("\n","")
        inputFileName="root://cmsxrootd.fnal.gov/"+inputFileName.replace("\n","")

    f=ROOT.TFile.Open(inputFileName)

    if not f.IsZombie():
        events=Events(inputFileName)
    else:
        print "Can't Open File: "+inputFileName
        continue

    for event in events:

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

        genElectrons = []

        for particle in particles:
            if isTCP == True: 
                if abs(particle.pdgId())==11 and particle.isDirectHardProcessTauDecayProductFinalState(): #TCP
                    genElectrons+=[particle]
            else: 
                if abs(particle.pdgId())==11 and particle.isHardProcess(): #DY
                    genElectrons+=[particle]
                

        selected_genElectrons = []
        for electron in genElectrons:
            if electron.pt()>7 and abs(electron.eta())<2.5: 
                selected_genElectrons+=[electron]
                h['gen_ePt'].Fill(electron.pt(), genweight)

        selected_genElectrons.sort(key=lambda x: x.pt(), reverse=True)
#        if len(selected_genElectrons) > 0: h['gen_ePt'].Fill(selected_genElectrons[0].pt(), genweight)

        selected_recoElectrons = []   
        selected_electrons = [[],[],[],[]]
        selected_IsoElectrons = [[],[],[],[]]
        for electron in electrons: 
            h['NElectron'].Fill(0.5, genweight)
            if electron.pt()>7 and abs(electron.eta())<2.5:
                selected_recoElectrons+=[electron]
                E_c = electron.superCluster().energy()
                h['NElectron'].Fill(1.5, genweight)
                h['ePt'].Fill(electron.pt(), genweight)

                if isCheckVar==True:
                    if electron.isEB():
                        h['EB_sigmaIetaIeta'].Fill(electron.full5x5_sigmaIetaIeta(), genweight)
                        h['EB_dEtaSeed'].Fill(dEtaInSeed(electron), genweight)
                        h['EB_dPhiIn'].Fill(electron.deltaPhiSuperClusterTrackAtVtx(), genweight)
                        h['EB_hOverE'].Fill(electron.hadronicOverEm(), genweight)
                        h['EB_relIso'].Fill(GsfEleEffAreaPFIsoCut(electron, rho), genweight)
                        h['EB_eInvpInv'].Fill(GsfEleEInverseMinusPInverse(electron), genweight)
                        h['EB_missHits'].Fill(GsfEleMissingHitsCut(electron), genweight)   
                        h['EB_dz'].Fill(electron.gsfTrack().dz(pv), genweight)             
                        h['EB_dxdy'].Fill(electron.gsfTrack().dxy(pv), genweight)
                    if electron.isEE():
                        h['EE_sigmaIetaIeta'].Fill(electron.full5x5_sigmaIetaIeta(), genweight)
                        h['EE_dEtaSeed'].Fill(dEtaInSeed(electron), genweight)
                        h['EE_dPhiIn'].Fill(electron.deltaPhiSuperClusterTrackAtVtx(), genweight)
                        h['EE_hOverE'].Fill(electron.hadronicOverEm(), genweight)
                        h['EE_relIso'].Fill(GsfEleEffAreaPFIsoCut(electron, rho), genweight)
                        h['EE_eInvpInv'].Fill(GsfEleEInverseMinusPInverse(electron), genweight)
                        h['EE_missHits'].Fill(GsfEleMissingHitsCut(electron), genweight)   
                        h['EE_dz'].Fill(electron.gsfTrack().dz(pv), genweight)             
                        h['EE_dxdy'].Fill(electron.gsfTrack().dxy(pv), genweight)

                if electron.isEB():
                    loc = "barrel"
                    for wp in list(workingPoint[loc]):
                        if wp == "veto":
                            workingPoint[loc][wp]['hOvE'] = 0.356
                            workingPoint[loc][wp]['relIso'] = 0.175
                            n = 0.5
                            w = 0
                        if wp == "loose": 
                            workingPoint[loc][wp]['hOvE'] = 0.298
                            workingPoint[loc][wp]['relIso'] = 0.0944
                            n = 1.5
                            w = 1
                        if wp == "medium": 
                            workingPoint[loc][wp]['hOvE'] = 0.253
                            workingPoint[loc][wp]['relIso'] = 0.0695
                            n = 2.5
                            w = 2
                        if wp == "tight": 
                            workingPoint[loc][wp]['hOvE'] = 0.0414
                            workingPoint[loc][wp]['relIso'] = 0.0588
                            n = 3.5
                            w = 3

                        if electron.full5x5_sigmaIetaIeta() < workingPoint[loc][wp]['sigmaIEIE'] \
                        and abs(dEtaInSeed(electron)) < workingPoint[loc][wp]['dEtaSeed'] \
                        and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < workingPoint[loc][wp]['dPhi'] \
                        and GsfEleEInverseMinusPInverse(electron) < workingPoint[loc][wp]['invEinvP'] \
                        and GsfEleMissingHitsCut(electron) <= workingPoint[loc][wp]['miss'] \
                        and electron.hadronicOverEm() < workingPoint[loc][wp]['hOvE'] \
                        and electron.passConversionVeto() \
                        and abs(electron.gsfTrack().dz(pv))<0.1 \
                        and abs(electron.gsfTrack().dxy(pv))<0.05:
                            selected_electrons[w]+=[electron]
                            h['NElectronID'].Fill(n, genweight)

                if electron.isEE():
                    loc = "endcap"
                    for wp in list(workingPoint[loc]):
                        if wp == "veto":
                            workingPoint[loc][wp]['hOvE'] = 0.211
                            workingPoint[loc][wp]['relIso'] = 0.159
                            n = 0.5
                            w = 0
                        if wp == "loose": 
                            workingPoint[loc][wp]['hOvE'] = 0.101
                            workingPoint[loc][wp]['relIso'] = 0.107
                            n = 1.5
                            w = 1
                        if wp == "medium": 
                            workingPoint[loc][wp]['hOvE'] = 0.0878
                            workingPoint[loc][wp]['relIso'] = 0.0821
                            n = 2.5
                            w = 2
                        if wp == "tight": 
                            workingPoint[loc][wp]['hOvE'] = 0.0641
                            workingPoint[loc][wp]['relIso'] = 0.0571
                            n = 3.5
                            w = 3

                        if electron.full5x5_sigmaIetaIeta() < workingPoint[loc][wp]['sigmaIEIE'] \
                        and abs(dEtaInSeed(electron)) < workingPoint[loc][wp]['dEtaSeed'] \
                        and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < workingPoint[loc][wp]['dPhi'] \
                        and GsfEleEInverseMinusPInverse(electron) < workingPoint[loc][wp]['invEinvP'] \
                        and GsfEleMissingHitsCut(electron) <= workingPoint[loc][wp]['miss'] \
                        and electron.hadronicOverEm() < workingPoint[loc][wp]['hOvE'] \
                        and electron.passConversionVeto() \
                        and abs(electron.gsfTrack().dz(pv))<0.2 \
                        and abs(electron.gsfTrack().dxy(pv))<0.1:
                            selected_electrons[w]+=[electron]
                            h['NElectronID'].Fill(n, genweight)

        for i in range(4):
            selected_electrons[i].sort(key=lambda x: x.pt(), reverse=True)
            selected_IsoElectrons[i].sort(key=lambda x: x.pt(), reverse=True)

        if len(selected_electrons[0])>0: h['ePt_veto'].Fill(selected_electrons[0][0].pt(), genweight)
        if len(selected_electrons[1])>0: h['ePt_loose'].Fill(selected_electrons[1][0].pt(), genweight)
        if len(selected_electrons[2])>0: h['ePt_medium'].Fill(selected_electrons[2][0].pt(), genweight)
        if len(selected_electrons[3])>0: h['ePt_tight'].Fill(selected_electrons[3][0].pt(), genweight)

        selected_recoElectrons.sort(key=lambda x: x.pt(), reverse=True)
#        if len(selected_recoElectrons)>0: h['ePt'].Fill(electron.pt(), genweight)

        if len(selected_recoElectrons)>0 and len(selected_genElectrons)>0:
            lE = selected_recoElectrons[0]
            lGenE = selected_genElectrons[0]
            if lE.charge() == lGenE.charge():
                genE = ROOT.TLorentzVector()
                genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                e = ROOT.TLorentzVector()
                e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                if genE.DeltaR(e) < 0.1:
                    h['matched_ePt'].Fill(e.Pt(), genweight)
                    h['matched_gen_ePt'].Fill(genE.Pt(), genweight)
            if len(selected_recoElectrons)>1 and len(selected_genElectrons)>1:
                lE = selected_recoElectrons[1]
                lGenE = selected_genElectrons[1]
                if lE.charge() == lGenE.charge():
                    genE = ROOT.TLorentzVector()
                    genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                    e = ROOT.TLorentzVector()
                    e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                    if genE.DeltaR(e) < 0.1:
                        h['matched_ePt'].Fill(e.Pt(), genweight)
                        h['matched_gen_ePt'].Fill(genE.Pt(), genweight)

                    # h['NElectron_matched'].Fill(0.5, genweight)
                    # h['NElectron_matched'].Fill(1.5, genweight)
                    # h['NElectron_matched'].Fill(2.5, genweight)
                    # h['NElectron_matched'].Fill(3.5, genweight)

        if len(selected_electrons[0])>0 and len(selected_genElectrons)>0:
            lE = selected_electrons[0][0]
            lGenE = selected_genElectrons[0]
            if lE.charge() == lGenE.charge():
                genE = ROOT.TLorentzVector()
                genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                e = ROOT.TLorentzVector()
                e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                if genE.DeltaR(e) < 0.1:
                    h['matched_ePt_veto'].Fill(e.Pt(), genweight)
                    h['matched_gen_ePt_veto'].Fill(genE.Pt(), genweight)
                    h['NElectronID_matched'].Fill(0.5, genweight)
            if len(selected_electrons[0])>1 and len(selected_genElectrons)>1:
                lE = selected_electrons[0][1]
                lGenE = selected_genElectrons[1]
                if lE.charge() == lGenE.charge():
                    genE = ROOT.TLorentzVector()
                    genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                    e = ROOT.TLorentzVector()
                    e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                    if genE.DeltaR(e) < 0.1:
                        h['matched_ePt_veto'].Fill(e.Pt(), genweight)
                        h['matched_gen_ePt_veto'].Fill(genE.Pt(), genweight)
                        h['NElectronID_matched'].Fill(0.5, genweight)

        if len(selected_electrons[1])>0 and len(selected_genElectrons)>0:
            lE = selected_electrons[1][0]
            lGenE = selected_genElectrons[0]
            if lE.charge() == lGenE.charge():
                genE = ROOT.TLorentzVector()
                genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                e = ROOT.TLorentzVector()
                e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                if genE.DeltaR(e) < 0.1:
                    h['matched_ePt_loose'].Fill(e.Pt(), genweight)
                    h['matched_gen_ePt_loose'].Fill(genE.Pt(), genweight)
                    h['NElectronID_matched'].Fill(1.5, genweight)
            if len(selected_electrons[1])>1 and len(selected_genElectrons)>1:
                lE = selected_electrons[1][1]
                lGenE = selected_genElectrons[1]
                if lE.charge() == lGenE.charge():
                    genE = ROOT.TLorentzVector()
                    genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                    e = ROOT.TLorentzVector()
                    e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                    if genE.DeltaR(e) < 0.1:
                        h['matched_ePt_loose'].Fill(e.Pt(), genweight)
                        h['matched_gen_ePt_loose'].Fill(genE.Pt(), genweight)
                        h['NElectronID_matched'].Fill(1.5, genweight)

        if len(selected_electrons[2])>0 and len(selected_genElectrons)>0:
            lE = selected_electrons[2][0]
            lGenE = selected_genElectrons[0]
            if lE.charge() == lGenE.charge():
                genE = ROOT.TLorentzVector()
                genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                e = ROOT.TLorentzVector()
                e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                if genE.DeltaR(e) < 0.1:
                    h['matched_ePt_medium'].Fill(e.Pt(), genweight)
                    h['matched_gen_ePt_medium'].Fill(genE.Pt(), genweight)
                    h['NElectronID_matched'].Fill(2.5, genweight)
            if len(selected_electrons[2])>1 and len(selected_genElectrons)>1:
                lE = selected_electrons[2][1]
                lGenE = selected_genElectrons[1]
                if lE.charge() == lGenE.charge():
                    genE = ROOT.TLorentzVector()
                    genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                    e = ROOT.TLorentzVector()
                    e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                    if genE.DeltaR(e) < 0.1:
                        h['matched_ePt_medium'].Fill(e.Pt(), genweight)
                        h['matched_gen_ePt_medium'].Fill(genE.Pt(), genweight)
                        h['NElectronID_matched'].Fill(2.5, genweight)

        if len(selected_electrons[3])>0 and len(selected_genElectrons)>0:
            lE = selected_electrons[3][0]
            lGenE = selected_genElectrons[0]
            if lE.charge() == lGenE.charge():
                genE = ROOT.TLorentzVector()
                genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                e = ROOT.TLorentzVector()
                e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                if genE.DeltaR(e) < 0.1:
                    h['matched_ePt_tight'].Fill(e.Pt(), genweight)
                    h['matched_gen_ePt_tight'].Fill(genE.Pt(), genweight)
                    h['NElectronID_matched'].Fill(3.5, genweight)
            if len(selected_electrons[3])>1 and len(selected_genElectrons)>1:
                lE = selected_electrons[3][1]
                lGenE = selected_genElectrons[1]
                if lE.charge() == lGenE.charge():
                    genE = ROOT.TLorentzVector()
                    genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                    e = ROOT.TLorentzVector()
                    e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                    if genE.DeltaR(e) < 0.1:
                        h['matched_ePt_tight'].Fill(e.Pt(), genweight)
                        h['matched_gen_ePt_tight'].Fill(genE.Pt(), genweight)
                        h['NElectronID_matched'].Fill(3.5, genweight)



out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
