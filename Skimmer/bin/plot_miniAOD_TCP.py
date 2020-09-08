import ROOT,sys, math
from DataFormats.FWLite import Events, Handle
import numpy as np

#from looseElectron import *

jobs=np.linspace(100, 1, 100)
#jobs=[1]

if len(sys.argv)>1:
    mass=sys.argv[1]
else:
    mass='10'
    print "Use default signal mass: 10 GeV."


# book histograms
h={}

h['NEvent'] = ROOT.TH1F ("hNEvent", "", 1, 0, 2)

h['Mu_Pt'] = ROOT.TH1F ("Mu_Pt", "", 500, 0, 500)
#h['Mu_Pt_selected'] = ROOT.TH1F ("Mu_Pt_selected", "", 500, 0, 500)

h['E_Pt'] = ROOT.TH1F ("E_Pt", "", 500, 0, 500)
#h['E_Pt_selected'] = ROOT.TH1F ("E_Pt_selected", "", 500, 0, 500)

h['J_Pt'] = ROOT.TH1F ("J_Pt", "", 2000, 0, 2000)
#h['J_Pt_selected'] = ROOT.TH1F ("J_Pt_selected", "", 2000, 0, 2000)

h['Tau_Pt'] = ROOT.TH1F ("Tau_Pt", "", 500, 0, 500)
#h['Tau_Pt_selected'] = ROOT.TH1F ("Tau_Pt_selected", "", 500, 0, 500)
#h['Tau_Pt_selected_dz'] = ROOT.TH1F ("Tau_Pt_selected_dz", "", 500, 0, 500)
#h['Tau_Pt_selected_ip'] = ROOT.TH1F ("Tau_Pt_selected_ip", "", 500, 0, 500)

h['eTau_Pt'] = ROOT.TH1F ("eTau_Pt", "", 500, 0, 500)
#h['eTau_Pt_selected'] = ROOT.TH1F ("eTau_Pt_selected", "", 500, 0, 500)
#h['eTau_Pt_selected_dz'] = ROOT.TH1F ("eTau_Pt_selected_dz", "", 500, 0, 500)
#h['eTau_Pt_selected_ip'] = ROOT.TH1F ("eTau_Pt_selected_ip", "", 500, 0, 500)

h['mTau_Pt'] = ROOT.TH1F ("mTau_Pt", "", 500, 0, 500)
#h['mTau_Pt_selected'] = ROOT.TH1F ("mTau_Pt_selected", "", 500, 0, 500)

h['M_ETau'] = ROOT.TH1F ("M_ETau", "", 100, 0, 100)
h['M_EeTau'] = ROOT.TH1F ("M_EeTau", "", 100, 0, 100)
h['M_MTau'] = ROOT.TH1F ("M_MTau", "", 100, 0, 100)
h['M_MmTau'] = ROOT.TH1F ("M_MmTau", "", 100, 0, 100)
#h['M_EeTau_ip'] = ROOT.TH1F ("M_EeTau_ip", "", 100, 0, 100)

h['NElectrons'] = ROOT.TH1F ("hNElectrons","",4,0,4)


handleMuon = Handle ('vector<pat::Muon>')
labelMuonSelected = ('selectedPATMuons', 'myMuons', 'PAT')

handleElectron = Handle ('vector<pat::Electron>')
labelElectronSelected = ('selectedPATElectrons', 'myElectrons', 'PAT')

handleJet = Handle ("vector<pat::Jet>")
labelJetSelected = ("selectedPATJets", "", "PAT")

handleTaus = Handle ('vector<pat::Tau>')
labelTausSelected = ('selectedPATTaus')

handleElectronCleanedTaus = Handle ('vector<pat::Tau>')
labelElectronCleanedTausSelected = ('selectedPATTausElectronCleaned')

handleMuonCleanedTaus = Handle ('vector<pat::Tau>')
labelMuonCleanedTausSelected = ('selectedPATTausMuonCleaned')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults')

prefix="root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer17DR94Premix/"

out=ROOT.TFile("ALP_m10_w1_htjmin400_RunIISummer17DR94Premix_hists_Slimmed.root", "recreate")

#prefix = ""

for job in jobs:

    #filename = "ALP_m10_w1_htjmin400_RunIISummer17DR94Premix_MINIAODSIM_Slimmed_noID_100.root"
    filename = "ALP_m10_w1_htjmin400_RunIISummer17DR94Premix_MINIAODSIM_Slimmed_"+str(int(job))+".root"
        
    print prefix+filename
    events=Events(prefix+filename)
    
    ntot=0 
    # loop over events
    for event in events:
        
        ntot+=1
        #print "Begin processing Event", ntot
    
        h['NEvent'].Fill(1)
    
        event.getByLabel(labelMuonSelected, handleMuon)
        muons_selected=handleMuon.product()
    
        event.getByLabel(labelElectronSelected, handleElectron)
        electrons_selected=handleElectron.product()
    
        h['NElectrons'].Fill(0.5, electrons_selected.size())
    
        event.getByLabel(labelJetSelected, handleJet)
        jets_selected=handleJet.product()
        
        event.getByLabel(labelTausSelected, handleTaus)
        taus_selected=handleTaus.product()
        
        event.getByLabel(labelElectronCleanedTausSelected, handleElectronCleanedTaus)
        etaus_selected=handleElectronCleanedTaus.product()
        
        event.getByLabel(labelMuonCleanedTausSelected, handleMuonCleanedTaus)
        mtaus_selected=handleMuonCleanedTaus.product()
    
        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()
    
        event.getByLabel(labelHLT, handleHLT)
        triggerResults=handleHLT.product()
        names = event.object().triggerNames(triggerResults)
    
        genMuons=[]
        genTaus=[]
        genElectrons=[]
        genNutaus=[]
        for particle in particles: 
            if abs(particle.pdgId())==15 and particle.isHardProcess():
                genTaus+=[particle]     
            if abs(particle.pdgId())==11 and particle.isDirectHardProcessTauDecayProductFinalState(): 
                genElectrons+=[particle]
            if abs(particle.pdgId())==13 and particle.isDirectHardProcessTauDecayProductFinalState(): 
                genMuons+=[particle]
            if abs(particle.pdgId())==16 and particle.isDirectHardProcessTauDecayProductFinalState(): 
                genNutaus+=[particle]
    
        #Muon selection
        selected_muons=[]
        if muons_selected.size()>0:
            for muon in muons_selected:
                selected_muons+=[muon]
            selected_muons.sort(key=lambda x: x.pt(), reverse=True)
            if len(selected_muons) > 0:
                h['Mu_Pt'].Fill(selected_muons[0].pt())
    
        #electron selection
        selected_electrons=[]
        if electrons_selected.size()>0:
            for electron in electrons_selected:
                selected_electrons+=[electron]
            selected_electrons.sort(key=lambda x: x.pt(), reverse=True)
            h['NElectrons'].Fill(3.5, len(selected_electrons))
            if len(selected_electrons) > 0:
                h['E_Pt'].Fill(selected_electrons[0].pt())
    
        # jet selection
        selected_jets=[]
        if jets_selected.size()>0:
            for jet in jets_selected:
                selected_jets+=[jet]
            selected_jets.sort(key=lambda x: x.pt(), reverse=True)
            if len(selected_jets) > 0:
                h['J_Pt'].Fill(selected_jets[0].pt())
    
        #tau selection
        selected_taus=[]
        if taus_selected.size()>0:
            for tau in taus_selected:
                if tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
                    #h['Tau_Pt'].Fill(tau.pt())
                    selected_taus+=[tau]
            selected_taus.sort(key=lambda x: x.pt(), reverse=True)
            if len(selected_taus) > 0:
                h['Tau_Pt'].Fill(selected_taus[0].pt())
                
        selected_etaus=[]
        if etaus_selected.size()>0:
            for tau in etaus_selected:
                if tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
                    #h['eTau_Pt'].Fill(tau.pt())
                    selected_etaus+=[tau]
            selected_etaus.sort(key=lambda x: x.pt(), reverse=True)
            if len(selected_etaus) > 0:
                h['eTau_Pt'].Fill(selected_etaus[0].pt())
                
        selected_mtaus=[]
        if mtaus_selected.size()>0:
            for tau in mtaus_selected:
                if tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
                    #h['mTau_Pt'].Fill(tau.pt())
                    selected_mtaus+=[tau]
            selected_mtaus.sort(key=lambda x: x.pt(), reverse=True)
            if len(selected_mtaus) > 0:
                h['mTau_Pt'].Fill(selected_mtaus[0].pt())
    
        if len(selected_etaus)>0 and len(selected_electrons)>0 and len(selected_jets)>0:
            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())
    
            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_etaus[0].pt(), selected_etaus[0].eta(), selected_etaus[0].phi(), selected_etaus[0].mass()) 
                
            if len(selected_jets)==1:
                j0=ROOT.TLorentzVector()
                j0.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j=j0
            else: 
                j1=ROOT.TLorentzVector()
                j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j2=ROOT.TLorentzVector()
                j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())
                if tau.DeltaR(j1)>0.3: 
                    j=j1
                else:
                    j=j2 
                    
            if e.DeltaR(tau)<0.4 and e.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8:
                h['M_EeTau'].Fill((e+tau).M())
    
        if len(selected_taus)>0 and len(selected_electrons)>0 and len(selected_jets)>0:
            e=ROOT.TLorentzVector()
            e.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass())
    
            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass()) 
                
            if len(selected_jets)==1:
                j0=ROOT.TLorentzVector()
                j0.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j=j0
            else: 
                j1=ROOT.TLorentzVector()
                j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j2=ROOT.TLorentzVector()
                j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())
                if tau.DeltaR(j1)>0.3: 
                    j=j1
                else:
                    j=j2 
    
            if e.DeltaR(tau)<0.4 and e.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8:
                h['M_ETau'].Fill((e+tau).M())
    
        if len(selected_taus)>0 and len(selected_muons)>0 and len(selected_jets)>0:
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())
    
            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass()) 
                
            if len(selected_jets)==1:
                j0=ROOT.TLorentzVector()
                j0.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j=j0
            else: 
                j1=ROOT.TLorentzVector()
                j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j2=ROOT.TLorentzVector()
                j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())
                if tau.DeltaR(j1)>0.3: 
                    j=j1
                else:
                    j=j2 
    
            if mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8:
                h['M_MTau'].Fill((mu+tau).M())
    
                
        if len(selected_mtaus)>0 and len(selected_muons)>0 and len(selected_jets)>0:
            mu=ROOT.TLorentzVector()
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())
    
            tau=ROOT.TLorentzVector()
            tau.SetPtEtaPhiM(selected_mtaus[0].pt(), selected_mtaus[0].eta(), selected_mtaus[0].phi(), selected_mtaus[0].mass()) 
                
            if len(selected_jets)==1:
                j0=ROOT.TLorentzVector()
                j0.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j=j0
            else: 
                j1=ROOT.TLorentzVector()
                j1.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass())
                j2=ROOT.TLorentzVector()
                j2.SetPtEtaPhiM(selected_jets[1].pt(), selected_jets[1].eta(), selected_jets[1].phi(), selected_jets[1].mass())
                if tau.DeltaR(j1)>0.3: 
                    j=j1
                else:
                    j=j2 
    
            if mu.DeltaR(tau)<0.4 and mu.DeltaR(j)>0.8 and tau.DeltaR(j)>0.8:
                h['M_MmTau'].Fill((mu+tau).M())
        
out.cd()
for key in h.keys():
    h[key].Write()
