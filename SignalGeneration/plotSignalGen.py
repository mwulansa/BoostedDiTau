import ROOT,sys
from DataFormats.FWLite import Events, Handle
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

jobs=np.linspace(100, 1, 100)
#jobs=[90]

masses=[10, 30, 50]

handleGenJet = Handle ('vector<reco::GenJet>')
labelGenJet = ('ak4GenJets')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('genParticles')

prefix="root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer16DR80Premix/"

for mass in masses:
    out=ROOT.TFile("h_plotSignalGen_m"+str(int(mass))+".root",'recreate')
    
    hJet1Pt = ROOT.TH1F ("hJet1Pt", "leading jet Pt;P_{t};N_{events}", 150, 0, 1500)
    hDiTauM = ROOT.TH1F ("hDiTauM ", "di-tau mass;M_{#tau#tau};N_{events}", 100, 0, 100)

    hEMu_M = ROOT.TH1F ("hEMu_M", "e - #mu mass;M_{#mu#mu};N_{events}", 100, 0, 100)
    hEMu_MuPt = ROOT.TH1F ("hEMu_MuPt", "muon Pt;P_{t};N_{events}", 100, 0, 1000)
    hEMu_EPt = ROOT.TH1F ("hEMu_EPt", "electron Pt;P_{t};N_{events}", 100, 0, 1000)
    hEMu_dR = ROOT.TH1F ("hEMu_dR", "e - #mu delta R;#delta R;N_{events}", 20, 0, 1)
    
    hMuMu_M = ROOT.TH1F ("hMuMu_M", "#mu - #mu mass;M_{#mu#mu};N_{events}", 100, 0, 100)
    hMuMu_Mu1Pt = ROOT.TH1F ("hMuMu_Mu1Pt", "leading muon Pt;P_{t};N_{events}", 100, 0, 1000)
    hMuMu_Mu2Pt = ROOT.TH1F ("hMuMu_Mu2Pt", "sub-leading muon Pt;P_{t};N_{events}", 100, 0, 1000)
    hMuMu_dR = ROOT.TH1F ("hMuMu_dR", "#mu - #mu delta R;#delta R;N_{events}", 20, 0, 1)

    hTauTau_M = ROOT.TH1F ("hTauTau_M", "#tau - #tau mass;M_{#tau#tau};N_{events}", 100, 0, 100)
    hTauTau_Tau1Pt = ROOT.TH1F ("hTauTau_Tau1Pt", "leading tau Pt;P_{t};N_{events}", 100, 0, 1000)
    hTauTau_Tau2Pt = ROOT.TH1F ("hTauTau_Tau2Pt", "sub-leading tau Pt;P_{t};N_{events}", 100, 0, 1000)
    hTauTau_dR = ROOT.TH1F ("hTauTau_dR", "#tau - #tau delta R;#delta R;N_{events}", 20, 0, 1)

    hETau_M = ROOT.TH1F ("hETau_M", "e - #tau mass;M_{e#tau};N_{events}", 100, 0, 100)
    hETau_TauPt = ROOT.TH1F ("hETau_TauPt", "tau Pt;P_{t};N_{events}", 100, 0, 1000)
    hETau_EPt = ROOT.TH1F ("hETau_EPt", "electron Pt;P_{t};N_{events}", 100, 0, 1000)
    hETau_dR = ROOT.TH1F ("hETau_dR", "e - #tau delta R;#delta R;N_{events}", 20, 0, 1)

    hMuTau_M = ROOT.TH1F ("hMuTau_M", "#mu - #tau mass;M_{#mu#tau};N_{events}", 100, 0, 100)
    hMuTau_TauPt = ROOT.TH1F ("hMuTau_TauPt", "tau Pt;P_{t};N_{events}", 100, 0, 1000)
    hMuTau_MuPt = ROOT.TH1F ("hMuTau_MuPt", "muon Pt;P_{t};N_{events}", 100, 0, 1000)
    hMuTau_dR = ROOT.TH1F ("hMuTau_dR", "#mu - #tau delta R;#delta R;N_{events}", 20, 0, 1)
    
    for job in jobs:
    
        events=Events(prefix+"ALP_m"+str(int(mass))+"_w1_htjmin400_RunIISummer16DR80Premix_AODSIM_"+str(int(job))+".root")
        nevt=0
        for event in events:
            nevt+=1
        
            event.getByLabel(labelGenJet, handleGenJet)
            #jets=handleGenJet.product()
        
            jets=[]
            for jet in handleGenJet.product():
                jets+=[jet]
                
            jets.sort(key=lambda x: x.pt(), reverse=True)
        
            jet1=ROOT.TLorentzVector(jets[0].px(), jets[0].py(), jets[0].pz(), jets[0].energy())
        
            hJet1Pt.Fill(jet1.Pt())
        
            event.getByLabel(labelGenParticle, handleGenParticle)
            particles=handleGenParticle.product()
        
            #for particle in particles:
            #    print "pdgId: "+str(particle.pdgId())+" pt: "+str(particle.pt())+" status: "+str(particle.status())
            #    if particle.pdgId() == 9999:
            #        for i in range(particle.numberOfDaughters()):
            #            print "Daughter: "+str(i)+" pdgId: "+str(particle.daughter(i).pdgId())+" pt: "+str(particle.daughter(i).pt())+" status: "+str(particle.daughter(i).status())
                #if particle.pdgId()==15:
                #    print particle.mother().pdgId()
            #sys.exit()
        
            taus_m=[]
            taus_p=[]
            for particle in particles:
                if particle.pdgId()==15 and particle.mother().pdgId()==9999:
                    taus_m+=[particle]
                if particle.pdgId()==-15 and particle.mother().pdgId()==9999:
                    taus_p+=[particle]
                    
            #taus_m.sort(key=lambda x: x.pt(), reverse=True)
            #taus_p.sort(key=lambda x: x.pt(), reverse=True)
        
            if len(taus_p)!=1 or len(taus_m)!=1:
                print "Something is wrong (debug 1)..."
            
            tau1=ROOT.TLorentzVector()
            tau1.SetPtEtaPhiM(taus_m[0].pt(), taus_m[0].eta(), taus_m[0].phi(), taus_m[0].mass())
            tau2=ROOT.TLorentzVector()
            tau2.SetPtEtaPhiM(taus_p[0].pt(), taus_p[0].eta(), taus_p[0].phi(), taus_p[0].mass())
        
            hDiTauM.Fill((tau1+tau2).M())
        
            leptons1=[]
            leptons2=[]
            for i in range(taus_m[0].numberOfDaughters()):
                if taus_m[0].daughter(i).pdgId()==11 or taus_m[0].daughter(i).pdgId()==13:
                    leptons1+=[taus_m[0].daughter(i)]
        
            for i in range(taus_p[0].numberOfDaughters()):
                if taus_p[0].daughter(i).pdgId()==-11 or taus_p[0].daughter(i).pdgId()==-13:
                    leptons2+=[taus_p[0].daughter(i)]
    
            if len(leptons1)>1 or len(leptons2)>1:
                print "Something is wrong (debug 2)..."
             
            if len(leptons1)==len(leptons2)==0:
                tau1=ROOT.TLorentzVector()
                tau1.SetPtEtaPhiM(taus_m[0].pt(), taus_m[0].eta(), taus_m[0].phi(), taus_m[0].mass())
                tau2=ROOT.TLorentzVector()
                tau2.SetPtEtaPhiM(taus_p[0].pt(), taus_p[0].eta(), taus_p[0].phi(), taus_p[0].mass())
                hTauTau_M.Fill((tau1+tau2).M())
                hTauTau_dR.Fill(tau1.DeltaR(tau2))
                hTauTau_Tau1Pt.Fill(tau1.Pt()) if tau1.Pt()>tau2.Pt() else hTauTau_Tau1Pt.Fill(tau2.Pt()) 
                hTauTau_Tau2Pt.Fill(tau2.Pt()) if tau1.Pt()>tau2.Pt() else hTauTau_Tau2Pt.Fill(tau1.Pt())
    
            elif len(leptons1)==len(leptons2)==1:
            
                if leptons1[0].charge() == leptons2[0].charge():
                    print "Something is wrong (debug 3)..."
        
                #if abs(leptons1[0].pdgId())==abs(leptons2[0].pdgId())==11: continue
        
                lepton1=ROOT.TLorentzVector()
                lepton1.SetPtEtaPhiM(leptons1[0].pt(), leptons1[0].eta(), leptons1[0].phi(), leptons1[0].mass())
                lepton2=ROOT.TLorentzVector()
                lepton2.SetPtEtaPhiM(leptons2[0].pt(), leptons2[0].eta(), leptons2[0].phi(), leptons2[0].mass())
        
                if abs(leptons1[0].pdgId())==abs(leptons2[0].pdgId())==13:
                    hMuMu_Mu1Pt.Fill(lepton1.Pt()) if lepton1.Pt()>lepton2.Pt() else hMuMu_Mu1Pt.Fill(lepton2.Pt())
                    hMuMu_Mu2Pt.Fill(lepton2.Pt()) if lepton1.Pt()>lepton2.Pt() else hMuMu_Mu2Pt.Fill(lepton1.Pt())
                    hMuMu_M.Fill((lepton1+lepton2).M())
                    hMuMu_dR.Fill(lepton1.DeltaR(lepton2))
                elif abs(leptons1[0].pdgId())==abs(leptons2[0].pdgId())==11: continue
                else:
                    hEMu_MuPt.Fill(lepton1.Pt()) if leptons1[0].pdgId()==13 else hEMu_MuPt.Fill(lepton2.Pt())
                    hEMu_EPt.Fill(lepton1.Pt()) if leptons1[0].pdgId()==11 else hEMu_EPt.Fill(lepton2.Pt())
                    hEMu_M.Fill((lepton1+lepton2).M())
                    hEMu_dR.Fill(lepton1.DeltaR(lepton2))
    
            else:
    
                if len(leptons1)==0:
                    tau=ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(taus_m[0].pt(), taus_m[0].eta(), taus_m[0].phi(), taus_m[0].mass())
    
                    if leptons2[0].pdgId()==-11:
                        electron=ROOT.TLorentzVector()
                        electron.SetPtEtaPhiM(leptons2[0].pt(), leptons2[0].eta(), leptons2[0].phi(), leptons2[0].mass())
                        hETau_M.Fill((tau+electron).M())
                        hETau_dR.Fill(tau.DeltaR(electron))
                        hETau_EPt.Fill(electron.Pt())
                        hETau_TauPt.Fill(tau.Pt())
                    else:
                        muon=ROOT.TLorentzVector()
                        muon.SetPtEtaPhiM(leptons2[0].pt(), leptons2[0].eta(), leptons2[0].phi(), leptons2[0].mass())
                        hMuTau_M.Fill((tau+muon).M())
                        hMuTau_dR.Fill(tau.DeltaR(muon))
                        hMuTau_MuPt.Fill(muon.Pt())
                        hMuTau_TauPt.Fill(tau.Pt())
                else:
                    tau=ROOT.TLorentzVector()
                    tau.SetPtEtaPhiM(taus_p[0].pt(), taus_p[0].eta(), taus_p[0].phi(), taus_p[0].mass())
                    if leptons1[0].pdgId()==11:
                        electron=ROOT.TLorentzVector()
                        electron.SetPtEtaPhiM(leptons1[0].pt(), leptons1[0].eta(), leptons1[0].phi(), leptons1[0].mass())
                        hETau_M.Fill((tau+electron).M())
                        hETau_dR.Fill(tau.DeltaR(electron))
                        hETau_EPt.Fill(electron.Pt())
                        hETau_TauPt.Fill(tau.Pt())
                    else:
                        muon=ROOT.TLorentzVector()
                        muon.SetPtEtaPhiM(leptons1[0].pt(), leptons1[0].eta(), leptons1[0].phi(), leptons1[0].mass())
                        hMuTau_M.Fill((tau+muon).M())
                        hMuTau_dR.Fill(tau.DeltaR(muon))
                        hMuTau_MuPt.Fill(muon.Pt())
                        hMuTau_TauPt.Fill(tau.Pt())

    out.cd()
    hJet1Pt.Write()
    hDiTauM.Write()
    hEMu_M.Write()
    hEMu_MuPt.Write()
    hEMu_EPt.Write()
    hEMu_dR.Write()
    hMuMu_M.Write()
    hMuMu_Mu1Pt.Write()
    hMuMu_Mu2Pt.Write()
    hMuMu_dR.Write()
    hTauTau_M.Write()
    hTauTau_Tau1Pt.Write()
    hTauTau_Tau2Pt.Write()
    hTauTau_dR.Write()
    hETau_M.Write()
    hETau_TauPt.Write()
    hETau_EPt.Write()
    hETau_dR.Write()
    hMuTau_M.Write()
    hMuTau_TauPt.Write()
    hMuTau_MuPt.Write()
    hMuTau_dR.Write()
    out.Close()

    
    

    
    
    
