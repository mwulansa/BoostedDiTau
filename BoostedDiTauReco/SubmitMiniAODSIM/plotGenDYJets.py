import ROOT,sys,os
from DataFormats.FWLite import Events, Handle
from looseElectron import *
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

#inputFileListDir="./filelists/"+Sample+"/"
inputFileListName=sys.argv[1]

#inputFileList=inputFileListDir+inputFileListName
inputFileList=inputFileListName

if len(sys.argv)>2:
    outputFileDir=sys.argv[2]
else:
    outputFileDir = "./plots/"
outputFileName = outputFileDir+"h_Gen_"+inputFileListName.split("/")[-1].replace(".txt",".root")
print(outputFileName)

handleGenJet = Handle ('vector<reco::GenJet>')
labelGenJet = ('slimmedGenJets')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

out=ROOT.TFile.Open(outputFileName,'recreate')

hNEvent = ROOT.TH1F ("hNEvent","Number of Events;;N_{events}", 2, 0, 2)
hJet1Pt = ROOT.TH1F ("hJet1Pt", "leading jet Pt;P_{t};N_{events}", 2000, 0, 2000)
    
hMuMu_M = ROOT.TH1F ("hMuMu_M", "#mu - #mu mass;M_{#mu#mu};N_{events}", 2000, 0, 400)
hMuMu_Mu1Pt = ROOT.TH1F ("hMuMu_Mu1Pt", "leading muon Pt;P_{t};N_{events}", 500, 0, 500)
hMuMu_Mu2Pt = ROOT.TH1F ("hMuMu_Mu2Pt", "sub-leading muon Pt;P_{t};N_{events}", 500, 0, 500)
hMuMu_dR = ROOT.TH1F ("hMuMu_dR", "#mu - #mu delta R;#delta R;N_{events}", 200, 0, 10)
hMuMu_dRjlMu = ROOT.TH1F ("hMuMu_dRjlMu", "#mu - #mu delta R;#delta R;N_{events}", 200, 0, 10)
hMuMu_dRjtMu = ROOT.TH1F ("hMuMu_dRjtMu", "#mu - #mu delta R;#delta R;N_{events}", 200, 0, 10)

hMuMuSelected_M = ROOT.TH1F ("hMuMuSelected_M", "#mu - #mu mass;M_{#mu#mu};N_{events}", 2000, 0, 400)
hMuMuSelected_Mu1Pt = ROOT.TH1F ("hMuMuSelected_Mu1Pt", "leading muon Pt;P_{t};N_{events}", 50, 0, 500)
hMuMuSelected_Mu2Pt = ROOT.TH1F ("hMuMuSelected_Mu2Pt", "sub-leading muon Pt;P_{t};N_{events}", 50, 0, 500)
hMuMuSelected_dR = ROOT.TH1F ("hMuMuSelected_dR", "#mu - #mu delta R;#delta R;N_{events}", 200, 0, 10)
hMuMuSelected_dRjlMu = ROOT.TH1F ("hMuMuSelected_dRjlMu", "#mu - #mu delta R;#delta R;N_{events}", 200, 0, 10)
hMuMuSelected_dRjtMu = ROOT.TH1F ("hMuMuSelected_dRjtMu", "#mu - #mu delta R;#delta R;N_{events}", 200, 0, 10)

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    print(inputFileName.replace("\n",""))
    inputFileName=inputFileName.replace("\n","")
    f=ROOT.TFile.Open(inputFileName)

    if not f.IsZombie():
        events=Events(inputFileName)
    else:
        print("Can't Open File: "+inputFileName)
        continue

    for event in events:
        
        event.getByLabel(labelGenJet, handleGenJet)
        jets=handleGenJet.product()

        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()

        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()

        hNEvent.Fill(0.5, 1)
        hNEvent.Fill(1.5, genweight)

        #print genweight, hNEvent.GetBinContent(2)

        genjets=[]
        for jet in jets:
            genjets+=[jet]
            
        genjets.sort(key=lambda x: x.pt(), reverse=True)

        if len(genjets)>0:
            hJet1Pt.Fill(genjets[0].pt(), genweight)

            j=ROOT.TLorentzVector()
            j.SetPtEtaPhiM(genjets[0].pt(), genjets[0].eta(), genjets[0].phi(), genjets[0].mass())

            muons=[]
            mup=[]
            mum=[]
            for particle in particles:
                if particle.pdgId()==13 and (particle.isHardProcess() or particle.isDirectHardProcessTauDecayProductFinalState()):
                    mup+=[particle]
                if particle.pdgId()==-13 and (particle.isHardProcess() or particle.isDirectHardProcessTauDecayProductFinalState()):
                    mum+=[particle]

            if len(mup)>1 or len(mum)>1:
                print("Gen Particle Selection is Wrong!")
        
            if len(mup)==1 and len(mum)==1:
                muons=[mup[0],mum[0]]

            muons.sort(key=lambda x: x.pt(), reverse=True)

            if len(muons)>0:
                mu1=ROOT.TLorentzVector()
                mu1.SetPtEtaPhiM(muons[0].pt(), muons[0].eta(), muons[0].phi(), muons[0].mass())
                mu2=ROOT.TLorentzVector()
                mu2.SetPtEtaPhiM(muons[1].pt(), muons[1].eta(), muons[1].phi(), muons[1].mass())

                hMuMu_M.Fill((mu1+mu2).M(), genweight)
                hMuMu_Mu1Pt.Fill(mu1.Pt(), genweight)
                hMuMu_Mu2Pt.Fill(mu2.Pt(), genweight)
                hMuMu_dR.Fill(mu1.DeltaR(mu2), genweight)
                hMuMu_dRjlMu.Fill(j.DeltaR(mu1), genweight)
                hMuMu_dRjtMu.Fill(j.DeltaR(mu2), genweight)

                selected_muons=[]
                for muon in muons:
                    #print muon.pt(), abs(muon.eta())
                    if muon.pt()>3 and abs(muon.eta())<2.4:
                        selected_muons+=[muon]
                        
                #print len(selected_muons)

                if len(selected_muons)>0:
                    mu1=ROOT.TLorentzVector()
                    mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass())

                    hMuMuSelected_Mu1Pt.Fill(mu1.Pt(), genweight)

                    if len(selected_muons)>1:
                        mu2=ROOT.TLorentzVector()
                        mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass())

                        hMuMuSelected_M.Fill((mu1+mu2).M(), genweight)
                        hMuMuSelected_Mu2Pt.Fill(mu2.Pt(), genweight)
                        hMuMuSelected_dR.Fill(mu1.DeltaR(mu2), genweight)
                        hMuMuSelected_dRjlMu.Fill(j.DeltaR(mu1), genweight)
                        hMuMuSelected_dRjtMu.Fill(j.DeltaR(mu2), genweight)

out.cd()
hNEvent.Write()
hJet1Pt.Write()

hMuMu_M.Write()
hMuMu_Mu1Pt.Write()
hMuMu_Mu2Pt.Write()
hMuMu_dR.Write()
hMuMu_dRjlMu.Write()
hMuMu_dRjtMu.Write()

hMuMuSelected_M.Write()
hMuMuSelected_Mu1Pt.Write()
hMuMuSelected_Mu2Pt.Write()
hMuMuSelected_dR.Write()
hMuMuSelected_dRjlMu.Write()
hMuMuSelected_dRjtMu.Write()
out.Close()                   
                          

        

        
            

