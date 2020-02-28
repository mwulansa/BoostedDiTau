import ROOT, sys, math, gc
from DataFormats.FWLite import Events, Handle
import numpy as np

handleTaus = Handle ('vector<pat::Tau>')
labelTaus = ('NewTauIDsEmbedded')
    
handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlineSlimmedPrimaryVertices')

jobs=np.linspace(100,51,50)
#jobs = np.linspace(50, 1, 50)

print jobs

#out=ROOT.TFile("./hists/h_plot_tauIds.root",'recreate')

prefix="root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/events/ALP/RunIISummer17DR94Premix/"

for job in jobs:

    out=ROOT.TFile("./hists/h_plot_tauIds_"+str(int(job))+".root",'recreate')
    filename = 'ALP_m50_w1_htjmin400_RunIISummer17DR94Premix_reMINIAODSIM_Cleaned_'+str(int(job))+'.root'
    print prefix+filename

    events=Events(prefix+filename)

    h={}
    h['tauPt'] = ROOT.TH1F ("h_tauPt", "", 500, 0, 500)
    h['tauPt_new'] = ROOT.TH1F ("h_tauPt_new", "", 500, 0, 500)
    h['NEvent'] = ROOT.TH1F ("hNEvent","",1,0,2)
    h["NTaus"] = ROOT.TH1F ("hNTaus","",4,0,4)
    h["NTaus_new"] = ROOT.TH1F ("hNTaus_new","",4,0,4)
    
    for event in events:

        #print h['NEvent']
        
        h['NEvent'].Fill(1)
    
        event.getByLabel(labelTaus, handleTaus)
        taus=handleTaus.product()
    
        h["NTaus"].Fill(.5, taus.size())
        h["NTaus_new"].Fill(.5, taus.size())
    
        event.getByLabel(labelVertex, handleVertex)
        vertex=handleVertex.product()
        pv = vertex[0].position()
    
        selected_taus=[]
        selected_taus_new=[]
        if taus.size()>0:
            for tau in taus:
                if tau.pt()>10 and abs(tau.eta())<2.3:
                    h["NTaus"].Fill(1.5)
                    if tau.tauID("decayModeFinding")>0.5 \
                    and tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") > -0.5 \
                    and tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"):
                        h["NTaus"].Fill(2.5)
                        selected_taus+=[tau]
                        dxy = abs(tau.leadChargedHadrCand().get().dxy(pv))
                        dz = abs(tau.leadChargedHadrCand().get().dz(pv))
                        if dxy <0.2 and dz < 0.5:
                            h["NTaus"].Fill(3.5)
                    if tau.tauID("decayModeFinding")>0.5 \
                    and tau.tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017") > -0.5 \
                    and tau.tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"):
                        h["NTaus_new"].Fill(2.5)
                        selected_taus_new+=[tau]
                        dxy = abs(tau.leadChargedHadrCand().get().dxy(pv))
                        dz = abs(tau.leadChargedHadrCand().get().dz(pv))
                        if dxy <0.2 and dz < 0.5:
                            h["NTaus_new"].Fill(3.5)
        
        selected_taus.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_taus) > 0:
                h['tauPt'].Fill(selected_taus[0].pt())
    
        selected_taus_new.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_taus_new) > 0:
                h['tauPt_new'].Fill(selected_taus_new[0].pt())

    out.cd()

    for key in h.keys():
        h[key].Write()

    out.Close()
