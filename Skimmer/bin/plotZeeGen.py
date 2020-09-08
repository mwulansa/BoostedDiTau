import ROOT, sys, os
from DataFormats.FWLite import Events, Handle
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

#jobs=np.linspace(100, 1, 100)
#jobs=[100]

#masses=[10, 30, 50]
#masses=[10]
#masses=[50]

handleGenJet = Handle ('vector<reco::GenJet>')
labelGenJet = ('slimmedGenJets')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

h={}
h['N_Events'] = ROOT.TH1F ('N_Events', '', 2, 0, 2)
h['M_Z'] = ROOT.TH1F ('M_Z', '', 100, 0, 1000)
h['ePt'] = ROOT.TH1F ('ePt', '', 500, 0, 500)


prefix = "root://xrootd.unl.edu/"

mass = "50_120"

#for mass in masses:
out=ROOT.TFile("h_plotZeeGen_m_"+str(mass)+"_94X.root",'recreate')

searchString = "/ZToEE*"+mass+"/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"

os.system('dasgoclient --query "'+searchString+'"')
query = 'dasgoclient --query "file dataset='+searchString+'"'
print query
files=os.popen(query).read().split()
print files
#sys.exit()

#files = ["/store/mc/RunIIFall17MiniAODv2/ZToEE_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/12E59949-4C44-E811-A054-0CC47A4C8E56.root"]
nevt=0
for fil in files:
    print prefix+fil
    events=Events(prefix+fil)

    for event in events:
        nevt+=1
    
        #if nevt > 100: sys.exit()
        #print "Start Processing Event:", nevt
    
        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()
    
        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()
    
        h['N_Events'].Fill(0.5, 1)
        h['N_Events'].Fill(1.5, genweight)
    
        Zs=[]
        for particle in particles:
            if particle.pdgId() == 23 and particle.isHardProcess():
                h['M_Z'].Fill(particle.mass(), genweight)
                Zs+=[particle]
    
        electrons = []
        for particle in particles:
            if abs(particle.pdgId()) == 11 and particle.isHardProcess():
                electrons+=[particle]
    
        if not len(electrons) == 2:
            print "ERROR!", len(electrons)
    
        selected_electrons = []
        for electron in electrons:
            if electron.pt()>7 and abs(electron.eta())<2.5:
                selected_electrons+=[electron]
    
        selected_electrons.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_electrons)>0:
            h['ePt'].Fill(selected_electrons[0].pt(), genweight)

out.cd()
for key in h.keys():
    h[key].Write()
out.Close()
