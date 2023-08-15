import os, sys, math
import ROOT
import argparse

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

if "-r" in opts:
    region = sys.argv[2]

#regions = ['MuMu', 'ETau', 'ETau_dRcut', 'ETau_dRcut_Metcut', 'MuTau', 'EMu']

parser = argparse.ArgumentParser(description="Normalize by cross-section")
parser.add_argument("-v", "--version", type=str, help="Version for doHadd")
parser.add_argument("-f", "--filename", type=str, help="filename from condor output")
parser.add_argument("-a", "--all", action="store_true", help="go through all samples")
parser.add_argument("-r", "--region", type=str, help="regions to plot")
parser.add_argument("--histo", type=str, nargs="+")
args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

pi = math.pi

c_t = 1
m_t = 1.77686 #Tau Mass
fa = 1000 #dimensional param

m_a = [10,20,30,40,50,60,70,80,90,100]#pseudoScalar Mass
xsec_gt = [1.0256092e-07, 1.930e-07 , 5.021e-07]

#Br = [0, 6.8 , 9.7 , 5.7 , 6.2 , 3.0 , 3.0 , 9.7 , 0.88 , 1.9 , 1.8 , 2.1 , 2.9] #so Model 1 is M[_][1]
Br=[6.7,3.4,2.5,1.9,1.5,1.2,.95,.79,.66,.57]#Branching ratio of TCP to tau tau from theory paper model 1
xsec_M = []

for i in range (9):
    Gamma_tt = ((c_t**2*m_a[i]/(8*pi))*m_t**2*math.sqrt(1-4*m_t**2/m_a[i]**2))/fa**2
    xsec_g_TCP = 1/ Gamma_tt#cross section of gg fusion to tcp
    xsec_M_list = xsec_g_TCP * Br[i]#Cross Section of gg fusion to TCP to tau Tau
    xsec_M.append(xsec_M_list)
    #print xsec_M[j][i]

###Rescale each
xsec_TCP10 = xsec_M[0]
xsec_TCP30 = xsec_M[2]
xsec_TCP50 = xsec_M[4]

print("xsec_TCP10", xsec_TCP10)
print("xsec_TCP30", xsec_TCP30)
print("xsec_TCP50", xsec_TCP50)

#----------------------------------

TCPxsec = {
    'm10':{'':xsec_M[0]*1.0256092e-07},
    'm30_HT-100to400':{'':xsec_M[2]*8.395e-06},
    'm30_HT-400toInf':{'':xsec_M[2]*1.930e-07},
    'm50_HT-100to400':{'':xsec_M[4]*1.233e-05},
    'm50_HT-400toInf':{'':xsec_M[4]*2.938e-07},
    'm30':{'':xsec_M[2]*8.887e-05+xsec_M[2]*8.395e-06+xsec_M[2]*1.930e-07},
    'm50':{'':xsec_M[4]*9.217e-05+xsec_M[4]*1.233e-05+xsec_M[4]*2.938e-07},
}

xsecs={
#     'TCP_m10':{'':xsec_TCP10},
     'TCP_Ntuple_m12_HT-0to100':{'':xsec_M[2]*8.887e-05},
     'TCP_Ntuple_m12_HT-100to400':{'':3.635E-06*10000000},
     'TCP_Ntuple_m12_HT-400toInf':{'':6.490e-08*10000000},
     'TCP_Ntuple_m30_HT-0to100':{'':xsec_M[2]*8.887e-05},
     'TCP_Ntuple_m30_HT-100to400':{'':xsec_M[2]*8.395e-06},
     'TCP_Ntuple_m30_HT-400toInf':{'':xsec_M[2]*1.930e-07},
     'TCP_Ntuple_m50_HT-0to100':{'':xsec_M[4]*9.217e-05},
     'TCP_Ntuple_m50_HT-100to400':{'':xsec_M[4]*1.233e-05},
     'TCP_Ntuple_m50_HT-400toInf':{'':xsec_M[4]*2.938e-07},
#     'DYJetsToLL_lowMassDY' : {
#         'M-10to50':15890.0},
     'DYJetsToLL_flat':{
#         'M-10to50':15890.0,                                                                                                                               
         'M-50':5398.0},
     'DYStitch':{
         'M-10to50':15890.0,
         'M-50_HT-70to100':146.5,
         'M-50_HT-100to200':160.7,
         'M-50_HT-200to400':48.63,
         'M-50_HT-400to600':6.993,
         'M-50_HT-600to800':1.761,
         'M-50_HT-800to1200':0.8021,
         'M-50_HT-1200to2500':0.1937,
         'M-50_HT-2500toInf':0.003514},
     'DYJetsToLL':{
#         'M-10to50':15890.0,
#         'M-10to50':15890.0*(1/20),
#         'M-50':5398.0},
         'M-50_HT-70to100':146.5,
         'M-50_HT-100to200':160.7, 
         'M-50_HT-200to400':48.63, 
         'M-50_HT-400to600':6.993, 
         'M-50_HT-600to800':1.761, 
         'M-50_HT-800to1200':0.8021, 
         'M-50_HT-1200to2500':0.1937, 
         'M-50_HT-2500toInf':0.003514},
     'DYJetsToLL_M-4to50':{
         'HT-70to100':321.2,
#         'HT-70to100':145.5,
         'HT-100to200':204.0,
         'HT-200to400':54.39,
         'HT-400to600':5.697,
         'HT-600toInf':1.85},
# #    'DYJetsToQQ':{
# #        'HT180':1208},
     # 'TTJets':{
     #      'TuneCP5':750.5},
     # 'TTTo2L2Nu':{
     #     'TuneCP5':88.29},
     # 'TTToSemiLeptonic':{
     #     'TuneCP5':365.34},
     # 'TTToHadronic':{
     #     'TuneCP5':377.96},
     'TT':{
         'TTTo2L2Nu_TuneCP5':88.2497,
         'TTToSemiLeptonic_TuneCP5':365.30899,
         'TTToHadronic_TuneCP5':377.9517},
     'ST':{
         's-channel':3.549,
         't-channel_antitop':26.2278,
         't-channel_top':44.07048,
         'tW_antitop':35.6,
         'tW_top':35.6},
         # 's_channel':3.549,
         # 't_channel_antitop':69.09,
         # 't_channel_top':115.3,
         # 'tW_antitop':34.97,
         # 'tW_top':34.91},
     'Diboson':{
         'WW':75.95,
         'WZ':27.59,
         'ZZ':12.17},
      'WJetsToLNu_flat':{
          'TuneCP5':52940.0},
      'WJetsToLNu':{
          'HT-70To100':1264.0,
          'HT-100To200':1343.0,
          'HT-200To400':359.6,
          'HT-400To600':48.85,
          'HT-600To800':12.05,
          'HT-800To1200':5.501,
          'HT-1200To2500':1.329},
          # 'HT-2500ToInf':0.03216},
          # 'HT-70to100':1264.0,
          # 'HT-100to200':1256.0,
          # 'HT-200to400':335.5,
          # 'HT-400to600':45.25,
          # 'HT-600to800':10.97,
          # 'HT-800to1200':4.933,
          # 'HT-1200to2500':1.16,
          # 'HT-2500toInf':0.008001},
      'QCD':{
          'HT50to100':185300000.0,
          'HT100to200':23590000.0,
          'HT200to300':1551000.0,
          'HT300to500':323400.0,
          'HT500to700':30140.0,
          'HT700to1000':6344.0,
          'HT1000to1500':1092.0,
          'HT1500to2000':99.76,
          'HT2000toInf':20.35}
}


xsecsQCD = {
    'HT-50to100':{'':185300000.0},
    'HT-100to200':{'':23590000.0},
    'HT-200to300':{'':1551000.0},
    'HT-300to500':{'':323400.0},
    'HT-500to700':{'':30140.0},
    'HT-700to1000':{'':6344.0},
    'HT-1000to1500':{'':1092.0},
    'HT-1500to2000':{'':99.76},
    'HT-2000toInf':{'':20.35}
}

print("TCPxsec", TCPxsec)

if args.filename :
    fileName = args.filename

def weightBackgroundHists(hists, files, version, study, var, Sample):
#    for sample in list(xsecs):
    for sample in Sample:
        for mass in list(xsecs[sample]):
#        for mass in list(xsecsQCD[sample]):
            if gen==True:
                filename="h_Gen_"+sample+"_"+mass+"_"+version+".root"
            else:
                filename = "h_"+fileName+"_"+sample+"_"+mass+"_"+version+".root"
                if "-al" in opts:
                    filename="h_debugMuTau_HighHT_Inclusive_Altered_"+sample+"_"+mass+"_"+version+".root"
#                    filename="h_debugMuTau_HighHT_FullyLeptonic_Inclusive_Altered_"+sample+"_"+mass+"_"+version+".root"
                    if "-2d" in opts:
                        filename="h_debugMuTau_HighHT_plot2DforTau_Inclusive_Altered_"+sample+"_"+mass+"_"+version+".root"
                if "-nm" in opts:
#                    filename="h_debugMuTau_HighHT_FullyLeptonic_Inclusive_"+sample+"_"+mass+"_"+version+".root"
                    filename="h_debugMuTau_HighHT_Inclusive_"+sample+"_"+mass+"_"+version+".root"
                    if "-2d" in opts:
                        filename="h_debugMuTau_HighHT_plot2DforTau_Inclusive_"+sample+"_"+mass+"_"+version+".root"
                    if "-mc" in opts:
                        filename="h_debugMuTau_HighHT_Inclusive_MCOnly_"+sample+"_"+mass+"_"+version+".root"
#            histname="h"+reco+"_"+var
            histname=var 
            print(filename, histname)
            fil=ROOT.TFile(filename, 'r')
            if fil.IsZombie(): continue
            files+=[fil]
            hist=fil.Get(histname) 
            nEvt=fil.Get("NEvents").GetBinContent(2)
            xsec=xsecs[sample][mass] 
#            xsec=xsecsQCD[sample][mass] 
            if list(xsecs[sample]).index(mass)==0:  
#            if list(xsecsQCD[sample]).index(mass)==0:  
                h = hist.Clone(sample+"_"+histname)
#                print(h)
                if scaling == 'xsection':
                    h.Scale((xsec/nEvt)*L)
                else: 
                    scale = 1/(h.Integral())
                    h.Scale(scale)
                hists[sample] = h 
            else: 
                if scaling == 'xsection': 
                    hist.Scale((xsec/nEvt)*L)
                else:
                    hist.Scale(scale)
                hists[sample].Add(hist)


hists = {}
files = []

#MuTauDebug - v11
#FullyLeptonic - v2
#JetHT - v2
#HTTrig - v2
#Inclusive - v4
#AlteredID - v6
#Nominal - v6
#2D = v7,v1

if "-nm" in opts:
    version = 'v13'
    iteration = 'v13'
    study ='Norminal'
    if "-2d" in opts:
        version = 'v5'
        iteration = 'v5'
        study ='2DNominal'
    if "-mc" in opts:
        version = 'v1'
        iteration = 'v1'
        study ='MCOnly'


if "-al" in opts:
    version = 'v12'
    iteration = 'v12'
    study ='AlteredID'
    if "-2d" in opts:
        version = 'v5'
        iteration = 'v5'
        study ='2DAltered'


L = 41480.0
#L = 1
#-------------

gen = False
scaling = 'xsection' #'xsection' else: 'unity'

# version = 'v5'
# study = 'applyBjetSF'
# iteration = 'v5'

if args.version :
    version = args.version
    iteration = args.version
    study = args.filename

#-------------

Sample = ['TTJets']

if "--tcp" in opts:
    Sample = ['TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf', 'TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf']
if "-l" in opts:
    Sample = ['TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf']
if "-h" in opts:
    Sample = ['TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf']
if args.all :
#    Sample = ['TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf','TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf','WJetsToLNu','DYJetsToLL','DYJetsToLL_M-4to50','Diboson','ST','QCD']
#    Sample = ['TCP_Ntuple_m12_HT-100to400', 'TCP_Ntuple_m12_HT-400toInf', 'TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf']
#    Sample = ['TCP_Ntuple_m12_HT-100to400', 'TCP_Ntuple_m12_HT-400toInf']
    Sample = ['TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf','TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf', 'TCP_Ntuple_m12_HT-100to400', 'TCP_Ntuple_m12_HT-400toInf']
#    Sample = ['TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf','TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf','DYJetsToLL','DYJetsToLL_M-4to50','QCD','WJetsToLNu','TT','Diboson','ST']
if "--tt" in opts:
    Sample = ['TTTo2L2Nu', 'TTToSemiLeptonic', 'TTToHadronic']
if "-b" in opts:
    Sample = ['TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf','TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf','QCD','WJetsToLNu','DYJetsToLL','TTJets','DYJetsToLL_M-4to50','Diboson']
if "-ht" in opts:
    Sample = ['QCD','WJetsToLNu','DYJetsToLL','DYJetsToLL_M-4to50']
if "-d" in opts:
    Sample = ['DYJetsToLL','DYJetsToLL_M-4to50']
if "-q" in opts:
    Sample = ['HT-50to100','HT-100to200','HT-200to300','HT-300to500','HT-500to700','HT-700to1000','HT-1000to1500','HT-1500to2000','HT-2000toInf']


#VARIABLE = ["TauPtJetPt","TauPtJet2Pt","TauPtMuonPt","TauPtdRl","TauPtdRj2tau","TauPtdRjtau","TauPtdRjmu","MuonPtdRl","TauPtdRgenMu","MuonPtdRgenMu"]
#VARIABLE = ['TauPtdRjmu','TauPtdRjtau','TauPtdRl','MuonPtdRl','TauPtMuonPt','TauPtJetPt','TauPtJet2Pt','TauPtdRj2tau','DimuonMass','MuonPtMuon2Pt','TauPtdRl2','TauPtdRgenMu','MuonPtdRgenMu']
#VARIABLE = ['TauPt','Nj','JetPt', 'MuonPt', 'MetPt','Mass','Mt']
#VARIABLE = ['MetPt']
#VARIABLE = ['Mass', 'Lepton1Pt', 'Lepton2Pt', 'JetPt', 'MetPt', 'Nj','dRl','dRj', 'dPhil', 'dPhi','Mtl1', 'Mtl2','Mtl','cosMl1','cosMl2','cosl','Count']
#VARIABLE = ['Mass', 'Lepton1Pt', 'Lepton2Pt', 'JetPt', 'MetPt', 'Mt','Nj','dRl','dRj', 'dPhil', 'dPhi','Count']
#VARIABLE = ['Lepton1Pt', 'Lepton2Pt', 'LeadingJetPt', 'Count', 'HT','dRl','dRj','Mass','MetPt']
VARIABLE = ['Lepton1Pt', 'Lepton2Pt', 'LeadingJetPt', 'Mass']
#VARIABLE = ["BFlavour_JetPt", "BFlavour_JetEta","CFlavour_JetPt", "CFlavour_JetEta", "LFlavour_JetPt", "LFlavour_JetEta", "BFlavour_BTagged_JetPt", "BFlavour_BTagged_JetEta","CFlavour_BTagged_JetPt", "CFlavour_BTagged_JetEta", "LFlavour_BTagged_JetPt", "LFlavour_BTagged_JetEta"]
#VARIABLE = ['MuonPt_SingleMuon','ElectronPt_SingleMuon','MetPt_SingleMuon','dRl_SingleMuon']
#VARIABLE = ['MuonPt_Both','ElectronPt_Both','dRl_Both','JetPt_Both','MetPt_Both','dRj_Both']
#VARIABLE = ['MuonPt_Both','ElectronPt_Both','dRl_Both','JetPt_Both','dRj_Both']
#VARIABLE = ['Mass', 'cosMl1','cosMl2','Count']
#VARIABLE = ['Mass']
#VARIABLE = ['TauPt', 'TauPt0','TauPt1','TauPt10','Nj','JetPt', 'MuonPt']
#VARIABLE = ['TauPtMass','TauPt0Mass','TauPt1Mass','TauPt10Mass','NJetMass']
#VARIABLE = ['Mt','MetPt']

REGION = []

histlist = []

#for region in REGION:
if "-r" in opts:
    for variable in VARIABLE:
        histlist.append(args.region+"_"+variable)
#        histlist.append(region+"_"+variable+"_loosedR")

if args.histo is not None:
    histlist = args.histo

for var in histlist:
    weightBackgroundHists(hists, files, version, study, var, Sample)

    if scaling == 'xsection':
#         out = ROOT.TFile("h_"+reco+"_"+var+".root",'recreate')
         out = ROOT.TFile("h_"+var+"_"+study+"_"+iteration+".root",'recreate')
         print(out)
    else: 
         out = ROOT.TFile("h_"+study+"_"+var+"_"+version+"_unity.root",'recreate')

    out.cd()

    for name in Sample:
#        hists[name].Rebin(10)                                                                                                           


        if name=='DYJetsToLL':
            hists[name].SetFillColor(ROOT.kRed-6)
        elif name=='DYJetsToLL_M-4to50':
            hists[name].SetFillColor(ROOT.kRed-9)
#        elif name=='DYJetsToLL_flat':
#            hists[name].SetFillColor(ROOT.kRed-6)
        # elif name == 'TCP_m10':
        #     hists[name].SetFillStyle(3335)
        #     hists[name].SetFillColor(ROOT.kBlue)
        # elif name == 'TCP_m30':
        #     hists[name].SetFillStyle(3335)
        #     hists[name].SetFillColor(ROOT.kGreen)
        # elif name == 'TCP_m50':
        #     hists[name].SetFillStyle(3335)
        #     hists[name].SetFillColor(ROOT.kBlack)
        elif name=='TTJets':
            hists[name].SetFillColor(ROOT.kOrange-4)
        elif name=='TTTo2L2Nu':
            hists[name].SetFillColor(ROOT.kOrange-4)
        elif name=='TTToSemiLeptonic':
            hists[name].SetFillColor(ROOT.kOrange-5)
        elif name=='TTToHadronic':
            hists[name].SetFillColor(ROOT.kOrange-6)
        elif name=='ST':
            hists[name].SetFillColor(ROOT.kMagenta-9)
        elif name == 'Diboson':
            hists[name].SetFillColor(ROOT.kAzure+7)
        elif name == 'QCD':
            hists[name].SetFillColor(ROOT.kGray)
        elif name=='WJetsToLNu':
            hists[name].SetFillColor(ROOT.kGreen-6)
#        elif name=='WJetsToLNu_flat':
#            hists[name].SetFillColor(ROOT.kGreen-6)

        hists[name].Write()
        print(name)

    out.Close()

print(VARIABLE)
print("Histlist")
print(*histlist, sep=", ")
print("Regions")
print(*REGION, sep=", ")
