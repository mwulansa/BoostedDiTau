import os, sys, math
import ROOT

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]

#histlist = sys.argv[1]

# Cross section calc----------------

# pi = math.pi

# c_t = 1
# m_t = 1.77686
# fa = 1000

# m_a = [10,30,50]
# #xsec_gt = [1.0256092e-07, 3.133914e-07, 5.072265e-07]

# xsec_gt = [1.0256092e-07, 1.930e-07 , 5.021e-07]

# Br = [[0,0,0], [6.7, 2.5, 1.5] , [9.7, 3.6, 2.1] , [5.7, 2.2, 1.4] , [6.2, 1.6, 0.79] , [3.0, 1.1, 0.66] , [3.0, 1.1, 0.66] , [9.7, 3.6, 2.1] , [0.88, 0.48, 0.43] , [1.9, 0.42, 0.19] , [1.8, 0.41, 0.19] , [2.1, 0.94, 0.66] , [2.9, 1.3, 0.85]] #so Model 1 is M[_][1]

# xsec_M = [[],[],[]]

# # for j in range (3):
# #     for i in range (13):
# #         Gamma_tt = ((c_t**2*m_a[j]/(8*pi))*m_t**2*math.sqrt(1-4*m_t**2/m_a[j]**2))/fa**2
# #         xsec_g_TCP = xsec_gt[j] / Gamma_tt
# #         xsec_M_list = xsec_g_TCP * Br[i][j]
# #         xsec_M[j].append(xsec_M_list)
# #         print xsec_M[j][i]
# #     print "Gamma_tt = ", Gamma_tt

# for i in range (9):
#     Gamma_tt = ((c_t**2*m_a[i]/(8*pi))*m_t**2*math.sqrt(1-4*m_t**2/m_a[i]**2))/fa**2
#     xsec_g_TCP = 1/ Gamma_tt#cross section of gg fusion to tcp
#     xsec_M_list = xsec_g_TCP * Br[i]#Cross Section of gg fusion to TCP to tau Tau
#     xsec_M.append(xsec_M_list)
#     #print xsec_M[j][i]
# ###Rescale each
# xsec_TCP10 = xsec_M[0]*41.53e3
# xsec_TCP30 = xsec_M[2]*41.53e3
# xsec_TCP50 = xsec_M[4]*41.53e3

# xsec_TCP10 = xsec_M[0][1]
# xsec_TCP30 = xsec_M[1][1]
# xsec_TCP50 = xsec_M[2][1]


pi = math.pi

c_t = 1
m_t = 1.77686 #Tau Mass
fa = 1000 #dimensional param

m_a = [10,20,30,40,50,60,70,80,90,100]#pseudoScalar Mass
xsec_gt = [1.0256092e-07, 1.930e-07 , 5.021e-07]

#Br = [0, 6.8 , 9.7 , 5.7 , 6.2 , 3.0 , 3.0 , 9.7 , 0.88 , 1.9 , 1.8 , 2.1 , 2.9] #so Model 1 is M[_][1]
Br=[6.7,3.4,2.5,1.9,1.5,1.2,.95,.79,.66,.57]#Branching ration of TCP to tau tau from theory paper model 1
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


#----------------------------------

xsecs={
#     'TCP_m10':{'':xsec_TCP10},
#     'TCP_m30':{
#     'TCP_m30_HT-0to100':{'':xsec_M[2]*8.887e-05},
     'TCP_Ntuple_m30_HT-100to400':{'':xsec_M[2]*8.395e-06},
     'TCP_Ntuple_m30_HT-400toInf':{'':xsec_M[2]*1.930e-07},
#     'TCP_Ntuple_m50_HT-0to100':{'':xsec_M[4]*9.217e-05},
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
#         'HT-70to100':321.2,
         'HT-70to100':145.5,
         'HT-100to200':204.0,
         'HT-200to400':54.39,
         'HT-400to600':5.697,
         'HT-600toInf':1.85},
# #    'DYJetsToQQ':{
# #        'HT180':1208},
     'TTJets':{
         'TuneCP5':750.5},
#     'ST':{
#         's-channel_4f_leptonDecays':3.549,
#         't-channel_antitop_4f_inclusiveDecays':26.38,
#         't-channel_top_4f_inclusiveDecays':44.33,
#         'tW_antitop_5f_inclusiveDecays':34.97,
#         'tW_top_5f_inclusiveDecays':34.91},
     'Diboson':{
         'WW':75.95,
         'WZ':27.59,
         'ZZ':12.17},
      'WJetsToLNu_flat':{
          'TuneCP5':52940.0},
      'WJetsToLNu':{
          'HT-70to100':1264.0,
          'HT-100to200':1256.0,
          'HT-200to400':335.5,
          'HT-400to600':45.25,
          'HT-600to800':10.97,
          'HT-800to1200':4.933,
          'HT-1200to2500':1.16,
          'HT-2500toInf':0.008001},
      'QCD':{
          'HT-50to100':185300000.0,
          'HT-100to200':23590000.0,
          'HT-200to300':1551000.0,
          'HT-300to500':323400.0,
          'HT-500to700':30140.0,
          'HT-700to1000':6344.0,
          'HT-1000to1500':1092.0,
          'HT-1500to2000':99.76,
          'HT-2000toInf':20.35}
# #    'QCD-10X':{
# #        'Pt-15to7000':1.34e9}
#     'WJetsToQQ':{
#         'HT-600ToInf':99.65},
# #        'HT180':3105.0},
#     'ZJetsToQQ':{
#         'HT600toInf':581.9}
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

def weightBackgroundHists(hists, files, version, study, var, Sample):
#    for sample in list(xsecs):
    for sample in Sample:
        for mass in list(xsecs[sample]):
#        for mass in list(xsecsQCD[sample]):
            if gen==True:
                filename="h_Gen_"+sample+"_"+mass+"_"+version+".root"
            else:
#                filename="h_TriggerStudy_EMu_"+sample+".root"          
#                filename="h_GenBaseline_"+sample+"_"+mass+"_"+version+".root" 
#                filename="h_TriggerStudy_"+sample+"_"+mass+"_"+version+".root"          
#                filename="h_ModeSelectionStudy_"+sample+"_"+mass+"_"+version+".root"
#                filename="h_ModeSelectionStudy_wIso_"+sample+"_"+mass+"_"+version+".root"
#                filename="h_studyMetRegion_"+sample+"_"+mass+"_"+version+".root"
#                filename="h_studyNbjetRegion_"+sample+"_"+mass+"_"+version+".root"
#                filename="h_studyWJetsRegion_JetTrig_"+sample+"_"+mass+"_"+version+".root"
#                filename="h_studyWJetsRegion_MuonTrig_"+sample+"_"+mass+"_"+version+".root"
##                filename="h_debugMuTau_HighHT_FullyLeptonic_"+sample+"_"+mass+"_"+version+".root"
##                filename="h_debugMuTau_HighHT_FullyLeptonic_JetHT_"+sample+"_"+mass+"_"+version+".root"
##                filename="h_debugMuTau_HighHT_FullyLeptonic_HTTrig_"+sample+"_"+mass+"_"+version+".root"
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
#                filename="h_debugMuTau_HighHT_"+sample+"_"+mass+"_"+version+".root"
#                filename="h_TriggerStudy_TightPtcut_"+sample+"_"+mass+"_"+version+".root"          
#                filename="h_TriggerStudy_HTMHTFirst_LooseJetPt_"+sample+"_"+mass+"_"+version+".root"          
#                filename="h_TriggerStudy_HTMHTFirst_TightJetPt_"+sample+"_"+mass+"_"+version+".root"          
#                filename="h_Baseline_wNIso_"+sample+"_"+mass+"_"+version+".root"          
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
    version = 'v8'
    iteration = 'v8'
    study ='Nominal'
    if "-2d" in opts:
        version = 'v5'
        iteration = 'v5'
        study ='2DNominal'

if "-al" in opts:
    version = 'v8'
    iteration = 'v8'
    study ='AlteredID'
    if "-2d" in opts:
        version = 'v5'
        iteration = 'v5'
        study ='2DAltered'

#study ='GenDY'
#study = 'mT'
#study = 'FullyLeptonic'

L = 41480.0
#L = 1
#-------------

gen = False
scaling = 'xsection' #'xsection' else: 'unity'

#-------------

#Sample = ['QCD', 'TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf', 'TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf']
#Sample = ['TTJets', 'TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf', 'TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf']
#Sample = ['WJetsToLNu', 'TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf', 'TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf','QCD','DYJetsToLL']
#Sample = ['QCD', 'TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf', 'TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf']

Sample = ['TTJets']

if "-t" in opts:
    Sample = ['TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf', 'TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf']
if "-l" in opts:
    Sample = ['TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf']
if "-h" in opts:
    Sample = ['TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf']
if "-a" in opts:
    Sample = ['TCP_Ntuple_m50_HT-100to400', 'TCP_Ntuple_m50_HT-400toInf','TCP_Ntuple_m30_HT-100to400', 'TCP_Ntuple_m30_HT-400toInf','QCD','WJetsToLNu','DYJetsToLL','TTJets','DYJetsToLL_M-4to50','Diboson']
if "-b" in opts:
    Sample = ['QCD','WJetsToLNu','DYJetsToLL','TTJets','DYJetsToLL_M-4to50','Diboson']
if "-ht" in opts:
    Sample = ['QCD','WJetsToLNu','DYJetsToLL','DYJetsToLL_M-4to50']
if "-d" in opts:
    Sample = ['DYJetsToLL','DYJetsToLL_M-4to50']
if "-q" in opts:
    Sample = ['HT-50to100','HT-100to200','HT-200to300','HT-300to500','HT-500to700','HT-700to1000','HT-1000to1500','HT-1500to2000','HT-2000toInf']


#VARIABLE = ["TauPtJetPt","TauPtJet2Pt","TauPtMuonPt","TauPtdRl","TauPtdRj2tau","TauPtdRjtau","TauPtdRjmu","MuonPtdRl","TauPtdRgenMu","MuonPtdRgenMu"]
#VARIABLE = ['TauPtdRjmu','TauPtdRjtau','TauPtdRl','MuonPtdRl','TauPtMuonPt','TauPtJetPt','TauPtJet2Pt','TauPtdRj2tau','DimuonMass','MuonPtMuon2Pt','TauPtdRl2','TauPtdRgenMu','MuonPtdRgenMu']
#VARIABLE = ['TauPt','Nj','JetPt', 'MuonPt', 'MetPt']
VARIABLE = ['MetPt']
#VARIABLE = ['TauPt', 'TauPt0','TauPt1','TauPt10','Nj','JetPt', 'MuonPt']
#VARIABLE = ['TauPtMass','TauPt0Mass','TauPt1Mass','TauPt10Mass','NJetMass']

if "-2dmass" in opts:
    REGION = ['hMuTau_SR_highMET_lowMt_dRcut']

if "-dr3" in opts:
    REGION = ["hMuTau_SS_dRcut_highMET_lowMt"]
    if "-2d" in opts:
        REGION = ["hMuTau_SS_highMET_lowMt_dRcut","hMuTau_SS_highMET_lowMt"]
        if "-dRl" in opts:
            REGION = ["hMuTau_SS_highMET_lowMt_dRcutl"]

if "-sr" in opts:
    REGION = ["hMuTau_SR_dRcut_highMET_lowMt"]
    if "-2d" in opts:
        REGION = ["hMuTau_SR_highMET_lowMt_dRcut","hMuTau_SR_highMET_lowMt"]

if "-dr1" in opts:
#    REGION = ['hMuTau_highMt_dRcut_highMET', 'hMuTau_highMt_highMET']
    REGION = ['hMuTau_highMt_dRcut']

if "-dr2" in opts:
    REGION = ["hMuTau_lowMET_dRcut_lowMt","hMuTau_lowMET_lowMt"]
#    REGION = ["hMuTau_lowMET_dRcut_lowMt"]
    if "-2d" in opts:
        REGION = ["hMuTau_OS_lowMET_lowMt_dRcut","hMuTau_OS_lowMET_lowMt"]
        if "-dRl" in opts:
            REGION = ["hMuTau_OS_lowMET_lowMt_dRcutl"]
        if "-dRjt" in opts:
            REGION = ["hMuTau_OS_lowMET_lowMt_dRcutjt"]
        if "-dRjm" in opts:
            REGION = ["hMuTau_OS_lowMET_lowMt_dRcutjm"]

if "-dr4" in opts:
    REGION = ['hMuTau_lowMET_dRcut_highMt','hMuTau_lowMET_highMt']

if "-dr5" in opts:
    REGION = ['hMuTau_SS_dRcut_highMET_highMt']

if "-dr6" in opts:
#    REGION = ["hMuTau_SS_lowMET_dRcut_lowMt","hMuTau_SS_lowMET_lowMt"]
    REGION = ["hMuTau_SS_lowMET_dRcut_lowMt"]
    if "-2d" in opts:
        REGION = ["hMuTau_SS_lowMET_lowMt_dRcut","hMuTau_SS_lowMET_lowMt"]
        if "-dRl" in opts:
            REGION = ["hMuTau_SS_lowMET_lowMt_dRcutl"]
        if "-dRjt" in opts:
            REGION = ["hMuTau_SS_lowMET_lowMt_dRcutjt"]
        if "-dRjm" in opts:
            REGION = ["hMuTau_SS_lowMET_lowMt_dRcutjm"]

if "-dr7" in opts:
    REGION = ['hMuTau_SS_lowMET_dRcut_highMt','hMuTau_SS_lowMET_highMt']

histlist = []

for region in REGION:
    histlist.append(region)
    if "-2d" not in opts: histlist.append(region)
    for variable in VARIABLE:
        histlist.append(region+"_"+variable)
#        histlist.append(region+"_"+variable+"_loosedR")

# histlist = ['hMuTau_SS','hMuTau_SS_TauPt','hMuTau_SS_Nj']
# histlist = ['hMuMu_lowMET','hMuMu_lowMET_Muon1Pt','hMuMu_lowMET_Muon2Pt','hMuMu_lowMET_JetPt','hMuMu_lowMET_MetPt','hMuMu_lowMET_Nj','hMuTau_lowMET','hMuTau_lowMET_TauPt','hMuTau_lowMET_JetPt','hMuTau_lowMET_Nj','hMuTau_lowMET_MuonPt']

# #-------Tau Pt
# histlist = ['hMuTau_SS_TauPt','hMuTau_SS_lowMET_TauPt','hMuTau_SS_lowMET_highMt_TauPt','hMuTau_SS_lowMET_lowMt_TauPt','hMuTau_SS_lowMET_dRcut_TauPt','hMuTau_SS_lowMET_dRcut_highMt_TauPt','hMuTau_SS_lowMET_dRcut_lowMt_TauPt','hMuTau_SS_dRcut_TauPt','hMuTau_SS_dRcut_highMET_TauPt','hMuTau_SS_dRcut_highMET_highMt_TauPt','hMuTau_SS_dRcut_highMET_lowMt_TauPt']

# #----------------SS
# if "-ss" in opts:
#     histlist = ['hMuTau_SS_TauPtJetPt','hMuTau_SS_TauPtJet2Pt','hMuTau_SS_TauPtMuonPt']

# #----------------DR1
# if "-dr1" in opts:
#     histlist = ['hMuTau_highMt_dRcut_highMET','hMuTau_highMt_dRcut_highMET_TauPt','hMuTau_highMt_dRcut_highMET_JetPt','hMuTau_highMt_dRcut_highMET_Nj','hMuTau_highMt_highMET','hMuTau_highMt_highMET_TauPt','hMuTau_highMt_highMET_JetPt','hMuTau_highMt_highMET_Nj','hMuTau_highMt_highMET_dRl','hMuTau_highMt_highMET_dRj','hMuTau_highMt_dRcut_highMET_TauPt0','hMuTau_highMt_dRcut_highMET_TauPt10','hMuTau_highMt_dRcut_highMET_TauPt1']
# #----------------DR2
# if "-dr2" in opts:
#     histlist = ['hMuTau_lowMET_dRcut_lowMt','hMuTau_lowMET_dRcut_lowMt_TauPt','hMuTau_lowMET_dRcut_lowMt_JetPt','hMuTau_lowMET_dRcut_lowMt_Nj','hMuTau_lowMET_lowMt','hMuTau_lowMET_lowMt_TauPt','hMuTau_lowMET_lowMt_JetPt','hMuTau_lowMET_lowMt_Nj','hMuTau_lowMET_lowMt_dRl','hMuTau_lowMET_lowMt_dRj','hMuTau_lowMET_dRcut_lowMt_MuonPt','hMuTau_lowMET_dRcut_lowMt_TauPt0','hMuTau_lowMET_dRcut_lowMt_TauPt1','hMuTau_lowMET_dRcut_lowMt_TauPt10']
# #----------------DR3
# if "-dr3" in opts:
#     histlist = ['hMuTau_SS_dRcut_highMET_lowMt','hMuTau_SS_dRcut_highMET_lowMt_TauPt','hMuTau_SS_dRcut_highMET_lowMt_JetPt','hMuTau_SS_dRcut_highMET_lowMt_Nj','hMuTau_SS_dRcut_MetPt','hMuTau_SS_dRcut_highMET_Mt','hMuTau_SS_dRcut_highMET_lowMt_TauPt1','hMuTau_SS_dRcut_highMET_lowMt_TauPt10','hMuTau_SS_dRcut_highMET_lowMt_TauPt0','hMuTau_SS_dRcut_highMET_lowMt_MuonPt','hMuTau_SS_dRcut_highMET_lowMt_MetPt']
#     if "-2d" in opts:
#         histlist = ['hMuTau_SS_highMET_lowMt_TauPtJetPt','hMuTau_SS_highMET_lowMt_TauPtJet2Pt','hMuTau_SS_highMET_lowMt_TauPtMuonPt', 'hMuTau_SS_highMET_lowMt_TauPtdRl','hMuTau_SS_highMET_lowMt_TauPtdRj2tau','hMuTau_SS_highMET_lowMt_TauPtdRjtau','hMuTau_SS_highMET_lowMt_TauPtdRjmu','hMuTau_SS_highMET_lowMt_dRcut_TauPtJetPt','hMuTau_SS_highMET_lowMt_dRcut_TauPtJet2Pt','hMuTau_SS_highMET_lowMt_dRcut_TauPtMuonPt', 'hMuTau_SS_highMET_lowMt_dRcut_TauPtdRl','hMuTau_SS_highMET_lowMt_dRcut_TauPtdRj2tau','hMuTau_SS_highMET_lowMt_dRcut_TauPtdRjtau','hMuTau_SS_highMET_lowMt_dRcut_TauPtdRjmu','hMuTau_SS_highMET_lowMt_dRcut_TauPtdRgenMu','hMuTau_SS_highMET_lowMt_TauPtdRgenMu','hMuTau_SS_highMET_lowMt_dRcut_MuonPtdRgenMu','hMuTau_SS_highMET_lowMt_MuonPtdRgenMu']
#         if "-dRl" in opts:
#             histlist = 'hMuTau_SS_highMET_lowMt_dRcutl_TauPtJetPt','hMuTau_SS_highMET_lowMt_dRcutl_TauPtJet2Pt','hMuTau_SS_highMET_lowMt_dRcutl_TauPtMuonPt', 'hMuTau_SS_highMET_lowMt_dRcutl_TauPtdRl','hMuTau_SS_highMET_lowMt_dRcutl_TauPtdRj2tau','hMuTau_SS_highMET_lowMt_dRcutl_TauPtdRjtau','hMuTau_SS_highMET_lowMt_dRcutl_TauPtdRjmu'
#         if "-dRjt" in opts:
#             histlist = 'hMuTau_SS_highMET_lowMt_dRcutjt_TauPtJetPt','hMuTau_SS_highMET_lowMt_dRcutjt_TauPtJet2Pt','hMuTau_SS_highMET_lowMt_dRcutjt_TauPtMuonPt', 'hMuTau_SS_highMET_lowMt_dRcutjt_TauPtdRl','hMuTau_SS_highMET_lowMt_dRcutjt_TauPtdRj2tau','hMuTau_SS_highMET_lowMt_dRcutjt_TauPtdRjtau','hMuTau_SS_highMET_lowMt_dRcutjt_TauPtdRjmu'
#         if "-dRjm" in opts:
#             histlist = 'hMuTau_SS_highMET_lowMt_dRcutjm_TauPtJetPt','hMuTau_SS_highMET_lowMt_dRcutjm_TauPtJet2Pt','hMuTau_SS_highMET_lowMt_dRcutjm_TauPtMuonPt', 'hMuTau_SS_highMET_lowMt_dRcutjm_TauPtdRl','hMuTau_SS_highMET_lowMt_dRcutjm_TauPtdRj2tau','hMuTau_SS_highMET_lowMt_dRcutjm_TauPtdRjtau','hMuTau_SS_highMET_lowMt_dRcutjm_TauPtdRjmu'
# #----------------AR/SR
# if "-sr" in opts:
#     histlist = ['hMuTau_SR_dRcut_highMET_lowMt','hMuTau_SR_dRcut_highMET_lowMt_TauPt','hMuTau_SR_dRcut_highMET_lowMt_JetPt','hMuTau_SR_dRcut_highMET_lowMt_Nj','hMuTau_SR_dRcut_highMET_lowMt_MuonPt','hMuTau_SR_dRcut_highMET_lowMt_MetPt','hMuTau_SR_dRcut_highMET_lowMt_Mt','hMuTau_SR_dRcut_highMET_lowMt_TauPt0','hMuTau_SR_dRcut_highMET_lowMt_TauPt10','hMuTau_SR_dRcut_highMET_lowMt_TauPt1','hMuTau_SR_dRcut_highMET_lowMt_MetPt']
#     if "-2d" in opts:
#         histlist = ['hMuTau_SR_highMET_lowMt_dRcut_TauPtJetPt','hMuTau_SR_highMET_lowMt_dRcut_TauPtMuonPt','hMuTau_SR_highMET_lowMt_dRcut_TauPtdRjmu','hMuTau_SR_highMET_lowMt_dRcut_TauPtdRl','hMuTau_SR_highMET_lowMt_dRcut_TauPtdRjtau','hMuTau_SR_highMET_lowMt_TauPtJetPt','hMuTau_SR_highMET_lowMt_TauPtMuonPt','hMuTau_SR_highMET_lowMt_TauPtdRjmu','hMuTau_SR_highMET_lowMt_TauPtdRl','hMuTau_SR_highMET_lowMt_TauPtdRjtau','hMuTau_SR_highMET_lowMt_dRcut_TauPtJet2Pt','hMuTau_SR_highMET_lowMt_TauPtJet2Pt']
# #----------------DR4
# if "-dr4" in opts:
#     histlist = ['hMuTau_lowMET_dRcut_highMt','hMuTau_lowMET_dRcut_highMt_TauPt','hMuTau_lowMET_dRcut_highMt_JetPt','hMuTau_lowMET_dRcut_highMt_Nj','hMuTau_lowMET_highMt','hMuTau_lowMET_highMt_TauPt','hMuTau_lowMET_highMt_JetPt','hMuTau_lowMET_highMt_Nj','hMuTau_lowMET_highMt_dRl','hMuTau_lowMET_highMt_dRj','hMuTau_lowMET_dRcut_highMt_TauPt0','hMuTau_lowMET_dRcut_highMt_TauPt10','hMuTau_lowMET_dRcut_highMt_TauPt1']
# #----------------DR5
# if "-dr5" in opts:
#     histlist = ['hMuTau_SS_dRcut_highMET_highMt','hMuTau_SS_dRcut_highMET_highMt_TauPt','hMuTau_SS_dRcut_highMET_highMt_JetPt','hMuTau_SS_dRcut_highMET_highMt_Nj','hMuTau_SS_dRcut_highMET_highMt_TauPt1','hMuTau_SS_dRcut_highMET_highMt_TauPt0','hMuTau_SS_dRcut_highMET_highMt_TauPt10']
# #----------------DR6
# if "-dr6" in opts:
#     histlist = ['hMuTau_SS_lowMET_dRcut_lowMt','hMuTau_SS_lowMET_dRcut_lowMt_TauPt','hMuTau_SS_lowMET_dRcut_lowMt_JetPt','hMuTau_SS_lowMET_dRcut_lowMt_Nj','hMuTau_SS_lowMET_lowMt','hMuTau_SS_lowMET_lowMt_TauPt','hMuTau_SS_lowMET_lowMt_JetPt','hMuTau_SS_lowMET_lowMt_Nj','hMuTau_SS_lowMET_lowMt_dRj','hMuTau_SS_lowMET_lowMt_dRl','hMuTau_SS_lowMET_dRcut_lowMt_TauPt0','hMuTau_SS_lowMET_dRcut_lowMt_TauPt10','hMuTau_SS_lowMET_dRcut_lowMt_TauPt1','hMuTau_SS_lowMET_dRcut_lowMt_MuonPt','hMuTau_SS_lowMET_lowMt_MuonPt','hMuTau_SS_lowMET_dRcut_lowMt_MetPt']
# #    histlist = ['hMuTau_SS_lowMET_dRcut_lowMt_TauPtJetPt','hMuTau_SS_lowMET_dRcut_lowMt_TauPtMuonPt','hMuTau_SS_lowMET_dRcut_lowMt_TauPtMetPt','hMuTau_SS_lowMET_dRcut_lowMt_TauPtMuMuMass','hMuTau_SS_lowMET_dRcut_lowMt_TauPtdRj','hMuTau_SS_lowMET_dRcut_lowMt_TauPtdRj2']
# #----------------DR7
# if "-dr7" in opts:
#     histlist = ['hMuTau_SS_lowMET_dRcut_highMt','hMuTau_SS_lowMET_dRcut_highMt_TauPt','hMuTau_SS_lowMET_dRcut_highMt_JetPt','hMuTau_SS_lowMET_dRcut_highMt_Nj','hMuTau_SS_lowMET_highMt','hMuTau_SS_lowMET_highMt_TauPt','hMuTau_SS_lowMET_highMt_JetPt','hMuTau_SS_lowMET_highMt_Nj','hMuTau_SS_lowMET_highMt_dRl','hMuTau_SS_lowMET_highMt_dRj','hMuTau_SS_lowMET_highMt_MuonPt','hMuTau_SS_lowMET_dRcut_highMt_MuonPt']


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
        # elif name=='ST':
        #     hists[name].SetFillColor(ROOT.kMagenta-9)
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
print(histlist)
