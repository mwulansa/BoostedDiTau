import subprocess
import sys,string,math,os
#import ConfigParser
import glob
import numpy as np
#from makeFileLists import *

isCopy = True
#Sample = sys.argv[1]

#Sample = 'DYJetsToLL_M-50'
#Sample = 'WJetsToLNu_HT-70to100'
#Sample = 'QCD_HT-200to300'
#Sample = 'DYJetsToLL_M-50_HT-70to100'
#Sample = 'DYJetsToLL_M-50_HT-100to200'
#Sample = 'DYJetsToLL_M-50_HT-200to400'
#Sample = 'DYJetsToLL_M-4to50_HT-400to600'
#Sample = 'DYJetsToLL_M-50_HT-600to800'
#Sample = 'DYJetsToLL_M-50_HT-800to1200'
#Sample = 'DYJetsToLL_M-50_HT-1200to2500'
#Sample = 'DYJetsToLL_M-50_HT-2500toInf'
#Sample = 'DYJetsToLL_M-10to50'
#Sample = 'TTJets'
#Sample = 'WW'
#Sample = 'WZ'
#Sample = 'ZZ'
#Sample = 'WJetsToLNu'
#Sample = 'SingleMuon'
#Sample = 'JetHT'

if len(sys.argv)>1:
    Sample = sys.argv[1]

version = 'v7'

#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_Baseline_EMu_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_TriggerStudy_TauTau_wNIso_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_TriggerStudy_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_TriggerStudy_TauTau_nIso_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_TriggerStudy_TightPtcut_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_studyVarIndependence_EMu_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_Baseline_wNIso_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_GenBaseline_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_studyVarIndependence_EMu_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_ModeSelectionStudy_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_ModeSelectionStudy_wIso_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_studyMetRegion_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_studyNbjetRegion_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_studyWJetsRegion_JetTrig_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_studyWJetsRegion_MuonTrig_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_studyWJetsRegion_"+Sample).read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_SingleMuEventCount_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_FullyLeptonic_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_FullyLeptonic_JetHT_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_FullyLeptonic_HTTrig_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_FullyLeptonic_Inclusive_Altered_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_FullyLeptonic_Inclusive_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_Inclusive_"+Sample+"_").read().split()
outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_Inclusive_Altered_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_FullyLeptonic_Data_Altered_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_FullyLeptonic_Data_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_Data_Altered_"+Sample+"_").read().split()
#outputfiles=os.popen("eos root://cmseos.fnal.gov ls /store/user/mwulansa/UL2017/ | grep h_debugMuTau_HighHT_Data_"+Sample+"_").read().split()

plotDir = "./output/"+version

if isCopy:
    for fil in outputfiles:
        print(fil)
        os.system("xrdcp root://cmseos.fnal.gov//store/user/mwulansa/UL2017/"+fil+" "+plotDir+"/"+fil)

#searchString = 'h_Baseline_EMu_'+Sample+'_*'
#searchString = 'h_TriggerStudy_TauTau_wNIso_'+Sample+'_*'
#searchString = 'h_TriggerStudy_'+Sample+'_*'
#searchString = 'h_TriggerStudy_TauTau_nIso_'+Sample+'_*'
#searchString = 'h_TriggerStudy_TightPtcut_'+Sample+'_*'
#searchString = 'h_GenBaseline_'+Sample+'_*'
#searchString = 'h_ModeSelectionStudy_'+Sample+'_*'
#searchString = 'h_ModeSelectionStudy_wIso_'+Sample+'_*'
#searchString = 'h_Baseline_wNIso_'+Sample+'_*'
#searchString = 'h_studyVarIndependence_EMu_'+Sample+'_*'
#searchString = 'h_studyMetRegion_'+Sample+'_*'
#searchString = 'h_studyNbjetRegion_'+Sample+'_*'
#searchString = 'h_studyWJetsRegion_JetTrig_'+Sample+'_*'
#searchString = 'h_studyWJetsRegion_MuonTrig_'+Sample+'_*'
#searchString = 'h_studyWJetsRegion_'+Sample+'_*'
#searchString = 'h_SingleMuEventCount_'+Sample+'_*'
#searchString = 'h_debugMuTau_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_FullyLeptonic_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_FullyLeptonic_JetHT_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_FullyLeptonic_HTTrig_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_Inclusive_'+Sample+'_*'
searchString = 'h_debugMuTau_HighHT_Inclusive_Altered_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_FullyLeptonic_Inclusive_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_FullyLeptonic_Data_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_FullyLeptonic_Inclusive_Altered_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_FullyLeptonic_Data_Altered_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_Data_Altered_'+Sample+'_*'
#searchString = 'h_debugMuTau_HighHT_Data_'+Sample+'_*'

print('ls '+plotDir+'/'+searchString)
os.system('ls '+plotDir+'/'+searchString)
if os.path.exists(searchString.replace("*",version+".root")):
    os.remove(searchString.replace("*",version+".root"))
print('hadd '+searchString.replace("*",version+".root")+' '+plotDir+'/'+searchString)
os.system('hadd '+searchString.replace("*",version+".root")+' '+plotDir+'/'+searchString)
