import os, sys, random
import numpy as np

version = "v4"
#sampleType = 'signal'
sampleType = 'data'

if sampleType == 'signal':
    SAMPLES = ["TCP_m50_HT-0to100", "TCP_m50_HT-100to400", "TCP_m50_HT-400toInf", "TCP_m30_HT-0to100", "TCP_m30_HT-100to400", "TCP_m30_HT-400toInf"]

if sampleType == 'background' :
#    SAMPLES = ["DYJetsToLL_M-50", "DYJetsToLL_M-10to50", "TTJets", "WW", "WZ", "ZZ", "WJetsToLNu"]
#    SAMPLES = ["DYJetsToLL_M-50_HT-70to100","DYJetsToLL_M-50_HT-100to200","DYJetsToLL_M-50_HT-200to400","DYJetsToLL_M-50_HT-400to600","DYJetsToLL_M-50_HT-600to800","DYJetsToLL_M-50_HT-800to1200","DYJetsToLL_M-50_HT-1200to2500","DYJetsToLL_M-50_HT-2500toInf"]

#    SAMPLES = ["WJetsToLNu_HT-70to100","WJetsToLNu_HT-100to200","WJetsToLNu_HT-200to400","WJetsToLNu_HT-400to600","WJetsToLNu_HT-600to800","WJetsToLNu_HT-800to1200","WJetsToLNu_HT-1200to2500","WJetsToLNu_HT-2500toInf"]

    SAMPLES = ["QCD_HT-50to100","QCD_HT-100to200","QCD_HT-200to300","QCD_HT-300to500","QCD_HT-500to700","QCD_HT-700to1000","QCD_HT-1000to1500","QCD_HT-1500to2000"]

if sampleType == 'data' :
    data = ['SingleMuon']
    SAMPLES = [ d+"_2017"+p for d in data for p in ['B','C','D','E','F']]
 
DATASETS = {
    "TCP_m50_HT-0to100" : {"/TCP_m50_HT-0to100.txt"},
    "TCP_m50_HT-100to400" : {"/TCP_m50_HT-100to400.txt"},
    "TCP_m50_HT-400toInf" : {"/TCP_m50_HT-400toInf.txt"},
    "TCP_m30_HT-0to100" : {"/TCP_m30_HT-0to100.txt"},
    "TCP_m30_HT-100to400" : {"/TCP_m30_HT-100to400.txt"},
    "TCP_m30_HT-400toInf" : {"/TCP_m30_HT-400toInf.txt"},
    "DYJetsToLL_M-50_HT-70to100" : {"/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "DYJetsToLL_M-50_HT-100to200" : {"/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "DYJetsToLL_M-50_HT-200to400" : {"/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "DYJetsToLL_M-50_HT-400to600" : {"/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "DYJetsToLL_M-50_HT-600to800" : {"/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "DYJetsToLL_M-50_HT-800to1200" : {"/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "DYJetsToLL_M-50_HT-1200to2500" : {"/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "DYJetsToLL_M-50_HT-2500toInf" : {"/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "TTJets" : {"/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM"},
    "WW" : {"/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "WZ" : {"/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "ZZ" : {"/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "WJetsToLNu_HT-70to100" : {"/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "WJetsToLNu_HT-100to200" : {"/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "WJetsToLNu_HT-200to400" : {"/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "WJetsToLNu_HT-400to600" : {"/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "WJetsToLNu_HT-600to800" : {"/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "WJetsToLNu_HT-800to1200" : {"/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM"},
    "WJetsToLNu_HT-1200to2500" : {"/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "WJetsToLNu_HT-2500toInf" : {"/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM"},
    "QCD_HT-50to100" : {"/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "QCD_HT-100to200" : {"/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "QCD_HT-200to300" : {"/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "QCD_HT-300to500" : {"/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "QCD_HT-500to700" : {"/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "QCD_HT-700to1000" : {"/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "QCD_HT-1000to1500" : {"/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"},
    "QCD_HT-1500to2000" : {"/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM"}
}

#print(list(DATASETS))

for sample in SAMPLES:
    fname = "crabConfig_"+sample+"_UL17.py"
    print(fname)
    f = open(fname, "w")
    if sampleType == 'data':
        f.writelines("""
from CRABClient.UserUtilities import config                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
config = config()                                                                                                                                                                                                                                                                       
config.General.requestName = 'Ntuple_"""+sample+"""-UL2017_MiniAODv2-v1_"""+version+"""'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True

config.JobType.psetName = '../rerunTauRecoOnMiniAOD_WithClean_Custom_Data.py'
#config.JobType.numCores = 4                                                                                                                                                                                                                                                             
#config.JobType.maxMemoryMB = 10000                                                                                                                                                                                                                                                      
config.Data.inputDataset = '/"""+sample.split("_")[0]+"""/"""+sample.split("_")[1]+"""-UL2017_MiniAODv2-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
#config.Data.splitting = 'LumiBased'                                                                                                                                                                                                                                                     
#config.Data.unitsPerJob = 50                                                                                                                                                                                                                                                            
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'

config.Data.outLFNDirBase = '/store/user/mwulansa/TCPNtuple/'
config.Data.publication = False
config.Data.outputDatasetTag = 'Ntuple_"""+sample+"""-UL2017_MiniAODv2-v1_"""+version+"""'

config.Site.storageSite = 'T3_US_FNALLPC'

config.Site.ignoreGlobalBlacklist = True
""")

    if sampleType != 'data':

        for dataset in list(DATASETS[sample]):
            print(dataset)
            f = open(fname, "w")
            f.writelines("""
    from CRABClient.UserUtilities import config

    config = config()
    """)
            if sampleType != 'data':
                f.writelines("""
    config.General.requestName = 'Ntuple_"""+sample+"""_Summer20UL17_"""+version+"""'

    config.General.workArea = 'crab_projects'
    config.General.transferOutputs = True
    config.General.transferLogs = False

    config.JobType.pluginName = 'Analysis'
    config.JobType.allowUndistributedCMSSW = True
    #config.JobType.numCores = 4                                                                                                               
    #config.JobType.maxMemoryMB = 10000                                                        
    #config.JobType.maxJobRuntimeMin = 2000
    """)
            if sampleType == 'signal':
                f.writelines("""
    config.JobType.psetName = '../rerunTauRecoOnMiniAOD_WithClean_Custom.py'

    config.Data.userInputFiles = open('.."""+dataset+"""').readlines()
    config.Data.outputPrimaryDataset = 'TCP_m50_UL17'
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 100
    config.Data.outLFNDirBase = '/store/user/mwulansa/TCPNtuple'
    config.Data.publication = False

    config.Data.outputDatasetTag = 'NTuple_"""+sample+"""_Summer20UL17_"""+version+"""'

    config.Site.storageSite = 'T3_US_FNALLPC'
    config.Site.blacklist = ['T3_KR_KNU', 'T3_FR_IPNL', 'T2_TR_METU', 'T2_TW_NCHC', 'T2_BE_IIHE', 'T3_US_Baylor','T2_US_Purdue']
    """)
            if sampleType == 'background':
                f.writelines("""
    config.JobType.psetName = '../rerunTauRecoOnMiniAOD_WithClean_Custom_Backgrounds.py'

    config.Data.inputDataset = '"""+dataset+"""'
    config.Data.inputDBS = 'global'
    config.Data.splitting = 'Automatic'
    #config.Data.splitting = 'FileBased'                                                                                                                         
    #config.Data.unitsPerJob = 2

    config.Data.outLFNDirBase = '/store/user/mwulansa/TCPNtuple/'
    config.Data.publication = False
    config.Data.outputDatasetTag = 'Ntuple_"""+sample+"""_Summer20UL17_"""+version+"""'  

    config.Site.storageSite = 'T3_US_FNALLPC'
    #config.Site.blacklist = ['T3_KR_KNU', 'T3_FR_IPNL', 'T2_TR_METU', 'T2_TW_NCHC', 'T2_BE_IIHE', 'T3_US_Baylor','T2_US_Purdue']
    """)


