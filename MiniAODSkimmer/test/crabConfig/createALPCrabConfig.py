import os, sys

mass = [ str(x) for x in range(10,75,5) ]
HT = ['100to400', '400toInf']
version = 'v2'

for m in mass:
    for ht in HT:
        fname = "crabConfig_ALP_M-"+m+"_HT-"+ht+".py"
        print(fname)
        f = open(fname, "w")
        f.writelines("""
from CRABClient.UserUtilities import config

config = config()
config.General.requestName = 'Ntuple_ALP_M-"""+m+"""_HT-"""+ht+"""_Summer20UL17_miniAODv2_"""+version+"""'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True

config.JobType.psetName = '../rerunTauRecoOnMiniAOD_WithClean_Custom.py'

config.Data.inputDataset = '/AToTauTau_ALP_M-"""+m+"""_HT-"""+ht+"""_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'

config.Data.outLFNDirBase = '/store/user/mwulansa/TCPNtuple/'
config.Data.publication = False
config.Data.outputDatasetTag = 'Ntuple_ALP_M-"""+m+"""_HT-"""+ht+"""_Summer20UL17_miniAODv2_"""+version+"""'

config.Site.storageSite = 'T3_US_FNALLPC'

""")
        f.close()
