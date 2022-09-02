
from CRABClient.UserUtilities import config

config = config()
config.General.requestName = 'Ntuple_SingleMuon_Run2017B-UL2017_MiniAODv2-v1_v3'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True

config.JobType.psetName = '../rerunTauRecoOnMiniAOD_WithClean_Custom_Data.py'
#config.JobType.numCores = 4
#config.JobType.maxMemoryMB = 10000

config.Data.inputDataset = '/SingleMuon/Run2017B-UL2017_MiniAODv2-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
#config.Data.splitting = 'LumiBased'                                                                                                                         
#config.Data.unitsPerJob = 50
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'

config.Data.outLFNDirBase = '/store/user/mwulansa/TCPNtuple/'
config.Data.publication = False
config.Data.outputDatasetTag = 'Ntuple_SingleMuon_Run2017B-UL2017_MiniAODv2-v1_v3'  

config.Site.storageSite = 'T3_US_FNALLPC'

config.Site.ignoreGlobalBlacklist = True

#config.Site.blacklist = ['T2_US_Nebraska','T2_IT_Pisa']
