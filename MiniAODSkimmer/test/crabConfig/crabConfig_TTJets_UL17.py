
from CRABClient.UserUtilities import config

config = config()
config.General.requestName = 'Ntuple_TTJets_Summer20UL17_v13'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 4                                                                                                               
#config.JobType.maxMemoryMB = 10000                                                        
#config.JobType.maxJobRuntimeMin = 2000

config.JobType.psetName = '../rerunTauRecoOnMiniAOD_WithClean_Custom_Backgrounds.py'

config.Data.inputDataset = '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
#config.Data.splitting = 'FileBased'                                                                                                                         
#config.Data.unitsPerJob = 1

config.Data.outLFNDirBase = '/store/user/mwulansa/TCPNtuple/'
config.Data.publication = False
config.Data.outputDatasetTag = 'Ntuple_TTJets_Summer20UL17_v13'  

config.Site.storageSite = 'T3_US_FNALLPC'

#config.Site.blacklist = ['T3_KR_KNU', 'T3_FR_IPNL', 'T2_TR_METU', 'T2_TW_NCHC', 'T2_BE_IIHE', 'T3_US_Baylor','T2_US_Purdue']
