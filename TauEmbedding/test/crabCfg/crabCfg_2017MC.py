from CRABClient.UserUtilities import config
config = config()

whichSample = 'DYJetsToLL_M-50'

if whichSample == 'DYJetsToLL_M-50':
    whichDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'

tag = whichSample+'_RunIISummer20UL17MiniAODv2'
    
######## change this your crab working directory name for each sample ######
config.General.requestName = 'tauEmbedding_mumuSelection_2017MC_analysis_'+tag
########################################################

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../runMuonSelection_cfg.py'
#config.JobType.inputFiles = ['../PileupHistogram-goldenJSON-13tev-2016-69200ub.root', '../PileupMC.root']
#config.JobType.numCores = 4
#config.JobType.pyCfgParams = ["numThreads=4"]
config.JobType.maxMemoryMB = 2000
config.JobType.maxJobRuntimeMin = 2000
config.JobType.allowUndistributedCMSSW = True

######### change the input dataset name ############
config.Data.inputDataset = whichDataset
#####################################################

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10

######### verify your destination directory #########
config.Data.outLFNDirBase = '/store/user/zhangj/TauEmbedding'
#####################################################

config.Data.publication = True

######## output dataset name, which will be published and used by the subsequent step, change this for different samples #######
config.Data.outputDatasetTag = 'tauEmbedding_mumuSelection_2017MC_analysis_'+tag
#############################################################################################

#config.Site.whitelist = ['T2_US_Florida']
#config.Site.whitelist = ['T3_US_FNALLPC','T2_US_Florida']

config.Site.storageSite = 'T2_US_Florida'
#config.Data.ignoreLocality = True
#config.Site.ignoreGlobalBlacklist = True

