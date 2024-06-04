from CRABClient.UserUtilities import config
config = config()

whichSample = 'DYJetsToLL_M-50'

if whichSample == 'DYJetsToLL_M-50':
    whichDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/zhangj-tauEmbedding_lhe_2017MC_analysis_DYJetsToLL_M-50_RunIISummer20UL17MiniAODv2-c8ed305529a83de0e071de8e0d861b62/USER'

tag = whichSample+'_RunIISummer20UL17MiniAODv2'
    
######## change this your crab working directory name for each sample ######
config.General.requestName = 'tauEmbedding_gen_sim_2017MC_analysis_'+tag
########################################################

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../gen_sim.py'
#config.JobType.inputFiles = ['../PileupHistogram-goldenJSON-13tev-2016-69200ub.root', '../PileupMC.root']
config.JobType.numCores = 1
config.JobType.pyCfgParams = ["numThreads=4"]
config.JobType.maxMemoryMB = 2000
config.JobType.maxJobRuntimeMin = 2000
config.JobType.allowUndistributedCMSSW = True

######### change the input dataset name ############
config.Data.inputDataset = whichDataset
#config.Data.inputDBS = 'phys03'
#####################################################

#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

######### verify your destination directory #########
config.Data.outLFNDirBase = '/store/user/zhangj/TauEmbedding'
#####################################################

config.Data.publication = True

######## output dataset name, which will be published and used by the subsequent step, change this for different samples #######
config.Data.outputDatasetTag = 'tauEmbedding_gen_sim_2017MC_analysis_'+tag
#############################################################################################

#config.Site.whitelist = ['T2_US_Florida']
#config.Site.whitelist = ['T3_US_FNALLPC','T2_US_Florida']

config.Site.storageSite = 'T2_US_Florida'
#config.Data.ignoreLocality = True
#config.Site.ignoreGlobalBlacklist = True

