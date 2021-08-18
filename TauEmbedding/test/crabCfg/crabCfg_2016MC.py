from CRABClient.UserUtilities import config
config = config()

#whichSample = 'DYJetsToLL_M-1to5_HT-70to100'
#whichSample = 'DYJetsToLL_M-1to5_HT-100to200'
#whichSample = 'DYJetsToLL_M-1to5_HT-200to400'
#whichSample = 'DYJetsToLL_M-1to5_HT-400to600'
#whichSample = 'DYJetsToLL_M-1to5_HT-600toInf'

#whichSample = 'DYJetsToLL_M-5to50_HT-70to100'
#whichSample = 'DYJetsToLL_M-5to50_HT-100to200'
#whichSample = 'DYJetsToLL_M-5to50_HT-200to400'
#whichSample = 'DYJetsToLL_M-5to50_HT-400to600'
#whichSample = 'DYJetsToLL_M-5to50_HT-600toInf'

#whichSample = 'DYJetsToLL_M-50_ext1'
#whichSample = 'DYJetsToLL_M-50_ext2'

#whichSample = 'DYJetsToLLNLO_M-10to50'
#whichSample = 'DYJetsToLLNLO_M-50'

#whichSample = 'TTJets'
#whichSample = 'TTJetsNLO'
#whichSample = 'TTTo2L2Nu'

#whichSample = 'WJetsToLNu'
#whichSample = 'WJetsToLNuNLO'

#whichSample = 'WW'
#whichSample = 'ZZ'
#whichSample = 'WZ'

#whichSample = 'QCDMuEnriched_Pt-15to20'
#whichSample = 'QCDMuEnriched_Pt-20to30'
#whichSample = 'QCDMuEnriched_Pt-30to50'
#whichSample = 'QCDMuEnriched_Pt-50to80'
#whichSample = 'QCDMuEnriched_Pt-80to120'
#whichSample = 'QCDMuEnriched_Pt-120to170'
#whichSample = 'QCDMuEnriched_Pt-170to300'
#whichSample = 'QCDMuEnriched_Pt-300to470'
#whichSample = 'QCDMuEnriched_Pt-470to600'
#whichSample = 'QCDMuEnriched_Pt-600to800'
#whichSample = 'QCDMuEnriched_Pt-800to1000'
whichSample = 'QCDMuEnriched_Pt-1000toInf'

#whichSample = 'QCD_Pt-15to30'
#whichSample = 'QCD_Pt-30to50'
#whichSample = 'QCD_Pt-50to80'
#whichSample = 'QCD_Pt-80to120'
#whichSample = 'QCD_Pt-120to170'
#whichSample = 'QCD_Pt-170to300'
#whichSample = 'QCD_Pt-300to470'
#whichSample = 'QCD_Pt-470to600'
#whichSample = 'QCD_Pt-600to800'
#whichSample = 'QCD_Pt-800to1000'
#whichSample = 'QCD_Pt-1000to1400'
#whichSample = 'QCD_Pt-1400to1800'
#whichSample = 'QCD_Pt-1800to2400'
#whichSample = 'QCD_Pt-2400to3200'
#whichSample = 'QCD_Pt-3200toInf'

if whichSample == 'DYJetsToLL_M-1to5_HT-70to100':
    whichDataset = '/DYJetsToLL_M-1to5_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
elif whichSample == 'DYJetsToLL_M-1to5_HT-100to200':
    whichDataset = '/DYJetsToLL_M-T1to5_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
elif whichSample == 'DYJetsToLL_M-1to5_HT-200to400':
    whichDataset = '/DYJetsToLL_M-1to5_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
elif whichSample == 'DYJetsToLL_M-1to5_HT-400to600':
    whichDataset = '/DYJetsToLL_M-1to5_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
elif whichSample == 'DYJetsToLL_M-1to5_HT-600toInf':
    whichDataset = '/DYJetsToLL_M-1to5_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
    
elif whichSample == 'DYJetsToLL_M-5to50_HT-70to100':
    whichDataset = '/DYJetsToLL_M-5to50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'DYJetsToLL_M-5to50_HT-100to200':
    whichDataset = '/DYJetsToLL_M-5to50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'DYJetsToLL_M-5to50_HT-200to400':
    whichDataset = '/DYJetsToLL_M-5to50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'DYJetsToLL_M-5to50_HT-400to600':
    whichDataset = '/DYJetsToLL_M-5to50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'DYJetsToLL_M-5to50_HT-600toInf':
    whichDataset = '/DYJetsToLL_M-5to50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
    
elif whichSample == 'DYJetsToLL_M-50_ext1':
    whichDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'DYJetsToLL_M-50_ext2':
    whichDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'

elif whichSample == 'WJetsToLNu':
    whichDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'

elif whichSample == 'WW':
    whichDataset = '/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'

elif whichSample == 'WZ':
    whichDataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'

elif whichSample == 'ZZ':
    whichDataset = '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
    
elif whichSample == 'TTJets':
    whichDataset = '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'

elif whichSample == 'DYJetsToLLNLO_M-50':
    whichDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM'
elif whichSample == 'DYJetsToLLNLO_M-10to50':
    whichDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'

elif whichSample == 'TTJetsNLO':
    whichDataset = '/TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'

elif whichSample == 'TTTo2L2Nu':
    whichDataset = '/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'

elif whichSample == 'WJetsToLNuNLO':
    whichDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM'

elif whichSample == 'QCDMuEnriched_Pt-15to20':
    whichDataset = '/QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-20to30':
    whichDataset = '/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-30to50':
    whichDataset = '/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-50to80':
    whichDataset = '/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-80to120':
    whichDataset = '/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-120to170':
    whichDataset = '/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-170to300':
    whichDataset = '/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-300to470':
    whichDataset = '/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-470to600':
    whichDataset = '/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-600to800':
    whichDataset = '/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-800to1000':
    whichDataset = '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'
elif whichSample == 'QCDMuEnriched_Pt-1000toInf':
    whichDataset = '/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'

elif whichSample == 'QCD_Pt-15to30':
    whichDataset = '/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-30to50':
    whichDataset = '/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-50to80':
    whichDataset = '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-80to120':
    whichDataset = '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-120to170':
    whichDataset = '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-170to300':
    whichDataset = '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-300to470':
    whichDataset = '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-470to600':
    whichDataset = '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-600to800':
    whichDataset = '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-800to1000':
    whichDataset = '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-1000to1400':
    whichDataset = '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-1400to1800':
    whichDataset = '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-1800to2400':
    whichDataset = '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-2400to3200':
    whichDataset = '/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
elif whichSample == 'QCD_Pt-3200toInf':
    whichDataset = '/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
    

    

tag = whichSample+'_RunIISummer16MiniAODv3'
    
######## change this your crab working directory name for each sample ######
config.General.requestName = 'tauEmbedding_mumuSelection_2016MC_analysis_'+tag
########################################################

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../runMumuAnalyzer_cfg.py'
config.JobType.inputFiles = ['../PileupHistogram-goldenJSON-13tev-2016-69200ub.root', '../PileupMC.root']
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
config.Data.outLFNDirBase = '/store/user/zhangj/'
#####################################################

config.Data.publication = True

######## output dataset name, which will be published and used by the subsequent step, change this for different samples #######
config.Data.outputDatasetTag = 'tauEmbedding_mumuSelection_2016MC_analysis_'+tag
#############################################################################################

#config.Site.whitelist = ['T2_US_Florida']
config.Site.whitelist = ['T3_US_FNALLPC','T2_US_Florida']

config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.ignoreLocality = True
config.Site.ignoreGlobalBlacklist = True

