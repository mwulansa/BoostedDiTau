from CRABClient.UserUtilities import config
config = config()

#whichSample = 'DYJetsToLL_M-50'

#whichSample = 'DYJetsToLL_M-50_HT-600to800'
#whichSample = 'DYJetsToLL_M-50_HT-800to1200'
#whichSample = 'DYJetsToLL_M-50_HT-1200to2500'
whichSample = 'DYJetsToLL_M-50_HT-2500toInf'

#whichSample = 'DYJetsToLL_M-4to50_HT-100to200'
#whichSample = 'DYJetsToLL_M-4to50_HT-200to400'
#whichSample = 'DYJetsToLL_M-4to50_HT-400to600'
#whichSample = 'DYJetsToLL_M-4to50_HT-600toInf'

#whichSample = 'DYJetsToLL_M-1to4_HT-100to200'
#whichSample = 'DYJetsToLL_M-1to4_HT-200to400'
#whichSample = 'DYJetsToLL_M-1to4_HT-400to600'
#whichSample = 'DYJetsToLL_M-1to4_HT-600toInf'

#whichSample = 'QCD_Pt_300to470'
#whichSample = 'QCD_Pt_470to600'
#whichSample = 'QCD_Pt_600to800'
#whichSample = 'QCD_Pt_800to1000'
#whichSample = 'QCD_Pt_1000to1400'
#whichSample = 'QCD_Pt_1400to1800'
#whichSample = 'QCD_Pt_1800to2400'
#whichSample = 'QCD_Pt_2400to3200'
#whichSample = 'QCD_Pt_3200toInf'

#whichSample = 'WJetsToLNu_HT-600To800'
#whichSample = 'WJetsToLNu_HT-800To1200'
#whichSample = 'WJetsToLNu_HT-1200To2500'
#whichSample = 'WJetsToLNu_HT-2500ToInf'

#whichSample = 'WZ'
#whichSample = 'ZZ'
#whichSample = 'WW'


if whichSample == 'DYJetsToLL_M-50':
    whichDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'

elif whichSample == 'DYJetsToLL_M-50_HT-600to800':
    whichDataset = '/DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'DYJetsToLL_M-50_HT-800to1200':
    whichDataset = '/DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'DYJetsToLL_M-50_HT-1200to2500':
    whichDataset = '/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'DYJetsToLL_M-50_HT-2500toInf':
    whichDataset = '/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
    
elif whichSample == 'DYJetsToLL_M-4to50_HT-100to200':
    whichDataset = '/DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'DYJetsToLL_M-4to50_HT-200to400':
    whichDataset = '/DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM'
elif whichSample == 'DYJetsToLL_M-4to50_HT-400to600':
    whichDataset = '/DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'DYJetsToLL_M-4to50_HT-600toInf':
    whichDataset = '/DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11_ext1-v1/AODSIM'
    
elif whichSample == 'DYJetsToLL_M-1to4_HT-100to200':
    whichDataset = '/DYJetsToLL_M-1to4_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'DYJetsToLL_M-1to50_4-200to400':
    whichDataset = '/DYJetsToLL_M-1to4_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'DYJetsToLL_M-1to4_HT-400to600':
    whichDataset = '/DYJetsToLL_M-1to4_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'DYJetsToLL_M-1to4_HT-600toInf':
    whichDataset = '/DYJetsToLL_M-1to4_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'


elif whichSample == 'QCD_Pt_300to470':
    whichDataset = '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'QCD_Pt_470to600':
    whichDataset = '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'QCD_Pt_600to800':
    whichDataset = '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'QCD_Pt_800to1000':
    whichDataset = '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'QCD_Pt_1000to1400':
    whichDataset = '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'QCD_Pt_1400to1800':
    whichDataset = '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'QCD_Pt_1800to2400':
    whichDataset = '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM'
elif whichSample == 'QCD_Pt_2400to3200':
    whichDataset = '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'QCD_Pt_3200toInf':
    whichDataset = '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'


elif whichSample == 'WJetsToLNu_HT-600To800':
    whichDataset = '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'WJetsToLNu_HT-800To1200':
    whichDataset = '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'WJetsToLNu_HT-1200To2500':
    whichDataset = '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'WJetsToLNu_HT-2500ToInf':
    whichDataset = '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v3/AODSIM'

    
    
elif whichSample == 'WJetsToLNu':
    whichDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'
elif whichSample == 'TTJets':
    whichDataset = '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'

    
elif whichSample == 'WW':
    whichDataset = '/WW_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'WZ':
    whichDataset = '/WZ_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
elif whichSample == 'ZZ':
    whichDataset = '/ZZ_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'


######## change this your crab working directory name for each sample ######
config.General.requestName = whichSample+'_RunIIFall17MiniAODv2_v1'
########################################################

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../run_miniAOD_cfg.py'
config.JobType.numCores = 4
config.JobType.pyCfgParams = ["numThreads=4"]
config.JobType.maxMemoryMB = 4000
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
config.Data.outputDatasetTag = whichSample+'_RunIIFall17MiniAODv2_v1'
#############################################################################################

config.Site.whitelist = ['T2_US_Florida']

config.Site.storageSite = 'T2_US_Florida'
config.Data.ignoreLocality = True
config.Site.ignoreGlobalBlacklist = True
