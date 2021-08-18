from CRABClient.UserUtilities import config
config = config()


#whichDataset = '/DoubleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD'
#whichDataset = '/DoubleMuon/Run2016C-17Jul2018-v1/MINIAOD'
#whichDataset = '/DoubleMuon/Run2016D-17Jul2018-v1/MINIAOD'
#whichDataset = '/DoubleMuon/Run2016E-17Jul2018-v1/MINIAOD'
whichDataset = '/DoubleMuon/Run2016F-17Jul2018-v1/MINIAOD'
#whichDataset = '/DoubleMuon/Run2016G-17Jul2018-v1/MINIAOD'
#whichDataset = '/DoubleMuon/Run2016H-17Jul2018-v1/MINIAOD'

#whichDataset = '/SingleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD'
#whichDataset = '/SingleMuon/Run2016C-17Jul2018-v1/MINIAOD'
#whichDataset = '/SingleMuon/Run2016D-17Jul2018-v1/MINIAOD'
#whichDataset = '/SingleMuon/Run2016E-17Jul2018-v1/MINIAOD'
#whichDataset = '/SingleMuon/Run2016F-17Jul2018-v1/MINIAOD'
#whichDataset = '/SingleMuon/Run2016G-17Jul2018-v1/MINIAOD'
#whichDataset = '/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD'

#whichDataset = '/JetHT/Run2016B-17Jul2018_ver2-v2/MINIAOD'
#whichDataset = '/JetHT/Run2016C-17Jul2018-v1/MINIAOD'
#whichDataset = '/JetHT/Run2016D-17Jul2018-v1/MINIAOD'
#whichDataset = '/JetHT/Run2016E-17Jul2018-v1/MINIAOD'
#whichDataset = '/JetHT/Run2016F-17Jul2018-v1/MINIAOD'
#whichDataset = '/JetHT/Run2016G-17Jul2018-v1/MINIAOD'
#whichDataset = '/JetHT/Run2016H-17Jul2018-v1/MINIAOD'

tag = whichDataset.split('/')[1]+'_'+whichDataset.split('/')[2]

config.General.requestName = 'tauEmbedding_mumuSelection_2016Data_analysis_'+tag
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../runMumuAnalyzer_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = whichDataset
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 200
config.Data.lumiMask = 'Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'
config.Data.runRange = '271036-284044'
#config.Data.runRange = '273158-274106'
config.Data.publication = False
config.Data.outputDatasetTag = 'tauEmbedding_mumuSelection_2016Data_analysis_'+tag

#config.Site.whitelist = ['T3_US_FNALLPC','T2_US_Florida']
config.Site.storageSite = 'T3_US_FNALLPC'
