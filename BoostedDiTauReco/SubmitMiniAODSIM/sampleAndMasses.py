import sys

#Sample = 'TCP'
#Sample = 'DYJetsToLL_m
#Sample = 'DYJetsToLL'
#Sample = 'TTJets'
#Sample = 'ST'
#Sample = 'Diboson'
#Sample = 'QCD'
#Sample = 'DYJetsToQQ'
#Sample = 'ZJetsToQQ'
#Sample = 'WJetsToQQ'

#Sample = 'DYJetsToLL_M-1to4_HT-600toInf'
#Sample = 'DYJetsToLL_M-4to50_HT-600toInf'

#Sample = "TCP"
Sample = 'DYJetsToLL'
#Sample = 'QCD_Pt'
#Sample = "Diboson"
#Sample = "WJetsToLNu"



isHad = False
isCopy = True
version = "v2"
isGen = False

if Sample == 'TCP':
    masses=['m10', 'm30', 'm50']
    #prefix="root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCP/OutputMiniAODSIM/"
    prefix = "root://cmseos.fnal.gov//store/user/zhangj/events/ALP/RunIISummer17DR94Premix/"


elif Sample == 'DYJetsToLL':
    prefix="root://cmsxrootd.fnal.gov/"
    preSearchString = "/DYJetsToLL_REPLACEME_TuneCP5_13TeV-madgraphMLM-pythia8/zhangj-DYJetsToLL_REPLACEME_RunIIFall17MiniAODv2_v1-547050f89421e80c43b97ce4c8917c0a/USER"
    #TSample = 'DYJetsToLL'
    masses = ["M-1to4_HT-600toInf", "M-4to50_HT-600toInf", "M-50_HT-600to800", "M-50_HT-800to1200", "M-50_HT-1200to2500"]

elif Sample == 'QCD_Pt':
    prefix="root://cmsxrootd.fnal.gov/"
    preSearchString = "/QCD_Pt_REPLACEME_TuneCP5_13TeV_pythia8/zhangj-QCD_Pt_REPLACEME_RunIIFall17MiniAODv2_v1-547050f89421e80c43b97ce4c8917c0a/USER"
    #TSample = 'DYJetsToLL'
    masses = ["300to470", "470to600", "600to800", "800to1000", "1000to1400", "1400to1800", "1800to2400", "2400to3200", "3200toInf"]

elif Sample == 'WJetsToLNu':
    prefix="root://cmsxrootd.fnal.gov/"
    preSearchString = "/WJetsToLNu_REPLACEME_TuneCP5_13TeV-madgraphMLM-pythia8/zhangj-WJetsToLNu_REPLACEME_RunIIFall17MiniAODv2_v1-547050f89421e80c43b97ce4c8917c0a/USER"
    #TSample = 'DYJetsToLL'
    masses = ["HT-600To800", "HT-800To1200", "HT-1200To2500", "HT-2500ToInf"]
    
    
##elif Sample == 'DYJetsToLL_M-1to4_HT-600toInf':
##    #SampleText = 'DYJetsToLL_M-4to50_HT-400to600'
##    prefix="root://cmsxrootd.fnal.gov/"
##    preSearchString = "/DYJetsToLL_M-1to4_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8/zhangj-DYJetsToLL_M-1to4_HT-600toInf_RunIIFall17MiniAODv2_v1-547050f89421e80c43b97ce4c8917c0a/USER"
##    TSample = 'DYJetsToLL'
##    TMass = "M-1to4_HT-600toInf"
##
##elif Sample == "DYJetsToLL_M-4to50_HT-600toInf":
##    prefix="root://cmsxrootd.fnal.gov/"
##    preSearchString = "/DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8/zhangj-DYJetsToLL_M-4to50_HT-600toInf_RunIIFall17MiniAODv2_v1-547050f89421e80c43b97ce4c8917c0a/USER"
##    TSample = 'DYJetsToLL'
##    TMass = "M-4to50_HT-600toInf"

elif Sample == 'DYJetsToLL_mini':
    SampleText = 'DYJetsToLL'
    masses=['M-1to5_HT-70to100', 'M-1to5_HT-100to200', 'M-1to5_HT-200to400', 'M-1to5_HT-400to600', 'M-1to5_HT-600toInf']
    preSearchString="/"+SampleText+"_REPLACEME_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X*/MINIAODSIM"
    prefix="root://xrootd.unl.edu/"
    
#elif Sample == 'DYJetsToLL':
#    masses=['M-1to4_HT-100to200', 'M-1to4_HT-200to400', 'M-1to4_HT-400to600', 'M-1to4_HT-600toIn#f', 'M-4to50_HT-100to200', 'M-4to50_HT-400to600','M-50_HT-70to100', 'M-50_HT-100to200', 'M-50_HT-200to400', 'M-50_HT-400to600', 'M-50_HT-600to800', 'M-50_HT-800to1200', 'M-50_HT-1200to2500', 'M-50_HT-2500toInf']
#    preSearchString="/"+Sample+"_REPLACEME_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix#*94X*v11-v*/AODSIM"
#    prefix="root://xrootd.unl.edu/"

# elif Sample == 'DYJetsToLL':
#     masses=['M-4to50_HT-200to400']
#     preSearchString="/"+Sample+"_REPLACEME_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix*94X*v11-v2*/AODSIM"
#     prefix="root://xrootd.unl.edu/"

# elif Sample == 'DYJetsToLL':
#     masses=['M-4to50_HT-600toInf']
#     preSearchString="/"+Sample+"_REPLACEME_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix*94X*v11_ext1-v1*/AODSIM"
#     prefix="root://xrootd.unl.edu/"
    
elif Sample == 'TTJets':
    masses=["DiLept"]
    preSearchString="/"+Sample+"_REPLACEME_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X*/AODSIM"
    prefix="root://xrootd.unl.edu/"

elif Sample == 'ST':
    masses=['s-channel_4f_leptonDecays_mtop1755','t-channel_antitop_4f_inclusiveDecays', 't-channel_top_4f_inclusiveDecays', 'tW_top_5f_inclusiveDecays','tW_antitop_5f_inclusiveDecays']
    preSearchString="/"+Sample+"_REPLACEME_*TuneCP5_PSweights_13TeV*/RunIIFallDRPremix*94X*v11*/AODSIM"
    prefix="root://xrootd.unl.edu/"

elif Sample == 'Diboson':
    masses=['WZ', 'WW', 'ZZ'] #ZZ is on tape
    preSearchString="/REPLACEME_TuneCP5_13TeV-pythia8/zhangj-REPLACEME_RunIIFall17MiniAODv2_v1-547050f89421e80c43b97ce4c8917c0a/USER"
    #preSearchString="/REPLACEME_TuneCP5_13TeV-pythia8*/RunIIFall17DRPremix*94X*v11*/AODSIM"
    prefix="root://cmsxrootd.fnal.gov/"

elif Sample == 'QCD':
    masses=['HT50to100','HT100to200','HT200to300','HT300to500','HT500to700','HT700to1000','HT1000to1500','HT1500to2000','HT2000toInf']
    preSearchString="/"+Sample+"_REPLACEME_TuneCP5*/RunIIFall17DRPremix*_94X*v11*/AODSIM"
    prefix="root://xrootd.unl.edu/"

elif Sample == 'DYJetsToQQ':
    masses = ['HT180']
    preSearchString = "/"+Sample+"_REPLACEME_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"
    prefix = "root://xrootd.unl.edu/"

elif Sample == 'ZJetsToQQ':
    masses = ['HT600toInf']
    preSearchString="/"+Sample+"_REPLACEME_13TeV-madgraph*/RunIISummer16MiniAODv2-PUMoriond17_80X*/MINIAODSIM"
    prefix = "root://xrootd.unl.edu/"

elif Sample == 'WJetsToQQ':
    masses = ['HT-600ToInf']
    preSearchString="/"+Sample+"_REPLACEME_TuneCUETP8M1*/RunIISummer16MiniAODv2*/MINIAODSIM"
    prefix = "root://xrootc.unl.edu/"

else:
    print "Please Specify Sample Name!"
    sys.exit()

    
