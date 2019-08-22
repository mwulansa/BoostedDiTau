#Sample='TCP'
#Sample='DYJetsToLL'
#Sample='DYJetsToLLNLO'
#Sample='TTJets'
Sample='ST'
#Sample="Diboson"

isCopy=True
version="v5"
isGen=False

if Sample == 'TCP':
    masses=['m10', 'm30', 'm50']
    prefix="root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCP/OutputMiniAODSIM/"
    
elif Sample == 'DYJetsToLL':
    masses=['M-5to50_HT-70to100', 'M-5to50_HT-100to200', 'M-5to50_HT-200to400', 'M-5to50_HT-400to600', 'M-5to50_HT-600toInf', 'M-50_HT-70to100', 'M-50_HT-100to200', 'M-50_HT-200to400', 'M-50_HT-400to600', 'M-50_HT-600to800', 'M-50_HT-800to1200', 'M-50_HT-1200to2500', 'M-50_HT-2500toInf']
    preSearchString="/"+Sample+"_REPLACEME_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X*/MINIAODSIM"
    prefix="root://xrootd.unl.edu/"

elif Sample == 'DYJetsToLLNLO':
    SampleText = 'DYJetsToLL'
    masses=['M-10to50', 'M-50']
    preSearchString="/"+SampleText+"_REPLACEME_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6*v1/MINIAODSIM"
    prefix="root://xrootd.unl.edu/"
    
elif Sample == 'TTJets':
    masses=["Dilept"]
    preSearchString="/"+Sample+"_REPLACEME*/RunIISummer16MiniAODv2-PUMoriond17_80X*/MINIAODSIM"
    prefix="root://xrootd.unl.edu/"

elif Sample == 'ST':
    masses=['s-channel_4f_leptonDecays','t-channel_antitop_4f_inclusiveDecays', 't-channel_top_4f_inclusiveDecays', 'tW_top_5f_inclusiveDecays','tW_antitop_5f_inclusiveDecays']
    preSearchString="/"+Sample+"_REPLACEME_*TuneCUETP8M1*/RunIISummer16MiniAODv2-PUMoriond17_80X*/MINIAODSIM"
    prefix="root://xrootd.unl.edu/"

elif Sample == 'Diboson':
    masses=['WZ', 'WW', 'ZZ']
    preSearchString="/REPLACEME_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X*/MINIAODSIM"
    prefix="root://xrootd.unl.edu/"

else:
    print "Please Specify Sample Name!"
    sys.exit()

    
