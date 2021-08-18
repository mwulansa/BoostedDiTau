import os, sys


massbins = ['15to20', '20to30', '30to50', '50to80', '80to120', '120to170', '170to300', '300to470', '470to600', '600to800', '800to1000', '1000toInf']
#massbins = ['15to30', '30to50', '50to80', '80to120', '120to170', '170to300', '300to470', '470to600', '600to800', '800to1000', '1000to1400', '1400to1800', '1800to2400', '2400to3200', '3200toInf']

for mass in massbins:

    path = os.popen('ls -d /eos/uscms/store/user/zhangj/QCD_Pt-{}_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/tauEmbedding_mumuSelection_2016MC_analysis_QCDMuEnriched_Pt-{}_RunIISummer16MiniAODv3/210628*/*/'.format(mass, mass)).read()

    #path = os.popen('ls -d /eos/uscms/store/user/zhangj/QCD_Pt_{}_TuneCUETP8M1_13TeV_pythia8/tauEmbedding_mumuSelection_2016MC_analysis_QCD_Pt-{}_RunIISummer16MiniAODv3/21062*/*/'.format(mass, mass)).read()

    path = path.replace('/eos/uscms','').replace('\n', '')

    #print repr(path)

    command = 'hadd ntuple_QCDMuEnriched_Pt-{}_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3.root `xrdfs root://cmseos.fnal.gov ls -u {} | grep .root`'.format(mass, path)

    print command

    os.system(command)
