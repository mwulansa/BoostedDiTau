import os, sys, random
#import numpy as np

version = "v5"
 
SAMPLES = {
    'ALP_M-10_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-10_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-10_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-10_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-15_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-15_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-15_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-15_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-20_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-20_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-20_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-20_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-25_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-25_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-25_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-25_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-30_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-30_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-30_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-30_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-35_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-35_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-35_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-35_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-40_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-40_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-40_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-40_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-45_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-45_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-45_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-45_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-50_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-50_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-50_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-50_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-55_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-55_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-55_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-55_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-60_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-60_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-60_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-60_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-65_HT-100to400_preVFPUL16':'/AToTauTau_ALP_M-65_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'ALP_M-65_HT-400toInf_preVFPUL16':'/AToTauTau_ALP_M-65_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-4to50_HT-70to100_preVFPUL16':'/DYJetsToLL_M-4to50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-4to50_HT-100to200_preVFPUL16':'/DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v3/MINIAODSIM',
    'DYJetsToLL_M-4to50_HT-200to400_preVFPUL16':'/DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-4to50_HT-400to600_preVFPUL16':'/DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v3/MINIAODSIM',
    'DYJetsToLL_M-4to50_HT-600toInf_preVFPUL16':'/DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-70to100_preVFPUL16':'/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-100to200_preVFPUL16':'/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-200to400_preVFPUL16':'/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-400to600_preVFPUL16':'/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-600to800_preVFPUL16':'/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-800to1200_preVFPUL16':'/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-1200to2500_preVFPUL16':'/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-2500toInf_preVFPUL16':'/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'WJetsToLNu_HT-70to100_preVFPUL16':'/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11_ext1-v3/MINIAODSIM',
    'WJetsToLNu_HT-100to200_preVFPUL16':'/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11_ext1-v3/MINIAODSIM',
    'WJetsToLNu_HT-200to400_preVFPUL16':'/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11_ext1-v3/MINIAODSIM',
    'WJetsToLNu_HT-400to600_preVFPUL16':'/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11_ext2-v3/MINIAODSIM',
    'WJetsToLNu_HT-600to800_preVFPUL16':'/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11_ext2-v3/MINIAODSIM',
    'WJetsToLNu_HT-800to1200_preVFPUL16':'/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11_ext2-v3/MINIAODSIM',
    'WJetsToLNu_HT-1200to2500_preVFPUL16':'/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11_ext2-v3/MINIAODSIM',
    'WJetsToLNu_HT-2500toInf_preVFPUL16':'/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11_ext2-v3/MINIAODSIM',
    'TTTo2L2Nu_preVFPUL16':'/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM',
    'TTToHadronic_preVFPUL16':'/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM',
    'TTToSemiLeptonic_preVFPUL16':'/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM',
    'ST_s_preVFPUL16':'/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM',
    'ST_t_antitop_preVFPUL16':'/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v3/MINIAODSIM',
    'ST_t_top_preVFPUL16':'/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v3/MINIAODSIM',
    'ST_tW_antitop_preVFPUL16':'/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM',
    'ST_tW_top_preVFPUL16':'/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM',
    'WW_preVFPUL16':'/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM',
    'WZ_preVFPUL16':'/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM',
    'ZZ_preVFPUL16':'/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM',
    'UpsilonToMuMu_pth80To200_preVFPUL16':'/UpsilonToMuMu_pth80To200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'UpsilonToMuMu_pth200To400_preVFPUL16':'/UpsilonToMuMu_pth200To400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'UpsilonToMuMu_pth400ToInf_preVFPUL16':'/UpsilonToMuMu_pth400ToInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'UpsilonToTauTau_pth80To200_preVFPUL16':'/UpsilonToTauTau_pth80To200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'UpsilonToTauTau_pth200To400_preVFPUL16':'/UpsilonToTauTau_pth200To400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
    'UpsilonToTauTau_pth400ToInf_preVFPUL16':'/UpsilonToTauTau_pth400ToInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM',
}

SAMPLES2 = {
    'ALP_M-10_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-10_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-10_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-10_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-15_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-15_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-15_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-15_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-20_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-20_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-20_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-20_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-25_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-25_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-25_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-25_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-30_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-30_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-30_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-30_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-35_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-35_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-35_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-35_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-40_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-40_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-40_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-40_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-45_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-45_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-45_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-45_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-50_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-50_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-50_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-50_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-55_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-55_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-55_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-55_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-60_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-60_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-60_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-60_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-65_HT-100to400_postVFPUL16':'/AToTauTau_ALP_M-65_HT-100to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'ALP_M-65_HT-400toInf_postVFPUL16':'/AToTauTau_ALP_M-65_HT-400toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-4to50_HT-70to100_postVFPUL16':'/DYJetsToLL_M-4to50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-4to50_HT-100to200_postVFPUL16':'/DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-4to50_HT-200to400_postVFPUL16':'/DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
    'DYJetsToLL_M-4to50_HT-400to600_postVFPUL16':'/DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-4to50_HT-600toInf_postVFPUL16':'/DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-70to100_postVFPUL16':'/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-100to200_postVFPUL16':'/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-200to400_postVFPUL16':'/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-400to600_postVFPUL16':'/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-600to800_postVFPUL16':'/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-800to1200_postVFPUL16':'/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-1200to2500_postVFPUL16':'/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'DYJetsToLL_M-50_HT-2500toInf_postVFPUL16':'/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'WJetsToLNu_HT-70to100_postVFPUL16':'/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17_ext1-v3/MINIAODSIM',
    'WJetsToLNu_HT-100to200_postVFPUL16':'/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17_ext1-v3/MINIAODSIM',
    'WJetsToLNu_HT-200to400_postVFPUL16':'/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17_ext1-v3/MINIAODSIM',
    'WJetsToLNu_HT-400to600_postVFPUL16':'/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17_ext2-v3/MINIAODSIM',
    'WJetsToLNu_HT-600to800_postVFPUL16':'/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17_ext2-v3/MINIAODSIM',
    'WJetsToLNu_HT-800to1200_postVFPUL16':'/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17_ext2-v3/MINIAODSIM',
    'WJetsToLNu_HT-1200to2500_postVFPUL16':'/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17_ext2-v3/MINIAODSIM',
    'WJetsToLNu_HT-2500toInf_postVFPUL16':'/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17_ext2-v3/MINIAODSIM',
    'TTTo2L2Nu_postVFPUL16':'/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
    'TTToHadronic_postVFPUL16':'/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
    'TTToSemiLeptonic_postVFPUL16':'/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
    'ST_s_postVFPUL16':'/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
    'ST_t_antitop_postVFPUL16':'/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v3/MINIAODSIM',
    'ST_t_top_postVFPUL16':'/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v3/MINIAODSIM',
    'ST_tW_antitop_postVFPUL16':'/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
    'ST_tW_top_postVFPUL16':'/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
    'WW_postVFPUL16':'/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
    'WZ_postVFPUL16':'/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
    'ZZ_postVFPUL16':'/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
    'UpsilonToMuMu_pth80To200_postVFPUL16':'/UpsilonToMuMu_pth80To200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'UpsilonToMuMu_pth200To400_postVFPUL16':'/UpsilonToMuMu_pth200To400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'UpsilonToMuMu_pth400ToInf_postVFPUL16':'/UpsilonToMuMu_pth400ToInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'UpsilonToTauTau_pth80To200_postVFPUL16':'/UpsilonToTauTau_pth80To200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'UpsilonToTauTau_pth200To400_postVFPUL16':'/UpsilonToTauTau_pth200To400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
    'UpsilonToTauTau_pth400ToInf_postVFPUL16':'/UpsilonToTauTau_pth400ToInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM',
}


for sample in SAMPLES2.keys():
    fname = "crabConfig_"+sample+".py"
    print(fname)
    f = open(fname, "w")

    f.writelines("""
from CRABClient.UserUtilities import config

config = config()
config.General.requestName = 'Ntuple_"""+sample+"""_miniAODv2_"""+version+"""'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True

config.JobType.psetName = '../rerunTauRecoOnMiniAOD_WithClean_Custom_Backgrounds.py'

config.Data.inputDataset = '"""+SAMPLES2[sample]+"""'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'                                                                                                                         
config.Data.unitsPerJob = 1

config.Data.outLFNDirBase = '/store/user/zhangj/TCPNtuple/'
config.Data.publication = False
config.Data.outputDatasetTag = 'Ntuple_"""+sample+"""_miniAODv2_"""+version+"""'

config.Site.storageSite = 'T2_US_Florida'
    """)


