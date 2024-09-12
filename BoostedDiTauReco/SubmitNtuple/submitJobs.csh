#!/bin/tcsh

setenv SAMPLE $1
setenv YEAR $2
setenv MODE $3

setenv CMSSW_BASE /uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_12_1_1

setenv OutputPrefix root://cmseos.fnal.gov//store/user/zhangj/UL2017/


setenv offset 0

foreach Mass (`ls filelists/$SAMPLE/UL$YEAR/`)
    setenv MASS $Mass
    setenv NQueue `ls filelists/$SAMPLE/UL${YEAR}/${MASS} | wc -l`
    # setenv NQueue 68
    # setenv NQueue 1
    echo $MASS
    setenv filelist ./filelists/${SAMPLE}/UL${YEAR}/${MASS}/${MASS}
    setenv logFile ${SAMPLE}_UL${YEAR}_${MASS}
    echo $filelist
    condor_submit condor.jdl                                                                                                        
end

# setenv MASS ALP_Ntuple_m_30_htj_400toInf_UL2018
# setenv MASS Ntuple_SingleMuon_Run2017E-UL2017_MiniAODv2-v1
# setenv MASS Ntuple_SingleMuon_Run2018A-UL2018_MiniAODv2
# setenv MASS Ntuple_WJetsToLNu_HT-2500toInf_Summer20UL18
# setenv MASS Ntuple_DYJetsToLL_M-50_HT-2500toInf_Summer20UL18

# setenv NQueue `ls filelists/SingleMuon/UL2017/${MASS} | wc -l`
# setenv NQueue `ls filelists/SingleMuon/UL2018/${MASS} | wc -l`
# setenv NQueue `ls singleMuonToResubmit/${MASS}* | wc -l`
# setenv NQueue 1
# echo ${NQueue}
# echo $MASS
# setenv filelist ./filelists/${SAMPLE}/UL${YEAR}/${MASS}/${MASS}
# setenv logFile ${SAMPLE}_UL${YEAR}_${MASS}
# echo $filelist
# condor_submit condor.jdl
