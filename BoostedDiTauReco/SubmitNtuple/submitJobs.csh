#!/bin/tcsh

setenv SAMPLE $1
setenv MODE $2

#setenv SAMPLE DYJetsToLL_M-50
#setenv SAMPLE DYJetsToLL_M-50_HTBinned
#setenv SAMPLE DYJetsToLL_M-4to50_HTBinned
#setenv SAMPLE WJetsToLNu_HTBinned
#setenv SAMPLE QCD_HTBinned
#setenv SAMPLE DYJetsToLL_M-10to50
#setenv SAMPLE TTJets
#setenv SAMPLE WJetsToLNu
#setenv SAMPLE WW
#setenv SAMPLE WZ
#setenv SAMPLE ZZ
#setenv SAMPLE JetHT
#setenv SAMPLE SingleMuon

setenv CMSSW_BASE /uscms/home/jingyu/nobackup/TCP/boostedDiTauReco/CMSSW_12_1_1

setenv OutputPrefix root://cmseos.fnal.gov//store/user/zhangj/UL2017/

set isFlat=True

if ($SAMPLE == "SingleMuon") then
    foreach Mass (`ls filelists/$SAMPLE`)
	setenv MASS $Mass
	#setenv HT DYJetsToLL_M-4to50
	setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
	echo Number of Jobs: $NQueue
	echo ./filelists/$SAMPLE/$MASS/${MASS}_Process.txt $OutputPrefix
	condor_submit condorData.jdl
    end
    set isFlat=False
endif

if ($SAMPLE == "DYJetsToLL_M-4to50") then
    foreach Mass (`ls filelists/$SAMPLE`)
	setenv MASS $Mass
	setenv HT DYJetsToLL_M-4to50
	setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
	echo Number of Jobs: $NQueue
	echo ./filelists/$SAMPLE/$MASS/${HT}_${MASS}_Process.txt $OutputPrefix
	condor_submit condorHT.jdl
    end
    set isFlat=False
endif

if ($SAMPLE == "DY1jToLL_M-1to10") then
    foreach Mass (`ls filelists/$SAMPLE`)
        setenv MASS $Mass
        setenv HT DY1jToLL_M-1to10
        setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
        echo Number of Jobs: $NQueue
        echo ./filelists/$SAMPLE/$MASS/${HT}_${MASS}_Process.txt $OutputPrefix
        condor_submit condorHT.jdl
    end
    set isFlat=False
endif

if ($SAMPLE == "DYJetsToLL_M-50") then
    foreach Mass (`ls filelists/$SAMPLE`)
	setenv MASS $Mass
	setenv HT DYJetsToLL_M-50
	setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
	echo Number of Jobs: $NQueue
	echo ./filelists/$SAMPLE/$MASS/${HT}_${MASS}_Process.txt $OutputPrefix
	condor_submit condorHT.jdl
    end
    set isFlat=False
endif

if ($SAMPLE == "WJetsToLNu") then
    foreach Mass (`ls filelists/$SAMPLE`)
	setenv MASS $Mass
	setenv HT WJetsToLNu
	setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
	echo Number of Jobs: $NQueue
	echo ./filelists/$SAMPLE/$MASS/${HT}_${MASS}_Process.txt $OutputPrefix
	condor_submit condorHT.jdl
    end
    set isFlat=False
endif

if ($SAMPLE == "QCD") then
    foreach Mass (`ls filelists/$SAMPLE`)
	setenv MASS $Mass
	setenv HT QCD
	setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
	echo Number of Jobs: $NQueue
	echo ./filelists/$SAMPLE/$MASS/${HT}_${MASS}_Process.txt $OutputPrefix
	condor_submit condorHT.jdl
    end
    set isFlat=False
endif

if ($SAMPLE == "YMuMu") then
    foreach Mass (`ls filelists/$SAMPLE`)
	setenv MASS $Mass
	setenv HT YMuMu
	setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
	echo Number of Jobs: $NQueue
	echo ./filelists/$SAMPLE/$MASS/${MASS}_Process.txt $OutputPrefix
	condor_submit condorHT.jdl
    end
    set isFlat=False
endif

if ($isFlat == "True") then 
    setenv NQueue `ls filelists/$SAMPLE | wc -l`
    echo MODE: $MODE
    echo $NQueue
    echo ./filelists/$SAMPLE/${SAMPLE}_Process.txt $OutputPrefix
    condor_submit condor.jdl
endif 
