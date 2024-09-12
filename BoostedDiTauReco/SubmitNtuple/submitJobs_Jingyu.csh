#!/bin/tcsh

setenv SAMPLE $1
setenv YEAR $2
setenv MODE $3

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

if ($SAMPLE == "SingleMuon") then
    foreach Mass (`ls filelists/$SAMPLE`)
	setenv MASS $Mass
        setenv HT SingleMuon
        setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
        echo Number of Jobs: $NQueue
        echo ./filelists/$SAMPLE/$MASS/${HT}_${MASS}_Process.txt $OutputPrefix
        condor_submit condorHT.jdl
    end
endif

if ($SAMPLE == "MuonEG") then
    foreach Mass (`ls filelists/$SAMPLE`)
	setenv MASS $Mass
        setenv HT MuonEG
        setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
        echo Number of Jobs: $NQueue
        echo ./filelists/$SAMPLE/$MASS/${HT}_${MASS}_Process.txt $OutputPrefix
        condor_submit condorHT.jdl
    end
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
