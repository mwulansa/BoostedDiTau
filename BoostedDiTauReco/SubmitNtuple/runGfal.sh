#!/bin/tcsh

#setenv SAMPLE DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
#setenv SAMPLENAME DYJetsToLL

setenv SAMPLE TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8
setenv SAMPLENAME TTJets

@ j = -1

foreach part (`gfal-ls root://cmsio7.rc.ufl.edu:1094//store/user/zhangj/TCPNtuple/$SAMPLE`)
    setenv PART $part
    if ($PART =~ *_v4) then
	echo $PART
	@ i = 0
	setenv DATE `gfal-ls root://cmsio7.rc.ufl.edu:1094//store/user/zhangj/TCPNtuple/$SAMPLE/$PART`
	echo $DATE
	foreach num (`gfal-ls root://cmsio7.rc.ufl.edu:1094//store/user/zhangj/TCPNtuple/$SAMPLE/$PART/$DATE`)
	    setenv NUM $num
	    echo $NUM
	    setenv namestr root://cmsio7.rc.ufl.edu:1094//store/user/zhangj/TCPNtuple/$SAMPLE/$PART/$DATE/$NUM/
	    foreach file (`gfal-ls $namestr`)
		if ($i % 100 == 0) then
		    @ j += 1
		    if (! -d filelists/$SAMPLENAME) then
			mkdir -p filelists/$SAMPLENAME
		    endif
		    echo $namestr$file > filelists/$SAMPLENAME/${SAMPLENAME}_${j}.txt
		    @ i += 1
		else
		    echo $namestr$file >> filelists/$SAMPLENAME/${SAMPLENAME}_${j}.txt
		    @ i += 1
		endif
	    end
	endif
    end
end
