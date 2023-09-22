#!/bin/tcsh

foreach sample (`gfal-ls root://cmsio7.rc.ufl.edu:1094//store/user/zhangj/SingleMuon/`)
    setenv SAMPLE $sample
    if ($SAMPLE =~ *WJetsToLNu_HT-200To400*) then
	rm sampleList/${SAMPLE}.txt
	echo $SAMPLE
	echo sampleList/$SAMPLE.txt
	foreach part (`gfal-ls root://cmsio7.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple/$SAMPLE`)
	    setenv PART $part
	    if ($PART =~ *_v31) then
		echo $PART
		setenv DATE `gfal-ls root://cmsio7.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple/$SAMPLE/$PART`
		foreach num (`gfal-ls root://cmsio7.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple/$SAMPLE/$PART/$DATE`)
		    setenv NUM $num
		    echo $NUM
		    set namestr = "root://cmsio7.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple/$SAMPLE/$PART/$DATE/$NUM/"
		    echo $namestr
		    gfal-ls root://cmsio7.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple/$SAMPLE/$PART/$DATE/$NUM > ${SAMPLE}_${NUM}.txt
		    awk -v prefix="$namestr" '{print prefix $0}' ${SAMPLE}_${NUM}.txt >> sampleList/${SAMPLE}.txt
		end
	    endif
	end
    endif
end
