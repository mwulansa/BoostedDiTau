#!/bin/tcsh

#Change the username and sample accordingly
foreach sample (`gfal-ls root://cmsio2.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple`)
    setenv SAMPLE $sample
    if ($SAMPLE =~ *DYJetsToLL*HT*) then
	foreach part (`gfal-ls root://cmsio2.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple/$SAMPLE`)
	    setenv PART $part
	    if ($PART =~ *_v2024) then
		echo $PART
		rm sampleList/${PART}.txt
		setenv DATE `gfal-ls root://cmsio2.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple/$SAMPLE/$PART`
		foreach num (`gfal-ls root://cmsio2.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple/$SAMPLE/$PART/$DATE`)
		    setenv NUM $num
		    echo $NUM
		    set namestr = "root://cmsio2.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple/$SAMPLE/$PART/$DATE/$NUM/"
		    echo $namestr
		    gfal-ls root://cmsio2.rc.ufl.edu:1094//store/user/mwulansa/TCPNtuple/$SAMPLE/$PART/$DATE/$NUM > sampleList/${PART}_${NUM}.txt
		    echo sampleList/${PART}.txt
		    awk -v prefix="$namestr" '{print prefix $0}' sampleList/${PART}_${NUM}.txt >> sampleList/${PART}.txt
		end
	    endif
	end
    endif
end
