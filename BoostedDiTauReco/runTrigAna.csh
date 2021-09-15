#!/bin/tcsh

foreach mass (30)
    foreach ht (0to100 100to400 400toInf)
	set i = 0
	foreach file (`ls filelists/TCP/${mass}/*${ht}*`)
	    #echo "python studyTrigEffSig.py $file 1 >& log1$i &"
	    #echo "python studyTrigEffSig.py $file 2 >& log2$i &"
	    echo "python studyChannelOverlapping.py $file >& log$i &"
	    @ i++
	end
    end
end
