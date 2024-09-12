#!/bin/tcsh 

setenv SAMPLE $1
setenv VER $2
setenv FNAME $3

foreach Mass (`ls filelists/$SAMPLE`)
    setenv MASS `echo $Mass | sed 's/-UL2017//' | sed 's/_UL17//'`
    setenv PROCESS `echo $SAMPLE | sed 's/_HTBinned//'`
    echo ${PROCESS}_${Mass}
    if ($PROCESS == "SingleMuon") then
	echo python3 doHadd.py ${MASS} ${VER} ${FNAME}
	python3 doHadd.py ${MASS} ${VER} ${FNAME}
    else
	echo python3 doHadd.py ${PROCESS}_${MASS} ${VER} ${FNAME}
	python3 doHadd.py ${PROCESS}_${MASS} ${VER} ${FNAME}
    endif
end
