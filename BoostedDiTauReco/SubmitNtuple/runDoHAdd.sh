#!/bin/tcsh 

setenv SAMPLE $1
setenv VER $2
setenv FNAME $3

foreach Mass (`ls filelists/$SAMPLE`)
    setenv MASS $Mass
    setenv PROCESS `echo $SAMPLE | sed 's/_HTBinned//'`
    echo ${PROCESS}_${Mass}
    echo python3 doHadd.py ${PROCESS}_${Mass} ${VER} ${FNAME}
    python3 doHadd.py ${PROCESS}_${Mass} ${VER} ${FNAME}
end
