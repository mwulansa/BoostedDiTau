#!/bin/tcsh

setenv CMSSW_BASE /uscms/home/mwulansa/nobackup/TCPNtuple/CMSSW_12_1_0_pre3

setenv OutputPrefix root://cmseos.fnal.gov//store/user/mwulansa/UL2017/

# setenv NQueue `ls filelists/$SAMPLE | wc -l`

#./makeTarBall.csh

foreach Mass (`ls filelists/AToTauTau`)
    setenv NQueue 1
    setenv MASS $Mass
    setenv MODE -s
    echo MODE: $MODE
    echo $NQueue
    echo ./filelists/AToTauTau/$MASS/${MASS}_Process.txt $OutputPrefix
    condor_submit condor_signal.jdl
end
