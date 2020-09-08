#!/bin/tcsh

setenv MAKETARBALL 1

#setenv SAMPLE DYJetsToLL
#setenv SAMPLE Diboson
#setenv SAMPLE TTJets

setenv CMSSW_BASE /uscms/home/mwulansa/nobackup/TCP/boostedDiTauReco/CMSSW_9_4_13

if ($MAKETARBALL == 1) then

    cd $CMSSW_BASE/src

#    tar -zcvf ../../CMSSW.tgz ../../CMSSW_9_4_13/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"

    tar -zcvf ../../CMSSW9413.tgz ../../CMSSW_9_4_13/ --exclude="*.root" --exclude="*.pdf" --exclude="*.gif" --exclude=.git --exclude="*.log" --exclude="*stderr" --exclude="*stdout"

    eosrm /eos/uscms/store/user/mwulansa/DIS/TCPAnalysis/CMSSW9413.tgz

    xrdcp ../../CMSSW9413.tgz root://cmseos.fnal.gov//store/user/mwulansa/DIS/TCPAnalysis/CMSSW9413.tgz

    cd $CMSSW_BASE/src/BoostedDiTau/Skimmer/bin/SubmitReMINIAOD/

endif

condor_submit condor.jdl

# foreach Mass (`ls filelists/$SAMPLE`)
#     setenv MASS $Mass
#     setenv NQueue `ls filelists/$SAMPLE/$MASS | wc -l`
#     echo $NQueue
#     echo ./configs/${SAMPLE}_${MASS}_Process.py
#     condor_submit condor.jdl
# end

# foreach file (configs/${SAMPLE}_WZ*)
#     setenv CFG $file
#     echo $CFG
#     condor_submit condor.jdl
# end
