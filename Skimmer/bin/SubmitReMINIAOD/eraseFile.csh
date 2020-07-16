#!bin/tcsh

setenv SAMPLE DYJetsToLL

foreach Mass (`ls filelists/$SAMPLE`)
    foreach JOB (`seq 1 45`)
	eosrm store/user/mwulansa/DIS/TCPAnalysis/Backgrounds/RunIIFall17DR94Premix/slimmed_${SAMPLE}_${MASS}_${JOB}.root
    end
end
