#!/bin/tcsh

#setenv SAMPLE DYJetsToLL
#setenv SAMPLE Diboson
#setenv SAMPLE TTJets
setenv SAMPLE QCD

foreach Mass (`ls filelists/$SAMPLE`)
    setenv MASS $Mass
    foreach f (filelists/$SAMPLE/$MASS/*)
	echo $f
	python createConfigs.py ${f} ${SAMPLE} ${MASS}
    end
end
