#!/bin/tcsh 

setenv VER $1
setenv FNAME $2

python3 doHadd.py TTTo2L2Nu ${VER} ${FNAME}
python3 doHadd.py TTToSemiLeptonic ${VER} ${FNAME}
python3 doHadd.py TTToHadronic ${VER} ${FNAME}
python3 doHadd.py WW ${VER} ${FNAME}
python3 doHadd.py WZ ${VER} ${FNAME}
python3 doHadd.py ZZ ${VER} ${FNAME}
python3 doHadd.py ST_tW_top ${VER} ${FNAME}
python3 doHadd.py ST_tW_antitop ${VER} ${FNAME}
python3 doHadd.py ST_t-channel_top ${VER} ${FNAME}
python3 doHadd.py ST_t-channel_antitop ${VER} ${FNAME}
python3 doHadd.py ST_s-channel ${VER} ${FNAME}
python3 doHadd.py TTJets ${VER} ${FNAME}
python3 doHadd.py DYJetsToLL ${VER} ${FNAME}
hadd -f h_studySVFit_DYJetsToLL_${VER}.root output/${VER}/h_${FNAME}_DYJetsToLL_0.root output/${VER}/h_${FNAME}_DYJetsToLL_1.root output/${VER}/h_${FNAME}_DYJetsToLL_2.root output/${VER}/h_${FNAME}_DYJetsToLL_3.root
cp h_${FNAME}_DYJetsToLL_${VER}.root h_${FNAME}_DYJetsFlat_DYJetsToLL_${VER}.root
