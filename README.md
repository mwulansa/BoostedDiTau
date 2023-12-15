# BoostedDiTau

Master branch based on CMSSW_12_1_X, can be used for 2017 (and 2016 and 2018) UL analysis.

### TCP Ntuples

To setup the package:

```
cmsrel CMSSW_12_1_0_pre3
cd CMSSW_12_1_0_pre3/src
git cms-init
cmsenv
git cms-addpkg PhysicsTools/PatAlgos
rm PhysicsTools/PatAlgos/plugins/PATTauSlimmer.cc 
cp /uscms_data/d3/rhabibul/MiniAODCleaner/CMSSW_10_6_20/src/PhysicsTools/PatAlgos/plugins/PATTauSlimmer.cc PhysicsTools/PatAlgos/plugins/PATTauSlimmer.cc
git clone git@github.com:jingyucms/BoostedDiTau.git
```

To produce the n-tuple:
```
cd BoostedDiTau/MiniAODSkimmer/test/
cmsRun rerunTauRecoOnMiniAOD_WithClean_Custom.py
```

To read the n-tuple:
```
cd BoostedDiTau/MiniAODSkimmer/test/
python3 readTCPNtuples.py
```

Misc:
1) Need to change the configuration to have the correct year for pu-reweighting in GenAnalyzer. 
