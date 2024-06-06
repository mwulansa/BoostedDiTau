#ifndef __TCPAnalysis_TrigObjectInfoDS_H__
#define __TCPAnalysis_TrigObjectInfoDS_H__

#include <vector>

struct TrigObjectInfo {
  float pt, eta, phi, mass;
  int isIsoEle,  isEleJet, isEleLeg, isJetLeg, isEle, isSingleJet, isJetHT, isPhoton, isMu, isIsoMu, isMuonEGmu, isMuonEGe, isMuonEG, isMuTrig, isPhoton175, isMuonEGnoDZ, isMuonEGnoDZmu, isMuonEGnoDZe;
};

typedef class std::vector<TrigObjectInfo> TrigObjectInfoDS;

#endif
