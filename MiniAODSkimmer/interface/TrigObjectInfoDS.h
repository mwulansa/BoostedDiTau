#ifndef __TCPAnalysis_TrigObjectInfoDS_H__
#define __TCPAnalysis_TrigObjectInfoDS_H__

#include <vector>

struct TrigObjectInfo {
  float pt, eta, phi, mass;
  int isIsoEle,  isEleJet, isEleLeg, isJetLeg, isEle, isSingleJet, isJetHT, isPhoton, isPhoton175, isMu, isIsoMu, isMuonEGmu, isMuonEGnoDZmu, isMuonEGe, isMuonEGnoDZe, isMuonEGnoDZ, isMuonEG, isMuTrig;
};

typedef class std::vector<TrigObjectInfo> TrigObjectInfoDS;

#endif
