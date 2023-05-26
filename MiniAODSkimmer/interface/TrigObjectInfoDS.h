#ifndef __TCPAnalysis_TrigObjectInfoDS_H__
#define __TCPAnalysis_TrigObjectInfoDS_H__

#include <vector>

struct TrigObjectInfo {
  float pt, eta, phi, mass;
  int isIsoEle,  isEleJet,  isEle,  isSingleJet,  isJetHT,  isPhoton;
  
};

typedef class std::vector<TrigObjectInfo> TrigObjectInfoDS;

#endif
