#ifndef __TCPAnalysis_JetInfoDS_H__
#define __TCPAnalysis_JetInfoDS_H__

#include <vector>

struct JetInfo {
  float pt, eta, phi, mass;
  float ptuncor; 
  float csvv2; 
  int id; //jetID 1, jetID + lept-veto 2 

  bool operator<(const JetInfo& e) const { return pt < e.pt; }
  
};

typedef class std::vector<JetInfo> JetInfoDS;

#endif
