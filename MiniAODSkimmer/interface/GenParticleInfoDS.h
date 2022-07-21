#ifndef __TCPAnalysis_GenParticleInfoDS_H__
#define __TCPAnalysis_GenParticleInfoDS_H__

#include <vector>

struct GenParticleInfo {
  float pt, eta, phi, mass;
  int pdgid;
  bool isfinalstate;

  bool operator<(const GenParticleInfo& p) const { return pt < p.pt; }
  
};

typedef class std::vector<GenParticleInfo> GenParticleInfoDS;

#endif
