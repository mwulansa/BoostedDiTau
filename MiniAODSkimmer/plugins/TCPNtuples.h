#ifndef MiniAODSkimmer_MiniAODCleaner_TCPNtuples_h
#define MiniAODSkimmer_MiniAODCleaner_TCPNtuples_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "BoostedDiTau/MiniAODSkimmer/interface/JetInfoDS.h"
#include "BoostedDiTau/MiniAODSkimmer/interface/MuonInfoDS.h"
#include "BoostedDiTau/MiniAODSkimmer/interface/ElectronInfoDS.h"
#include "BoostedDiTau/MiniAODSkimmer/interface/TauInfoDS.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

using namespace edm;
using namespace std;

class TCPNtuples : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit TCPNtuples(const edm::ParameterSet&);
  ~TCPNtuples() override {}
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void Reset() {
    event_ = -9999;
    run_ = -9999;
    lumiblock_ = -9999;
    met_ = -9999.;
    metphi_ = -9999;
    weight1_ = -9999;
    weight2_ = -9999;
    jetInfoData->clear();
    muonInfoData->clear();
    electronInfoData->clear();
    tauInfoDataUnCleaned->clear();
    tauInfoDataECleaned->clear();
    tauInfoDataMCleaned->clear();
    tauInfoDataBoosted->clear();
  }

  JetInfoDS* jetInfoData;
  MuonInfoDS* muonInfoData;
  ElectronInfoDS* electronInfoData;
  TauInfoDS* tauInfoDataUnCleaned;
  TauInfoDS* tauInfoDataECleaned;
  TauInfoDS* tauInfoDataMCleaned;
  TauInfoDS* tauInfoDataBoosted;

private:

  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;

  TTree *tree;
  edm::EDGetTokenT< std::vector<pat::MET> > MET_;
  edm::EDGetTokenT< std::vector<pat::Jet> > Jets_;
  edm::EDGetTokenT< std::vector<pat::Muon> > Muons_;
  edm::EDGetTokenT< std::vector<pat::Electron> > Electrons_;
  edm::EDGetTokenT< std::vector<reco::Vertex> > Vertices_;
  edm::EDGetTokenT<double> rhoTag_;
  EffectiveAreas effectiveAreas_;
  edm::EDGetTokenT< std::vector<pat::Tau> > TausUnCleaned_;
  edm::EDGetTokenT< std::vector<pat::Tau> > TausECleaned_;
  edm::EDGetTokenT< std::vector<pat::Tau> > TausMCleaned_;
  edm::EDGetTokenT< std::vector<pat::Tau> > TausBoosted_;
  
  int event_;
  int run_;
  int lumiblock_;
  float weight1_;  // gen weight
  float weight2_;  // pu weight
  float met_;
  float metphi_;

  void fillTauInfoDS(const std::vector<pat::Tau>& TauCollection, int whichColl);
  float deltaR(float phi1, float phi2, float eta1, float eta2);
};

#endif
