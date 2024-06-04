#ifndef TauAnalysis_ClassicSVfit_fastMTTNtuples_h
#define TauAnalysis_ClassicSVfit_fastMTTNtuples_h

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
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

using namespace edm;
using namespace std;

class fastMTTNtuples : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit fastMTTNtuples(const edm::ParameterSet&);
  ~fastMTTNtuples() override {}
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void Reset() {
    event_ = -9999;
    run_ = -9999;
    lumiblock_ = -9999;
    evtInfo_.mu_pt = -9999.;
    evtInfo_.mu_eta = -9999.;
    evtInfo_.mu_phi = -9999.;
    evtInfo_.mu_m = -9999.;
    evtInfo_.tau_pt = -9999.;
    evtInfo_.tau_eta = -9999.;
    evtInfo_.tau_phi = -9999.;
    evtInfo_.tau_m = -9999.;
    evtInfo_.pt = -9999.;
    evtInfo_.phi = -9999.;
    evtInfo_.ptUncor = -9999.;
    evtInfo_.phiUncor = -9999.;
    evtInfo_.ptJECUp = -9999.;
    evtInfo_.phiJECUp = -9999.;
    evtInfo_.ptJERUp = -9999.;
    evtInfo_.phiJERUp = -9999.;
    evtInfo_.ptUncUp = -9999.;
    evtInfo_.phiUncUp = -9999.;
    evtInfo_.ptJECDown = -9999.;
    evtInfo_.phiJECDown = -9999.;
    evtInfo_.ptJERDown = -9999.;
    evtInfo_.phiJERDown = -9999.;
    evtInfo_.ptUncDown = -9999.;
    evtInfo_.phiUncDown = -9999.;
    evtInfo_.covXX = -9999.;
    evtInfo_.covYY = -9999.;
    evtInfo_.covXY = -9999.;
  }

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
  edm::EDGetTokenT< std::vector<reco::Vertex> > Vertices_;
  edm::EDGetTokenT< std::vector<pat::Tau> > TausMCleaned_;
  edm::EDGetTokenT< std::vector<pat::Tau> > TausBoosted_;
  
  int event_;
  int run_;
  int lumiblock_;

  struct EvtInfo {
    EvtInfo() {
      mu_pt = tau_pt = 0.;
      mu_eta = tau_eta = 0.;
      mu_phi = tau_phi = 0.;
      mu_m = tau_m = 0.;
      pt = phi = 0.;
      ptUncor = phiUncor = 0.;
      ptJECUp = phiJECUp = 0.;
      ptJERUp = phiJERUp = 0.;
      ptUncUp = phiUncUp = 0.;
      ptJECDown = phiJECDown = 0.;
      ptJERDown = phiJERDown = 0.;
      ptUncDown = phiUncDown = 0.;
      covXX = covYY = covXY = 0.;
    }
    float mu_pt, tau_pt;
    float mu_eta, tau_eta;
    float mu_phi, tau_phi;
    float mu_m, tau_m;
    float pt, phi;
    float ptUncor, phiUncor;
    float ptJECUp, phiJECUp;
    float ptJERUp, phiJERUp;
    float ptUncUp, phiUncUp;
    float ptJECDown, phiJECDown;
    float ptJERDown, phiJERDown;
    float ptUncDown, phiUncDown;
    float covXX, covYY, covXY;
  };
  
  EvtInfo evtInfo_;
  float deltaR(float phi1, float phi2, float eta1, float eta2);
};

#endif
