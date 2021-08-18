// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "TH1D.h"
#include "TTree.h"

using namespace edm;
using namespace std;

class LumiAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:
  explicit LumiAnalyzer(const edm::ParameterSet&);
  ~LumiAnalyzer() override {}
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
 
  struct lumiInfo {
    lumiInfo() {
      event = 0;
      count = weight = 0.;
    }
    int event;
    float count, weight;
  };

  // --------------member data ---------------
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  
  TTree *tree;
  TH1D *Nevts;

  lumiInfo lumiInfo_;
};

LumiAnalyzer::LumiAnalyzer(const edm::ParameterSet& iConfig) :
  genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))) {
  usesResource(TFileService::kSharedResource);
}

void LumiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void LumiAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  fs->mkdir( "lumi" );
  
  Nevts = fs->make<TH1D>("Nevts", "", 2, 0, 2);
  tree = fs->make<TTree>("lumiTree", "");
  tree->Branch("lumiInfo", &lumiInfo_, "event/I:count/F:weight/F");
}

void LumiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  int Event = iEvent.id().event();

  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(genEventInfoToken_, genEventInfo);

  float genWeight = genEventInfo->weight();

  Nevts->Fill(0.5, 1);
  Nevts->Fill(1.5, genWeight);
  lumiInfo_.event = Event;
  lumiInfo_.count = 1.;
  lumiInfo_.weight = genWeight;

  tree->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(LumiAnalyzer);
  
  
