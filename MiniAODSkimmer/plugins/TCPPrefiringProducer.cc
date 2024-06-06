#include <memory>
#include <iostream>
#include <regex>

// user include files                                                                                                                                                                                                                  
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TTree.h"

using namespace edm;
using namespace std;

class TCPPrefiring : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit TCPPrefiring(const edm::ParameterSet&);
  ~TCPPrefiring() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;

  edm::EDGetTokenT< float > prefweight_token;
  edm::EDGetTokenT< float > prefweightup_token;
  edm::EDGetTokenT< float > prefweightdown_token;

  TTree *tree;

  int event_;
  float prefiringweight_;
  float prefiringweightup_;
  float prefiringweightdown_;
};

TCPPrefiring::TCPPrefiring(const edm::ParameterSet& iConfig):
  prefweight_token(consumes< float >(iConfig.getParameter<edm::InputTag>("PrefiringWeight"))),
  prefweightup_token(consumes< float >(iConfig.getParameter<edm::InputTag>("PrefiringWeightUp"))),
  prefweightdown_token(consumes< float >(iConfig.getParameter<edm::InputTag>("PrefiringWeightDown"))){
  usesResource(TFileService::kSharedResource);
}

void TCPPrefiring::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TCPPrefiring::beginJob(){

  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("prefiringTree", "");

  tree->Branch("event", &event_, "event/I");
  tree->Branch("prefiringWeight", &prefiringweight_, "prefiringWeight/F");
  tree->Branch("prefiringWeightUp", &prefiringweightup_, "prefiringWeightUp/F");
  tree->Branch("prefiringWeightDown", &prefiringweightdown_, "prefiringWeightDown/F");

}

void TCPPrefiring::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  int Event = iEvent.id().event();

  edm::Handle< float > theprefweight;
  iEvent.getByToken(prefweight_token, theprefweight ) ;
  float _prefiringweight =(*theprefweight);

  edm::Handle< float > theprefweightup;
  iEvent.getByToken(prefweightup_token, theprefweightup ) ;
  float _prefiringweightup =(*theprefweightup);

  edm::Handle< float > theprefweightdown;
  iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
  float _prefiringweightdown =(*theprefweightdown);

  event_ = Event;
  prefiringweight_ = _prefiringweight;
  prefiringweightup_ = _prefiringweightup;
  prefiringweightdown_ = _prefiringweightdown;

  tree->Fill();
}

DEFINE_FWK_MODULE(TCPPrefiring);
