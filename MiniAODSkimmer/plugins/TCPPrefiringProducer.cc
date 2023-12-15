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

  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;

  TTree *tree;

  int event_;
  double prefiringweight_;
  double prefiringweightup_;
  double prefiringweightdown_;
};

TCPPrefiring::TCPPrefiring(const edm::ParameterSet& iConfig):
  prefweight_token(consumes< double >(iConfig.getParameter<edm::InputTag>("PrefiringWeight"))),
  prefweightup_token(consumes< double >(iConfig.getParameter<edm::InputTag>("PrefiringWeightUp"))),
  prefweightdown_token(consumes< double >(iConfig.getParameter<edm::InputTag>("PrefiringWeightDown"))){
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
  tree->Branch("prefiringWeight", &prefiringweight_, "prefiringWeight/D");
  tree->Branch("prefiringWeightUp", &prefiringweightup_, "prefiringWeightUp/D");
  tree->Branch("prefiringWeightDown", &prefiringweightdown_, "prefiringWeightDown/D");

}

void TCPPrefiring::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  int Event = iEvent.id().event();

  edm::Handle< double > theprefweight;
  iEvent.getByToken(prefweight_token, theprefweight ) ;
  double _prefiringweight =(*theprefweight);

  edm::Handle< double > theprefweightup;
  iEvent.getByToken(prefweightup_token, theprefweightup ) ;
  double _prefiringweightup =(*theprefweightup);

  edm::Handle< double > theprefweightdown;
  iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
  double _prefiringweightdown =(*theprefweightdown);

  event_ = Event;
  prefiringweight_ = _prefiringweight;
  prefiringweightup_ = _prefiringweightup;
  prefiringweightdown_ = _prefiringweightdown;

  tree->Fill();
}

DEFINE_FWK_MODULE(TCPPrefiring);
