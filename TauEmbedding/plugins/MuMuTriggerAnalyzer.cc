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

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace edm;
using namespace std;

class MuMuTriggerAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit MuMuTriggerAnalyzer(const edm::ParameterSet&);
  ~MuMuTriggerAnalyzer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;


  edm::EDGetTokenT< edm::TriggerResults > TriggerResults_;
  edm::EDGetTokenT< edm::View<reco::CompositeCandidate> > ZmumuCandidates_;
  
  TTree *tree;
  int event_;
  int isDoubleMu_;
  int isSingleMu_;
  int isJetHT_;
};

MuMuTriggerAnalyzer::MuMuTriggerAnalyzer(const edm::ParameterSet& iConfig) :
  TriggerResults_(consumes< edm::TriggerResults >(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  ZmumuCandidates_(consumes< edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("ZmumuCandidatesCollection"))) {
  usesResource(TFileService::kSharedResource);
}

void MuMuTriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void MuMuTriggerAnalyzer::beginJob() {
  
  edm::Service<TFileService> fs;
  
  tree = fs->make<TTree>("genTree", "");
  tree->Branch("event", &event_, "event/I");
  tree->Branch("isDoubleMu", &isDoubleMu_, "isDoubleMu/I");
  tree->Branch("isSingleMu", &isSingleMu_, "isSingleMu/I");
  tree->Branch("isJetHT", &isJetHT_, "isJetHT/I");
}

void MuMuTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  int Event = iEvent.id().event();
  event_ = Event;

  edm::Handle< edm::View<reco::CompositeCandidate> > ZmumuCandidatesHandle;
  iEvent.getByToken(ZmumuCandidates_, ZmumuCandidatesHandle);
  edm::View<reco::CompositeCandidate> ZmumuCandidates = *ZmumuCandidatesHandle;

  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByToken(TriggerResults_, triggerResultsHandle);
  TriggerResults triggerResults = *triggerResultsHandle;
  auto & names = iEvent.triggerNames(*triggerResultsHandle);

  isDoubleMu_ = 0;
  isSingleMu_ = 0;
  isJetHT_ = 0;
  
  if ( triggerResults.accept(names.triggerIndex("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7")) || triggerResults.accept(names.triggerIndex("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6")) ) {
    isDoubleMu_ = 1;
  }

  if ( triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v4")) || triggerResults.accept(names.triggerIndex("HLT_IsoTkMu24_v4")) ) {
    isSingleMu_ = 1;
  }

  if ( triggerResults.accept(names.triggerIndex("HLT_PFJet450_v9")) || triggerResults.accept(names.triggerIndex("HLT_PFHT900_v6")) ) {
    isJetHT_ = 1;
  }
  
  tree->Fill();
}

DEFINE_FWK_MODULE(MuMuTriggerAnalyzer);
