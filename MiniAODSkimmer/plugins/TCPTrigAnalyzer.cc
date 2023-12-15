// system include files
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

class TCPTrigNtuples : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit TCPTrigNtuples(const edm::ParameterSet&);
  ~TCPTrigNtuples() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;

  //  edm::EDGetTokenT< edm::TriggerResults > triggerBits_;
  edm::EDGetTokenT< edm::TriggerResults > TriggerResults_;
  
  TTree *tree;
  int event_;
  bool isSingleJet450_;
  bool isSingleJet500_;
  bool isHT_;
  bool isHTMHT_;
  bool isIsoMu_;
  bool isIsoMuTau_;
  bool isMu_;
  bool isIsoEle_;
  bool isEle_;
  bool isEleJet_;
  bool isPhoton200_;
  bool isEleTau_;
  bool isDoubleMu_;
  bool isDoubleIsoEG_;
  bool isDoubleEG_;
  bool isMuonEG_;
  bool isDoubleTauMedium_;
  bool isDoubleTauTight_;
  bool isSingleTauMET_;
};

TCPTrigNtuples::TCPTrigNtuples(const edm::ParameterSet& iConfig) :
  //triggerBits_(consumes< edm::TriggerResults >(iConfig.getParameter<edm::InputTag>("bits"))),
  TriggerResults_(consumes< edm::TriggerResults >(iConfig.getParameter<edm::InputTag>("TriggerResults"))) {
  usesResource(TFileService::kSharedResource);
  //std::cout << "debug0" << "\n";
}

void TCPTrigNtuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TCPTrigNtuples::beginJob() {
  
  edm::Service<TFileService> fs;
  fs->mkdir( "trigger" );
  
  tree = fs->make<TTree>("triggerTree", "");
  
  tree->Branch("event", &event_, "event/I");
  
  tree->Branch("isSingleJet450", &isSingleJet450_, "isSingleJet450/O");
  tree->Branch("isSingleJet500", &isSingleJet500_, "isSingleJet500/O");
  tree->Branch("isHT", &isHT_, "isHT/O");
  tree->Branch("isHTMHT", &isHTMHT_, "isHTMHT/O");
  tree->Branch("isIsoMu", &isIsoMu_, "isIsoMu/O");
  tree->Branch("isIsoMuTau", &isIsoMuTau_, "isIsoMuTau/O");
  tree->Branch("isMu", &isMu_, "isMu/O");
  tree->Branch("isIsoEle", &isIsoEle_, "isIsoEle/O");
  tree->Branch("isEle", &isEle_, "isEle/O");
  tree->Branch("isPhoton200", &isPhoton200_, "isPhoton200/O");
  tree->Branch("isEleJet", &isEleJet_, "isEleJet/O");
  tree->Branch("isEleTau", &isEleTau_, "isEleTau/O");
  tree->Branch("isDoubleMu", &isDoubleMu_, "isDoubleMu/O");
  tree->Branch("isDoubleIsoEG", &isDoubleIsoEG_, "isDoubleIsoEG/O");
  tree->Branch("isDoubleEG", &isDoubleMu_, "isDoubleEG/O");
  tree->Branch("isMuonEG", &isMuonEG_, "isMuonEG/O");
  tree->Branch("isDoubleTauMedium", &isDoubleTauMedium_, "isDoubleTauMedium/O");
  tree->Branch("isDoubleTauTight", &isDoubleTauTight_, "isDoubleTauTight/O");
  tree->Branch("isSingleTauMET", &isSingleTauMET_, "isSingleTauMET/O");
}

void TCPTrigNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  int Event = iEvent.id().event();
  event_ = Event;

  //  edm::Handle<edm::TriggerResults> triggerBits;
  //  iEvent.getByToekn(triggerBits_, triggerBits);

  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByToken(TriggerResults_, triggerResultsHandle);

  TriggerResults triggerResults = *triggerResultsHandle;
  auto & names = iEvent.triggerNames(*triggerResultsHandle);

  //std::cout << "debug3" << "\n";

  isSingleJet450_ = 0;
  isSingleJet500_ = 0;
  isHT_ = 0;
  isHTMHT_ = 0;
  isIsoMu_ = 0;
  isIsoMuTau_ = 0;
  isMu_ = 0;
  isIsoEle_ = 0;
  isEle_ = 0;    // for 2018 data, non-iso electron trigger is available for all runs
  isEleJet_ = 0;
  isPhoton200_ = 0;
  isEleTau_ = 0;
  isDoubleMu_ = 0;
  isDoubleIsoEG_ = 0;
  isDoubleEG_ = 0;
  isMuonEG_ = 0;
  isDoubleTauMedium_ = 0;
  isDoubleTauTight_ = 0;
  isSingleTauMET_ = 0;

  //std::cout << "debug4" << "\n";

    for (unsigned i = 0, n = triggerResults.size(); i < n; ++i){

    std::string Trigger;
    Trigger = names.triggerName(i);

    std::regex SingleJet450("(HLT_PFJet450_v)(.*)");
    std::regex SingleJet500("(HLT_PFJet500_v)(.*)");
    std::regex HT("(HLT_PFHT1050_v)(.*)");
    std::regex HTMHT("(HLT_PFHT500_PFMET100_PFMHT100_IDTight_v)(.*)");
    std::regex IsoMu("(HLT_IsoMu27_v)(.*)");
    std::regex IsoMuTau("(HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v)(.*)");
    std::regex Mu("(HLT_Mu50_v)(.*)");
    std::regex IsoEle35("(HLT_Ele35_WPTight_Gsf_v)(.*)");
    std::regex IsoEle32("(HLT_Ele32_WPTight_Gsf_L1DoubleEG_v)(.*)");
    std::regex EleTau("(HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v)(.*)");
    std::regex DoubleMu("(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v)(.*)");
    std::regex DoubleIsoEG("(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v)(.*)");
    std::regex DoubleEG("(HLT_DoubleEle33_CaloIdL_MW_v)(.*)");
    std::regex Muon8EG("(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v)(.*)");
    std::regex Muon23EG("(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v)(.*)");
    std::regex DoubleTauMedium("(HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v)(.*)");
    std::regex DoubleTauTight35("(HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v)(.*)");
    std::regex DoubleTauTight40("(HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v)(.*)");
    std::regex SingleTauMET("(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v)(.*)");
    std::regex Ele("(HLT_Ele115_CaloIdVT_GsfTrkIdT_v)(.*)");
    std::regex EleJet("(HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v)(.*)");
    std::regex Photon200("(HLT_Photon200_v)(.*)");

    if ( std::regex_match(Trigger, Ele) && (triggerResults.accept(i) == 1) ) isEle_ = 1;

    if ( std::regex_match(Trigger, EleJet) && (triggerResults.accept(i) == 1) ) isEleJet_ = 1;

    if ( std::regex_match(Trigger, Photon200) && (triggerResults.accept(i) == 1) ) isPhoton200_ = 1;

    if ( std::regex_match(Trigger, SingleJet450) && (triggerResults.accept(i) == 1) ) isSingleJet450_ = 1;

    if ( std::regex_match(Trigger, SingleJet500) && (triggerResults.accept(i) == 1) ) isSingleJet500_ = 1;

    if ( std::regex_match(Trigger, HT) && (triggerResults.accept(i) == 1) ) isHT_ = 1;

    if ( std::regex_match(Trigger, HTMHT) && (triggerResults.accept(i) == 1) ) isHTMHT_ = 1;

    if ( std::regex_match(Trigger, IsoMu) && (triggerResults.accept(i) == 1) ) isIsoMu_ = 1;

    if ( std::regex_match(Trigger, IsoMuTau) && (triggerResults.accept(i) == 1) ) isIsoMuTau_ = 1;

    if ( std::regex_match(Trigger, Mu) && (triggerResults.accept(i) == 1) ) isMu_ = 1;

    if ( ( std::regex_match(Trigger, IsoEle35) && (triggerResults.accept(i) == 1)) || 
	 ( std::regex_match(Trigger, IsoEle32) && (triggerResults.accept(i) == 1)) ) isIsoEle_ = 1;

    if ( std::regex_match(Trigger, EleTau) && (triggerResults.accept(i) == 1) ) isEleTau_ = 1;

    if ( std::regex_match(Trigger, DoubleMu) && (triggerResults.accept(i) == 1) ) isDoubleMu_ = 1;

    if ( std::regex_match(Trigger, DoubleIsoEG) && (triggerResults.accept(i) == 1) ) isDoubleIsoEG_ = 1;

    if ( std::regex_match(Trigger, DoubleEG) && (triggerResults.accept(i) == 1) ) isDoubleEG_ = 1;

    if ( ( std::regex_match(Trigger, Muon8EG) && (triggerResults.accept(i) == 1)) ||
	 ( std::regex_match(Trigger, Muon23EG) && (triggerResults.accept(i) == 1)) ) isMuonEG_ = 1; 

    if ( std::regex_match(Trigger, DoubleTauMedium) && (triggerResults.accept(i) == 1) ) isDoubleTauMedium_ = 1;

    if ( ( std::regex_match(Trigger, DoubleTauTight35) && (triggerResults.accept(i) == 1)) ||
	 ( std::regex_match(Trigger, DoubleTauTight40) && (triggerResults.accept(i) == 1)) ) isDoubleTauTight_ = 1;

    if ( std::regex_match (Trigger, SingleTauMET) && (triggerResults.accept(i) == 1) ) isSingleTauMET_ = 1;

    }

  
  /*

  if ( triggerResults.accept(names.triggerIndex("HLT_PFJet450_v17")) ) isSingleJet_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_PFHT1050_v14")) ) isHT_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v8")) ) isHTMHT_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_IsoMu27_v13")) ) isIsoMu_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v8")) ) isIsoMuTau_ = 1;
 
  if ( triggerResults.accept(names.triggerIndex("HLT_Mu50_v11")) ) isMu_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_Ele35_WPTight_Gsf_v7")) ||
       triggerResults.accept(names.triggerIndex("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7")) ) isIsoEle_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v9")) ) isEleTau_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v3")) ) isDoubleMu_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v17")) ) isDoubleIsoEG_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_DoubleEle33_CaloIdL_MW_v15")) ) isDoubleEG_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10")) ||
       triggerResults.accept(names.triggerIndex("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12")) ) isMuonEG_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v8")) ) isDoubleTauMedium_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v8")) ||
       triggerResults.accept(names.triggerIndex("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v8")) ) isDoubleTauTight_ = 1;

  if ( triggerResults.accept(names.triggerIndex("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8")) ) isSingleTauMET_ = 1;
  */
  tree->Fill();
}

DEFINE_FWK_MODULE(TCPTrigNtuples);
