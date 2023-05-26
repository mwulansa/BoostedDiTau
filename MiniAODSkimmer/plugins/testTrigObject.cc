// system include files
#include <memory>
#include <iostream>
#include <regex>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "BoostedDiTau/MiniAODSkimmer/interface/TrigObjectInfoDS.h"

#include "TTree.h"

using namespace edm;
using namespace std;

class TCPTrigObjectAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit TCPTrigObjectAnalyzer (const edm::ParameterSet&);
  ~TCPTrigObjectAnalyzer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
  edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescales_;

  TTree *tree;
  int event_;
  TrigObjectInfoDS *trigObjectInfoData;
};

TCPTrigObjectAnalyzer::TCPTrigObjectAnalyzer(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")))
{
  usesResource(TFileService::kSharedResource);
}

void TCPTrigObjectAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TCPTrigObjectAnalyzer::beginJob(){

  edm::Service<TFileService> fs;
  fs->mkdir( "triggerObject" );

  tree = fs->make<TTree>("TriggerObjectTree", "");

  tree->Branch("event", &event_, "event/I");

  trigObjectInfoData = new TrigObjectInfoDS();

  tree->Branch("TriggerObjects", "TrigObjectInfoDS", &trigObjectInfoData);
}


void TCPTrigObjectAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  trigObjectInfoData->clear();

  int Event = iEvent.id().event();
  event_ = Event;

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
    obj.unpackPathNames(names);
    std::vector pathNamesAll = obj.pathNames(false);
    std::vector pathNamesLast = obj.pathNames(true);
    // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
    // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
    // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
    //  std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
    TrigObjectInfo trigObj;
    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {

      if ( pathNamesAll[h].find("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v") == std::string::npos &&
	   pathNamesAll[h].find("HLT_Ele35_WPTight_Gsf_v") == std::string::npos && 
	   pathNamesAll[h].find("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v") == std::string::npos &&  
	   pathNamesAll[h].find("HLT_Ele115_CaloIdVT_GsfTrkIdT_v") == std::string::npos &&
	   pathNamesAll[h].find("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v") == std::string::npos &&
	   pathNamesAll[h].find("HLT_PFHT1050_v") == std::string::npos &&
	   pathNamesAll[h].find("HLT_PFJet500_v") == std::string::npos &&
           pathNamesAll[h].find("HLT_Photon200_v") == std::string::npos
	   ) continue;

      trigObj.pt = obj.pt();
      trigObj.eta = obj.eta();
      trigObj.mass = obj.mass();
      trigObj.phi = obj.phi();

      if (pathNamesAll[h].find("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v") != std::string::npos) trigObj.isEleJet = 1; // 1 isEleJet
      
      bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
      bool isLF = obj.hasPathName( pathNamesAll[h], true, false );
      if ( !isBoth && !isLF ) continue;

      if ( pathNamesAll[h].find("HLT_Ele35_WPTight_Gsf_v") != std::string::npos || pathNamesAll[h].find("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v") != std::string::npos ) trigObj.isIsoEle = 1; // 0 isIsoEle                                                                                                                       
      
      if ( pathNamesAll[h].find("HLT_Ele115_CaloIdVT_GsfTrkIdT_v") != std::string::npos ) trigObj.isEle = 1; // 2 isEle

      if ( pathNamesAll[h].find("HLT_PFJet500") != std::string::npos ) trigObj.isSingleJet = 1;
      if ( pathNamesAll[h].find("HLT_PFHT1050") != std::string::npos ) trigObj.isJetHT = 1;
      if ( pathNamesAll[h].find("HLT_Photon200") != std::string::npos ) trigObj.isPhoton = 1;

      trigObjectInfoData->push_back(trigObj);

    }
  }
  tree->Fill();
}


//define this as a plug-in
DEFINE_FWK_MODULE(TCPTrigObjectAnalyzer);
