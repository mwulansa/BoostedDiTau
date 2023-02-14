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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "BoostedDiTau/MiniAODSkimmer/interface/GenParticleInfoDS.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace edm;
using namespace std;

class GenAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit GenAnalyzer(const edm::ParameterSet&);
  ~GenAnalyzer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void Reset() {
    genJetInfo_.pt = -9999;
    genJetInfo_.eta = -9999;
    genJetInfo_.phi = -9999;
    genJetInfo_.mass = -9999;
    genParticleInfoData->clear();
  }

  GenParticleInfoDS* genParticleInfoData; 

private:

  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;

  struct genJetInfo {
    genJetInfo() {
      pt = eta = phi = mass = 0.;
    }
    float pt, eta, phi, mass;
  };

  edm::EDGetTokenT< std::vector<reco::GenParticle> > GenParticles_;
  edm::EDGetTokenT< std::vector<reco::GenJet> > GenJets_;
  edm::EDGetTokenT< GenEventInfoProduct > genEventInfoToken_;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  
  edm::LumiReWeighting LumiWeights_;
  
  TTree *tree;
  int event_;
  float genWeight_;
  float puWeight_;
  genJetInfo genJetInfo_;

  std::string puDataFileName_;
  std::string puMCFileName_;
};

GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig) :
  GenParticles_(consumes< std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("GenParticleCollection"))),
  GenJets_(consumes< std::vector<reco::GenJet> > (iConfig.getParameter<edm::InputTag>("GenJetCollection"))),
  genEventInfoToken_(consumes< GenEventInfoProduct >(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
  pileupSummaryToken_(consumes< std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfo"))),
  puDataFileName_((iConfig.getParameter<edm::FileInPath>("puDataFileName")).fullPath()),
  puMCFileName_((iConfig.getParameter<edm::FileInPath>("puMCFileName")).fullPath()){
  
  usesResource(TFileService::kSharedResource);
  //std::cout << "debug0 " << puDataFileName_ << ' ' << puMCFileName_ << '\n';

  LumiWeights_ = edm::LumiReWeighting(puMCFileName_, puDataFileName_, "input_Event/N_TrueInteractions", "pileup");
  
}

void GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void GenAnalyzer::beginJob() {
  
  edm::Service<TFileService> fs;
  
  tree = fs->make<TTree>("genTree", "");
  tree->Branch("event", &event_, "event/I");
  tree->Branch("genWeight", &genWeight_, "genWeight/F");
  tree->Branch("puWeight", &puWeight_, "puWeight/F");
  tree->Branch("genJetInfo", &genJetInfo_, "pt/F:eta/F:phi/F:mass/F");
  
  genParticleInfoData = new GenParticleInfoDS();
  tree->Branch("GenParticleInfo", "GenParticleInfoDS", &genParticleInfoData);
}

void GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  Reset();

  int Event = iEvent.id().event();

  edm::Handle< std::vector<reco::GenParticle> > GenParticleHandle;
  iEvent.getByToken(GenParticles_, GenParticleHandle);
  auto GenParticles = *GenParticleHandle;

  edm::Handle< std::vector<reco::GenJet> > GenJetHandle;
  iEvent.getByToken(GenJets_, GenJetHandle);
  auto GenJets = *GenJetHandle;

  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(genEventInfoToken_, genEventInfo);

  float genWeight = genEventInfo->weight();

  event_ = Event;
  genWeight_ = genWeight;

  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByToken(pileupSummaryToken_, PupInfo);

  //std::cout << "debug1" << '\n';
  std::vector<PileupSummaryInfo>::const_iterator PVI;

  float Tnpv = -1;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    
    int BX = PVI->getBunchCrossing();
    
    if(BX == 0) { 
      Tnpv = PVI->getTrueNumInteractions();
      continue;
    } 
  }
  puWeight_ = LumiWeights_.weight( Tnpv );


  if (GenJets.size() > 0) {
    genJetInfo_.pt = GenJets[0].pt();
    genJetInfo_.eta = GenJets[0].eta();
    genJetInfo_.phi = GenJets[0].phi();
    genJetInfo_.mass = GenJets[0].mass();
  } 
  
  for (unsigned int i = 0; i < GenParticles.size(); ++i) {
    auto genParticle = GenParticles[i];
    GenParticleInfo g;
    g.ishardprocess = genParticle.isHardProcess();
    g.isdirecthardprocesstaudecayproductfinalstate = genParticle.isDirectHardProcessTauDecayProductFinalState();
    g.pt = genParticle.pt();
    g.eta = genParticle.eta();
    g.phi = genParticle.phi();
    g.mass = genParticle.mass();
    g.pdgid = genParticle.pdgId();
    if ((genParticle.isHardProcess() and abs(genParticle.pdgId()) == 15) or genParticle.isDirectHardProcessTauDecayProductFinalState()) {
      g.isfinalstate = true;
    }
    if (abs(genParticle.pdgId()) == 15 and genParticle.isHardProcess()) g.isfinalstate = false;
    genParticleInfoData->push_back(g);
  }
  
  tree->Fill();
}

DEFINE_FWK_MODULE(GenAnalyzer);
