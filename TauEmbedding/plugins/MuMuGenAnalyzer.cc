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

#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace edm;
using namespace std;

class MuMuGenAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit MuMuGenAnalyzer(const edm::ParameterSet&);
  ~MuMuGenAnalyzer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;

  struct genJetInfo {
    genJetInfo() {
      pt = eta = phi = 0.;
    }
    float pt, eta, phi;
  };

  struct genParticleInfo {
    genParticleInfo() {
      pdgID = 0;
      pt = eta = phi = mass = 0.;
    }
    int pdgID;
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
  genParticleInfo genLmInfo_;
  genParticleInfo genLpInfo_;
  float mLL_;
  float dRLL_;

  std::string puDataFileName_;
  std::string puMCFileName_;
};

MuMuGenAnalyzer::MuMuGenAnalyzer(const edm::ParameterSet& iConfig) :
  GenParticles_(consumes< std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("GenParticleCollection"))),
  GenJets_(consumes< std::vector<reco::GenJet> > (iConfig.getParameter<edm::InputTag>("GenJetCollection"))),
  genEventInfoToken_(consumes< GenEventInfoProduct >(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
  pileupSummaryToken_(consumes< std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfo"))),
  puDataFileName_(iConfig.getParameter<std::string>("puDataFileName")),
  puMCFileName_(iConfig.getParameter<std::string>("puMCFileName")){
  
  usesResource(TFileService::kSharedResource);
  //std::cout << "debug0 " << puDataFileName_ << ' ' << puMCFileName_ << '\n';

  LumiWeights_ = edm::LumiReWeighting(puMCFileName_, puDataFileName_, "input_Event/N_TrueInteractions", "pileup");
  
}

void MuMuGenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void MuMuGenAnalyzer::beginJob() {
  
  edm::Service<TFileService> fs;
  
  tree = fs->make<TTree>("genTree", "");
  tree->Branch("event", &event_, "event/I");
  tree->Branch("genWeight", &genWeight_, "genWeight/F");
  tree->Branch("puWeight", &puWeight_, "puWeight/F");
  tree->Branch("genJetInfo", &genJetInfo_, "pt/F:eta/F:phi/F");
  tree->Branch("genLmInfo", &genLmInfo_, "pdgID/I:pt/F:eta/F:phi/F:mass/F");
  tree->Branch("genLpInfo", &genLpInfo_, "pdgID/I:pt/F:eta/F:phi/F:mass/F");
  tree->Branch("mLL", &mLL_, "mLL/F");
  tree->Branch("dRLL", &dRLL_, "dRLL/F");
}

void MuMuGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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
  //puWeight_ = -999.;
  //std::cout << "debug2 " << Tnpv << '\n';


  if (GenJets.size() > 0) {
    genJetInfo_.pt = GenJets[0].pt();
    genJetInfo_.eta = GenJets[0].eta();
    genJetInfo_.phi = GenJets[0].phi();
  } else {
    genJetInfo_.pt = -1;
    genJetInfo_.eta = -999;
    genJetInfo_.phi = -999;
  }

  genLmInfo_.pdgID = 0;
  genLmInfo_.pt = -1;
  genLmInfo_.eta = -999;
  genLmInfo_.phi = -999;
  genLmInfo_.mass = -1;
  genLpInfo_.pdgID = 0;
  genLpInfo_.pt = -1;
  genLpInfo_.eta = -999;
  genLpInfo_.phi = -999;
  genLpInfo_.mass = -1;
  
  for (unsigned int i = 0; i < GenParticles.size(); ++i) {
    auto GenParticle = GenParticles[i];
    if (GenParticle.isHardProcess()) {
      if (GenParticle.pdgId() == 11 || GenParticle.pdgId() == 13 || GenParticle.pdgId() == 15) {
	genLmInfo_.pdgID = GenParticle.pdgId();
	genLmInfo_.pt = GenParticle.pt();
	genLmInfo_.eta = GenParticle.eta();
	genLmInfo_.phi = GenParticle.phi();
	genLmInfo_.mass = GenParticle.mass();
      }
      if (GenParticle.pdgId() == -11 || GenParticle.pdgId() == -13 || GenParticle.pdgId() == -15) {
	genLpInfo_.pdgID = GenParticle.pdgId();
	genLpInfo_.pt = GenParticle.pt();
	genLpInfo_.eta = GenParticle.eta();
	genLpInfo_.phi = GenParticle.phi();
	genLpInfo_.mass = GenParticle.mass();
      }
    }
  }
  
  mLL_ = -1;
  dRLL_ = -1;
  if (genLpInfo_.pdgID !=0 and genLmInfo_.pdgID !=0) {
    TLorentzVector lm;
    TLorentzVector lp;
    lm.SetPtEtaPhiM(genLmInfo_.pt,
		    genLmInfo_.eta,
		    genLmInfo_.phi,
		    genLmInfo_.mass);
    lp.SetPtEtaPhiM(genLpInfo_.pt,
		    genLpInfo_.eta,
		    genLpInfo_.phi,
		    genLpInfo_.mass);
    mLL_ = (lm+lp).M();
    dRLL_ = lm.DeltaR(lp);
  }
  tree->Fill();
}

DEFINE_FWK_MODULE(MuMuGenAnalyzer);
