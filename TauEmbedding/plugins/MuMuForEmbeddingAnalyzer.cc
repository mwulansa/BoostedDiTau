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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

using namespace edm;
using namespace std;

class MuMuForEmbeddingAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit MuMuForEmbeddingAnalyzer(const edm::ParameterSet&);
  ~MuMuForEmbeddingAnalyzer() override {}
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;

  float deltaR(float phi1, float phi2, float eta1, float eta2);

  struct jetInfo {
    jetInfo() {
      pt = eta = phi = 0.;
    }
    float pt, eta, phi;
  };

  struct mumuInfo {
    mumuInfo() {
      mass = dR = pt = eta = phi = 0.;
    }
    float mass, dR, pt, eta, phi;
  };

  struct globalInfo {
    globalInfo() {
      met = metphi = 0.;
      Nb = 0;
    }
    float met, metphi;
    int Nb;
  };

  // --------------member data ---------------
  edm::EDGetTokenT< edm::View<reco::CompositeCandidate> > ZmumuCandidates_;
  edm::EDGetTokenT< std::vector<pat::Jet> > Jets_;
  edm::EDGetTokenT< std::vector<pat::MET> > MET_;
  
  double ZMass = 91.0;
  TH1D *m_mumu_z;
  TH1D *m_mumu_boosted;
  TTree *tree;

  int event_;
  jetInfo jetInfo_;
  mumuInfo mumuInfo_;
  mumuInfo mumuInfoBoosted_;
  globalInfo globalInfo_;
};

MuMuForEmbeddingAnalyzer::MuMuForEmbeddingAnalyzer(const edm::ParameterSet& iConfig) :
  ZmumuCandidates_(consumes< edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("ZmumuCandidatesCollection"))),
  Jets_(consumes< vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("JetCollection"))),
  MET_(consumes< vector<pat::MET> > (iConfig.getParameter<edm::InputTag>("METCollection"))){
  usesResource(TFileService::kSharedResource);
}

void MuMuForEmbeddingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void MuMuForEmbeddingAnalyzer::beginJob() {
  
  edm::Service<TFileService> fs;
  fs->mkdir( "analysis" );
  
  m_mumu_z = fs->make<TH1D>("h_m_mumu_z", "", 1000, 0, 1000);
  m_mumu_boosted = fs->make<TH1D>("h_m_mumu_boosted", "", 1000, 0, 1000);
  tree = fs->make<TTree>("analysisTree", "");
  tree->Branch("event", &event_, "event/I");
  tree->Branch("jetInfo", &jetInfo_, "pt/F:eta/F:phi/F");
  tree->Branch("mumuInfo", &mumuInfo_, "mass/F:dR/F:pt/F:eta/F:phi/F");
  tree->Branch("mumuInfoBoosted", &mumuInfoBoosted_, "mass/F:dR/F:pt/F:eta/F:phi/F");
  tree->Branch("globalInfo", &globalInfo_, "met/F:metphi/F:Nb/I");
}

void MuMuForEmbeddingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  int Event = iEvent.id().event();
  
  edm::Handle< edm::View<reco::CompositeCandidate> > ZmumuCandidatesHandle;
  iEvent.getByToken(ZmumuCandidates_, ZmumuCandidatesHandle);
  edm::View<reco::CompositeCandidate> ZmumuCandidates = *ZmumuCandidatesHandle;

  edm::Handle< std::vector<pat::Jet> > JetsHandle;
  iEvent.getByToken(Jets_, JetsHandle);
  auto Jets = *JetsHandle;

  edm::Handle< std::vector<pat::MET> > METHandle;
  iEvent.getByToken(MET_, METHandle);
  auto MET = *METHandle;

  event_ = Event;
  
  int Nb = 0;
  if (Jets.size() > 0) {
    float leadingPt = -1;
    int theLeadingJet = -1;
    for (unsigned int i = 0; i < Jets.size(); ++i) {
      auto pfjet = Jets[i];
      float NHF  = pfjet.neutralHadronEnergyFraction();
      float NEMF = pfjet.neutralEmEnergyFraction();
      float CHF  = pfjet.chargedHadronEnergyFraction();
      // float MUF  = pfjet.muonEnergyFraction();
      float CEMF = pfjet.chargedEmEnergyFraction();
      auto NumConst = pfjet.chargedMultiplicity()+pfjet.neutralMultiplicity();
      // auto NumNeutralParticles =pfjet.neutralMultiplicity();
      auto CHM = pfjet.chargedMultiplicity();
      bool isLoose = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(pfjet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(pfjet.eta())>2.4) && abs(pfjet.eta())<=2.7;
      if (pfjet.pt() > leadingPt and isLoose) {
	leadingPt = pfjet.pt();
	theLeadingJet = i;
      }
      if (isLoose and pfjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.9535) {
	Nb += 1;
      }
    }
    if (theLeadingJet >= 0){
      jetInfo_.pt = Jets[theLeadingJet].pt();
      jetInfo_.eta = Jets[theLeadingJet].eta();
      jetInfo_.phi = Jets[theLeadingJet].phi();
    } else {
      jetInfo_.pt = -1;
      jetInfo_.eta = -999;
      jetInfo_.phi = -999;
    }
  } else {
    jetInfo_.pt = -1;
    jetInfo_.eta = -999;
    jetInfo_.phi = -999;
  }
  
  globalInfo_.Nb = Nb;
  globalInfo_.met = MET[0].pt();
  globalInfo_.metphi = MET[0].phi();
  
  const reco::CompositeCandidate* chosenCand = nullptr;
  float dR = 9999;
  float chosenCandDR(0.);
  const reco::CompositeCandidate* chosenZCand = nullptr;
  float massDifference = -1.0;
  float chosenZCandDR(0.);
  for (edm::View<reco::CompositeCandidate>::const_iterator iCand = ZmumuCandidates.begin(); iCand != ZmumuCandidates.end(); ++iCand) {

    float phi1 = iCand->daughter(0)->phi();
    float phi2 = iCand->daughter(1)->phi();
    float eta1 = iCand->daughter(0)->eta();
    float eta2 = iCand->daughter(1)->eta();
    
    float dr = deltaR(phi1, phi2, eta1, eta2);
    
    if (dr < dR) {
      dR = dr;
      chosenCand = &(*iCand);
      chosenCandDR = dr;
    }
    if (std::abs(ZMass - iCand->mass()) < massDifference || massDifference < 0) {
      massDifference = std::abs(ZMass - iCand->mass());
      chosenZCand = &(*iCand);
      chosenZCandDR = dr;
    }
  }
  
  m_mumu_boosted->Fill(chosenCand->mass());
  m_mumu_z->Fill(chosenZCand->mass());

  mumuInfoBoosted_.mass = chosenCand->mass();
  mumuInfo_.mass = chosenZCand->mass();
  mumuInfoBoosted_.dR = chosenCandDR;
  mumuInfo_.dR = chosenZCandDR;
  mumuInfoBoosted_.pt = chosenCand->pt();
  mumuInfo_.pt = chosenZCand->pt();
  mumuInfoBoosted_.eta = chosenCand->eta();
  mumuInfo_.eta = chosenZCand->eta();
  mumuInfoBoosted_.phi = chosenCand->phi();
  mumuInfo_.phi = chosenZCand->phi();

  tree->Fill();
}

float MuMuForEmbeddingAnalyzer::deltaR(float phi1, float phi2, float eta1, float eta2) {
  
  const float dphi = reco::deltaPhi(phi1, phi2);
  const float deta = eta1 - eta2;
  const float dr = std::sqrt(deta*deta + dphi*dphi);
  return dr;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuForEmbeddingAnalyzer);
  
