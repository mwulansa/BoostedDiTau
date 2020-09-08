#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "TLorentzVector.h"


#include "TMath.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/RefToBaseVector.h"

using namespace edm;
using namespace std;
//
// class declaration
//

class PATMuonBaseLineSelection : public edm::stream::EDFilter<> {
public:
  explicit PATMuonBaseLineSelection(const edm::ParameterSet&);
  ~PATMuonBaseLineSelection();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginStream(edm::StreamID) override;
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;


  edm::EDGetTokenT<pat::MuonCollection> muonSrc_;
  edm::EDGetTokenT<reco::VertexCollection> vtx_;
};

PATMuonBaseLineSelection::PATMuonBaseLineSelection(const edm::ParameterSet& iConfig):
  muonSrc_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  vtx_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex")))
{
  produces<pat::MuonCollection>( "myMuons" );
}


PATMuonBaseLineSelection::~PATMuonBaseLineSelection(){}

// ------------ method called on each new Event  ------------
bool PATMuonBaseLineSelection::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonSrc_,muons);
  unique_ptr<pat::MuonCollection> passedmuons(new pat::MuonCollection);
  
  Handle<reco::VertexCollection> Vertex;
  iEvent.getByToken(vtx_,Vertex);
  const reco::Vertex& pv=*Vertex->begin();
  math::XYZPoint p = pv.position();

  for(pat::MuonCollection::const_iterator imu = muons->begin(); imu !=muons->end(); ++imu){

    if ( (imu->pt() > 3.0) &&
	 (abs(imu->eta()) < 2.4) &&
	 imu->isLooseMuon() ) {
      double dxy = abs(imu->innerTrack()->dxy(p));
      double dz = abs(imu->innerTrack()->dz(p));
      //std::cout << imu->pt() << " | " << dxy << " | " << dz << "\n";
      if ( (dxy < 0.2) && (dz < 0.5) ) {
	passedmuons->push_back(*imu);
      }
    }
  }

  iEvent.put(move(passedmuons), "myMuons");
  return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
PATMuonBaseLineSelection::beginStream(edm::StreamID)
{}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PATMuonBaseLineSelection::endStream() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PATMuonBaseLineSelection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PATMuonBaseLineSelection);
