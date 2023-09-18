// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//
class EmbeddingPFCandProducer : public edm::stream::EDProducer<> {
   public:
      explicit EmbeddingPFCandProducer(const edm::ParameterSet&);
      ~EmbeddingPFCandProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::Muon>> muonsSrc_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonsEmbedSrc_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonsSrcEmbedding_;
  edm::EDGetTokenT<pat::ElectronCollection> electronsSrc_;
  edm::EDGetTokenT<pat::ElectronCollection> electronsSrcEmbedding_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedCandSrc_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedCandSrcEmbedding_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTrackSrc_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTrackSrcEmbedding_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> offlineSlimmedPrimaryVerticesSrc_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> offlineSlimmedPrimaryVerticesSrcEmbedding_;
  
  edm::ParameterSet* cfg_;

};

EmbeddingPFCandProducer::EmbeddingPFCandProducer(const edm::ParameterSet& iConfig):
  muonsSrc_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonsSrc"))),
  muonsEmbedSrc_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonsEmbedSrc"))),
  muonsSrcEmbedding_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonsSrcEmbedding"))),
  electronsSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronsSrc"))),
  electronsSrcEmbedding_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronsSrcEmbedding"))),
  packedCandSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedCandSrc"))),
  packedCandSrcEmbedding_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedCandSrcEmbedding"))),
  lostTrackSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostTrackSrc"))),
  lostTrackSrcEmbedding_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostTrackSrcEmbedding"))),
  offlineSlimmedPrimaryVerticesSrc_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offlineSlimmedPrimaryVerticesSrc"))),
  offlineSlimmedPrimaryVerticesSrcEmbedding_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offlineSlimmedPrimaryVerticesSrcEmbedding"))){
  //register your products
  cfg_ = const_cast<edm::ParameterSet*>(&iConfig);

  produces<std::vector<pat::Muon>>("slimmedMuonsEmbedded");
  produces<pat::ElectronCollection>("slimmedElectronsEmbedded");
  produces<pat::PackedCandidateCollection >("packedPFCandidatesEmbedded");
  produces<pat::PackedCandidateCollection >("lostTracksEmbedded");
  produces<std::vector<reco::Vertex> >("offlineSlimmedPrimaryVerticesEmbedded");
}


EmbeddingPFCandProducer::~EmbeddingPFCandProducer() {}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
EmbeddingPFCandProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   edm::Handle< edm::View<pat::Muon> > muonHandle;
   iEvent.getByToken(muonsSrc_, muonHandle);
   edm::View<pat::Muon> coll_muons = *muonHandle;

   edm::Handle< std::vector<pat::Muon> > muons;
   iEvent.getByToken(muonsEmbedSrc_, muons);

   edm::Handle< std::vector<pat::Muon> > muonsEmbedding;
   iEvent.getByToken(muonsSrcEmbedding_, muonsEmbedding);

   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronsSrc_, electrons);

   edm::Handle<pat::ElectronCollection> electronsEmbedding;
   iEvent.getByToken(electronsSrcEmbedding_, electronsEmbedding);

   edm::Handle<pat::PackedCandidateCollection> packedCands;
   iEvent.getByToken(packedCandSrc_, packedCands);

   edm::Handle<pat::PackedCandidateCollection> packedCandsEmbedding;
   iEvent.getByToken(packedCandSrcEmbedding_, packedCandsEmbedding);

   edm::Handle<pat::PackedCandidateCollection> lostTracks;
   iEvent.getByToken(lostTrackSrc_, lostTracks);

   edm::Handle<pat::PackedCandidateCollection> lostTracksEmbedding;
   iEvent.getByToken(lostTrackSrcEmbedding_, lostTracksEmbedding);

   edm::Handle<std::vector<reco::Vertex>> offlineSlimmedPrimaryVertices;
   iEvent.getByToken(offlineSlimmedPrimaryVerticesSrc_, offlineSlimmedPrimaryVertices);

   edm::Handle<std::vector<reco::Vertex>> offlineSlimmedPrimaryVerticesEmbedding;
   iEvent.getByToken(offlineSlimmedPrimaryVerticesSrcEmbedding_, offlineSlimmedPrimaryVerticesEmbedding);


   std::unique_ptr< std::vector<pat::Muon> > slimmedMuonsEmbedded(new std::vector<pat::Muon>);
   std::unique_ptr< pat::ElectronCollection > slimmedElectronsEmbedded(new pat::ElectronCollection);
   std::unique_ptr<pat::PackedCandidateCollection> packedCandsEmbedded(new pat::PackedCandidateCollection);
   std::unique_ptr<pat::PackedCandidateCollection> lostTracksEmbedded(new pat::PackedCandidateCollection);
   std::unique_ptr<std::vector<reco::Vertex>> offlineSlimmedPrimaryVerticesEmbedded(new std::vector<reco::Vertex>);

   
   std::vector<reco::CandidatePtr> mSourceCandPtrs;
   for (edm::View<pat::Muon>::const_iterator muon = coll_muons.begin(); muon!= coll_muons.end();  ++muon) {
     for( unsigned int i=0; i < muon->numberOfSourceCandidatePtrs(); ++i){
       mSourceCandPtrs.push_back(muon->sourceCandidatePtr(i));
     }
   }

   //std::cout << mSourceCandPtrs.size() << "\n";

   for( size_t i = 0; i < packedCands->size(); ++i) {
     reco::CandidatePtr ptr2PF(packedCands,i);
     if (std::find(mSourceCandPtrs.begin(),mSourceCandPtrs.end(),ptr2PF) == mSourceCandPtrs.end()) {
       packedCandsEmbedded->push_back((*packedCands)[i]);
       packedCandsEmbedded->back().setCovarianceVersion(1);
     }
   }

//   std::vector<reco::PFCandidateRef> mSourceRefs;
//   for (edm::View<pat::Muon>::const_iterator muon = coll_muons.begin(); muon!= coll_muons.end();  ++muon) {
//     mSourceRefs.push_back(muon->pfCandidateRef());
//   }
  
   for( size_t i = 0; i < muons->size(); ++i) {
     //reco::CandidatePtr mptr(muons,i);
     if (std::find(mSourceCandPtrs.begin(),mSourceCandPtrs.end(),(*muons)[i].sourceCandidatePtr(0)) == mSourceCandPtrs.end()) slimmedMuonsEmbedded -> push_back((*muons)[i]);
   }
   std::cout << "Debug: " << muons->size() << " | " << slimmedMuonsEmbedded->size() << " | " <<  muonsEmbedding->size() << "\n";

   for( size_t i = 0; i < packedCandsEmbedding->size(); ++i) {
     reco::CandidatePtr ptr2PF(packedCandsEmbedding,i);
     packedCandsEmbedded->push_back((*packedCandsEmbedding)[i]);
     packedCandsEmbedded->back().setCovarianceVersion(1);
   }

   for( size_t i = 0; i < muonsEmbedding->size(); ++i) {
     slimmedMuonsEmbedded -> push_back((*muonsEmbedding)[i]);
   }
   
   for( size_t i = 0; i < lostTracks->size(); ++i) {
     lostTracksEmbedded->push_back((*lostTracks)[i]);
     //std::cout << (*lostTracks)[i].covarianceVersion() << "\n";
     lostTracksEmbedded->back().setCovarianceVersion(1);
   }
   
   for( size_t i = 0; i < lostTracksEmbedding->size(); ++i) {
     lostTracksEmbedded->push_back((*lostTracksEmbedding)[i]);
     lostTracksEmbedded->back().setCovarianceVersion(1);
   }
   
   for ( size_t i = 0; i < offlineSlimmedPrimaryVerticesEmbedding->size(); ++i) {
     offlineSlimmedPrimaryVerticesEmbedded -> push_back((*offlineSlimmedPrimaryVerticesEmbedding)[i]);
   }

   for ( size_t i = 0; i < offlineSlimmedPrimaryVertices->size(); ++i) {
     if (i == 0) continue;
     offlineSlimmedPrimaryVerticesEmbedded -> push_back((*offlineSlimmedPrimaryVertices)[i]);
   }

   for ( size_t i = 0; i < electronsEmbedding->size(); ++i) {
     slimmedElectronsEmbedded -> push_back((*electronsEmbedding)[i]);
   }

   for ( size_t i = 0; i < electrons->size(); ++i) {
     slimmedElectronsEmbedded -> push_back((*electrons)[i]);
   }
   
   iEvent.put(std::move(slimmedMuonsEmbedded),"slimmedMuonsEmbedded");
   iEvent.put(std::move(slimmedElectronsEmbedded),"slimmedElectronsEmbedded");
   iEvent.put(std::move(packedCandsEmbedded),"packedPFCandidatesEmbedded");
   iEvent.put(std::move(lostTracksEmbedded),"lostTracksEmbedded");
   iEvent.put(std::move(offlineSlimmedPrimaryVerticesEmbedded),"offlineSlimmedPrimaryVerticesEmbedded");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
EmbeddingPFCandProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
EmbeddingPFCandProducer::endStream() {
}

 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EmbeddingPFCandProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(EmbeddingPFCandProducer);
