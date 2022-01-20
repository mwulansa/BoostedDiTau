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
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedCandSrc_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedCandSrcEmbedding_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTrackSrc_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTrackSrcEmbedding_;
  
  edm::ParameterSet* cfg_;

};

EmbeddingPFCandProducer::EmbeddingPFCandProducer(const edm::ParameterSet& iConfig):
  muonsSrc_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonsSrc"))),
  packedCandSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedCandSrc"))),
  packedCandSrcEmbedding_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedCandSrcEmbedding"))),
  lostTrackSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostTrackSrc"))),
  lostTrackSrcEmbedding_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostTrackSrcEmbedding"))){
  //register your products
  cfg_ = const_cast<edm::ParameterSet*>(&iConfig);
  
  produces<pat::PackedCandidateCollection >("packedPFCandidatesEmbedded");
  produces<pat::PackedCandidateCollection >("lostTracksEmbedded");
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

   edm::Handle<pat::PackedCandidateCollection> packedCands;
   iEvent.getByToken(packedCandSrc_, packedCands);

   edm::Handle<pat::PackedCandidateCollection> packedCandsEmbedding;
   iEvent.getByToken(packedCandSrcEmbedding_, packedCandsEmbedding);

   edm::Handle<pat::PackedCandidateCollection> lostTracks;
   iEvent.getByToken(lostTrackSrc_, lostTracks);

   edm::Handle<pat::PackedCandidateCollection> lostTracksEmbedding;
   iEvent.getByToken(lostTrackSrcEmbedding_, lostTracksEmbedding);

   std::cout << lostTracks->size() << " | " << lostTracksEmbedding->size() << "\n";
  
   std::unique_ptr<pat::PackedCandidateCollection> packedCandsEmbedded(new pat::PackedCandidateCollection);

   std::unique_ptr<pat::PackedCandidateCollection> lostTracksEmbedded(new pat::PackedCandidateCollection);

   std::cout << coll_muons.size() << " | " << packedCands->size() << " | " << packedCandsEmbedding->size() << " | " << packedCandsEmbedded->size() << "\n";
   
   std::vector<reco::CandidatePtr> mSourceCandPtrs;
   for (edm::View<pat::Muon>::const_iterator muon = coll_muons.begin(); muon!= coll_muons.end();  ++muon) {
     for( unsigned int i=0; i < muon->numberOfSourceCandidatePtrs(); ++i){
       mSourceCandPtrs.push_back(muon->sourceCandidatePtr(i));
     }
   }

   //std::cout << (*packedCands)[0].covarianceVersion() << "|" << (*packedCands)[0].bestTrack()->covarianceVersion()<< "\n";
   for( size_t i = 0; i < packedCands->size(); ++i) {
     reco::CandidatePtr ptr2PF(packedCands,i);
     if (std::find(mSourceCandPtrs.begin(),mSourceCandPtrs.end(),ptr2PF) == mSourceCandPtrs.end()) {
       packedCandsEmbedded->push_back((*packedCands)[i]);
     }
   }

   //std::cout << (*packedCandsEmbedding)[0].covarianceVersion() << "\n";
   for( size_t i = 0; i < packedCandsEmbedding->size(); ++i) {
     reco::CandidatePtr ptr2PF(packedCandsEmbedding,i);
     packedCandsEmbedded->push_back((*packedCandsEmbedding)[i]);
     packedCandsEmbedded->back().setCovarianceVersion(0);
   }
   //std::cout << (*packedCandsEmbedding)[0].covarianceVersion() << "\n";

   std::cout << coll_muons.size() << " | " << packedCands->size() << " | " << packedCandsEmbedding->size() << " | " << packedCandsEmbedded->size() << "\n";

   for( size_t i = 0; i < lostTracks->size(); ++i) {
     lostTracksEmbedded->push_back((*lostTracks)[i]);
   }

   for( size_t i = 0; i < lostTracksEmbedding->size(); ++i) {
     lostTracksEmbedded->push_back((*lostTracksEmbedding)[i]);
     lostTracksEmbedded->back().setCovarianceVersion(0);
   }
   
//
//   //Get the PFCandidates being pointed to by pat::Electrons
//   std::vector<reco::CandidatePtr> mSourceCandPtrs;
//   
//   if (coll_muons.isValid()) {
//     for (pat::MuRefVector::const_iterator iElectron = electrons->begin(); iElectron != electrons->end(); ++iElectron)
//       {
//	 
//	 for( unsigned int i=0; i < (*iElectron)->numberOfSourceCandidatePtrs(); ++i) {
//	   eSourceCandPtrs.push_back((*iElectron)->sourceCandidatePtr(i));
//	 }                      
//       }
//   }
//   for( size_t i = 0; i < packedCands->size(); ++i) {
//     //bool ElectronFlag= false;
//     //if((*packedCands)[i].isElectron())
//     if((*packedCands)[i].pdgId()==11) {
//       reco::CandidatePtr ptr2PF(packedCands,i);
//       //std::cout<< " ====packed Candidate is an electron=== "<<std::endl;
//       if (std::find(eSourceCandPtrs.begin(),eSourceCandPtrs.end(),ptr2PF) != eSourceCandPtrs.end()) {}
//       else {
//	 packedCandsExcludingElectrons->push_back((*packedCands)[i]);
//       }
//     }
//     else {
//       packedCandsExcludingElectrons->push_back((*packedCands)[i]);
//     }    
//   }
   iEvent.put(std::move(packedCandsEmbedded),"packedPFCandidatesEmbedded");
   iEvent.put(std::move(lostTracksEmbedded),"lostTracksEmbedded");
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
