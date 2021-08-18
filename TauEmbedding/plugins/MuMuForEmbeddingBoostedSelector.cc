// -*- C++ -*-
//
// Package:    TauAnalysis/EmbeddingProducer
// Class:      MuMuForEmbeddingBoostedSelector
// 

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Math/interface/deltaPhi.h"

//
// class declaration
//

class MuMuForEmbeddingBoostedSelector : public edm::stream::EDProducer<> {
   public:
      explicit MuMuForEmbeddingBoostedSelector(const edm::ParameterSet&);
      ~MuMuForEmbeddingBoostedSelector() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginStream(edm::StreamID) override;
      void produce(edm::Event&, const edm::EventSetup&) override;
      void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<reco::CompositeCandidate>> BoostedmumuCandidates_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MuMuForEmbeddingBoostedSelector::MuMuForEmbeddingBoostedSelector(const edm::ParameterSet& iConfig) :
   BoostedmumuCandidates_(consumes< edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("ZmumuCandidatesCollection")))
{
  
   produces<edm::RefVector<pat::MuonCollection>>();
   
   //now do what ever other initialization is needed
}


MuMuForEmbeddingBoostedSelector::~MuMuForEmbeddingBoostedSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuMuForEmbeddingBoostedSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle< edm::View<reco::CompositeCandidate> > BoostedmumuCandidatesHandle;
   iEvent.getByToken(BoostedmumuCandidates_, BoostedmumuCandidatesHandle);
   edm::View<reco::CompositeCandidate> BoostedmumuCandidates = *BoostedmumuCandidatesHandle;
   
   const reco::CompositeCandidate* chosenCand = nullptr;
   double dR = 9999;
   for (edm::View<reco::CompositeCandidate>::const_iterator iCand = BoostedmumuCandidates.begin(); iCand != BoostedmumuCandidates.end(); ++iCand) {

     const float dphi = reco::deltaPhi(iCand->daughter(0)->phi(), iCand->daughter(1)->phi());
     const float deta = iCand->daughter(0)->eta() - iCand->daughter(1)->eta();
       
     const float dr = std::sqrt(deta*deta + dphi*dphi);
     
     if (dr < dR) {
       dR = dr;
       chosenCand = &(*iCand);
     }
     
   }
   std::unique_ptr<edm::RefVector<pat::MuonCollection>> prod(new edm::RefVector<pat::MuonCollection>());
   prod->reserve(2);
   prod->push_back(chosenCand->daughter(0)->masterClone().castTo<pat::MuonRef>());
   prod->push_back(chosenCand->daughter(1)->masterClone().castTo<pat::MuonRef>());
   iEvent.put(std::move(prod));
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
MuMuForEmbeddingBoostedSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MuMuForEmbeddingBoostedSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MuMuForEmbeddingSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MuMuForEmbeddingSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MuMuForEmbeddingSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MuMuForEmbeddingSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuMuForEmbeddingBoostedSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuForEmbeddingBoostedSelector);
