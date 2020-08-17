#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "PhysicsTools/SelectorUtils/interface/CutApplicatorBase.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "TLorentzVector.h"


#include "TMath.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "PhysicsTools/SelectorUtils/interface/CutApplicatorWithEventContentBase.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/RefToBaseVector.h"

using namespace edm;
using namespace std;
//
// class declaration
//

class PATElectronBaseLineSelection : public edm::stream::EDFilter<> {
public:
  explicit PATElectronBaseLineSelection(const edm::ParameterSet&);
  ~PATElectronBaseLineSelection();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  float dEtaInSeed(pat::ElectronCollection::const_iterator ele);
  float GsfEleEInverseMinusPInverse(pat::ElectronCollection::const_iterator ele);
  int GsfEleMissingHitsCut(pat::ElectronCollection::const_iterator ele);
  bool GsfEleConversionVetoCut(pat::ElectronCollection::const_iterator ele,edm::Event& iEvent);

private:
  bool isDebug_ ;
  
  
private:
  virtual void beginStream(edm::StreamID) override;
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;

  // ----------member data ---------------------------
  float dEtaInSeedCut;
  float GsfEleEInverseMinusPInverseCut;
  
  edm::EDGetTokenT<pat::ElectronCollection> electronSrc_;
  edm::EDGetTokenT<double> rho_;
  edm::EDGetTokenT<reco::ConversionCollection> convs_;
  edm::EDGetTokenT<reco::BeamSpot> thebs_;
  edm::EDGetTokenT<reco::VertexCollection> vtx_;
};

PATElectronBaseLineSelection::PATElectronBaseLineSelection(const edm::ParameterSet& iConfig):
  electronSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  rho_(consumes<double>(iConfig.getParameter<edm::InputTag>("Rho"))),
  convs_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conv"))),
  thebs_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BM"))),
  vtx_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex")))
  //isDebug_ (iConfig.getUntrackedParameter<bool>("isDebug", false))
{

  isDebug_ = iConfig.getUntrackedParameter<bool>("isDebug", false);
  produces<pat::ElectronCollection>( "myElectrons" );
}


PATElectronBaseLineSelection::~PATElectronBaseLineSelection(){}

float PATElectronBaseLineSelection::dEtaInSeed (pat::ElectronCollection::const_iterator ele){
  return ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull() ?
    ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
}
float PATElectronBaseLineSelection::GsfEleEInverseMinusPInverse (pat::ElectronCollection::const_iterator ele)
{
  const float ecal_energy_inverse = 1.0/ele->ecalEnergy();
  const float eSCoverP = ele->eSuperClusterOverP();
  return std::abs(1.0 - eSCoverP)*ecal_energy_inverse;
}
int PATElectronBaseLineSelection::GsfEleMissingHitsCut(pat::ElectronCollection::const_iterator ele)
{

  constexpr reco::HitPattern::HitCategory missingHitType =
    reco::HitPattern::MISSING_INNER_HITS;
    const unsigned mHits =
      ele->gsfTrack()->hitPattern().numberOfAllHits(missingHitType);
    return mHits;
}


bool PATElectronBaseLineSelection::GsfEleConversionVetoCut(pat::ElectronCollection::const_iterator ele ,edm::Event& iEvent)
{

  edm::Handle<reco::ConversionCollection> convs;
  iEvent.getByToken(convs_,convs);
  edm::Handle<reco::BeamSpot> thebs;
  iEvent.getByToken(thebs_,thebs);
  if(thebs.isValid() && convs.isValid() ) {
    return !ConversionTools::hasMatchedConversion(*ele,convs,
						  thebs->position());
  } else {
    edm::LogWarning("GsfEleConversionVetoCut")
      << "Couldn't find a necessary collection, returning true!";
    return true;
  }
}

// ------------ method called on each new Event  ------------
bool PATElectronBaseLineSelection::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronSrc_,electrons);
  unique_ptr<pat::ElectronCollection> passedelectrons(new pat::ElectronCollection);
  
  Handle<reco::VertexCollection> Vertex;
  iEvent.getByToken(vtx_,Vertex);
  const reco::Vertex& pv=*Vertex->begin();
  math::XYZPoint p = pv.position();
  
  edm::Handle<double>_rhoHandle;
  iEvent.getByToken(rho_,_rhoHandle);

  for(pat::ElectronCollection::const_iterator iele = electrons->begin() ; iele !=electrons->end(); ++iele){
    dEtaInSeedCut =abs(dEtaInSeed(iele));
    GsfEleEInverseMinusPInverseCut = GsfEleEInverseMinusPInverse(iele);
    
    //-------------------The new H/E variables------------//
    double HoE=iele->hadronicOverEm();
    double E_c = iele->superCluster()->energy();
    double rho = _rhoHandle.isValid() ? (*_rhoHandle) : 0;
    double dxy = abs(iele->gsfTrack()->dxy(p));
    double dz = abs(iele->gsfTrack()->dz(p));
    //std::cout << rho << " | " << E_c << " | " << 0.05 + 1.16/E_c + 0.0324*rho/E_c << " | " << 0.0441 + 2.54/E_c + 0.183*rho/E_c << "\n";
    //if ((iele -> pt() > 7) && (abs(iele -> eta()) < 2.5) ) {
    //	 passedelectrons->push_back(*iele);
	
    if ((iele -> pt() > 7) && (abs(iele -> eta()) < 2.5)) {
      if( iele->isEB()) {
	if( (iele->full5x5_sigmaIetaIeta()<0.0112)  &&
	    //	    (HoE < (0.05 + 1.16/E_c + 0.0324*rho/E_c)) &&
	    (abs(iele->deltaPhiSuperClusterTrackAtVtx()) <0.0884) &&
	    (GsfEleEInverseMinusPInverseCut < 0.193) &&
	    (dEtaInSeedCut < 0.00377) &&
	    (GsfEleMissingHitsCut(iele) <= 1 ) &&
	    (GsfEleConversionVetoCut(iele,iEvent)) && 
	    (dxy < 0.05) &&
	    (dz < 0.1) ) {
	  passedelectrons->push_back(*iele);
	  //if (iele -> pt() > 10 && iele -> pt() < 11 && isDebug_) {
	  if (isDebug_) {
	      std::cout << "Debug: "
			<< iele->pt() << " | "
			<< dxy << " | "
			<< dz << " | "
			<< GsfEleConversionVetoCut(iele,iEvent) << " | "
			<< dEtaInSeedCut << "|"
			<< GsfEleEInverseMinusPInverseCut << " | "
			<< abs(iele->deltaPhiSuperClusterTrackAtVtx()) << " | "
			<< iele->full5x5_sigmaIetaIeta() << " | "
	      		<< HoE << " | "
	       		<< (0.05 + 1.16/E_c + 0.0324*rho/E_c) << " | "
			<< GsfEleConversionVetoCut(iele,iEvent) << "\n";
	  }
	}
      }
      if(iele->isEE()) {
	if( (iele->full5x5_sigmaIetaIeta() < 0.0425) &&
	    //	    (HoE < (0.0441 + 2.54/E_c + 0.183*rho/E_c)) &&
	    (abs(iele->deltaPhiSuperClusterTrackAtVtx()) < 0.169) &&
	    (GsfEleEInverseMinusPInverseCut < 0.111) &&
	    (dEtaInSeedCut <0.00674) &&
	    (GsfEleMissingHitsCut(iele) <= 1 ) &&
	    (GsfEleConversionVetoCut(iele,iEvent)) &&
	    (dxy < 0.1) &&
	    (dz < 0.2) ) {
	  passedelectrons->push_back(*iele);
	  //if (iele -> pt() > 10 && iele -> pt() < 11 && and isDebug_) {
	  if (isDebug_) {
	      std::cout << "Debug: "
			<< iele->pt() << " | "
			<< dxy << " | "
			<< dz << " | "
			<< GsfEleConversionVetoCut(iele,iEvent) << " | "
			<< dEtaInSeedCut << "|"
			<< GsfEleEInverseMinusPInverseCut << " | "
			<< abs(iele->deltaPhiSuperClusterTrackAtVtx()) << " | "
			<< iele->full5x5_sigmaIetaIeta() << " | "
	           	<< HoE << " | "
       			<< (0.05 + 1.16/E_c + 0.0324*rho/E_c) << " | "
			<< GsfEleConversionVetoCut(iele,iEvent) << "\n";
	  }
	}
      }
    }
  }

  iEvent.put(move(passedelectrons), "myElectrons");
  return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
PATElectronBaseLineSelection::beginStream(edm::StreamID)
{}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PATElectronBaseLineSelection::endStream() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PATElectronBaseLineSelection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PATElectronBaseLineSelection);
