// -*- C++ -*-
//
// Package:    ElectronWork/ElectronNtupler
// Class:      PhotonNtuplerMiniAOD
// 
/**\class PhotonNtuplerMiniAOD PhotonNtuplerMiniAOD.cc ElectronWork/PhotonNtuplerMiniAOD/plugins/PhotonNtuplerMiniAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"



//
// class declaration
//

class PhotonNtuplerMiniAOD : public edm::EDAnalyzer {
   public:
      explicit PhotonNtuplerMiniAOD(const edm::ParameterSet&);
      ~PhotonNtuplerMiniAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      enum PhotonMatchType {UNMATCHED = 0, 
			    MATCHED_FROM_GUDSCB,
			    MATCHED_FROM_PI0,
			    MATCHED_FROM_OTHER_SOURCES};
  
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      int matchToTruth(const pat::Photon &pho, 
            const edm::Handle<edm::View<reco::GenParticle>>  &prunedGenParticles);

       void findFirstNonPhotonMother(const reco::Candidate *particle,
				     int &ancestorPID, int &ancestorStatus);

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<edm::View<pat::Photon> > photonCollectionToken_;
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMapToken_; 
      edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolationToken_; 
      edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolationToken_; 
      edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolationToken_; 

  TTree *photonTree_;

  Int_t nPV_;        // number of reconsrtucted primary vertices
  Float_t rho_;      // the rho variable


  // all photon variables
  Int_t nPhotons_;

  std::vector<Float_t> pt_;
  std::vector<Float_t> eta_;
  std::vector<Float_t> phi_;

  std::vector<Float_t> full5x5_sigmaIetaIeta_;
  std::vector<Float_t> hOverE_;
  std::vector<Float_t> hasPixelSeed_;
  std::vector<Float_t> r9_;

  std::vector<Float_t> isoChargedHadrons_;
  std::vector<Float_t> isoNeutralHadrons_;
  std::vector<Float_t> isoPhotons_;

  std::vector<Float_t> isoChargedHadronsWithEA_;
  std::vector<Float_t> isoNeutralHadronsWithEA_;
  std::vector<Float_t> isoPhotonsWithEA_;

  std::vector<Int_t>   isTrue_;

};

//
// constants, enums and typedefs
//
// Effective areas for photons from Savvas's slides
// for phys14 PU20bx25, described here:
// https://indico.cern.ch/event/367861/contribution/3/material/slides/0.pdf
namespace EffectiveAreas {
  const int nEtaBins = 7;
  const float etaBinLimits[nEtaBins+1] = {
    0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};

  const float areaPhotons[nEtaBins] = {
    0.0894361 ,  0.0750406 ,  0.0422692 ,  0.0561305 ,  0.0881777 ,  0.114391  , 0.168394    
  };
  const float areaNeutralHadrons[nEtaBins] = {
    0.00491696 ,  0.0108215  ,  0.00185819 ,  0.0036691  ,  0.00619475 ,  0.0130409  ,  0.169902   
  };
  const float areaChargedHadrons[nEtaBins] = {
    0.00885603 ,  0.00615045 ,  0.00861475 ,  0.00412423 ,  0.0113191  ,  0.00846383 ,  0.0038625  
  };

  // const float areaPhotons[nEtaBins] = {
  //   0.0894, 0.0750, 0.0423, 0.0561, 0.0882, 0.1144, 0.1684
  // };
  // const float areaNeutralHadrons[nEtaBins] = {
  //   0.049, 0.0108, 0.0019, 0.0037, 0.0062, 0.0130, 0.1699
  // };
  // const float areaChargedHadrons[nEtaBins] = {
  //   0.0089, 0.0062, 0.0086, 0.0041, 0.0113, 0.0085, 0.0039
  // };
}

//
// static data member definitions
//

//
// constructors and destructor
//
PhotonNtuplerMiniAOD::PhotonNtuplerMiniAOD(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  photonCollectionToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"))),
  rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
  full5x5SigmaIEtaIEtaMapToken_(consumes <edm::ValueMap<float> >
				(iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap"))),
  phoChargedIsolationToken_(consumes <edm::ValueMap<float> >
			    (iConfig.getParameter<edm::InputTag>("phoChargedIsolation"))),
  phoNeutralHadronIsolationToken_(consumes <edm::ValueMap<float> >
				  (iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation"))),
  phoPhotonIsolationToken_(consumes <edm::ValueMap<float> >
			   (iConfig.getParameter<edm::InputTag>("phoPhotonIsolation")))
{

  edm::Service<TFileService> fs;
  photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");
  
  photonTree_->Branch("nPV"        ,  &nPV_     , "nPV/I");
  photonTree_->Branch("rho"        ,  &rho_ , "rho/F");

  photonTree_->Branch("nPho",  &nPhotons_ , "nPho/I");
  photonTree_->Branch("pt"  ,  &pt_    );
  photonTree_->Branch("eta" ,  &eta_ );
  photonTree_->Branch("phi" ,  &phi_ );

  photonTree_->Branch("hOverE",  &hOverE_);
  photonTree_->Branch("hasPixelSeed"           ,  &hasPixelSeed_);
  photonTree_->Branch("full5x5_sigmaIetaIeta"  , &full5x5_sigmaIetaIeta_);
  photonTree_->Branch("r9",  &r9_);
  photonTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_);
  photonTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_);
  photonTree_->Branch("isoPhotons"             , &isoPhotons_);

  photonTree_->Branch("isoChargedHadronsWithEA"      , &isoChargedHadronsWithEA_);
  photonTree_->Branch("isoNeutralHadronsWithEA"      , &isoNeutralHadronsWithEA_);
  photonTree_->Branch("isoPhotonsWithEA"             , &isoPhotonsWithEA_);

  photonTree_->Branch("isTrue"                 , &isTrue_);
 
}


PhotonNtuplerMiniAOD::~PhotonNtuplerMiniAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonNtuplerMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  // using namespace reco;
  
  // Get photon collection
  edm::Handle<edm::View<pat::Photon> > collection;
  iEvent.getByToken(photonCollectionToken_, collection);

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  //const reco::Vertex &pv = vertices->front();
  nPV_    = vertices->size();

  reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (reco::VertexCollection::const_iterator vtx = vertices->begin(); 
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    if (  !vtx->isFake()
	  //!(vtx->chi2()==0 && vtx->ndof()==0)  // This line is for AOD
        &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
        && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  if ( firstGoodVertex==vertices->end() )
    return; // skip event if there are no good PVs

  // Pruned particles are the one containing "important" stuff
  Handle<edm::View<reco::GenParticle> > prunedGenParticles;
  iEvent.getByToken(prunedGenToken_,prunedGenParticles);
  
  // Get rho
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  // Get the full5x5 sieie map
  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
  iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);

  // Get the isolation maps
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  iEvent.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  iEvent.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  iEvent.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);

  // Clear vectors
  nPhotons_ = 0;
  pt_.clear();
  eta_.clear();
  phi_.clear();
  full5x5_sigmaIetaIeta_.clear();
  hOverE_.clear();
  hasPixelSeed_.clear();
  r9_.clear();
  isoChargedHadrons_.clear();
  isoNeutralHadrons_.clear();
  isoPhotons_.clear();
  isoChargedHadronsWithEA_.clear();
  isoNeutralHadronsWithEA_.clear();
  isoPhotonsWithEA_.clear();
  isTrue_.clear();

  // Loop over photons
  
  // const auto& pho_refs = collection->refVector();
  
  // for( const auto& pho : pho_refs ) {

  for( View<pat::Photon>::const_iterator pho = collection->begin();
       pho != collection->end(); pho++){
    
    // Kinematics (nobody uses photons below 15 GeV, and they are not stored in miniAOD anyways)
    if( pho->pt() < 15 ) 
      continue;
    
    nPhotons_++;
    pt_  .push_back( pho->pt() );
    eta_ .push_back( pho->superCluster()->eta() );
    phi_ .push_back( pho->superCluster()->phi() );

    //const edm::Ptr<reco::Photon> phoPtr( pho );
    const edm::Ptr<pat::Photon> phoPtr( collection, pho - collection->begin() );

    full5x5_sigmaIetaIeta_ .push_back( (*full5x5SigmaIEtaIEtaMap)[ phoPtr ] );
    hOverE_                .push_back( pho->hadTowOverEm() );
    hasPixelSeed_          .push_back( pho->hasPixelSeed() );
    r9_                    .push_back( pho->userFloat("r9_NoZS"));

    isoChargedHadrons_ .push_back( (*phoChargedIsolationMap)[phoPtr] );
    isoNeutralHadrons_ .push_back( (*phoNeutralHadronIsolationMap)[phoPtr] );
    isoPhotons_        .push_back( (*phoPhotonIsolationMap)[phoPtr] );

    // Compute isolation with effective area correction for PU
    // Find eta bin first. If eta>2.5, the last eta bin is used.
    int etaBin = 0; 
    while ( etaBin < EffectiveAreas::nEtaBins-1 
	    && abs( pho->superCluster()->eta() ) > EffectiveAreas::etaBinLimits[etaBin+1] )
      { ++etaBin; };
    isoPhotonsWithEA_        .push_back( std::max( (float)0.0, (*phoPhotonIsolationMap)       [phoPtr] 
						  - rho_ * EffectiveAreas::areaPhotons[etaBin] ) );
    isoNeutralHadronsWithEA_ .push_back( std::max( (float)0.0, (*phoNeutralHadronIsolationMap)[phoPtr] 
						  - rho_ * EffectiveAreas::areaNeutralHadrons[etaBin] ) );
    isoChargedHadronsWithEA_ .push_back( std::max( (float)0.0, (*phoChargedIsolationMap)      [phoPtr] 
						  - rho_ * EffectiveAreas::areaChargedHadrons[etaBin] ) );

    isTrue_.push_back( matchToTruth(*pho, prunedGenParticles) );

   }
   
  // Save the info
  photonTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonNtuplerMiniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonNtuplerMiniAOD::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PhotonNtuplerMiniAOD::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PhotonNtuplerMiniAOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
PhotonNtuplerMiniAOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PhotonNtuplerMiniAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonNtuplerMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int PhotonNtuplerMiniAOD::matchToTruth(const pat::Photon &pho, 
				       const edm::Handle<edm::View<reco::GenParticle>>  
				       &prunedGenParticles)
{
  // 
  // Explicit loop and geometric matching method 
  //

  // Find the closest status 1 gen photon to the reco photon
  double dR = 999;
  const reco::Candidate *closestPhoton = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not photon or not status 1
    if( abs(particle->pdgId()) != 22 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( pho.p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestPhoton = particle;
    }
  }
  // See if the closest photon (if it exists) is close enough.
  // If not, no match found.
  if( !(closestPhoton != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // Find ID of the parent of the found generator level photon match
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);

  // Allowed parens: quarks pdgId 1-5, or a gluon 21
  // Note: in miniAOD, photons are sometimes traced to the incoming protons 
  //    (pid=2212, status=4)
  std::vector<int> allowedParents { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -21, 21, 2212 };
  if( !(std::find(allowedParents.begin(), 
		 allowedParents.end(), ancestorPID)
	!= allowedParents.end()) ){
    // So it is not from g, u, d, s, c, b. Check if it is from pi0 or not. 
    if( abs(ancestorPID) == 111 )
      return MATCHED_FROM_PI0;
    else
      return MATCHED_FROM_OTHER_SOURCES;
    
  }
  // if we matched an ancestor to a proton, but it is not an incoming
  // (status 4) proton, then we still consider this photon to be "from other sources"
  if( ancestorPID==2212 && ancestorStatus!=4 ){
    return MATCHED_FROM_OTHER_SOURCES;
  }

  return MATCHED_FROM_GUDSCB;
   
}

void PhotonNtuplerMiniAOD::findFirstNonPhotonMother(const reco::Candidate *particle,
						      int &ancestorPID, int &ancestorStatus){
  
  if( particle == 0 ){
    printf("PhotonNtuplerMiniAOD: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-photon parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 22 ){
    findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }
  
  return;
}


//define this as a plug-in
DEFINE_FWK_MODULE(PhotonNtuplerMiniAOD);
