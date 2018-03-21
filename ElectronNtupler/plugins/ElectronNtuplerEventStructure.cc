// -*- C++ -*-
//
// Package:    ElectronWork/ElectronNtupler
// Class:      ElectronNtuplerEventStructure
// 
/**\class ElectronNtuplerEventStructure ElectronNtuplerEventStructure.cc ElectronWork/ElectronNtuplerEventStructure/plugins/ElectronNtuplerEventStructure.cc

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

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"



//
// class declaration
//

class ElectronNtuplerEventStructure : public edm::EDAnalyzer {
   public:
      explicit ElectronNtuplerEventStructure(const edm::ParameterSet&);
      ~ElectronNtuplerEventStructure();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  enum ElectronMatchType {UNMATCHED = 0, 
			  TRUE_PROMPT_ELECTRON, 
			  TRUE_ELECTRON_FROM_TAU,
			  TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      bool checkAncestor(const reco::Candidate *gen, int ancestorPid);
      int matchToTruth(const pat::Electron &el, const edm::Handle<edm::View<reco::GenParticle>>  &prunedGenParticles);
      void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);
      int matchToTruthAlternative(const pat::Electron &el);
  
      void printAllZeroMothers(const reco::Candidate *particle);

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
      edm::EDGetTokenT<double> rhoToken_;

  TTree *electronTree_;

  // Vars for PVs
  Int_t pvNTracks_;

  // Vars for pile-up
  Int_t nPUTrue_;    // true pile-up
  Int_t nPU_;        // generated pile-up
  Int_t nPV_;        // number of reconsrtucted primary vertices
  Float_t rho_;      // the rho variable

  // all electron variables
  Int_t nElectrons_;

  std::vector<Float_t> pt_;
  std::vector<Float_t> etaSC_;
  std::vector<Float_t> phiSC_;
  std::vector<Float_t> dEtaIn_;
  std::vector<Float_t> dPhiIn_;
  std::vector<Float_t> hOverE_;
  std::vector<Float_t> sigmaIetaIeta_;
  std::vector<Float_t> full5x5_sigmaIetaIeta_;
  std::vector<Float_t> isoChargedHadrons_;
  std::vector<Float_t> isoNeutralHadrons_;
  std::vector<Float_t> isoPhotons_;
  std::vector<Float_t> isoChargedFromPU_;
  std::vector<Float_t> relIsoWithDBeta_;
  std::vector<Float_t> ooEmooP_;
  std::vector<Float_t> d0_;
  std::vector<Float_t> dz_;
  std::vector<Int_t>   expectedMissingInnerHits_;
  std::vector<Int_t>   passConversionVeto_;     
  std::vector<Int_t>   isTrueElectron_;
  std::vector<Int_t>   isTrueElectronAlternative_; 
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
ElectronNtuplerEventStructure::ElectronNtuplerEventStructure(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  pileupToken_(consumes<edm::View<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileup"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
{

  edm::Service<TFileService> fs;
  electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");
  
  electronTree_->Branch("pvNTracks"    ,  &pvNTracks_ , "pvNTracks/I");

  electronTree_->Branch("nPV"        ,  &nPV_     , "nPV/I");
  electronTree_->Branch("nPU"        ,  &nPU_     , "nPU/I");
  electronTree_->Branch("nPUTrue"    ,  &nPUTrue_ , "nPUTrue/I");
  electronTree_->Branch("rho"        ,  &rho_ , "rho/F");

  electronTree_->Branch("nEle"    ,  &nElectrons_ , "nEle/I");
  electronTree_->Branch("pt"    ,  &pt_    );
  electronTree_->Branch("etaSC" ,  &etaSC_ );
  electronTree_->Branch("phiSC" ,  &phiSC_ );
  electronTree_->Branch("dEtaIn",  &dEtaIn_);
  electronTree_->Branch("dPhiIn",  &dPhiIn_);
  electronTree_->Branch("hOverE",  &hOverE_);
  electronTree_->Branch("sigmaIetaIeta",         &sigmaIetaIeta_);
  electronTree_->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_);
  electronTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_);
  electronTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_);
  electronTree_->Branch("isoPhotons"             , &isoPhotons_);
  electronTree_->Branch("isoChargedFromPU"       , &isoChargedFromPU_);
  electronTree_->Branch("relIsoWithDBeta"        , &relIsoWithDBeta_);
  electronTree_->Branch("ooEmooP", &ooEmooP_);
  electronTree_->Branch("d0"     , &d0_);
  electronTree_->Branch("dz"     , &dz_);
  electronTree_->Branch("expectedMissingInnerHits", &expectedMissingInnerHits_);
  electronTree_->Branch("passConversionVeto", &passConversionVeto_);
  electronTree_->Branch("isTrueElectron"    , &isTrueElectron_);
  electronTree_->Branch("isTrueElectronAlternative"    , &isTrueElectronAlternative_);
 
}


ElectronNtuplerEventStructure::~ElectronNtuplerEventStructure()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronNtuplerEventStructure::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  // Pruned particles are the one containing "important" stuff
  Handle<edm::View<reco::GenParticle> > prunedGenParticles;
  iEvent.getByToken(prunedGenToken_,prunedGenParticles);
  
  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  Handle<edm::View<pat::PackedGenParticle> > packedGenParticles;
  iEvent.getByToken(packedGenToken_,packedGenParticles);

  // Get Pileup info
  Handle<edm::View<PileupSummaryInfo> > pileupHandle;
  iEvent.getByToken(pileupToken_, pileupHandle);
  for( auto & puInfoElement : *pileupHandle){
    if( puInfoElement.getBunchCrossing() == 0 ){
      nPU_    = puInfoElement.getPU_NumInteractions();
      nPUTrue_= puInfoElement.getTrueNumInteractions();
    }
  }
  // vector<PileupSummaryInfo>::const_iterator PVI;
  // for (PVI = pileupHandle->begin(); PVI != pileupHandle->end(); ++PVI) {
  //   if (PVI->getBunchCrossing() == 0) {
  //     nPU_     = PVI->getPU_NumInteractions();
  //     nPUTrue_ = PVI->getTrueNumInteractions();
  //   }
  // }

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  //const reco::Vertex &pv = vertices->front();
  nPV_    = vertices->size();

  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin(); 
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    if (  /*!vtx->isFake() &&*/ 
	!(vtx->chi2()==0 && vtx->ndof()==0) 
	&&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	&& fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  if ( firstGoodVertex==vertices->end() )
    return; // skip event if there are no good PVs

  // Seems always zero. Not stored in miniAOD...?
  pvNTracks_ = firstGoodVertex->nTracks();

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;


   // Get electron collection
   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);

   // Loop over electrons
   nElectrons_ = 0;
   pt_.clear();
   etaSC_.clear();
   phiSC_.clear();
   dEtaIn_.clear();
   dPhiIn_.clear();
   hOverE_.clear();
   sigmaIetaIeta_.clear();
   full5x5_sigmaIetaIeta_.clear();
   isoChargedHadrons_.clear();
   isoNeutralHadrons_.clear();
   isoPhotons_.clear();
   isoChargedFromPU_.clear();
   relIsoWithDBeta_.clear();
   ooEmooP_.clear();
   d0_.clear();
   dz_.clear();
   expectedMissingInnerHits_.clear();
   passConversionVeto_.clear();     
   isTrueElectron_.clear();
   isTrueElectronAlternative_.clear(); 

   for (const pat::Electron &el : *electrons) {

     // Kinematics
     if( el.pt() < 10 ) 
       continue;

     nElectrons_++;
     pt_.push_back( el.pt() );
     etaSC_.push_back( el.superCluster()->eta() );
     phiSC_.push_back( el.superCluster()->phi() );
     
     // ID and matching
     dEtaIn_.push_back( el.deltaEtaSuperClusterTrackAtVtx() );
     dPhiIn_.push_back( el.deltaPhiSuperClusterTrackAtVtx() );
     hOverE_.push_back( el.hcalOverEcal() );
     sigmaIetaIeta_.push_back( el.sigmaIetaIeta() );
     full5x5_sigmaIetaIeta_.push_back( el.full5x5_sigmaIetaIeta() );
     // |1/E-1/p| = |1/E - EoverPinner/E| is computed below
     // The if protects against ecalEnergy == inf or zero
     // (always the case for miniAOD for electrons <5 GeV)
     if( el.ecalEnergy() == 0 ){
       printf("Electron energy is zero!\n");
       ooEmooP_.push_back( 1e30 );
     }else if( !std::isfinite(el.ecalEnergy())){
       printf("Electron energy is not finite!\n");
       ooEmooP_.push_back( 1e30 );
     }else{
       ooEmooP_.push_back( fabs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() ) );
     }

     // Isolation
     GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
     // Compute isolation with delta beta correction for PU
     isoChargedHadrons_.push_back( pfIso.sumChargedHadronPt );
     isoNeutralHadrons_.push_back( pfIso.sumNeutralHadronEt );
     isoPhotons_.push_back( pfIso.sumPhotonEt );
     isoChargedFromPU_.push_back( pfIso.sumPUPt );
     float absiso = pfIso.sumChargedHadronPt 
       + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
     relIsoWithDBeta_.push_back( absiso/el.pt() );
     
     // Impact parameter
     d0_.push_back( (-1) * el.gsfTrack()->dxy(firstGoodVertex->position() ) );
     dz_.push_back( el.gsfTrack()->dz( firstGoodVertex->position() ) );
     // printf("     dz= %f\n", dz_);
     
     // Conversion rejection
     // pre-72X method below is commented out
     //  expectedMissingInnerHits_.push_back( el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() );
     // since 72X, the access of missing hits is this:
     expectedMissingInnerHits_.push_back(el.gsfTrack()->hitPattern()
					 .numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) );
     passConversionVeto_.push_back( el.passConversionVeto() );
     
     // Match to generator level truth
 
     // 
     // Explicit loop over gen candidates method
     //
     isTrueElectron_.push_back( matchToTruth( el, prunedGenParticles) );
     isTrueElectronAlternative_.push_back( matchToTruthAlternative( el ) );

     // Use this for debugging if needed, prints decay history.
     // Works with standard matching:
     //   printAllZeroMothers( el.genParticle() );

   }
   
   // Save this electron's info
   electronTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
ElectronNtuplerEventStructure::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronNtuplerEventStructure::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
ElectronNtuplerEventStructure::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ElectronNtuplerEventStructure::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ElectronNtuplerEventStructure::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ElectronNtuplerEventStructure::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronNtuplerEventStructure::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool ElectronNtuplerEventStructure::checkAncestor(const reco::Candidate *gen, int ancestorPid){

  // General sanity check
  if( gen == 0 ){
    printf("ElectronNtuplerEventStructure::checkAncestor: ERROR null particle is passed in, ignore it.\n");
    return false;
  }

  // If this is true, we found our target ancestor
  if( abs( gen->pdgId() ) == ancestorPid )
    return true;

  // Go deeper and check all mothers
  for(size_t i=0;i< gen->numberOfMothers();i++) {
    if ( checkAncestor( gen->mother(i), ancestorPid) )
      return true;
  }
  
  return false;
}

int ElectronNtuplerEventStructure::matchToTruth(const pat::Electron &el, 
				  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el.p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtuplerEventStructure: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }
  
  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void ElectronNtuplerEventStructure::findFirstNonElectronMother(const reco::Candidate *particle,
						 int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("ElectronNtuplerEventStructure: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

int ElectronNtuplerEventStructure::matchToTruthAlternative(const pat::Electron &el){

     //
     // genParticle method
     //
     int result = UNMATCHED;

     const reco::GenParticle * gen = el.genParticle();
     if( gen != 0 ){
       int pid = gen->pdgId();
       int status = gen->status();
       bool isFromZ   = checkAncestor(gen, 23);
       bool isFromW   = checkAncestor(gen, 24);
       bool isFromTau = checkAncestor(gen, 15);
       // Check if it is a true prompt electron
       if( abs( pid ) == 11 // this is electron
	   && (status == 1 || status == 22 || status == 23 ) // NOTE: Pythia8 status here 22/23 (for Pythia6 would be 3)
	   && (isFromZ || isFromW ) && !isFromTau // comes from Z or W+-, but not from tau
	   )
	 {
	   result = TRUE_PROMPT_ELECTRON;
	 } else if ( abs( pid ) == 11 
		     && (status == 1 || status == 22 || status == 23 ) 
		     && (isFromTau ) 
		     ) 
	 {
	   // This is a true electron, but it comes from tau
	   result = TRUE_ELECTRON_FROM_TAU;
	 } else if ( abs( pid ) == 11 )
	 {
	   // This is a true electron, but it comes from something else
	   const reco::Candidate *mom = el.mother(0);
	   int momPid = -999;
	   if ( mom != 0 )
	     momPid = mom->pdgId();
	   printf("pid= %d  status= %d isFromZ= %d isFromW= %d  isFromTau= %d  momPid= %d\n", 
		  pid,  status, isFromZ, isFromW, isFromTau, momPid);
	   result = TRUE_NON_PROMPT_ELECTRON;
	 } else {
	 printf("The reco electron has a truth match with pid= %d\n", pid);
       }
     }

     return result;
}

void ElectronNtuplerEventStructure::printAllZeroMothers(const reco::Candidate *particle){
  
  if( particle == 0 ){
    printf("ElectronNtuplerEventStructure::printAllZeroMothers: reached the top of the decay tree\n");
    return;
  }
  
  printf("ElectronNtuplerEventStructure::printAllZeroMothers: ancestor ID= %d, status= %d\n",
	 particle->pdgId(), particle->status() );

  printAllZeroMothers( particle->mother(0) );

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronNtuplerEventStructure);
