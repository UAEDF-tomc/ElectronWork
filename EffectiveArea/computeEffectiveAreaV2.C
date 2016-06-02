#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "TStyle.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TTree.h"
#include "TList.h"
#include "TString.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TLegend.h"

#include <vector>

enum EffectiveAreaType {
  EA_CHARGED=0,
  EA_PHOTON,
  EA_NEUTRAL_HADRON,
  EA_NEUTRAL_TOTAL };

const TString eaTypeString[4] = {
  "charged",
  "photon",
  "neutral_hadron",
  "neutral_total"};

const TString eaTypeAxisLabel[4] = {
  "ISO_{ch.had}",
  "ISO_{pho}",
  "ISO_{neu.had.}",
  "ISO_{pho}+ISO_{neu.had.}"};

// For this exercise, we use two MC samples: the signal sample
// and a background-reach sample. Both have the same ntuple structure.
//
// Signal sample: DYToLL
const TString fileNameSignal = 
  "~/DYJetsToLL_madgraph_80X.root";
//"/afs/cern.ch/user/i/ikrav/workspace/ntuples/Spring16/DYJetsToLL_madgraph_80X.root";
// Directory and tree name:
const TString treeName = "ntupler/ElectronTree";

const bool verbose = false;
const bool smallEventCount = false; // DEBUG

const bool useWeights = true;

const float boundaryBarrelEndcap = 1.479;

// Selection cuts
// Kinematics
const float ptCut = 20; 

const int nEtaBins = 7;
const float etaBinLimits[nEtaBins+1] = {0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};

const int rhoBinsPlots  = 35;
const float rhoMinPlots = 0;
const float rhoMaxPlots = 35;

const float rhoMinFit   = 3;
const float rhoMaxFit   = 25;

//
// Forward declarations
//
void drawAndFitEA(EffectiveAreaType eaType,
		  int etaBin, TProfile *hist, float &area, float &areaErr);
void drawIsoVsRho(int etaBin, TH2F *hist);

//
// Main program
//

void computeEffectiveAreaV2(EffectiveAreaType eaType = EA_NEUTRAL_TOTAL){

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 

  // General settings
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);

  // Book histograms
  TH2F *hIsoPhoNhVsRho[nEtaBins];
  TProfile *profIsoPhoNhVsRho[nEtaBins];
  TString hNameBase = "hIsoPhoNhVsRho";
  TString profNameBase = "profIsoPhoNhVsRho";
  for(int i=0; i<nEtaBins; i++){
    TString hName = hNameBase + TString::Format("_%d",i);
    TString profName = profNameBase + TString::Format("_%d",i);
    hIsoPhoNhVsRho[i] = new TH2F(hName,"",
				 rhoBinsPlots, rhoMinPlots, rhoMaxPlots, 
				 1000, 0, 50);
    hIsoPhoNhVsRho[i]->Sumw2();
    hIsoPhoNhVsRho[i]->GetXaxis()->SetTitle("rho");
    hIsoPhoNhVsRho[i]->GetYaxis()->SetTitle(eaTypeAxisLabel[eaType]);
    profIsoPhoNhVsRho[i] = new TProfile(profName,"",
					rhoBinsPlots, rhoMinPlots, rhoMaxPlots);
    profIsoPhoNhVsRho[i]->Sumw2();
    profIsoPhoNhVsRho[i]->GetXaxis()->SetTitle("rho");
    TString averageString = TString::Format("<%s>", 
					    eaTypeAxisLabel[eaType].Data());
    profIsoPhoNhVsRho[i]->GetYaxis()->SetTitle(averageString);
  }

  //
  // Open a file and find the tree with electron data
  //
  TFile *fileSignal     = new TFile(fileNameSignal);
  if( !fileSignal ){
    printf("Failed to open the input files, check\n   %s\n", 
	   fileNameSignal.Data());
    assert(0);
  }
  TTree *treeSignal     = (TTree*)fileSignal->Get(treeName);
  if( !treeSignal ){
    printf("Failed to find the tree %s\n", treeName.Data() );
    assert(0);
  }

  // 
  // Set up the branches of interest
  //
  // Declare variables
  //
  // Event-level variables:
  int nEle; // the number of reconstructed electrons in the event
  float genWeight;
  float rho;
  // Per-eletron variables
  // Kinematics
  std::vector <float> *elePt = 0;         // electron PT
  std::vector <float> *eleEtaSC = 0;      // supercluser eta
  std::vector <float> *elePhiSC = 0;      // supercluser phi
  // Variables for analysis
  std::vector <float> *isoChargedHadrons = 0;
  std::vector <float> *isoNeutralHadrons = 0;
  std::vector <float> *isoPhotons = 0;
  std::vector <int> *isTrue = 0;
  // Other vars  
  // Impact parameters
  std::vector <float> *eleD0 = 0;      // r-phi plane impact parameter
  std::vector <float> *eleDZ = 0;      // r-z plane impact parameter
  // Matching track-supercluster
  std::vector <float> *eleDEtaIn = 0;  // deltaEtaIn
  std::vector <float> *eleDPhiIn = 0;  // deltaPhiIn
  // Misc ID variables
  std::vector <float> *eleHoverE = 0;  // H/E  
  std::vector <float> *eleFull5x5SigmaIEtaIEta = 0;  
  std::vector <float> *eleOOEMOOP = 0; // |1/E - 1/p|
  // Conversion rejection
  std::vector <float> *eleExpectedMissingInnerHits = 0;
  std::vector <float> *elePassConversionVeto = 0;


  // Declare branches
  TBranch *b_nEle = 0;
  TBranch *b_genWeight = 0;
  TBranch *b_rho = 0;
  TBranch *b_elePt = 0;
  TBranch *b_eleEtaSC = 0;
  TBranch *b_elePhiSC = 0;
  TBranch *b_isoChargedHadrons = 0;
  TBranch *b_isoNeutralHadrons = 0;
  TBranch *b_isoPhotons = 0;
  TBranch *b_isTrue;
  // Other vars
  TBranch *b_eleD0 = 0;
  TBranch *b_eleDZ = 0;
  TBranch *b_eleDEtaIn = 0;
  TBranch *b_eleDPhiIn = 0;
  TBranch *b_eleHoverE = 0;
  TBranch *b_eleFull5x5SigmaIEtaIEta = 0;
  TBranch *b_eleOOEMOOP = 0;
  TBranch *b_eleExpectedMissingInnerHits = 0;
  TBranch *b_elePassConversionVeto = 0;


  // Connect variables and branches to the tree with the data
  treeSignal->SetBranchAddress("nEle", &nEle, &b_nEle);
  treeSignal->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
  treeSignal->SetBranchAddress("rho", &rho, &b_rho);
  treeSignal->SetBranchAddress("pt", &elePt, &b_elePt);
  treeSignal->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  treeSignal->SetBranchAddress("phiSC", &elePhiSC, &b_elePhiSC);
  treeSignal->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
  treeSignal->SetBranchAddress("isoNeutralHadrons", &isoNeutralHadrons, &b_isoNeutralHadrons);
  treeSignal->SetBranchAddress("isoPhotons",        &isoPhotons,        &b_isoPhotons);
  treeSignal->SetBranchAddress("isTrue",    &isTrue,    &b_isTrue);
  treeSignal->SetBranchAddress("d0",                &eleD0,             &b_eleD0);
  treeSignal->SetBranchAddress("dz",                &eleDZ,             &b_eleDZ);
  treeSignal->SetBranchAddress("dEtaIn",            &eleDEtaIn,         &b_eleDEtaIn);
  treeSignal->SetBranchAddress("dPhiIn",            &eleDPhiIn,         &b_eleDPhiIn);
  treeSignal->SetBranchAddress("hOverE",            &eleHoverE,         &b_eleHoverE);
  treeSignal->SetBranchAddress("full5x5_sigmaIetaIeta", &eleFull5x5SigmaIEtaIEta,
			       &b_eleFull5x5SigmaIEtaIEta);
  treeSignal->SetBranchAddress("ooEmooP",           &eleOOEMOOP,        &b_eleOOEMOOP);
  treeSignal->SetBranchAddress("expectedMissingInnerHits", &eleExpectedMissingInnerHits, 
			       &b_eleExpectedMissingInnerHits);
  treeSignal->SetBranchAddress("passConversionVeto",       &elePassConversionVeto,
			       &b_elePassConversionVeto);


  // 
  // Loop over events
  //
  UInt_t maxEvents = treeSignal->GetEntries();
  if( smallEventCount )
    maxEvents = 1000000;
  if(verbose)
    printf("Start loop over events, total events = %lld\n", 
	   treeSignal->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = treeSignal->LoadTree(ievent);
    
    // Load the value of the number of the electrons in the event    
    b_nEle->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nEle);
    
    // Get data for all electrons in this event, only vars of interest
    b_genWeight->GetEntry(tentry);
    b_rho->GetEntry(tentry);
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_elePhiSC->GetEntry(tentry);
    b_isoChargedHadrons->GetEntry(tentry);
    b_isoNeutralHadrons->GetEntry(tentry);
    b_isoPhotons->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);
    // Other vars
    b_eleD0->GetEntry(tentry);
    b_eleDZ->GetEntry(tentry);
    b_eleDEtaIn->GetEntry(tentry);
    b_eleDPhiIn->GetEntry(tentry);
    b_eleHoverE->GetEntry(tentry);
    b_eleFull5x5SigmaIEtaIEta->GetEntry(tentry);
    b_eleOOEMOOP->GetEntry(tentry);
    b_eleExpectedMissingInnerHits->GetEntry(tentry);
    b_elePassConversionVeto->GetEntry(tentry);

    // Nested loops over the electrons
    for(int iele = 0; iele < nEle; iele++){

      // Check kinematics:
      if( !(elePt->at(iele) > ptCut) )
	continue;

      // Check truth match
      if( isTrue->at(iele) != 1 ) continue;

      // Loose ID of 2012 (VETO WP)
      const bool useID = false;
      if( useID ){
	if( abs(eleEtaSC->at(iele)) <  boundaryBarrelEndcap ){
	  if( abs(eleDEtaIn->at(iele))>0.007 ) continue;
	  if( abs(eleDPhiIn->at(iele))>0.8 ) continue;
	  if( eleFull5x5SigmaIEtaIEta->at(iele) >0.01 ) continue;
	  if( eleHoverE->at(iele) > 0.15 ) continue;
	  if( abs(eleD0->at(iele)) > 0.04 ) continue;
	  if( abs(eleDZ->at(iele)) > 0.2 ) continue;
	}else{
	  if( abs(eleDEtaIn->at(iele))>0.001 ) continue;
	  if( abs(eleDPhiIn->at(iele))>0.7 ) continue;
	  if( eleFull5x5SigmaIEtaIEta->at(iele) >0.03 ) continue;
	  if( abs(eleD0->at(iele)) > 0.04 ) continue;
	  if( abs(eleDZ->at(iele)) > 0.2 ) continue;
	}
      } // end if loose ID


      // Find eta bin
      if( abs(eleEtaSC->at(iele))>etaBinLimits[nEtaBins] ) continue;
      int ieta = 0; 
      while ( ieta < nEtaBins-1 
	      && abs(eleEtaSC->at(iele)) > etaBinLimits[ieta+1] )
	{ ++ieta; };

      // Look up the isolation type we need
      double iso = 0;
      if( eaType == EA_CHARGED ){
	iso = isoChargedHadrons->at(iele);
      }else if ( eaType == EA_PHOTON ) {
	iso = isoPhotons->at(iele);
      }else if ( eaType == EA_NEUTRAL_HADRON ) {
	iso = isoNeutralHadrons->at(iele);
      }else if ( eaType == EA_NEUTRAL_TOTAL ) {
	iso = isoNeutralHadrons->at(iele) + isoPhotons->at(iele);
      }else{
	printf("Unknown isolation type requested, exiting.\n");
	assert(0);
      }

      float weight = genWeight;
      if( !useWeights )
	weight = 1;

      //if( iso>0 ){ // DEBUG
      hIsoPhoNhVsRho[ieta]->Fill( rho, iso, weight);
      profIsoPhoNhVsRho[ieta]->Fill( rho, iso, weight);
	//}

    } // end loop over the electrons

  } // end loop over events
  printf("\n");

  // 
  // Loop over eta bins and fit the effective areas
  //
  float effArea[nEtaBins];
  float effAreaErr[nEtaBins];
  for(int ieta=0; ieta<nEtaBins; ieta++){
    drawAndFitEA(eaType, ieta, profIsoPhoNhVsRho[ieta], 
		 effArea[ieta], effAreaErr[ieta]);
    drawIsoVsRho(ieta, hIsoPhoNhVsRho[ieta]);
  }

  // 
  // Print the result
  // 
  TString outFileName = TString::Format("figures/ea_%s_isolation_constants.txt",
					eaTypeString[eaType].Data());
  ofstream outFile;
  outFile.open(outFileName.Data());
  // Header for the table
  TString firstLine = 
    TString::Format("\nEffective areas for %s isolation:\n", 
		    eaTypeString[eaType].Data());
  outFile << firstLine.Data();
  // Start building a line with a C++ array of effective areas
  TString singleLine = TString::Format("const float ea_%s_iso[%d] = { ",
				       eaTypeString[eaType].Data(), nEtaBins);
  // Loop over eta bins
  for(int i=0; i<nEtaBins; i++){
    TString thisLine = 
      TString::Format("eta bin [%4.2f, %4.2f]: EA = %7.4f +- %7.4f\n",
		      etaBinLimits[i], etaBinLimits[i+1], 
		      effArea[i], effAreaErr[i]);
    outFile << thisLine.Data();
    singleLine += TString::Format("%7.4f", effArea[i]);
    if( i == nEtaBins-1 )
      singleLine += TString("};\n");
    else
      singleLine += TString(", ");
  }
  outFile << endl;
  outFile << singleLine.Data();
  outFile.close();
  TString typeCommand = TString::Format("cat %s\n", outFileName.Data());
  gSystem->Exec(typeCommand.Data());

}

void drawAndFitEA(EffectiveAreaType eaType, int etaBin, 
		  TProfile *hist, float &area, float &areaErr){

  //
  // Draw plots
  //
  // TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  // c1->cd();
  // hIsoPhoNhVsRho[0]->Draw("colz");

  TString canvasName = "EAfit_eta_";
  canvasName += etaBin;

  TCanvas *c2 = new TCanvas(canvasName,canvasName,10,10,600,600);
  c2->cd();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1);
  hist->GetYaxis()->SetRangeUser(0,8);

  hist->Draw("pe");

  TF1 *flin = new TF1("flin","pol1",rhoMinFit, rhoMaxFit);
  hist->Fit("flin","R");
  
  // Retrieve the slope of the polynomial
  area = flin->GetParameter(1);
  areaErr = flin->GetParError(1);
  
  c2->Update();
  TString outPlotFileName = TString::Format("figures/ea_%s_isolation_eta%d.png",
					    eaTypeString[eaType].Data(),
					    etaBin);
  c2->Print(outPlotFileName);

  return;
}

void drawIsoVsRho(int etaBin, TH2F *hist){

  TString canvasName = "isoVSrho_eta_";
  canvasName += etaBin;

  TCanvas *c1 = new TCanvas(canvasName,canvasName,10,10,600,600);
  c1->cd();
  hist->Draw("colz");
  c1->Update();

  return;
}

