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

//
// Signal sample: DYToLL
const TString fileNameSignal = 
  "/afs/cern.ch/user/r/rkamalie/workspace/public/DY_Spring15_Asympt50ns_24june2015.root";
// Directory and tree name:
const TString treeName = "ntupler/ElectronTree";

const bool verbose = false;
const bool smallEventCount = true; // DEBUG

// Selection cuts
// Kinematics
const float ptCut = 20; 

const float cutoffFraction = 0.90;

const int nEtaBins = 5;
const float etaBinLimits[nEtaBins+1] = {0.0, 0.8, 1.3, 2.0, 2.2, 2.5};

const int rhoBinsPlots  = 50;
const float rhoMinPlots = 0;
const float rhoMaxPlots = 50;

const int isoBinsPlots  = 1100;
const float isoMinPlots = -1;
const float isoMaxPlots = 10;

const float rhoMinFit   = 3;
const float rhoMaxFit   = 30;

//
// Forward declarations
//
void drawIsoVsRho(int etaBin, TH2F *hist);
void drawCutoffsAndFit(int etaBin, TH1F *hist, float &a, float &b);

//
// Main program
//

void computeIsolationCutoffs(EffectiveAreaType eaType = EA_NEUTRAL_TOTAL){

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 

  // General settings
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);

  // Book histograms
  TH2F *hIsoPhoNhVsRho[nEtaBins];
  TString hNameBase = "hIsoPhoNhVsRho";
  TH1F *hCutoffs[nEtaBins];
  TString hCutoffNameBase = "hCutoffs";
  for(int i=0; i<nEtaBins; i++){
    TString hName = hNameBase + TString::Format("_%d",i);
    hIsoPhoNhVsRho[i] = new TH2F(hName,"",
				 rhoBinsPlots, rhoMinPlots, rhoMaxPlots, 
				 isoBinsPlots, isoMinPlots, isoMaxPlots);
    hIsoPhoNhVsRho[i]->GetXaxis()->SetTitle("rho");
    hIsoPhoNhVsRho[i]->GetYaxis()->SetTitle("ISO_{pho}+ISO_{neu.had.}");

    TString hCutoffName = hCutoffNameBase + TString::Format("_%d",i);
    hCutoffs[i] = new TH1F(hCutoffName, "", rhoBinsPlots, rhoMinPlots, rhoMaxPlots);
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
  std::vector <float> *elePassConversionVeto = 0;


  // Declare branches
  TBranch *b_nEle = 0;
  TBranch *b_rho = 0;
  TBranch *b_elePt = 0;
  TBranch *b_eleEtaSC = 0;
  TBranch *b_elePhiSC = 0;
  TBranch *b_isoChargedHadrons = 0;
  TBranch *b_isoNeutralHadrons = 0;
  TBranch *b_isoPhotons = 0;
  TBranch *b_isTrue;
  // Other vars
  TBranch *b_elePassConversionVeto = 0;


  // Connect variables and branches to the tree with the data
  treeSignal->SetBranchAddress("nEle", &nEle, &b_nEle);
  treeSignal->SetBranchAddress("rho", &rho, &b_rho);
  treeSignal->SetBranchAddress("pt", &elePt, &b_elePt);
  treeSignal->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  treeSignal->SetBranchAddress("phiSC", &elePhiSC, &b_elePhiSC);
  treeSignal->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
  treeSignal->SetBranchAddress("isoNeutralHadrons", &isoNeutralHadrons, &b_isoNeutralHadrons);
  treeSignal->SetBranchAddress("isoPhotons",        &isoPhotons,        &b_isoPhotons);
  treeSignal->SetBranchAddress("isTrue",    &isTrue,    &b_isTrue);
  treeSignal->SetBranchAddress("passConversionVeto",       &elePassConversionVeto,
			       &b_elePassConversionVeto);


  // 
  // Loop over events
  //
  UInt_t maxEvents = treeSignal->GetEntries();
  if( smallEventCount )
    maxEvents = 2000000;
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
    b_rho->GetEntry(tentry);
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_elePhiSC->GetEntry(tentry);
    b_isoChargedHadrons->GetEntry(tentry);
    b_isoNeutralHadrons->GetEntry(tentry);
    b_isoPhotons->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);
    // Other vars
    b_elePassConversionVeto->GetEntry(tentry);

    // Nested loops over the electrons
    for(int iele = 0; iele < nEle; iele++){

      // Check kinematics:
      if( !(elePt->at(iele) > ptCut) )
	continue;

      // Check truth match
      if( isTrue->at(iele) != 1 ) continue;

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
      
      //if( iso>0 ){ // DEBUG
	hIsoPhoNhVsRho[ieta]->Fill( rho, iso);
	//}

    } // end loop over the electrons

  } // end loop over events
  printf("\n");

  float A[nEtaBins];
  float B[nEtaBins];
  for(int ieta=0; ieta<nEtaBins; ieta++){
    // Loop over rhos and find a cut-off for each rho for this eta range
    for(int iRho = 1; iRho <= rhoBinsPlots; iRho++){
      // First find the total number of electrons
      float total = 0;
      for(int iIso = 1; iIso<= isoBinsPlots+1; iIso++){ // Includes overflows!
	total += hIsoPhoNhVsRho[ieta]->GetBinContent(iRho, iIso);
      }
      // Second, find the cutoff
      float cutoff = 999;
      // No meaningful cutoff if statistics is low
      if( total > 10 ) {
	// Define max entry count within the cutoff
	float maxCount = cutoffFraction * total;
	float count = 0;
	for(int iIso = 1; iIso<= isoBinsPlots+1; iIso++){ // Includes overflows!
	  count += hIsoPhoNhVsRho[ieta]->GetBinContent(iRho, iIso);
	  if(count < maxCount ){
	    cutoff =  hIsoPhoNhVsRho[ieta]->GetYaxis()->GetBinUpEdge(iIso);
	  }else
	    break;
	}
      } // end if enough statistics
      hCutoffs[ieta]->SetBinContent(iRho, cutoff);
      
    } // end loop over rho

    drawIsoVsRho(ieta, hIsoPhoNhVsRho[ieta]);
    float a, b;
    drawCutoffsAndFit (ieta, hCutoffs[ieta], a, b);
    A[ieta] = a;
    B[ieta] = b;
  } //end loop over eta bins

  TString singleLineA = TString::Format("const float cutoff_A[%d] = { ",nEtaBins);
  TString singleLineB = TString::Format("const float cutoff_B[%d] = { ",nEtaBins);
  for(int ieta = 0; ieta<nEtaBins; ieta++){
    singleLineA += TString::Format("%7.4f", A[ieta]);
    singleLineB += TString::Format("%7.4f", B[ieta]);
    if( ieta == nEtaBins-1 ){
      singleLineA += TString("};\n");
      singleLineB += TString("};\n");
    }else{
      singleLineA += TString(", ");
      singleLineB += TString(", ");
    }
  }
  printf("\n%s", singleLineA.Data());
  printf("%s", singleLineB.Data());

  TFile *fout = new TFile("cutoffs.root","recreate");
  fout->cd();
  for(int ieta=0; ieta<nEtaBins; ieta++){
    hCutoffs[ieta]->Write();
  }
  fout->Close();
  delete fout;
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

void drawCutoffsAndFit(int etaBin, TH1F *hist, float &a, float &b){

  TString canvasName = "cutoffs_eta_";
  canvasName += etaBin;

  TCanvas *c1 = new TCanvas(canvasName,canvasName,10,10,600,600);
  c1->cd();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1);
  hist->GetYaxis()->SetRangeUser(-1,10);
  hist->Draw("P");

  TF1 *func = new TF1("func", "pol1",rhoMinFit, rhoMaxFit);
  hist->Fit("func","R");
  a = func->GetParameter(0);
  b = func->GetParameter(1);
  
  c1->Update();

		      
  return;
}

