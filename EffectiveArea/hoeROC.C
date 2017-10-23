#include "TTree.h"
#include "TBranch.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLatex.h"

// Constants and settings

const bool verbose = false;
const bool smallEventCount = false;

const TString finSigName = "~/DYJetsToLL_cutID_tuning_92X_v1.root";
const TString finBgName  = "~/TTJets_cutID_92X_v1.root";

const TString treeName = "ntupler/ElectronTree";

const double ptMin = 20;
const double ptMax = 1000;

// Calibrated constants
bool useScaling = true;

const bool selectBarrel = true;

// Barrel constants
const float constRho = 0.0368;      // 90% contour, barrel, H/E ~= constRho*rho/E
const float constE   = 1.64;        // mean-based, barrel, H/E ~= constE/E
// const float constRho = 0.0368;      // 90% contour, barrel, H/E ~= constRho*rho/E
// const float constE   = 3.09;        // 90% contour, barrel, H/E ~= constE/E
 // const float constRho = 0.0368;     // 90% contour, barrel, H/E ~= constRho*rho/E
 // const float constE   = 3.91;       // 95% contour, barrel, H/E ~= constE/E
//
// const float constRho = 0.032;    // 95% contour, barrel, H/E ~= constRho*rho/E
// const float constE = 3.91;       // 95% contour, barrel, H/E ~= constE/E

// Endcap constants
// const float constRho = 0.201;       // 90% contour, endcap, H/E ~= constRho*rho/E
// const float constE   = 2.65;        // mean-based, endcap, H/E ~= constE/E
// const float constRho = 0.237;    // 95% contour, endcap, H/E ~= constRho*rho/E
// const float constE = 7.01;       // 95% contour, endcap, H/E ~= constE/E

const int nCuts = 1000;
const double cutMinVal = 0.00;
const double cutMaxVal = 0.50;
double cutVal[nCuts];

double effSigNum[nCuts];
double effSigDen[nCuts];
double effSig   [nCuts];

double effBgNum[nCuts];
double effBgDen[nCuts];
double effBg   [nCuts];

// Forward declarations
void computeEfficiency(bool selectTrue, bool selectBarrel);

// Main function
void hoeROC(){

  // Define the array of cut values
  for(int icut=0; icut<nCuts; icut++){
    cutVal[icut] = cutMinVal + icut*(cutMaxVal - cutMinVal) / (nCuts - 1);
  }

  // Compute signal efficiency
  bool selectTrue = true;
  computeEfficiency(selectTrue, selectBarrel);

  // Compute background efficiency
  selectTrue = false;
  computeEfficiency(selectTrue, selectBarrel);
  
  // Plot efficiencies vs cut
  TGraph *effSigGraph = new TGraph(nCuts, cutVal, effSig);
  //effSigGraph->Draw("ALP");

  TGraph *effBgGraph = new TGraph(nCuts, cutVal, effBg);
  effBgGraph->Draw("ALP");

  // Define ROC
  double rejBg[nCuts];
  for(int i=0; i<nCuts; i++)
    rejBg[i] = 1 - effBg[i];
  TGraph *rocGraph = new TGraph(nCuts, effSig, rejBg);
  rocGraph->Draw("ALP");

  TFile fout("roc.root","recreate");
  fout.cd();
  effSigGraph->Write("effSigGraph");
  effBgGraph->Write("effBgGraph");
  rocGraph->Write("rocGraph");
  fout.Close();
}

void computeEfficiency(bool selectTrue, bool selectBarrel){

  double *effNum = effSigNum;
  double *effDen = effSigDen;
  double *eff    = effSig;

  TString finName = finSigName;

  if( !selectTrue ){
    effNum = effBgNum;
    effDen = effBgDen;
    eff    = effBg;

    finName = finBgName;
  }

  // Reset the container arrays
  for(int i=0; i<nCuts; i++){
    effNum[i] = 0;
    effDen[i] = 0;
  }

  // Read in the trees
  TFile *fin = new TFile(finName);
  TTree *tree = (TTree*)fin->Get(treeName);
  
  // 
  // Set up the branches of interest
  //
  // Declare variables
  //
  // Event-level variables:
  int nEle; // the number of reconstructed electrons in the event
  float rho;
  // Per-electron variables
  std::vector <float> *elePt = 0;         // electron PT
  std::vector <float> *eleESC = 0;        // supercluser energy
  std::vector <float> *eleEtaSC = 0;      // supercluser eta
  std::vector <float> *eleHoverE = 0;     // H/E  
  std::vector <int> *isTrue = 0;          // truth matching
 
  // Declare branches
  TBranch *b_nEle = 0;
  TBranch *b_rho = 0;
  TBranch *b_elePt = 0;
  TBranch *b_eleESC = 0;
  TBranch *b_eleEtaSC = 0;
  TBranch *b_isTrue;
  TBranch *b_eleHoverE = 0;

  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nEle",  &nEle,      &b_nEle);
  tree->SetBranchAddress("rho",   &rho,       &b_rho);
  tree->SetBranchAddress("pt",    &elePt,     &b_elePt);
  tree->SetBranchAddress("eSC",   &eleESC,    &b_eleESC);
  tree->SetBranchAddress("etaSC", &eleEtaSC,  &b_eleEtaSC);
  tree->SetBranchAddress("isTrue",&isTrue,    &b_isTrue);
  tree->SetBranchAddress("hOverE",&eleHoverE, &b_eleHoverE);

  // 
  // Loop over events
  //
  UInt_t maxEvents = tree->GetEntries();
  if( smallEventCount )
    maxEvents = 100000;
  if(verbose)
    printf("Start loop over events, total events = %lld\n", 
           tree->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);
    
    // Load the value of the number of the electrons in the event    
    b_nEle->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nEle);
    
    // Get data for all electrons in this event, only vars of interest
    b_rho->GetEntry(tentry);
    b_elePt->GetEntry(tentry);
    b_eleESC->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);
    b_eleHoverE->GetEntry(tentry);
  
    // Nested loops over the electrons
    for(int iele = 0; iele < nEle; iele++){

      // Apply selection
      
      // True or fake electron
      if( selectTrue ){
	// Keep only true electrons, as requested
	if( !( isTrue->at(iele) == 1 ) )
	  continue;
      } else {
	// Keep only background electrons, as requested
	if( !( isTrue->at(iele)==0 || isTrue->at(iele)==3 ) )
	  continue;
      }

      // Pt range
      if( elePt->at(iele) < ptMin  )
        continue;
      if( elePt->at(iele) > ptMax  )
        continue;

      // Eta cuts
      if( selectBarrel ){
	// Keep only barrel electrons
	if( !(abs(eleEtaSC->at(iele))<1.4442) )
	  continue;
      }else{
	// Keep only barrel electrons
	if( !(abs(eleEtaSC->at(iele))>1.4442 
	      && abs(eleEtaSC->at(iele))<2.5) )
	  continue;
      }

      // Quantity to cut on
      float hOverE = eleHoverE->at(iele);
      float E = eleESC->at(iele);
      float adjustedHoverE = hOverE - constRho * rho / E - constE / E;
      if( !useScaling )
	adjustedHoverE = hOverE;

      // Loop over H/E cut values
      for(int icut = 0; icut<nCuts; icut++){

	bool pass = false;
	effDen[icut] += 1;
	if( adjustedHoverE < cutVal[icut] ){
	  effNum[icut] += 1;
	}

      } // end loop over cuts

    } // end loop over electrons

  } // end loop over events

  // Compute the efficiency dependence on cuts
  for(int i=0; i<nCuts; i++){
    eff[i] = effNum[i] / effDen[i];
    printf(" cut= %f    %f    %f   eff= %f\n", cutVal[i], effNum[i],
	   effDen[i], eff[i]);
  }

}

