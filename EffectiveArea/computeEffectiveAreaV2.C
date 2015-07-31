#include <iostream>
#include <fstream>
#include "TRandom.h"
#include "TArrow.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1D.h"
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

// For this exercise, we use two MC samples: the signal sample
// and a background-reach sample. Both have the same ntuple structure.
//
// Signal sample: DYToLL
const TString fileNameSignal = 
  "/afs/cern.ch/user/r/rkamalie/workspace/public/DY_Run2Asympt25ns_miniAOD_20july2015.root";
  //"/afs/cern.ch/user/r/rkamalie/workspace/public/DY_Run2Asympt50ns_miniAOD_21july2015.root";
  //"/afs/cern.ch/user/r/rkamalie/workspace/public/DY_Spring15_Asympt50ns_24june2015.root";
  //"/home/hep/ikrav/work/ntuples/PHYS14/DYJetsToLL_PU20bx25_event_structure.root";
  // "/home/hep/ikrav/work/ntuples/PHYS14/TTJets_PU20bx25_event_structure.root";
// Directory and tree name:
const TString treeName = "ntupler/ElectronTree";

const bool verbose = false;
const bool smallEventCount = false; // DEBUG

const bool useWeights = false;

// Use either the mean of iso(rho) or the XX% cutoff contour
const bool useIsoMean = false;
const float cutoffFrac = 0.90;

const float boundaryBarrelEndcap = 1.479;

// Selection cuts
// Kinematics
const float ptCut = 20; 

const int nEtaBins = 7;
//const float etaBinLimits[nEtaBins+1] = {0.0, 0.8, 1.3, 2.0, 2.2, 2.5};
const float etaBinLimits[nEtaBins+1] = {0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};

const int rhoBinsPlots  = 40;
const float rhoMinPlots = 0;
const float rhoMaxPlots = 40;

// Generous binning for cutoff computation
const int isoBinsPlots  = 400;
// Start a bit above zero, so that the iso=zero falls into underflow bin
const float isoMinPlots = 0.001;
const float isoMaxPlots = 20;

const float rhoMinFit   = 3;
const float rhoMaxFit   = 30;

enum SleepMode {UNSEEN = 0,
		SHORT_SLEEP = 1,
		CR = 2};
const int sleepMode = SHORT_SLEEP;

TCanvas *cutoffsCanv;

//
// Forward declarations
//
void drawAndFitEA(EffectiveAreaType eaType,
		  int etaBin, TH1D *hist, float &area, float &areaErr);
void drawIsoVsRho(int etaBin, TH2F *hist);

TH1D * getCutoffContour( TH2F *hist2d, TH1F *hIsoRel, float cutoffFrac, int etaBin );

float getCutoffValue(TH1D *hist, float frac);

void getCutoffError(TH1D *hist, TH1F* hrel, 
		    float rhoMin, float rhoMax, float cutoffFrac, 
		    float &cutoff, float &cutoffErr, int etaBin);

void drawProgressBar(float progress);

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

  // Set up canvas for recurrent plots
  cutoffsCanv = new TCanvas("cutoffsCanv","cutoffsCanv",10,10,900,500);
  cutoffsCanv->Divide(2,1);

  // Book histograms
  TH2F *hIsoVsRho[nEtaBins];
  TH1F *hIsoRelToRho[nEtaBins];
  TProfile *profIsoVsRho[nEtaBins];
  TString hNameBase = "hIsoVsRho";
  TString profNameBase = "profIsoVsRho";
  for(int i=0; i<nEtaBins; i++){
    TString hName = hNameBase + TString::Format("_%d",i);
    TString profName = profNameBase + TString::Format("_%d",i);
    hIsoVsRho[i] = new TH2F(hName,"",
				 rhoBinsPlots, rhoMinPlots, rhoMaxPlots, 
				 isoBinsPlots, isoMinPlots, isoMaxPlots);
    hIsoVsRho[i]->Sumw2();
    hIsoVsRho[i]->GetXaxis()->SetTitle("rho");
    hIsoVsRho[i]->GetYaxis()->SetTitle("ISO_{pho}+ISO_{neu.had.}");
    //
    TString hName2 = TString::Format("hIsoRelToRho_%d",i);
    // Since the histogram below contains isolation/rho ratio,
    // increase the number of bins.
    hIsoRelToRho[i] = new TH1F(hName2, "", isoBinsPlots*10, isoMinPlots, isoMaxPlots);
    hIsoRelToRho[i]->Sumw2();
    //
    profIsoVsRho[i] = new TProfile(profName,"",
					rhoBinsPlots, rhoMinPlots, rhoMaxPlots);
    profIsoVsRho[i]->Sumw2();
    profIsoVsRho[i]->GetXaxis()->SetTitle("rho");
    profIsoVsRho[i]->GetYaxis()->SetTitle("<ISO_{pho}+ISO_{neu.had.}>");
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

    if( ievent%100000 == 0 || ievent == maxEvents-1){
      // printf("."); fflush(stdout);
      drawProgressBar( (1.0*ievent+1)/maxEvents);
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
    b_eleExpectedMissingInnerHits->GetEntry(tentry);
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

      float weight = genWeight;
      if( !useWeights )
	weight = 1;

      hIsoVsRho[ieta]->Fill( rho, iso, weight);
      profIsoVsRho[ieta]->Fill( rho, iso, weight);
      hIsoRelToRho[ieta]->Fill( iso/rho, weight);

    } // end loop over the electrons

  } // end loop over events
  printf("\n");

  // 
  // Loop over eta bins and fit the effective areas
  //
  float effArea[nEtaBins];
  float effAreaErr[nEtaBins];
  for(int ieta=0; ieta<nEtaBins; ieta++){
    // Prepare a histogram for effective area extraction
    // TProfile inherits from TH1D, so we can just do this assignment below.
    TH1D *hist = profIsoVsRho[ieta];
    if( !useIsoMean )
      hist = getCutoffContour( hIsoVsRho[ieta], hIsoRelToRho[ieta],
			       cutoffFrac, ieta );
    drawAndFitEA(eaType, ieta, hist,
		 effArea[ieta], effAreaErr[ieta]);
    drawIsoVsRho(ieta, hIsoVsRho[ieta]);
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
		  TH1D *hist, float &area, float &areaErr){

  //
  // Draw plots
  //
  // TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  // c1->cd();
  // hIsoVsRho[0]->Draw("colz");

  TString canvasName = "EAfit_eta_";
  canvasName += etaBin;

  TCanvas *c2 = new TCanvas(canvasName,canvasName,10,10,600,600);
  c2->cd();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1);
  hist->GetYaxis()->SetRangeUser(0,15);

  hist->Draw("pe");

  TF1 *flin = new TF1("flin","pol1",rhoMinFit, rhoMaxFit);
  flin->SetParameters(1,0.1);
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

  TString outPlotFileName = TString::Format("figures/isoVsRho_eta%d.png",
					    etaBin);
  c1->Print(outPlotFileName);

  return;
}

TH1D * getCutoffContour(TH2F *hist2d, TH1F *hIsoRel, float cutoffFrac, int etaBin){

  TString hname = TString::Format("hIsoCutoffs_%d", etaBin);
  TH1D *hCutoffs = new TH1D(hname,"",rhoBinsPlots, rhoMinPlots, rhoMaxPlots); 

  // Loop over rhos and find a cut-off for each rho for this eta range
  for(int iRho = 1; iRho <= rhoBinsPlots; iRho++){

    // Prepare a slice for given rho bin
    TH1D *slice = hist2d->ProjectionY("tmp",iRho,iRho);
    float cutoff = getCutoffValue( slice, cutoffFrac);
    float cutoffErr = 0.5; // some dummy value
    hCutoffs->SetBinContent(iRho, cutoff);
    if( !useWeights) {
      float rhoMin = hist2d->GetXaxis()->GetBinLowEdge(iRho);
      float rhoMax = hist2d->GetXaxis()->GetBinUpEdge(iRho);
      getCutoffError(slice, hIsoRel, rhoMin, rhoMax, cutoffFrac, cutoff, cutoffErr, etaBin);
    }else{
      // We can't do proper estimate of the cutoff errors with
      // genWeights that can be positive or negative with the current
      // method, so set the error to some dummy value, set above already.
    }      
    hCutoffs->SetBinError(iRho, cutoffErr);
    
  } // end loop over rho

  return hCutoffs;
}


float getCutoffValue(TH1D *hist, float frac){

  // First find the total number of electrons
  float total = 0;
  int nbins = hist->GetNbinsX();
  for(int ibin = 0; ibin<= nbins+1; ibin++){ 
    // This loop includes underflows and overflows!
    total += hist->GetBinContent(ibin);
  }
  
  // Second, find the cutoff
  float cutoff = 999;
  // No meaningful cutoff if statistics is low
  if( total < 5 )
    return cutoff;
  
  // Define max entry count within the cutoff
  float maxCount = frac * total;
  float count = 0;
  for(int ibin = 0; ibin<= nbins+1; ibin++){ 
    // This loop includes underflows and overflows!
    count += hist->GetBinContent(ibin);
    if(count < maxCount ){
      cutoff =  hist->GetXaxis()->GetBinUpEdge(ibin);
    }else
      break;
  }
  
  return cutoff;
}

// This function computes the error on the cutoff, and updates
// the cutoff itself if it is 999, based on toy simulation.
void getCutoffError(TH1D *hist, TH1F* hrel, 
		    float rhoMin, float rhoMax, float cutoffFrac, 
		    float &cutoff, float &cutoffErr, int etaBin){

  // The idea for each toy sample: we use the relative-to-rho
  // isolation distribution shape (distribution of Isolation/Rho ratio)
  // to generate the isolation distribution for each toy sample.
  // The ensemble is generated for a specific rho range [rhomin, rhomax],
  // and each toy event's isolation is a random number generated 
  // randomly from the (Iso/Rho) histogram, multiplied by rhoAve = (rhomax+rhomin)/2.
  // This obviously assumes that the shape of the non-zero part of the
  // isolation distribution scales with rho.
  // The zero component of the isolation distribution doesn't scale 
  // with rho the same way, so it is treated separately in toys.

  float rhoAve = (rhoMax+rhoMin)/2;

  // A histogram for a visual representation of the process
  TH1D *htoyGlobal = (TH1D*)hist->Clone("htoyGlobal");
  htoyGlobal->Reset();
  htoyGlobal->SetLineColor(kRed);

  // A histogram to hold the distribution of cutoffs across all toys
  TH1D *hcutoffs = new TH1D("hcutoffs","", isoBinsPlots*10, isoMinPlots, isoMaxPlots);
  hcutoffs->Reset();

  // Run an ensemble of toys
  const int ntoys = 100;
  int toySampleSize = hist->GetSumOfWeights();
  int toySampleZeroBin = hist->GetBinContent(0);
  for(int itoy = 0; itoy<ntoys; itoy++){

    // Generate toy sample size: the zero bin and rest    
    int thisToySampleSizeNonZero =  gRandom->Poisson(toySampleSize-toySampleZeroBin);
    int thisToySampleSizeZero = gRandom->Poisson(toySampleZeroBin);
    
    // Generate all toy events for this toy sample
    TH1D *htoy = (TH1D*)hist->Clone("htoy");
    htoy->Reset();
    for(int i=0; i<thisToySampleSizeNonZero; i++) {
      float thisVal = hrel->GetRandom()*rhoAve;
      htoy->Fill( thisVal );
      htoyGlobal->Fill(thisVal);
    }
    // Fill the underflow bin
    htoy->SetBinContent(0, thisToySampleSizeZero);

    float cutoffTmp = getCutoffValue(htoy, cutoffFrac);
    hcutoffs->Fill(cutoffTmp);

    delete htoy;
  } // end loop over toys
  htoyGlobal->Scale( hist->GetSumOfWeights()/htoyGlobal->GetSumOfWeights());

  cutoffsCanv->cd(1);
  hcutoffs->GetXaxis()->SetTitle("cutof position [GeV]");
  hcutoffs->Draw();

  cutoffErr = hcutoffs->GetRMS();
  if( hcutoffs->GetSumOfWeights() == 0 ){
    // If no entries within the limits of the histogram,
    // set the error to some large value.
    // This could happen for example if the size of the toy sample
    // is zero.
    cutoffErr = 999;    
  }
  if( cutoff == 999 ){
    printf(" low statistics in this sample, taking cutoff from toys\n");
    cutoff = hcutoffs->GetMean();
  }
  //printf("cutoff = %f +- %f\n", cutoff, cutoffErr);

  cutoffsCanv->cd(2);
  hist->Draw("pe");
  htoyGlobal->Draw("hist,same");

  TLegend *leg = new TLegend(0.4,0.7, 0.85, 0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hist,"MC distribution","l");
  leg->AddEntry(htoyGlobal, "Shape from toy MC","l");
  leg->Draw("same");

  float histmax = hist->GetMaximum();
  TArrow *arr2 = new TArrow(cutoff-cutoffErr, histmax*0.5,
			    cutoff+cutoffErr, histmax*0.5,
			    0.05, "|-|");
  TArrow *arr = new TArrow(cutoff, histmax*0.5, cutoff, 0);
  arr->SetLineColor(kBlue);
  arr->SetLineWidth(2);
  arr->Draw();

  arr2->SetLineColor(kBlue);
  arr2->SetLineWidth(2);
  arr2->Draw();
  
  if( sleepMode == UNSEEN ){
    // Do not even show anything
  }else{
    cutoffsCanv->Update();
    printf("rhoAve=%.1f   cutoff = %f +- %f\n", 
	   rhoAve, cutoff, cutoffErr);
    if( sleepMode == SHORT_SLEEP ){
      gSystem->Sleep(300);    
    }else if( sleepMode == CR ){
      getchar();
    }else{
      printf("Unknown sleep mode, die\n");
      assert(0);
    }
    // Save into a file
    TString outPlotFileName = TString::Format("figures/cutoffErrCanvas_rho%2.1f-%2.1f_eta%d.root",
					      rhoMin, rhoMax, 
					      etaBin);
    cutoffsCanv->SaveAs(outPlotFileName);

  }
  


  delete htoyGlobal;
  delete hcutoffs;
  delete arr;
  delete arr2;
			    
  return;
}

void drawProgressBar(float progress){

  const int barWidth = 70;
  
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
  
  if( progress >= 1.0 )
    std::cout << std::endl;

  return;
}
