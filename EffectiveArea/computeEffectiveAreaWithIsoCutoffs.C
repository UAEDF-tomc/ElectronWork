#include <iostream>
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
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
#include <cassert>
#include <stdio.h>
#include <stdlib.h>

enum EffectiveAreaType {
  EA_UNDEFINED=-1,
  EA_CHARGED=0,
  EA_PHOTON,
  EA_NEUTRAL_HADRON,
  EA_NEUTRAL_TOTAL };

const TString eaTypeString[4] = {
  "charged hadron",
  "photon",
  "neutral hadron",
  "neutral hadron and photon"};

enum MethodType {
  METHOD_UNDEFINED = -1,
  METHOD_TOY_MC,
  METHOD_EFF_CURVE};


//
// Signal sample: DYToLL
const TString fileNameSignal = 
  // "/afs/cern.ch/user/i/ikrav/workspace/ntuples/Spring16/DYJetsToLL_madgraph_80X.root";
  "/user/tomc/eleIdTuning/tuples/DYJetsToLL_cutID_tuning_94X_v1.root";
// Directory and tree name:
const TString treeName = "ntupler/ElectronTree";

const bool verbose = false;
const bool smallEventCount = false; // DEBUG

const float minEntries = 100;

const MethodType method = METHOD_EFF_CURVE;

// Selection cuts
// Kinematics
const float ptCut = 20; 

const float cutoffFraction = 0.90;

const int nEtaBins = 7;
// const float etaBinLimits[nEtaBins+1] = {0.0, 0.8, 1.3, 2.0, 2.2, 2.5};
//const float etaBinLimits[nEtaBins+1] = {0.0, 1.0, 1.479, 2.0, 2.2, 2.5};
const float etaBinLimits[nEtaBins+1] = 
  {0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};

const int rhoBinsPlots  = 50;
const float rhoMinPlots = 0;
const float rhoMaxPlots = 50;

const int isoBinsPlots  = 1100;
const float isoMinPlots = -1;
const float isoMaxPlots = 10;

// Limit the fit range to avoid low statistics edge effects
const float rhoMinFit   = 5;
const float rhoMaxFit   = 35;

// Global variables
TH2F *dummy; // A histogram for drawing graphs on top of it
EffectiveAreaType eaTypeGlobal = EA_UNDEFINED;

//
// Forward declarations
//
void drawIsoVsRho(int etaBin, TH2F *hist);
void computeCutoff(TH1D *hist, float &total, float &cutoff);
void computeCutoffAndErrorMethodToy(TH1D *hist, float &cutoff, float &cutoffErrPlus, 
                              float &cutoffErrMinus, TCanvas *canv);
void computeCutoffAndErrorMethodEff(TH1D *hist, float &cutoff, float &cutoffErrPlus,
                              float &cutoffErrMinus, TCanvas *canv);
void interpolate( float x1, float x2, float y1, float y2, float &x, float y);
void drawCutoffsAndFit(int etaBin, TH1F *hist, TGraphAsymmErrors *graph, float &a, float &b, float &bErr);

double box(double x);
Double_t isoShape(Double_t *x, Double_t *par);

//
// Main program
//

void computeEffectiveAraWithIsoCutoffs(EffectiveAreaType eaType = EA_NEUTRAL_TOTAL){
  system("mkdir -p effAreas/" +  eaTypeString[eaTypeGlobal] + "/");
  eaTypeGlobal = eaType;

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 

  // General settings
  gStyle->SetOptFit();
  gStyle->SetOptStat();

  // Book histograms
  TH2F *hIsoPhoNhVsRho[nEtaBins];
  TString hNameBase = "hIsoPhoNhVsRho";
  TH1F *hCutoffs[nEtaBins];
  TGraphAsymmErrors *hCutoffsGraph[nEtaBins];
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
    TString hCutoffGraphName = hCutoffName + TString("_graph");
    hCutoffsGraph[i] = new TGraphAsymmErrors(rhoBinsPlots);
    hCutoffsGraph[i]->SetName(hCutoffGraphName);
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
  float BErr[nEtaBins];
  float cutoff, cutoffErrPlus, cutoffErrMinus;
  TCanvas *canv = new TCanvas("slices","",10,10,600,600);
  //canv->SetLogy();
  for(int ieta=0; ieta<nEtaBins; ieta++){
    // Loop over rhos and find a cut-off for each rho for this eta range
    for(int iRho = 1; iRho <= rhoBinsPlots; iRho++){
      // Create a rho slice for the 2D histogram of the given eta bin
      TH1D *hSlice =  hIsoPhoNhVsRho[ieta]->ProjectionY("_py",iRho, iRho);
      hSlice->Rebin(5);
      if( method == METHOD_TOY_MC ){
        // Method with errors on the cutoff based on a toy ensemble
        computeCutoffAndErrorMethodToy(hSlice, cutoff, cutoffErrPlus, cutoffErrMinus, canv);
      }else if( method == METHOD_EFF_CURVE ){
        // Method with errors on the cutoff based on the efficiency curve analysis
        computeCutoffAndErrorMethodEff(hSlice, cutoff, cutoffErrPlus, cutoffErrMinus, canv);
      }else{
        printf("Unknown method, crashing\n");
        assert(0);
      }
      hCutoffs[ieta]->SetBinContent(iRho, cutoff);
      hCutoffs[ieta]->SetBinError(iRho, std::min(cutoffErrPlus, cutoffErrMinus));
      float x = hCutoffs[ieta]->GetBinCenter(iRho);
      float dx = 0.5 * hCutoffs[ieta]->GetBinWidth(iRho);
      hCutoffsGraph[ieta]->SetPoint(iRho, x, cutoff);
      hCutoffsGraph[ieta]->SetPointError(iRho, dx, dx, cutoffErrMinus, cutoffErrPlus);
      
      delete hSlice;
    } // end loop over rho

    drawIsoVsRho(ieta, hIsoPhoNhVsRho[ieta]);
    float a, b, bErr;
    drawCutoffsAndFit (ieta, hCutoffs[ieta], hCutoffsGraph[ieta], a, b, bErr);
    A[ieta] = a;
    B[ieta] = b;
    BErr[ieta] = bErr;
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

  // Print the effective area constants in CMSSW-like format
  FILE *f = fopen("effAreas/" +  eaTypeString[eaTypeGlobal] + "/effAreas.txt", "w");
  fprintf(f,"\nCMSSW-like printout of the effective areas\n");
  fprintf(f,"# |eta| min   |eta| max   effective area    error\n");
  for(int ieta = 0; ieta<nEtaBins; ieta++){
    fprintf(f,"%10.3f   %10.3f   %10.4f   %10.4f\n", 
           etaBinLimits[ieta], etaBinLimits[ieta+1],
           B[ieta], BErr[ieta]);
  }
  fprintf(f,"\n");
  fclose(f);

  TFile *fout = new TFile("effAreas/" +  eaTypeString[eaTypeGlobal] + "/cutoffs.root","recreate");
  fout->cd();
  for(int ieta=0; ieta<nEtaBins; ieta++){
    hCutoffs[ieta]->Write();
    hCutoffsGraph[ieta]->Write();
    hIsoPhoNhVsRho[ieta]->Write();
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

void computeCutoff(TH1D *hist, float &total, float &cutoff){

  total = 0;
  cutoff = 0;

  // First find the total number of electrons
  int histBins = hist->GetNbinsX();
  for(int iIso = 1; iIso<= histBins+1; iIso++){ // Includes overflows!
    total += hist->GetBinContent(iIso);
  }
  // printf("Total=%f\n", total);

  // Second, find the cutoff
  // Define max entry count within the cutoff
  float maxCount = cutoffFraction * total;
  float count = 0;
  for(int iIso = 1; iIso<= histBins+1; iIso++){ // Includes overflows!
    count += hist->GetBinContent(iIso);
    if(count < maxCount ){
      cutoff =  hist->GetXaxis()->GetBinUpEdge(iIso);
      // printf("   x= %f   count= %f   maxcount= %f\n", cutoff, count, maxCount);
    }else{
      // printf("   bailing out on i= %d \n", iIso);
      break;
    }
  }
  return;
}

void computeCutoffAndErrorMethodToy(TH1D *hist, float &cutoff, float &cutoffErrPlus, 
                                    float &cutoffErrMinus, TCanvas *canv){


  float total = 0;
  computeCutoff(hist, total, cutoff);
  
  const float largeCutoffError = 10;
  float cutoffErr = 0;
  if( total < 500 ) {
    // not enough data for reliable cutoff computations
    cutoff = 999;
    cutoffErr = largeCutoffError;
  }else{
    // sufficient statistics, proceed with error calculations
    
    // Compute the error
    TF1 *isofunc = new TF1("isofunc","isoShape",-0.1, 6, 4);
    isofunc->SetParLimits(0, 0, 1e8);
    isofunc->SetParLimits(1, 0, 1e8);
    isofunc->SetParLimits(2, 0.1, 10);
    isofunc->SetParLimits(3, 0.0, 10);
    isofunc->SetParameters(10,10,2,1);
    
    isofunc->SetNpx(1000);
    isofunc->SetLineColor(kRed);
    hist->Fit("isofunc","R");
    
    isofunc->SetRange( hist->GetXaxis()->GetBinLowEdge(1), 100);
                       // hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()));

    const int nPseudoExp = 10;
    // Create a histogram to hold pseudo eperiments so that it contains the long tail
    TH1D *tmpExperimentHist = new TH1D("tmpExperimentHist","",
                                       10000, hist->GetXaxis()->GetBinLowEdge(1),
                                       100);
    // TH1D *tmpExperimentHist = (TH1D*)hist->Clone("tmpExperimentHist");
    TH1D *tmpCutoffHist = (TH1D*)hist->Clone("tmpCutoffHist");
    tmpCutoffHist->Reset();
    float tmpTotal, tmpCutoff;
    for(int iexp = 0; iexp < nPseudoExp; iexp++){
      tmpExperimentHist->Reset();
      tmpExperimentHist->FillRandom("isofunc", hist->GetSumOfWeights() );
      // tmpExperimentHist->Print("all");
      computeCutoff( tmpExperimentHist, tmpTotal, tmpCutoff );
      tmpCutoffHist->Fill( tmpCutoff );
      // printf("pseudo %3d   cutpff= %f\n", iexp, tmpCutoff);
    }
    // tmpCutoffHist->Print();
    cutoffErr = tmpCutoffHist->GetRMS();
    if( cutoffErr == 0 )
      cutoffErr = largeCutoffError;
    printf("cutoff= %f   cutoff toy= %f   cutoffErr= %f\n", cutoff, tmpCutoffHist->GetMean(), cutoffErr);
    delete tmpExperimentHist;
    delete tmpCutoffHist;
    //delete isofunc; // deleting function causes crash for some reason
  } // end if enough statistics
  
  cutoffErrPlus = cutoffErr;
  cutoffErrMinus = cutoffErr;

  canv->cd();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1);
  hist->Draw("pe");
  canv->Update();
  // float xxx;
  // cin >> xxx;
}

void computeCutoffAndErrorMethodEff(TH1D *hist, float &cutoff, float &cutoffErrPlus,
                                    float &cutoffErrMinus, TCanvas *canv){

  int   nbins = hist->GetNbinsX();
  float xlow = hist->GetXaxis()->GetBinLowEdge(1);
  float xhigh = hist->GetXaxis()->GetBinUpEdge(nbins);

  // initialize to some large values
  cutoff = 2*xhigh;
  cutoffErrPlus = 2*xhigh;
  cutoffErrMinus = 2*xhigh;
  
  // Compute the total event counts
  float total = 0;
  float totalNonZero =0;
  float totalErr = 0;
  for(int i=1; i<=nbins+1; i++) { // include overflows
    total += hist->GetBinContent(i);
    totalErr += hist->GetBinError(i) * hist->GetBinError(i);
    if( ! (hist->GetXaxis()->GetBinLowEdge(i) <= 0) )
      totalNonZero += hist->GetBinContent(i);
  }
  totalErr = sqrt(totalErr);
  if( total < minEntries)
    return;

  TGraphErrors *grEff = new TGraphErrors(0);
  grEff->SetMarkerStyle(20);
  grEff->SetMarkerSize(0.5);
  grEff->SetFillColor(kMagenta);

  for(int i=1; i<=nbins; i++){
    
    // Find numerator for the efficiency
    float inRange = 0;
    float inRangeErr = 0;
    for(int irange = 1; irange <= i; irange++){
      inRange += hist->GetBinContent(irange);
      inRangeErr += hist->GetBinError(irange) * hist->GetBinError(irange);
    }
    inRangeErr = sqrt(inRangeErr);
    // Find the efficiency
    float eff = 0;
    float effErr = 1;
    if(total>0){
      eff = inRange/total;
      effErr = sqrt( eff*(1-eff)/totalNonZero );
    }
    // Save into graph
    float cutVal = hist->GetXaxis()->GetBinUpEdge(i);
    float cutValErr = hist->GetXaxis()->GetBinWidth(i)/2.0;
    int newPointIndex = grEff->GetN();
    grEff->SetPoint( newPointIndex, cutVal, eff);
    grEff->SetPointError( newPointIndex, cutValErr, effErr);
  } 

  // Next, find the efficiency cutoff value and its errors
  float cutoffLow = xlow;
  float cutoffHigh = xlow;
  cutoff = xlow;
  int npoints = grEff->GetN();
  Double_t *xArray      = grEff->GetX();
  Double_t *effArray    = grEff->GetY();
  Double_t *effErrArray = grEff->GetEY();
  // Central values
  for(int i=1; i < npoints; i++){ // Do not start with point zero
    if( effArray[i] > cutoffFraction){
      // Found the first efficiency value above threshold
      // Interpolate and find the cutoff value 
      interpolate( xArray[i-1], xArray[i], effArray[i-1],effArray[i],
                   cutoff, cutoffFraction);
      break;
    }
  }
  // Error up
  for(int i=1; i < npoints; i++){ // Do not start with point zero
    if( effArray[i] + effErrArray[i] > cutoffFraction){
      // Found the first efficiency value above threshold
      // Interpolate and find the cutoff value
      interpolate( xArray[i-1], xArray[i], 
                   effArray[i-1] + effErrArray[i-1],
                   effArray[i] + effErrArray[i],
                   cutoffLow, cutoffFraction);
      break;
    }
  }
  // Error down
  for(int i=1; i < npoints; i++){ // Do not start with point zero
    if( effArray[i] - effErrArray[i] > cutoffFraction){
      // Found the first efficiency value above threshold
      // Interpolate and find the cutoff value
      interpolate( xArray[i-1], xArray[i], 
                   effArray[i-1] - effErrArray[i-1],
                   effArray[i] - effErrArray[i],
                   cutoffHigh, cutoffFraction);
      break;
    }
  }
  cutoffErrPlus = cutoffHigh - cutoff;
  cutoffErrMinus = cutoff - cutoffLow;


  printf("Method eff curve: Cutoff is at %f + %f - %f\n", cutoff, 
         cutoffErrPlus,
         cutoffErrMinus);


  //
  // Draw the result
  //
  canv->cd();
  gStyle->SetOptStat(0);

  if(dummy != 0)
    delete dummy;
  dummy = new TH2F("dummy","",100, xlow, xhigh, 100, 0, 1.2);
  dummy->GetXaxis()->SetTitle("cut value");
  dummy->GetYaxis()->SetTitle("efficiency");
  dummy->Draw();

  grEff->Draw("E3,same");
  grEff->SetLineWidth(0); // suppress drawing error bars
  grEff->DrawClone("P,same");
  // Draw lines to show the cutoff, etc
  TLine *hline = new TLine(xlow, cutoffFraction, xhigh, cutoffFraction);
  hline->Draw("same");
  TLine *vline1 = new TLine(cutoffLow,  0, cutoffLow,  1.2);
  TLine *vline2 = new TLine(cutoff,     0, cutoff,     1.2);
  TLine *vline3 = new TLine(cutoffHigh, 0, cutoffHigh, 1.2);
  vline1->Draw("same");
  vline2->Draw("same");
  vline3->Draw("same");
  canv->Update();

}

void interpolate( float x1, float x2, float y1, float y2, float &x, float y){

  x = x1 + (y-y1)*(x2-x1)/(y2-y1);

  // printf("x1= %f  x2= %f  y1= %f  y2= %f  x= %f  y= %f\n",
  //     x1, x2, y1, y2, x, y);
  return;
}

void drawCutoffsAndFit(int etaBin, TH1F *hist, TGraphAsymmErrors *graph, float &a, 
                       float &b, float &bErr){

  printf("Start fitting\n");

  TString canvasName = "cutoffs_eta_";
  canvasName += etaBin;

  TCanvas *c1 = new TCanvas(canvasName,canvasName,10,10,600,600);
  c1->cd();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1);
  hist->GetYaxis()->SetRangeUser(-1,10);
  TString yaxisTitle = eaTypeString[eaTypeGlobal] + TString(" isolation");
  hist->GetXaxis()->SetTitle("rho");
  hist->GetYaxis()->SetTitle(yaxisTitle);
  hist->Draw("P");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->Draw("PE,same");

  TF1 *func = new TF1("func", "pol1",rhoMinFit, rhoMaxFit);
  graph->Fit("func","R");
  printf("Finished fitting\n");
  printf("The histogram is\n");
  a = func->GetParameter(0);
  b = func->GetParameter(1);
  bErr = func->GetParError(1);
  
  c1->Update();
  TString cFileName = TString("effAreas/" +  eaTypeString[eaTypeGlobal] + "/") + canvasName + TString(".png");
  c1->Print(cFileName);
                      
  return;
}


const double boxlow = -0.05;
const double boxhigh = 0.05;

double box(double x){

  double result = 0;
  if( x > boxlow && x < boxhigh )
    result = 1;

  return result;
}

Double_t isoShape(Double_t *x, Double_t *par) 
{

  Double_t xx = x[0];
  Double_t boxNorm        = par[0];
  Double_t landauNorm     = par[1];
  Double_t landauLocation = par[2];
  Double_t landauScale    = par[3];

  Double_t result = boxNorm*box(xx) 
    + landauNorm * TMath::Landau(xx, landauLocation, landauScale);
  
  return result;
}

// Compiled
int main(int argc, char *argv[]){
  gROOT->SetBatch();
  computeEffectiveAraWithIsoCutoffs(EA_CHARGED);
  computeEffectiveAraWithIsoCutoffs(EA_PHOTON);
  computeEffectiveAraWithIsoCutoffs(EA_NEUTRAL_HADRON);
  computeEffectiveAraWithIsoCutoffs(EA_NEUTRAL_TOTAL);
}
