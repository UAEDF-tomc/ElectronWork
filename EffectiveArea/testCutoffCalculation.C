#include "TGraphErrors.h"
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
#include "TLine.h"


const TString fname = "cutoffs.root";
const TString hname = "hIsoPhoNhVsRho_1";

const double effCutoff = 0.90;
const double cutoffFraction = effCutoff;
const int sliceIndex = 25;

void methodEffCurve(TH1D *hist);
void methodToy(TH1D *hist);
void interpolate( float x1, float x2, float y1, float y2, float &x, float y);
void computeCutoff(TH1D *hist, float &total, float &cutoff);
void computeCutoffAndError(TH1D *hist, float &cutoff, float &cutoffErr, TCanvas *canv);
double box(double x);
Double_t isoShape(Double_t *x, Double_t *par);


void testCutoffCalculation(){

  TFile *f = new TFile(fname);
  TH2F *hfull = (TH2F*)f->Get(hname);

  TH1D *hslice = hfull->ProjectionY("_py",sliceIndex, sliceIndex);
  hslice->SetMarkerStyle(20);
  hslice->SetMarkerSize(1);

  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  c1->Draw();

  hslice->Draw("PE");
  c1->Update();

  methodEffCurve(hslice);

  TH1D *hslice2 = (TH1D*)hslice->Clone("hslice2");
  hslice2->Rebin(4);
  methodToy(hslice2);
  

}

void methodEffCurve(TH1D *hist){

  int   nbins = hist->GetNbinsX();
  float xlow = hist->GetXaxis()->GetBinLowEdge(1);
  float xhigh = hist->GetXaxis()->GetBinUpEdge(nbins);
  
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
  printf("total= %f   totalNonZero= %f\n", total, totalNonZero);

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
  float xcutoffLow = xlow;
  float xcutoffHigh = xlow;
  float xcutoff = xlow;
  int npoints = grEff->GetN();
  Double_t *xArray      = grEff->GetX();
  Double_t *effArray    = grEff->GetY();
  Double_t *effErrArray = grEff->GetEY();
  // Central values
  for(int i=1; i < npoints; i++){ // Do not start with point zero
    if( effArray[i] > effCutoff){
      // Found the first efficiency value above threshold
      // Interpolate and find the cutoff value 
      interpolate( xArray[i-1], xArray[i], effArray[i-1],effArray[i],
		   xcutoff, effCutoff);
      break;
    }
  }
  // Error up
  for(int i=1; i < npoints; i++){ // Do not start with point zero
    if( effArray[i] + effErrArray[i] > effCutoff){
      // Found the first efficiency value above threshold
      // Interpolate and find the cutoff value
      interpolate( xArray[i-1], xArray[i], 
		   effArray[i-1] + effErrArray[i-1],
		   effArray[i] + effErrArray[i],
		   xcutoffLow, effCutoff);
      break;
    }
  }
  // Error down
  for(int i=1; i < npoints; i++){ // Do not start with point zero
    if( effArray[i] - effErrArray[i] > effCutoff){
      // Found the first efficiency value above threshold
      // Interpolate and find the cutoff value
      interpolate( xArray[i-1], xArray[i], 
		   effArray[i-1] - effErrArray[i-1],
		   effArray[i] - effErrArray[i],
		   xcutoffHigh, effCutoff);
      break;
    }
  }
  float xcutoffErrPlus = xcutoffHigh-xcutoff;
  float xcutoffErrMinus = xcutoff-xcutoffLow;
  // A copy to be used later for drawing
  TGraphErrors *grEffClone = (TGraphErrors*)grEff->Clone("grEffClone");


  printf("Method eff curve: Cutoff is at %f + %f - %f\n", xcutoff, 
	 xcutoffErrPlus,
	 xcutoffErrMinus);


  //
  // Draw the result
  //
  TCanvas *c2 = new TCanvas("c2","c2",100,10,600,600);
  gStyle->SetOptStat(0);
  c2->Draw();

  TH2F *dummy = new TH2F("dummy","",100, xlow, xhigh, 100, 0, 1.2);
  dummy->GetXaxis()->SetTitle("cut value");
  dummy->GetYaxis()->SetTitle("efficiency");
  dummy->GetYaxis()->SetTitleOffset(1.4);
  // dummy->GetXaxis()->SetRangeUser(xcutoff -10*xcutoffErrMinus,
  // 				  xcutoff +10*xcutoffErrPlus);
  dummy->Draw();

  grEff->Draw("E3,same");
  grEff->SetLineWidth(0); // suppress drawing error bars
  grEff->DrawClone("P,same");
  // Draw lines to show the cutoff, etc
  TLine *hline = new TLine(xlow, effCutoff, xhigh, effCutoff);
  hline->Draw("same");
  TLine *vline1 = new TLine(xcutoffLow,  0, xcutoffLow,  1.2);
  TLine *vline2 = new TLine(xcutoff,     0, xcutoff,     1.2);
  TLine *vline3 = new TLine(xcutoffHigh, 0, xcutoffHigh, 1.2);
  vline1->Draw("same");
  vline2->Draw("same");
  vline3->Draw("same");
  c2->Update();

  TCanvas *c3 = new TCanvas("c3","c3",200,10,600,600);
  c3->Draw();

  float x1 = xcutoff -5*xcutoffErrMinus;
  float x2 =  xcutoff +5*xcutoffErrPlus;
  float y1 = effCutoff-0.01;
  float y2 = effCutoff+0.01;
  TH2F *dummy2 = (TH2F*)dummy->Clone("dummy2");
  dummy2->Draw();
  dummy2->GetXaxis()->SetRangeUser(x1,x2);
  dummy2->GetYaxis()->SetRangeUser(y1,y2);
  dummy2->Draw();

  grEffClone->Draw("E3,same");
  grEffClone->SetLineWidth(0);
  grEffClone->Draw("P,same");

  // Draw lines to show the cutoff, etc
  TLine *hline2 = new TLine(x1, effCutoff, x2, effCutoff);
  TLine *vline21 = new TLine(xcutoffLow,  y1, xcutoffLow,  y2);
  TLine *vline22 = new TLine(xcutoff,     y1, xcutoff,     y2);
  TLine *vline23 = new TLine(xcutoffHigh, y1, xcutoffHigh, y2);
  hline2->Draw("same");
  vline21->Draw("same");
  vline22->Draw("same");
  vline23->Draw("same");
  c3->Update();

}

void interpolate( float x1, float x2, float y1, float y2, float &x, float y){

  x = x1 + (y-y1)*(x2-x1)/(y2-y1);

  // printf("x1= %f  x2= %f  y1= %f  y2= %f  x= %f  y= %f\n",
  // 	 x1, x2, y1, y2, x, y);
  return;
}

void methodToy(TH1D* hist){

  TCanvas *c4 = new TCanvas("c4","c4",500,10,600,600);
  float cutoffToy, cutoffErrToy;

  computeCutoffAndError(hist, cutoffToy, cutoffErrToy, c4);
  printf("Toy method: cutoff = %f +- %f\n", cutoffToy, cutoffErrToy);
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

void computeCutoffAndError(TH1D *hist, float &cutoff, float &cutoffErr, TCanvas *canv){


  float total = 0;
  computeCutoff(hist, total, cutoff);
  
  const float largeCutoffError = 10;
  if( total < 500 ) {
    // not enough data for reliable cutoff computations
    cutoff = 999;
    cutoffErr = largeCutoffError;
  }else{
    // sufficient statistics, proceed with error calculations
    
    // Compute the error
    TF1 *isofunc = new TF1("isofunc","isoShape2",-0.1, 6, 7);
    isofunc->SetParLimits(0, 0, 1e8);
    isofunc->SetParLimits(1, 0, 1e8);
    isofunc->SetParLimits(2, 0.1, 10);
    isofunc->SetParLimits(3, 0.0, 10);
    isofunc->SetParLimits(4, 0, 1e8);
    isofunc->SetParLimits(5, 0.2, 1.5);
    isofunc->SetParLimits(6, 0, 2);
    isofunc->SetParameters(10,10,2,1, 10, 0.9, 0.1);
    
    isofunc->SetNpx(1000);
    isofunc->SetLineColor(kRed);
    hist->Fit("isofunc","R");
    
    isofunc->SetRange( hist->GetXaxis()->GetBinLowEdge(1), 100);
		       // hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()));

    const int nPseudoExp = 100;
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
  
  canv->cd();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1);
  hist->Draw("pe");
  canv->Update();
  // float xxx;
  // cin >> xxx;
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

Double_t isoShape2(Double_t *x, Double_t *par) 
{

  Double_t xx = x[0];
  Double_t boxNorm        = par[0];
  Double_t landauNorm     = par[1];
  Double_t landauLocation = par[2];
  Double_t landauScale    = par[3];
  Double_t gausNorm       = par[4];
  Double_t gausMean       = par[5];
  Double_t gausSigma      = par[6];

  Double_t result = boxNorm*box(xx) 
    + landauNorm * TMath::Landau(xx, landauLocation, landauScale)
    + gausNorm * TMath::Gaus(xx, gausMean, gausSigma);

  return result;
}
