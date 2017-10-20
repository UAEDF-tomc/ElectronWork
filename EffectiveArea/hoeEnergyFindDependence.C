#include "TMath.h"
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
const TString finName    = "hoeVsE_EE.root";
const TString etaLabel   = "Endcap";
const TString hinName    = "hoeVsE";
const TString xAxisTitle = "E_{SC}";
const TString yAxisTitle = "H/E";

// Contours to work with
const float cutoffFraction = 0.95;

// Range of the variable to fit
const float fitMin = 20;
const float fitMax = 1000;

// Binning in the variable of interest Var. We assume that
// the input 2D histogram is finely binned in the variable Var.
// For energy, assuming that the incoming hist has limits
// from 0 to 1000 GeV with 1000 original bins total:
//    0 - 100    2 GeV bins, count 50
//  100 - 200    5 GeV bins, count 20
//  200 - 400   10 GeV bins, count 20
//  400 - 600   50 GeV bins, count 4
//  600 -1000  100 GeV bins, count 4
bool binningReady = false;
// Total number of expected original bins:
const int nVarOriginalBins = 1000; 
const int nVarBins = 98;
int varBinIndexLow[nVarBins];
int varBinIndexHigh[nVarBins];
double varBinLimits[nVarBins+1];
double varBinWeightedCenters[nVarBins];
double varBinRangePos[nVarBins+1];
double varBinRangeNeg[nVarBins+1];

// Forward declarations
void setupBinning(TH2F *hist);
void findCutoff(TH1D *hist, float &val, float &valErrPos, float &valErrNeg);
void interpolate( float x1, float x2, float y1, float y2, float &x, float y);
// Plot the main 2D histogram analyzed in this script
void plot2D(TH2F *hist);
// Plot the projection onto X-axis, the histogram of the
// varilable depdence on which we study
void plotVariable(TH2F *hist);

// Main method
void hoeEnergyFindDependence(){

  // Read in the 2D histogram to analyze
  TFile *fin = new TFile(finName);
  if( !fin ) {
    printf("File %s not found\n", finName.Data());
    assert(0);
  }
  TH2F *hist2D = (TH2F*)fin->Get(hinName);
  if( !hist2D ){
    printf("Histogram %s not found in file %s\n", hinName.Data(), finName.Data());
    assert(0);
  }

  setupBinning(hist2D);

  // Loop over all slices along X and determine the cut-off
  TGraphAsymmErrors *graph = new TGraphAsymmErrors();
  for(int ibin=1; ibin<= nVarBins; ibin++){

    if( !binningReady ){
      printf("Not properly set up\n");
      assert(0);
    }
    TString pname = TString::Format("_py%d\n",ibin);
    TH1D *hist1D = hist2D->ProjectionY(pname, varBinIndexLow[ibin-1], 
    				       varBinIndexHigh[ibin-1]);
    // Find the cutoff and the error on it
    float val, valErrPos, valErrNeg;
    findCutoff(hist1D, val, valErrPos, valErrNeg);
    // Save into the graph. Note that calling SetPoint(N) creates
    // a new point in the TGraph because the last point index is (N-1)
    int nCurrentPoints = graph->GetN();
    float xCurrent    = varBinWeightedCenters[nCurrentPoints];
    float xRangeLow   = varBinRangeNeg[nCurrentPoints];
    float xRangeHigh  = varBinRangePos[nCurrentPoints];
    graph->SetPoint(nCurrentPoints, xCurrent, val);
    graph->SetPointError(nCurrentPoints, xRangeLow, xRangeHigh, 
    			 valErrPos, valErrNeg);
    delete hist1D;
  }

  //
  // Fit the dependence
  //
  TF1 *func = new TF1("func","[0]+[1]/x",fitMin, fitMax);
  graph->Fit("func","R");
  float A  = func->GetParameter(0);
  float dA = func->GetParError(0);
  float B  = func->GetParameter(1);
  float dB = func->GetParError(1);

  //
  // Draw the dependence
  //

  TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);
  c1->cd();
  gStyle->SetOptStat(0);

  double xlow, xhigh, ylow, yhigh;
  xlow = 0;
  ylow = 0;
  xhigh = hist2D->GetXaxis()->GetBinUpEdge(hist2D->GetNbinsX());

  // Find Y range
  double xtmp, ytmp;
  graph->GetPoint(0, xtmp, yhigh);
  for(int jbin = 0; jbin < graph->GetN(); jbin++){
    graph->GetPoint(jbin, xtmp, ytmp);
    float ytmpErr = graph->GetErrorYhigh(jbin);
    double yPlusErrTmp = ytmp + ytmpErr;
    if( yhigh < yPlusErrTmp)
      yhigh = yPlusErrTmp;
  }
  // Increase vertical range by another 50%
  yhigh *= 1.5;
  printf("%f   %f   %f   %f\n", xlow, xhigh, ylow, yhigh);
  TH2F *dummy = new TH2F("dummy","",100, xlow, xhigh, 100, ylow, yhigh);
  dummy->GetXaxis()->SetTitle(xAxisTitle);
  dummy->GetYaxis()->SetTitle(yAxisTitle);
  dummy->GetYaxis()->SetTitleOffset(1.4);
  dummy->Draw();

  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.0);
  graph->Draw("P,same");

  // Labels
  TLatex *lat1 = new TLatex(0.15, 0.8, etaLabel);
  lat1->SetNDC(kTRUE);
  lat1->Draw();

  TString fracString = TString::Format("contours at %.0f %%", 100*cutoffFraction);
  TLatex *lat2 = new TLatex(0.15, 0.75, fracString);
  lat2->SetNDC(kTRUE);
  lat2->SetTextSize(0.03);
  lat2->Draw();

  TString slopeString = TString::Format("f(E)=A+B/E, B=%.2e #pm %.2e",B, dB);
  TLatex *lat3 = new TLatex(0.15, 0.7, slopeString);
  lat3->SetNDC(kTRUE);
  lat3->SetTextSize(0.03);
  lat3->Draw();

  // Other plots
  plot2D(hist2D);
  plotVariable(hist2D);

}

void setupBinning(TH2F *hist){

  if( nVarOriginalBins != hist->GetNbinsX() ){
    printf("Unexpected number of bins, exiting\n");
    assert(0);
  }

  // Original bins range 1 - 100
  int maxBin1 = 50;
  int step1 = 2;
  int startBin = 0;
  int startBinOrig = 0;
  for(int i=0; i<maxBin1; i += 1){
    varBinIndexLow[startBin + i]  = startBinOrig + step1*i + 1;
    varBinIndexHigh[startBin + i] = startBinOrig + step1*i + step1;
  }
  // Original bins range 100 - 200
  int maxBin2 = 20;
  int step2 = 5;
  startBin = maxBin1;
  startBinOrig = maxBin1*step1;
  for(int i=0; i<maxBin2; i += 1){
    varBinIndexLow [startBin+i]  = startBinOrig + step2*i + 1;
    varBinIndexHigh[startBin+i]  = startBinOrig + step2*i + step2;
  }
  // Original bins range 200 - 400
  int maxBin3 = 20;
  int step3   = 10;
  startBin = maxBin1 + maxBin2;
  startBinOrig = maxBin1*step1 + maxBin2*step2;
  for(int i=0; i<maxBin3; i += 1){
    varBinIndexLow [startBin+i]  = startBinOrig + step3*i + 1;
    varBinIndexHigh[startBin+i]  = startBinOrig+ step3*i + step3;
  }
  // Original bins range 400 - 600
  int maxBin4 = 4;
  int step4   = 50;
  startBin = maxBin1 + maxBin2 + maxBin3;
  startBinOrig = maxBin1*step1 + maxBin2*step2 + maxBin3*step3;
  for(int i=0; i<maxBin4; i += 1){
    varBinIndexLow [startBin+i]  = startBinOrig + step4*i + 1;
    varBinIndexHigh[startBin+i]  = startBinOrig + step4*i + step4;
  }
  // Original bins range 600 - 1000
  int maxBin5 = 4;
  int step5   = 100;
  startBin = maxBin1 + maxBin2 + maxBin3 + maxBin4;
  startBinOrig = maxBin1*step1 + maxBin2*step2 
    + maxBin3*step3 + maxBin4*step4;
  for(int i=0; i<maxBin5; i += 1){
    varBinIndexLow [startBin+i]  = startBinOrig + step5*i + 1;
    varBinIndexHigh[startBin+i]  = startBinOrig + step5*i + step5;
  }

  TAxis *varAxis = hist->GetXaxis();
  varBinLimits[0] = varAxis->GetBinLowEdge( varBinIndexLow[0] );
  for(int i=0; i<nVarBins; i++){
    varBinLimits[i+1] = varAxis->GetBinUpEdge( varBinIndexHigh[i] );
  }

  // Find weighted centers
  TH1D *varSpectrum = hist->ProjectionX();
  for(int iNew=0; iNew<nVarBins; iNew++){

    // Create a histogram that preserves only the bin range
    // of the original histogram that correpsonds to a given
    // bin in the new variable bin histogram
    TH1D *histTmp = (TH1D*)varSpectrum->Clone("histTmp");
    for(int iOrig=0; iOrig<= histTmp->GetNbinsX()+1; iOrig++){
      if(iOrig < varBinIndexLow[iNew] || iOrig > varBinIndexHigh[iNew] ){
	histTmp->SetBinContent(iOrig, 0);
	histTmp->SetBinError(iOrig, 0);
      }
    }

    // TH2F *histTmp = (TH2F*)hist->Clone("histTmp");
    // for(int iOrig=0; iOrig<= histTmp->GetNbinsX()+1; iOrig++){
    //   for(int jOrig=0; jOrig<=histTmp->GetNbinsY()+1; jOrig++){
    // 	if(iOrig < varBinIndexLow[iNew] || iOrig > varBinIndexHigh[iNew] ){
    // 	  histTmp->SetBinContent(iOrig, jOrig, 0);
    // 	  histTmp->SetBinError(iOrig, jOrig, 0);
    // 	}
    //   }
    // }

    // Find the mean of the projection of this histogram on X
    double mean = histTmp->GetMean();
    varBinWeightedCenters[iNew] = mean;
    varBinRangeNeg[iNew] = mean - varBinLimits[iNew];
    varBinRangePos[iNew] = varBinLimits[iNew+1] - mean;
    delete histTmp;
  }

  binningReady = true;
  return;
}


void findCutoff(TH1D *hist, float &val, float &valErrPos, float &valErrNeg){

  val = 0;
  valErrPos = 0;
  valErrNeg = 0;

  // If there are very few points, give up
  if( hist->GetEntries() < 10 ) 
    return;

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
  // printf("total= %f   totalNonZero= %f\n", total, totalNonZero);

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
    if( effArray[i] > cutoffFraction){
      // Found the first efficiency value above threshold
      // Interpolate and find the cutoff value 
      interpolate( xArray[i-1], xArray[i], effArray[i-1],effArray[i],
		   xcutoff, cutoffFraction);
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
		   xcutoffLow, cutoffFraction);
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
		   xcutoffHigh, cutoffFraction);
      break;
    }
  }
  float xcutoffErrPlus = xcutoffHigh-xcutoff;
  float xcutoffErrMinus = xcutoff-xcutoffLow;

  val       = xcutoff;
  valErrPos = xcutoffErrPlus;
  valErrNeg = xcutoffErrMinus;

  printf("  %f    %f     %f\n", val, valErrPos, valErrNeg);

  return;
}

void interpolate( float x1, float x2, float y1, float y2, float &x, float y){

  x = x1 + (y-y1)*(x2-x1)/(y2-y1);

  // printf("x1= %f  x2= %f  y1= %f  y2= %f  x= %f  y= %f\n",
  // 	 x1, x2, y1, y2, x, y);
  return;
}

void plot2D(TH2F *hist){

  TCanvas *c2 = new TCanvas("c2","c2",100,10, 800, 800);
  c2->cd();
  gStyle->SetOptStat(0);
  
  hist->GetXaxis()->SetTitle(xAxisTitle);
  hist->GetYaxis()->SetTitle(yAxisTitle);
  hist->GetYaxis()->SetTitleOffset(1.4);

  hist->Draw("colz");

  TLatex *lat = new TLatex(0.7, 0.8, etaLabel);
  lat->SetNDC(kTRUE);
  lat->Draw();

}

void plotVariable(TH2F *hist){

  TCanvas *c3 = new TCanvas("c3","c3",200,10, 800, 800);
  c3->cd();
  gStyle->SetOptStat(0);
  
  TH1D *hist1D = hist->ProjectionX();
  hist1D->GetXaxis()->SetTitle(xAxisTitle);
  
  hist1D->Draw("hist");

}

