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
const TString finName    = "hVsRho_EE_DY.root";
const TString etaLabel   = "Endcap";
const TString hinName    = "hVsRho";
const TString xAxisTitle = "rho";
const TString yAxisTitle = "HCAL energy [GeV]";

// Contours to work with
const float cutoffFraction = 0.90;

// Range of the variable to fit
const float fitMin = 5;
const float fitMax = 30;

// Forward declarations
void findCutoff(TH1D *hist, float &val, float &valErrPos, float &valErrNeg);
void interpolate( float x1, float x2, float y1, float y2, float &x, float y);
// Plot the main 2D histogram analyzed in this script
void plot2D(TH2F *hist);
// Plot the projection onto X-axis, the histogram of the
// varilable depdence on which we study
void plotVariable(TH2F *hist);

// Main method
void hoeRhoFindDependence(){

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

  // Loop over all slices along X and determine the cut-off
  TGraphAsymmErrors *graph = new TGraphAsymmErrors();
  for(int ibin=1; ibin<= hist2D->GetNbinsX(); ibin++){
    // One bin slice:
    TH1D *hist1D = hist2D->ProjectionY("_py", ibin, ibin);
    // Find the cutoff and the error on it
    float val, valErrPos, valErrNeg;
    findCutoff(hist1D, val, valErrPos, valErrNeg);
    // Save into the graph. Note that calling SetPoint(N) creates
    // a new point in the TGraph because the last point index is (N-1)
    int nCurrentPoints = graph->GetN();
    float xCurrent = hist2D->GetXaxis()->GetBinCenter(ibin);
    float xWidth   = hist2D->GetXaxis()->GetBinWidth(ibin);
    graph->SetPoint(nCurrentPoints, xCurrent, val);
    graph->SetPointError(nCurrentPoints, xWidth/2.0, xWidth/2.0, 
			 valErrPos, valErrNeg);
  }

  //
  // Fit the dependence
  //
  TF1 *func = new TF1("func","[0]+[1]*x",fitMin, fitMax);
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
    double yPlusErrTmp = ytmp + graph->GetErrorYhigh(jbin);
    if( yhigh < yPlusErrTmp)
      yhigh = yPlusErrTmp;
  }
  // Increase vertical range by another 50%
  yhigh *= 1.5;
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

  TString slopeString = TString::Format("slope %.2e #pm %.2e",B, dB);
  TLatex *lat3 = new TLatex(0.15, 0.7, slopeString);
  lat3->SetNDC(kTRUE);
  lat3->SetTextSize(0.03);
  lat3->Draw();

  // Other plots
  plot2D(hist2D);
  plotVariable(hist2D);

}

void findCutoff(TH1D *hist, float &val, float &valErrPos, float &valErrNeg){

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

