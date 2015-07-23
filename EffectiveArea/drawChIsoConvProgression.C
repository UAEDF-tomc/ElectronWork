#include "TGraphErrors.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
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

// Signal sample: DYToLL
const TString fileNameSignal = 
  "/afs/cern.ch/user/r/rkamalie/workspace/public/DY_Spring15_Asympt50ns_24june2015.root";
// Directory and tree name:
const TString treeName = "ntupler/ElectronTree";

const bool verbose = false;
const bool smallEventCount = false; // DEBUG

float ptCut = 20;

const int nRhoBins = 10;
const float rhoBinLimits[nRhoBins+1] = {
  0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};

// 
// Forward declarations
//
float getCutoff(TH1F *hist, float frac);

//
// Main program
//

void drawChIsoConvProgression(){

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 

  // General settings
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);

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

  float meanp[nRhoBins];
  float meanpErr[nRhoBins];
  TH1F *hIsoFromPU[nRhoBins];
  TH1F *hIso[nRhoBins];
  for(int i=0; i<nRhoBins; i++){
    meanp[i] = 0;
    meanpErr[i] = 0;
    TString hname = TString::Format("hIsoFromPU_%d", i);
    hIsoFromPU[i] = new TH1F(hname,hname,50,0,20);
    hIsoFromPU[i]->GetXaxis()->SetTitle("isolation, GeV");
    TString hname2 = TString::Format("hIso_%d", i);
    hIso[i] = new TH1F(hname2,hname2,100,0,20);
    hIso[i]->GetXaxis()->SetTitle("isolation, GeV");
  }

  int nEle; // the number of reconstructed electrons in the event
  float rho;
  std::vector <float> *elePt = 0;         // electron PT
  std::vector <float> *eleEtaSC = 0;      // supercluser eta
  std::vector <float> *isoChargedHadrons = 0;
  std::vector <float> *isoChargedFromPU = 0;
  std::vector <int> *isTrue = 0;

  // Declare branches
  TBranch *b_nEle = 0;
  TBranch *b_rho = 0;
  TBranch *b_elePt = 0;
  TBranch *b_eleEtaSC = 0;
  TBranch *b_isoChargedHadrons = 0;
  TBranch *b_isoChargedFromPU = 0;
  TBranch *b_isTrue = 0;

  treeSignal->SetBranchAddress("nEle", &nEle, &b_nEle);
  treeSignal->SetBranchAddress("rho", &rho, &b_rho);
  treeSignal->SetBranchAddress("pt", &elePt, &b_elePt);
  treeSignal->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  treeSignal->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
  treeSignal->SetBranchAddress("isoChargedFromPU", &isoChargedFromPU, &b_isoChargedFromPU);
  treeSignal->SetBranchAddress("isTrue",    &isTrue,    &b_isTrue);

  // 
  // Loop over events pass 1
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
    b_rho->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nEle);

    // Find rho bin 
    int rhoBin = nRhoBins-1; // set to max
    for(int ibin=0; ibin<nRhoBins; ibin++){
      if( rho >= rhoBinLimits[ibin] && rho<rhoBinLimits[ibin+1] ){
	rhoBin = ibin;
	break;
      }
    }	
  
    // Get data for all electrons in this event, only vars of interest
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_isoChargedHadrons->GetEntry(tentry);
    b_isoChargedFromPU->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);

    // Nested loops over the electrons
    for(int iele = 0; iele < nEle; iele++){

      // Check kinematics:
      if( !(elePt->at(iele) > ptCut) )
	continue;

      // Check truth match
      if( isTrue->at(iele) != 1 ) continue;

      // Keep only electrons of the central barrel for this study
      if( abs(eleEtaSC->at(iele))>0.8 ) continue;

      hIsoFromPU[rhoBin]->Fill( isoChargedFromPU->at(iele) );

    } // end loop over electrons

  } // end loop over events

  //
  // Loop over rho bins and determine all the Poisson parameters
  //
  // Define the Poisson function
  TF1 *fp[nRhoBins];
  float norm = 1.0;
  float mean = 1;
  for(int i=0; i<nRhoBins; i++){
    TString fname = TString::Format("fp_%d",i);
    fp[i] = new TF1(fname,"[0]*TMath::Poisson(x,[1])",0,20);
    fp[i]->SetLineColor(kRed);
    fp[i]->SetParameters( norm, mean);
  }

  float centerRho[nRhoBins];
  float deltaRho[nRhoBins];
  for(int i=0; i<nRhoBins; i++){
    TString cname = TString::Format("c1_%d",i);
    TCanvas *c1 = new TCanvas(cname,cname,10,10,600,600);
    hIsoFromPU[i]->Fit(fp[i]->GetName(),"R");
    c1->Update();
    meanp[i] = fp[i]->GetParameter(1);
    meanpErr[i] = fp[i]->GetParError(1);
    centerRho[i] = (rhoBinLimits[i+1] + rhoBinLimits[i])/2.0;
    deltaRho[i] = (rhoBinLimits[i+1] - rhoBinLimits[i])/2.0;
  }

  TCanvas *crho = new TCanvas("crho","crho",20,20,600,600);
  TGraphErrors *grMean = new TGraphErrors(nRhoBins, centerRho, meanp, deltaRho, meanpErr);
  grMean->SetTitle("");
  grMean->SetMarkerSize(1);
  grMean->SetMarkerStyle(20);
  grMean->Draw("APE");
  grMean->GetXaxis()->SetTitle("rho");
  grMean->GetYaxis()->SetTitle("Poisson mean from isolation fit");
  grMean->Draw("APE");

  // 
  // Loop over events pass 2
  //
  printf("\nPass 2\n");
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = treeSignal->LoadTree(ievent);
    
    // Load the value of the number of the electrons in the event    
    b_nEle->GetEntry(tentry);
    b_rho->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nEle);
  
    // Get data for all electrons in this event, only vars of interest
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_isoChargedHadrons->GetEntry(tentry);
    b_isoChargedFromPU->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);

    // Nested loops over the electrons
    for(int iele = 0; iele < nEle; iele++){

      // Check kinematics:
      if( !(elePt->at(iele) > ptCut) )
	continue;

      // Check truth match
      if( isTrue->at(iele) != 1 ) continue;

      // Keep only electrons of the central barrel for this study
      if( abs(eleEtaSC->at(iele))>0.8 ) continue;

      // This is just a simulation so for each electron regardless of true rho
      // we simulate all sorts of rhos
      for(int ibin=0; ibin<nRhoBins; ibin++){
	hIso[ibin]->Fill( isoChargedHadrons->at(iele) 
			 + fp[ibin]->GetRandom());
      }	

    } // end loop over electrons

  } // end loop over events

  // Find 95% cutoffs
  const float frac = 0.95;
  float cutoffs[nRhoBins];
  float means[nRhoBins];
  TCanvas *ccut = new TCanvas("ccut","ccut",10,10,600,600);
  hIso[0]->Draw("hist");
  TLegend *leg10 = new TLegend(0.30,0.6, 0.88, 0.88);
  leg10->SetFillStyle(0);
  leg10->SetBorderSize(0);
  for(int i=0; i<nRhoBins; i++){
    cutoffs[i] = getCutoff( hIso[i], frac);
    means[i] = hIso[i]->GetMean();
    hIso[i]->Scale(1.0/hIso[i]->GetSumOfWeights());
    hIso[i]->SetLineWidth(2);
    hIso[i]->SetLineColor(kOrange+10-i);
    hIso[i]->Draw("hist,same");
    TString legString = TString::Format("ISO@PU for %2.0f<#rho<%2.0f", 
					rhoBinLimits[i],rhoBinLimits[i+1]);
    if(i == nRhoBins-1)
      legString = TString::Format("ISO@PU for #rho>%2.0f", rhoBinLimits[i]);
    leg10->AddEntry(hIso[i], legString,"l");
  }
  leg10->Draw("same");

  TCanvas *ccut2 = new TCanvas("ccut2","ccut2",10,10,600,600);
  TH2F *dummy = new TH2F("dummy","",100,0,50,100,0,20);
  dummy->GetXaxis()->SetTitle("rho");
  dummy->GetYaxis()->SetTitle("isolation");
  dummy->Draw();

  TGraph *gCutoffs = new TGraph(nRhoBins, centerRho, cutoffs);
  gCutoffs->SetMarkerSize(1);
  gCutoffs->SetMarkerStyle(20);
  gCutoffs->Draw("P");

  TGraph *gMeans = new TGraph(nRhoBins, centerRho, means);
  gMeans->SetMarkerSize(1);
  gMeans->SetMarkerStyle(20);
  gMeans->SetMarkerColor(kBlue);
  gMeans->Draw("P");

  TLegend *leg1 = new TLegend(0.15, 0.7, 0.60, 0.88);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry(gMeans,"mean of ISO","p");
  leg1->AddEntry(gCutoffs,"95% cutoff of ISO","p");
  leg1->Draw("same");

}

float getCutoff(TH1F *hist, float frac){

  int n = hist->GetNbinsX();
  // Find the sum counting the overflow bin
  float total = 0;
  for(int i=1; i<=n+1; i++){
    total += hist->GetBinContent(i);
  }

  // Now find the bin of the cutoff
  float cutoff = 999;
  float maxCount = frac * total;
  float count = 0;
  for(int iBin = 1; iBin<=n; iBin++){ 
    count += hist->GetBinContent(iBin);
    if(count < maxCount ){
      cutoff =  hist->GetXaxis()->GetBinUpEdge(iBin);
    }else
      break;
  }

  return cutoff;
}

