#include "TGraph.h"
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

// Background sample: TT
const TString fileNameBg = 
  "/afs/cern.ch/user/r/rkamalie/workspace/public/TT_Spring15_Asympt50ns_24june2015.root";

const bool verbose = false;
const bool smallEventCount = false; // DEBUG

float ptCut = 20;

const float eaDefault = 0.0973;
const float eaCustom  = 0.1682;

// 
// Forward declarations
//
float getCutoff(TH1F *hist, float frac);

//
// Main program
//

void drawIsolationROC(){

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

  TFile *fileBg     = new TFile(fileNameBg);
  if( !fileBg ){
    printf("Failed to open the input files, check\n   %s\n", 
	   fileNameBg.Data());
    assert(0);
  }
  TTree *treeBg     = (TTree*)fileBg->Get(treeName);
  if( !treeBg ){
    printf("Failed to find the tree %s\n", treeName.Data() );
    assert(0);
  }

  float rho;
  int nEle; // the number of reconstructed electrons in the event
  std::vector <float> *elePt = 0;         // electron PT
  std::vector <float> *eleEtaSC = 0;      // supercluser eta
  std::vector <float> *isoChargedHadrons = 0;
  std::vector <float> *isoNeutralHadrons = 0;
  std::vector <float> *isoPhotons = 0;
  std::vector <int> *isTrue = 0;

  // Declare branches
  TBranch *b_rho = 0;
  TBranch *b_nEle = 0;
  TBranch *b_elePt = 0;
  TBranch *b_eleEtaSC = 0;
  TBranch *b_isoChargedHadrons = 0;
  TBranch *b_isoNeutralHadrons = 0;
  TBranch *b_isoPhotons = 0;
  TBranch *b_isTrue = 0;

  // Containers for efficiency related numbers
  const int nCuts = 10000;
  float totalSig = 0;
  float totalBg = 0;
  float passDefaultSig[nCuts];
  float passDefaultBg[nCuts];
  float passCustomSig[nCuts];
  float passCustomBg[nCuts];
  float cutVal[nCuts];
  const float cutMin = 0.00;
  const float cutMax = 0.40;
  const float delta = (cutMax-cutMin)/(nCuts-1);
  for(int i=0; i<nCuts; i++){
    passDefaultSig[i] = 0;
    passDefaultBg[i] = 0;
    passCustomSig[i] = 0;
    passCustomBg[i] = 0;
    cutVal[i] = cutMin + (i+1)*delta;
  }

  //========================================================
  // Work with signal 
  //========================================================

  treeSignal->SetBranchAddress("rho", &rho, &b_rho);
  treeSignal->SetBranchAddress("nEle", &nEle, &b_nEle);
  treeSignal->SetBranchAddress("pt", &elePt, &b_elePt);
  treeSignal->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  treeSignal->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
  treeSignal->SetBranchAddress("isoNeutralHadrons", &isoNeutralHadrons, &b_isoNeutralHadrons);
  treeSignal->SetBranchAddress("isoPhotons"       , &isoPhotons, &b_isoPhotons);
  treeSignal->SetBranchAddress("isTrue",    &isTrue,    &b_isTrue);

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
    b_rho->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nEle);
    
    // Get data for all electrons in this event, only vars of interest
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_isoChargedHadrons->GetEntry(tentry);
    b_isoNeutralHadrons->GetEntry(tentry);
    b_isoPhotons->GetEntry(tentry);
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
      
      float isoCh = isoChargedHadrons->at(iele);
      float isoNh = isoNeutralHadrons->at(iele);
      float isoPh = isoPhotons->at(iele);
      float pt = elePt->at(iele);
      //
      float isoDefault = (isoCh + std::max(0.0, (double)(isoNh + isoPh - rho*eaDefault)))/pt;
      float isoCustom = (isoCh + std::max(0.0, (double)(isoNh + isoPh - rho*eaCustom)))/pt;

      totalSig++;
      for(int iCut = 0; iCut<nCuts; iCut++){
	if( isoDefault < cutVal[iCut] )
	  passDefaultSig[iCut] += 1;
	if( isoCustom < cutVal[iCut] )
	  passCustomSig[iCut] += 1;
      }

    } // end loop over electrons

  } // end loop over events

  //========================================================
  // Work with background
  //========================================================

  treeBg->SetBranchAddress("rho", &rho, &b_rho);
  treeBg->SetBranchAddress("nEle", &nEle, &b_nEle);
  treeBg->SetBranchAddress("pt", &elePt, &b_elePt);
  treeBg->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  treeBg->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
  treeBg->SetBranchAddress("isoNeutralHadrons", &isoNeutralHadrons, &b_isoNeutralHadrons);
  treeBg->SetBranchAddress("isoPhotons"       , &isoPhotons, &b_isoPhotons);
  treeBg->SetBranchAddress("isTrue",    &isTrue,    &b_isTrue);

  // 
  // Loop over events
  //
  maxEvents = treeBg->GetEntries();
  if( smallEventCount )
    maxEvents = 1000000;
  if(verbose)
    printf("Start loop over events, total events = %lld\n", 
	   treeBg->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = treeBg->LoadTree(ievent);
    
    // Load the value of the number of the electrons in the event    
    b_nEle->GetEntry(tentry);
    b_rho->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nEle);
    
    // Get data for all electrons in this event, only vars of interest
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_isoChargedHadrons->GetEntry(tentry);
    b_isoNeutralHadrons->GetEntry(tentry);
    b_isoPhotons->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);

    // Nested loops over the electrons
    for(int iele = 0; iele < nEle; iele++){

      // Check kinematics:
      if( !(elePt->at(iele) > ptCut) )
	continue;

      // Check truth match
      if( !( isTrue->at(iele) ==0 || isTrue->at(iele)==3) ) continue;

      // Keep only electrons of the central barrel for this study
      if( abs(eleEtaSC->at(iele))>0.8 ) continue;
      
      float isoCh = isoChargedHadrons->at(iele);
      float isoNh = isoNeutralHadrons->at(iele);
      float isoPh = isoPhotons->at(iele);
      float pt = elePt->at(iele);
      //
      float isoDefault = (isoCh + std::max(0.0, (double)(isoNh + isoPh - rho*eaDefault)))/pt;
      float isoCustom = (isoCh + std::max(0.0, (double)(isoNh + isoPh - rho*eaCustom)))/pt;

      totalBg++;
      for(int iCut = 0; iCut<nCuts; iCut++){
	if( isoDefault < cutVal[iCut] )
	  passDefaultBg[iCut] += 1;
	if( isoCustom < cutVal[iCut] )
	  passCustomBg[iCut] += 1;
      }

    } // end loop over electrons

  } // end loop over events

  //========================================================
  // compute and draw ROCs
  //========================================================
  float sigEffDefault[nCuts];
  float sigEffCustom[nCuts];
  float bgRejDefault[nCuts];
  float bgRejCustom[nCuts];

  for(int i=0; i<nCuts; i++){
    sigEffDefault[i] = passDefaultSig[i]/totalSig;
    bgRejDefault[i]  = 1.0 - passDefaultBg[i]/totalBg;

    sigEffCustom[i] = passCustomSig[i]/totalSig;
    bgRejCustom[i]  = 1.0 - passCustomBg[i]/totalBg;
    if(i%1000==0)
      printf(" i=%d   cut=%f   sigEff=%f   bgRej=%f\n", i, cutVal[i], sigEffCustom[i], bgRejCustom[i]);
  }
  TGraph *rocDefault = new TGraph(nCuts, sigEffDefault, bgRejDefault);
  TGraph *rocCustom  = new TGraph(nCuts, sigEffCustom, bgRejCustom);

  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,600);

  TH2F *dummy = new TH2F("dummy","",100,0,1,100,0,1);
  dummy->GetXaxis()->SetTitle("signal efficiency");
  dummy->GetYaxis()->SetTitleOffset(1.4);
  dummy->GetYaxis()->SetTitle("background rejection");
  dummy->GetXaxis()->SetRangeUser(0.30, 1.0);
  dummy->GetYaxis()->SetRangeUser(0.90, 1.0);
  dummy->Draw();

  rocDefault->SetLineColor(kBlue);
  rocCustom->SetLineColor(kRed);

  rocDefault->SetLineWidth(2);
  rocCustom->SetLineWidth(2);

  rocDefault->Draw("L,same");
  rocCustom->Draw("L,same");

  TLegend *leg = new TLegend(0.2,0.2, 0.7, 0.4);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(rocDefault,"ISO with standard EA","l");
  leg->AddEntry(rocCustom,"ISO with EA from 90% cutoff","l");
  leg->Draw();
}

