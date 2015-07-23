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

// 
// Forward declarations
//
float getCutoff(TH1F *hist, float frac);

//
// Main program
//

void drawChIsoConvolution(){

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

  // Define histograms
  TH1F *hChIso = new TH1F("hChIso","",100,0,20);
  hChIso->GetXaxis()->SetTitle("charged hadron isolation, GeV");
  hChIso->SetLineWidth(2);
  TH1F *hChIsoConvPoisson = (TH1F*)hChIso->Clone("hChIsoConvPoisson");
  hChIsoConvPoisson->SetLineColor(kRed);
  TH1F *hChIsoConvGauss = (TH1F*)hChIso->Clone("hChIsoConvGauss");
  hChIsoConvGauss->SetLineColor(kBlue);

  // Define functions
  TF1 *fp = new TF1("fp","[0]*TMath::Poisson(x,[1])",0,20);
  // The Poisson mean of 7.5 is what we see for the pile-up distribution for 
  // rho range 35-40.
  float norm = 1.0;
  float mean = 7.5;
  fp->SetLineColor(kRed);
  fp->SetParameters( norm, mean);

  TF1 *fg = new TF1("fg","gaus(0) / sqrt(2*TMath::Pi())", 0, 20);
  float sigma = 0.5;
  fg->SetParameters( norm, mean, sigma);
  fg->SetLineColor(kBlue);
  fg->SetNpx(1000);


  int nEle; // the number of reconstructed electrons in the event
  std::vector <float> *elePt = 0;         // electron PT
  std::vector <float> *eleEtaSC = 0;      // supercluser eta
  std::vector <float> *isoChargedHadrons = 0;
  std::vector <int> *isTrue = 0;

  // Declare branches
  TBranch *b_nEle = 0;
  TBranch *b_elePt = 0;
  TBranch *b_eleEtaSC = 0;
  TBranch *b_isoChargedHadrons = 0;
  TBranch *b_isTrue = 0;

  treeSignal->SetBranchAddress("nEle", &nEle, &b_nEle);
  treeSignal->SetBranchAddress("pt", &elePt, &b_elePt);
  treeSignal->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  treeSignal->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
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
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nEle);
    
    // Get data for all electrons in this event, only vars of interest
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_isoChargedHadrons->GetEntry(tentry);
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

      float iso = isoChargedHadrons->at(iele);
      hChIso->Fill( iso );
      hChIsoConvPoisson->Fill( iso + fp->GetRandom() );
      hChIsoConvGauss  ->Fill( iso + fg->GetRandom() );

    } // end loop over electrons

  } // end loop over events

  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,600);
  hChIso->Draw("hist");
  norm = hChIso->GetSumOfWeights();
  // Update function parameters: now we know what the norm is
  fp->SetParameters(norm,mean);
  fg->SetParameters(norm,mean,sigma);
  // Draw the functions
  fp->Draw("same");
  fg->Draw("same");

  TLegend *leg1 = new TLegend(0.4,0.7,0.88,0.9);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry(hChIso,"charged hadrons isolation in MC", "l");
  leg1->AddEntry(fp, "relaistic ISO_{ch} from PU (35<#rho<40)","l");
  leg1->AddEntry(fg, "a narrow (unrealistic) ISO_{ch} from PU", "l");
  leg1->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",10,10,800,600);
  hChIsoConvGauss->Draw("hist");
  hChIsoConvPoisson->Draw("same,hist");

  float frac = 0.95;
  float cutoffGauss   = getCutoff(hChIsoConvGauss  , frac);
  float cutoffPoisson = getCutoff(hChIsoConvPoisson, frac);
  printf("cutG %f   cutP  %f\n", cutoffGauss, cutoffPoisson);

  float ymax = hChIsoConvGauss->GetMaximum();
  TArrow *arrowGauss = new TArrow( cutoffGauss, ymax*0.5, cutoffGauss, 0, 0.02,"|>");
  arrowGauss->SetLineColor(kBlue);
  arrowGauss->SetLineWidth(2);
  arrowGauss->Draw();

  TArrow *arrowPoisson = new TArrow( cutoffPoisson, ymax*0.4, cutoffPoisson, 0, 0.02,"|>");
  arrowPoisson->SetLineColor(kRed);
  arrowPoisson->SetLineWidth(2);
  arrowPoisson->Draw();

  float histMean = hChIsoConvGauss->GetMean();
  TLine *meanLine = new TLine(histMean, 0, histMean, ymax*0.1);
  meanLine->SetLineWidth(3);
  meanLine->SetLineColor(kGreen+2);
  meanLine->SetLineStyle(kDashed);
  meanLine->Draw();

  TLatex *textMean = new TLatex(histMean-1, ymax*0.12, "mean");
  textMean->SetTextColor(kGreen+2);
  textMean->SetTextSize(0.04);
  textMean->Draw();

  TLatex *textRed = new TLatex(cutoffPoisson-1, ymax*0.42, "95% iso-Red");
  textRed->SetTextColor(kRed);
  textRed->SetTextSize(0.04);
  textRed->Draw();

  TLatex *textBlue = new TLatex(cutoffGauss-1, ymax*0.52, "95% iso-Blue");
  textBlue->SetTextColor(kBlue);
  textBlue->SetTextSize(0.04);
  textBlue->Draw();

  TLegend *leg2 = new TLegend(0.5,0.7,0.88,0.9);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(hChIsoConvPoisson,"ISO_{ch}+ISO_{ch from PU} realistic", "l");
  leg2->AddEntry(hChIsoConvGauss,"ISO_{ch}+ISO_{ch from PU} narrow", "l");
  leg2->Draw();


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

