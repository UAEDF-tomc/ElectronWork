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
const bool smallEventCount = true; // DEBUG

//
// Main program
//

void drawIsoWithCutoffs(){

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

  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,600);
  c1->SetRightMargin(0.25);

  TH2F *hIsoVsRho = new TH2F("hIsoVsRho","",50,0,50, 55, -1, 10);
  hIsoVsRho->GetXaxis()->SetTitle("rho");
  hIsoVsRho->GetYaxis()->SetTitle("neutral hadron + photon isolation, GeV");



  treeSignal->Draw("isoNeutralHadrons+isoPhotons:rho>>hIsoVsRho",
		   "isTrue==1 && abs(etaSC)<0.8 && pt>20","colz",1000000);

  const int nc = 10;
  const float cutoff_A[nc] = {-0.1097, -0.0566, 0.0233, 0.0992, 0.2072,
			      0.3542, 0.5068, 0.7352, 1.1225, 2.1576};
  const float cutoff_B[nc] = {0.0820, 0.0889, 0.0957, 0.1038, 0.1116,
			      0.1197, 0.1325, 0.1481, 0.1708, 0.2054};

  TF1 *func[nc];
  for(int i=0; i<nc; i++){
    TString funname = TString::Format("func_%d",i);
    func[i] = new TF1(funname, "pol1", 3, 50);
    func[i]->SetParameters(cutoff_A[i], cutoff_B[i]);
    func[i]->SetLineWidth(2);
    func[i]->SetLineColor(kRed-3);
    func[i]->Draw("same");
  }

  // Slope from the mean of the iso distribution
  float slope = 0.0973;
  // Cut value (convert to comparable rel iso cut by dividing by <pt_ele>=40 GeV)
  //    Cut corresponding to rel iso of 0.00125 (flat eff)
  //float cut = 0.05;
  //    Cut corresponding to rel iso of 0.05
  //float cut = 2.0;
  //    Cut corresponding to the const term of the mean fit (~0.013 of rel iso cut)
  float cut = 0.520;

  TF1 *fcut = new TF1("fcut","pol1",3,50);
  fcut->SetParameters(cut, slope);
  fcut->SetLineWidth(3);
  fcut->SetLineColor(kGreen);
  fcut->Draw("same");

  // You can't get the palette without actually drawing the histogram
  // first! So we call update here.
  c1->Update();

  TPaletteAxis *palette = (TPaletteAxis*)hIsoVsRho->GetListOfFunctions()->FindObject("palette");  
  cout << palette->GetX1() << endl;
  palette->SetX1NDC(0.85);
  palette->SetX2NDC(0.90);
  palette->Draw();
  c1->Update();

  TLatex *txt[nc];
  const float xpos[nc] = {0.60, 0.76, 0.76, 0.76, 0.76,
			  0.76, 0.76, 0.76, 0.76, 0.76};
  const float ypos[nc] = {0.92, 0.87, 0.75, 0.69, 0.63,
			  0.59, 0.55, 0.51, 0.48, 0.45};
  const TString percentage[nc] = {"95%","90%","85%","80%","75%",
			   "70%","65%","60%","55%","50%"};
  for(int i=0; i<nc; i++){
    txt[i] = new TLatex(xpos[i],ypos[i], percentage[i]);
    txt[i]->SetTextColor(kRed-3);
    txt[i]->SetNDC(kTRUE);
    txt[i]->SetTextSize(0.03);
    txt[i]->Draw("same");
  }

}
