#include <algorithm>
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


// Fall17
// 90% contours
const float ea_neutral_total_iso[7] = {  0.1566,  0.1626,  0.1073,  0.0854,  0.1051,  0.1204,  0.1524};

// Spring16 
// 90% contours
//const float ea_neutral_total_iso[7] = {0.1703,  0.1715,  0.1213,  0.1230,  0.1635, 0.1937, 0.2393};
// Mean-based
//const float ea_neutral_total_iso[7] = {  0.0932,  0.0919,  0.0512,  0.0591,  0.0882,  0.1101,  0.1373};

//Spring15 50ns:
// With no gen level weights:
//const float ea_neutral_total_iso[5] = {  0.0973,  0.0954,  0.0632,  0.0727,  0.1337};
// With gen level weights:
//const float ea_neutral_total_iso[5] = {  0.0958,  0.0940,  0.0616,  0.0708,  0.1321};

// Tests:
//const float ea_neutral_total_iso[5] = {  0.0655,  0.0674,  0.0319,  0.0551,  0.1149};
//const float ea_neutral_total_iso[5] = {  0.0828,  0.0821,  0.0563,  0.0640,  0.1138};
//const float ea_neutral_total_iso[5] = {  0.0735,  0.0745,  0.0494,  0.0552,  0.1102};
//const float ea_neutral_total_iso[5] = {  0.0566,  0.0611,  0.0406,  0.0507,  0.1017};
//const float ea_neutral_total_iso[5] = {  0.1251,  0.1335,  0.0892,  0.1338,  0.2195};
//const float ea_neutral_total_iso[5] = {  0.1682,  0.0954,  0.0632,  0.0727,  0.1337};

// Phys14:
//const float ea_neutral_total_iso[5] = {  0.1041,  0.1038,  0.0590,  0.0825,  0.1551};

// For this exercise, we use two MC samples: the signal sample
// and a background-reach sample. Both have the same ntuple structure.
//
// Signal sample: DYToLL
const TString fileNameSignal =
  "/eos/user/i/ikrav/ElectronID/92X/DYJetsToLL_cutID_tuning_92X_v1.root";

// Directory and tree name:
const TString treeName = "ntupler/ElectronTree";

const bool verbose = false;
const bool smallEventCount = false;

const bool useWeights = true;

const float boundaryBarrelEndcap = 1.479; 

// Selection cuts
// Kinematics
const float ptCut = 20; 

const float isoCut = 0.100;

const int nEtaBins = 5;
const float etaBinLimits[nEtaBins+1] = {0.0, 0.8, 1.3, 2.0, 2.2, 2.5};

const int nNvtxBins = 15;
const float nVtxBinLimits[nNvtxBins+1] = {5, 10, 13, 15, 
					  16, 17, 18, 19, 20, 
					  22, 25, 30, 35, 40, 45, 50};

//
// Forward declarations
//
void computeEfficiency(TH1F *hnum, TH1F *hdenum, TH1F *heff);

//
// Main program
//

void validateCorrectedIsolationV2(bool forBarrel = true){

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 

  // General settings
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);

  // Book histograms

  // Efficiencies
  TH1F *numEffRaw   = new TH1F("numEffRaw"  ,"", nNvtxBins, nVtxBinLimits);
  TH1F *numEffDBeta = new TH1F("numEffDBeta","", nNvtxBins, nVtxBinLimits);
  TH1F *numEffEA    = new TH1F("numEffEA"   ,"", nNvtxBins, nVtxBinLimits);

  TH1F *denomEff    = new TH1F("denomEff"   ,"", nNvtxBins, nVtxBinLimits);

  TH1F *histEffRaw   = new TH1F("histEffRaw"  ,"", nNvtxBins, nVtxBinLimits);
  TH1F *histEffDBeta = new TH1F("histEffDBeta","", nNvtxBins, nVtxBinLimits);
  TH1F *histEffEA    = new TH1F("histEffEA"   ,"", nNvtxBins, nVtxBinLimits);

  // Fake rates
  TH1F *numFakeRaw   = new TH1F("numFakeRaw"  ,"", nNvtxBins, nVtxBinLimits);
  TH1F *numFakeDBeta = new TH1F("numFakeDBeta","", nNvtxBins, nVtxBinLimits);
  TH1F *numFakeEA    = new TH1F("numFakeEA"   ,"", nNvtxBins, nVtxBinLimits);

  TH1F *denomFake    = new TH1F("denomFake"   ,"", nNvtxBins, nVtxBinLimits);

  TH1F *histFakeRaw   = new TH1F("histFakeRaw"  ,"", nNvtxBins, nVtxBinLimits);
  TH1F *histFakeDBeta = new TH1F("histFakeDBeta","", nNvtxBins, nVtxBinLimits);
  TH1F *histFakeEA    = new TH1F("histFakeEA"   ,"", nNvtxBins, nVtxBinLimits);

  // Since we might use weighted events:
  numEffRaw    ->Sumw2(); 
  numEffDBeta  ->Sumw2(); 
  numEffEA     ->Sumw2(); 
  //
  denomEff     ->Sumw2(); 
  //
  // histEffRaw   ->Sumw2(); 
  // histEffDBeta ->Sumw2(); 
  // histEffEA    ->Sumw2(); 
  //
  numFakeRaw   ->Sumw2(); 
  numFakeDBeta ->Sumw2(); 
  numFakeEA    ->Sumw2(); 
  //
  denomFake    ->Sumw2(); 
  //
  // histFakeRaw  ->Sumw2();
  // histFakeDBeta->Sumw2();
  // histFakeEA   ->Sumw2();


  histEffRaw->GetYaxis()->SetRangeUser(0.7, 1.0);
  histEffRaw->GetXaxis()->SetTitle("Nvtx");
  histEffRaw->GetYaxis()->SetTitle("Efficiency");
  histEffRaw->GetYaxis()->SetTitleOffset(1.4);

  histFakeRaw->GetYaxis()->SetRangeUser(0.0, 0.5);
  histFakeRaw->GetXaxis()->SetTitle("Nvtx");
  histFakeRaw->GetYaxis()->SetTitle("Fake rate");
  histFakeRaw->GetYaxis()->SetTitleOffset(1.4);

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
  float genWeight;
  int nPV; //reconstructed primary vertices
  // Per-eletron variables
  // Kinematics
  std::vector <float> *elePt = 0;         // electron PT
  std::vector <float> *eleEtaSC = 0;      // supercluser eta
  std::vector <float> *elePhiSC = 0;      // supercluser phi
  // Variables for analysis
  std::vector <float> *isoChargedHadrons = 0;
  std::vector <float> *isoNeutralHadrons = 0;
  std::vector <float> *isoPhotons = 0;
  std::vector <float> *isoChargedFromPU = 0;
  std::vector <int> *isTrue = 0;
  // Other vars  
  // Impact parameters
  std::vector <float> *eleD0 = 0;      // r-phi plane impact parameter
  std::vector <float> *eleDZ = 0;      // r-z plane impact parameter
  // Matching track-supercluster
  std::vector <float> *eleDEtaIn = 0;  // deltaEtaIn
  std::vector <float> *eleDPhiIn = 0;  // deltaPhiIn
  // Misc ID variables
  std::vector <float> *eleHoverE = 0;  // H/E  
  std::vector <float> *eleFull5x5SigmaIEtaIEta = 0;  
  std::vector <float> *eleOOEMOOP = 0; // |1/E - 1/p|
  // Conversion rejection
  std::vector <float> *eleExpectedMissingInnerHits = 0;
  std::vector <float> *elePassConversionVeto = 0;


  // Declare branches
  TBranch *b_nEle = 0;
  TBranch *b_rho = 0;
  TBranch *b_genWeight = 0;
  TBranch *b_nPV = 0;
  TBranch *b_elePt = 0;
  TBranch *b_eleEtaSC = 0;
  TBranch *b_elePhiSC = 0;
  TBranch *b_isoChargedHadrons = 0;
  TBranch *b_isoNeutralHadrons = 0;
  TBranch *b_isoPhotons = 0;
  TBranch *b_isoChargedFromPU = 0;
  TBranch *b_isTrue;
  // Other vars
  TBranch *b_eleD0 = 0;
  TBranch *b_eleDZ = 0;
  TBranch *b_eleDEtaIn = 0;
  TBranch *b_eleDPhiIn = 0;
  TBranch *b_eleHoverE = 0;
  TBranch *b_eleFull5x5SigmaIEtaIEta = 0;
  TBranch *b_eleOOEMOOP = 0;
  TBranch *b_eleExpectedMissingInnerHits = 0;
  TBranch *b_elePassConversionVeto = 0;


  // Connect variables and branches to the tree with the data
  treeSignal->SetBranchAddress("nEle", &nEle, &b_nEle);
  treeSignal->SetBranchAddress("rho", &rho, &b_rho);
  treeSignal->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
  treeSignal->SetBranchAddress("nPV", &nPV, &b_nPV);
  treeSignal->SetBranchAddress("pt", &elePt, &b_elePt);
  treeSignal->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  treeSignal->SetBranchAddress("phiSC", &elePhiSC, &b_elePhiSC);
  treeSignal->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
  treeSignal->SetBranchAddress("isoNeutralHadrons", &isoNeutralHadrons, &b_isoNeutralHadrons);
  treeSignal->SetBranchAddress("isoPhotons",        &isoPhotons,        &b_isoPhotons);
  treeSignal->SetBranchAddress("isoChargedFromPU",  &isoChargedFromPU,  &b_isoChargedFromPU);
  treeSignal->SetBranchAddress("isTrue",    &isTrue,    &b_isTrue);
  treeSignal->SetBranchAddress("d0",                &eleD0,             &b_eleD0);
  treeSignal->SetBranchAddress("dz",                &eleDZ,             &b_eleDZ);
  treeSignal->SetBranchAddress("dEtaIn",            &eleDEtaIn,         &b_eleDEtaIn);
  treeSignal->SetBranchAddress("dPhiIn",            &eleDPhiIn,         &b_eleDPhiIn);
  treeSignal->SetBranchAddress("hOverE",            &eleHoverE,         &b_eleHoverE);
  treeSignal->SetBranchAddress("full5x5_sigmaIetaIeta", &eleFull5x5SigmaIEtaIEta,
			       &b_eleFull5x5SigmaIEtaIEta);
  treeSignal->SetBranchAddress("ooEmooP",           &eleOOEMOOP,        &b_eleOOEMOOP);
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
    b_genWeight->GetEntry(tentry);
    b_nPV->GetEntry(tentry);
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_elePhiSC->GetEntry(tentry);
    b_isoChargedHadrons->GetEntry(tentry);
    b_isoNeutralHadrons->GetEntry(tentry);
    b_isoPhotons->GetEntry(tentry);
    b_isoChargedFromPU->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);
    // Other vars
    b_eleD0->GetEntry(tentry);
    b_eleDZ->GetEntry(tentry);
    b_eleDEtaIn->GetEntry(tentry);
    b_eleDPhiIn->GetEntry(tentry);
    b_eleHoverE->GetEntry(tentry);
    b_eleFull5x5SigmaIEtaIEta->GetEntry(tentry);
    b_eleOOEMOOP->GetEntry(tentry);
    b_eleExpectedMissingInnerHits->GetEntry(tentry);
    b_elePassConversionVeto->GetEntry(tentry);

    // Nested loops over the electrons
    for(int iele = 0; iele < nEle; iele++){

      // Check kinematics:
      if( !(elePt->at(iele) > ptCut) )
	continue;

      // Loose ID of 2012 (VETO WP)
      const bool useID = false; 
      if( useID ){
	if( abs(eleEtaSC->at(iele)) <  boundaryBarrelEndcap ){
	  if( abs(eleDEtaIn->at(iele))>0.007 ) continue;
	  if( abs(eleDPhiIn->at(iele))>0.8 ) continue;
	  if( eleFull5x5SigmaIEtaIEta->at(iele) >0.01 ) continue;
	  if( eleHoverE->at(iele) > 0.15 ) continue;
	  if( abs(eleD0->at(iele)) > 0.04 ) continue;
	  if( abs(eleDZ->at(iele)) > 0.2 ) continue;
	}else{
	  if( abs(eleDEtaIn->at(iele))>0.001 ) continue;
	  if( abs(eleDPhiIn->at(iele))>0.7 ) continue;
	  if( eleFull5x5SigmaIEtaIEta->at(iele) >0.03 ) continue;
	  if( abs(eleD0->at(iele)) > 0.04 ) continue;
	  if( abs(eleDZ->at(iele)) > 0.2 ) continue;
	}
      } // end if loose ID
    
      // Choose barrel or endcap
      if(forBarrel){
	if( !(abs(eleEtaSC->at(iele))<boundaryBarrelEndcap) ) continue;
      }	else {
	if( !(abs(eleEtaSC->at(iele))>boundaryBarrelEndcap) ) continue;
      }
	  
      // Find eta bin
      if( abs(eleEtaSC->at(iele))>etaBinLimits[nEtaBins] ) continue;
      int ieta = 0; 
      while ( ieta < nEtaBins-1 
	      && abs(eleEtaSC->at(iele)) > etaBinLimits[ieta+1] )
	{ ++ieta; };

      // Compute relative combined corrected isolation
      double isoCh = isoChargedHadrons->at(iele);
      double isoNh = isoNeutralHadrons->at(iele);
      double isoPh = isoPhotons->at(iele);
      double isoChPU = isoChargedFromPU->at(iele);
      float pt = elePt->at(iele);
      float relRawIso = (isoCh + isoNh + isoPh)/pt;
      float relIsoWithDBeta = (isoCh + std::max(0.0, isoNh + isoPh - 0.5*isoChPU))/pt;
      float relIsoWithEA = (isoCh + std::max(0.0, isoNh + isoPh
					     - rho*ea_neutral_total_iso[ieta]) )/pt;

      float weight = genWeight;
      if( !useWeights)
	weight = 1;

      // Fill efficiencie related histograms
      int isTrueEle = isTrue->at(iele); 
      if( isTrueEle == 1 ){

	denomEff->Fill(nPV, weight);
	
	if( relRawIso < isoCut )
	  numEffRaw->Fill(nPV, weight);
	
	if( relIsoWithDBeta < isoCut )
	  numEffDBeta->Fill(nPV, weight);
	
	if( relIsoWithEA < isoCut )
	  numEffEA->Fill(nPV, weight);

      }else if ( isTrueEle == 0 || isTrueEle == 3 ){ 

	denomFake->Fill(nPV, weight);
	
	if( relRawIso < isoCut )
	  numFakeRaw->Fill(nPV, weight);
	
	if( relIsoWithDBeta < isoCut )
	  numFakeDBeta->Fill(nPV, weight);
	
	if( relIsoWithEA < isoCut )
	  numFakeEA->Fill(nPV, weight);

      }

    } // end loop over the electrons

  } // end loop over events
  printf("\n");

  //
  // Compute and draw efficiencies
  //
  
  computeEfficiency(numEffRaw  , denomEff, histEffRaw);
  computeEfficiency(numEffDBeta, denomEff, histEffDBeta);
  computeEfficiency(numEffEA   , denomEff, histEffEA);

  computeEfficiency(numFakeRaw  , denomFake, histFakeRaw);
  computeEfficiency(numFakeDBeta, denomFake, histFakeDBeta);
  computeEfficiency(numFakeEA   , denomFake, histFakeEA);

  //
  // Draw efficiencies
  //
  TString canvasName = "c1";
  TCanvas *c1 = new TCanvas(canvasName,canvasName,10,10,600,600);
  c1->cd();

  histEffRaw  ->SetMarkerStyle(24);
  histEffDBeta->SetMarkerStyle(20);
  histEffEA   ->SetMarkerStyle(20);

  histEffRaw  ->SetMarkerColor(kBlack);
  histEffDBeta->SetMarkerColor(kRed);
  histEffEA   ->SetMarkerColor(kBlue);

  histEffRaw->Draw("PE");
  histEffDBeta->Draw("pe,same");
  histEffEA->Draw("pe,same");

  TLegend *leg = new TLegend(0.2, 0.2, 0.5, 0.4);
  leg->AddEntry(histEffRaw  , "uncorrected iso", "pl");
  leg->AddEntry(histEffDBeta, "iso with #Delta#beta corr", "pl");
  leg->AddEntry(histEffEA, "iso with #rho*EA corr", "pl");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("same");

  //
  // Draw fake rates
  //
  TString canvasName2 = "c2";
  TCanvas *c2 = new TCanvas(canvasName2,canvasName2,10,10,600,600);
  c2->cd();

  histFakeRaw  ->SetMarkerStyle(24);
  histFakeDBeta->SetMarkerStyle(20);
  histFakeEA   ->SetMarkerStyle(20);

  histFakeRaw  ->SetMarkerColor(kBlack);
  histFakeDBeta->SetMarkerColor(kRed);
  histFakeEA   ->SetMarkerColor(kBlue);

  histFakeRaw->Draw("PE");
  histFakeDBeta->Draw("pe,same");
  histFakeEA->Draw("pe,same");

  TLegend *leg2 = new TLegend(0.2, 0.6, 0.5, 0.8);
  leg2->AddEntry(histFakeRaw  , "uncorrected iso", "pl");
  leg2->AddEntry(histFakeDBeta, "iso with #Delta#beta corr", "pl");
  leg2->AddEntry(histFakeEA, "iso with #rho*EA corr", "pl");
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->Draw("same");


}

void computeEfficiency(TH1F *hnum, TH1F *hdenum, TH1F *heff){

  heff->SetMarkerSize(1);

  for(int i=1; i<=nNvtxBins; i++){
    float a = hnum->GetBinContent(i);
    float b = hdenum->GetBinContent(i);
    float da = hnum->GetBinError(i);
    float db = hdenum->GetBinError(i);
    // Non-trivial manipulation is to handle the case
    // of weighted events. Convert denum into pass and fail,
    // that are independent.
    float pass = a;
    float dpass = da;
    float fail = b - a;
    float dfail = sqrt(db*db - da*da);
    //
    float eff = 0;
    float effErr = 0;
    if ( b!= 0 ){
      eff = a/b;
      // The formula below is from simple propagation of errors for
      // pass/(pass+fail) where the errors are dpass and dfail.
      effErr = sqrt( 
		    (fail*fail * dpass*dpass + pass*pass * dfail*dfail)
		    / (b*b*b*b) 
		     );
      //effErr = sqrt(eff*(1-eff)/b);
    }
    heff->SetBinContent(i,eff);
    heff->SetBinError(i,effErr);
  }

  return;
}

