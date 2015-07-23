{

  int maxEvents = 5000000;
  bool drawLine = false;
  float cutVal = 2.0;

  // Signal uncorrected
  TFile *fsig = new TFile("/afs/cern.ch/user/r/rkamalie/workspace/public/DY_Spring15_Asympt50ns_24june2015.root");
  TTree *treeSignal     = (TTree*)fsig->Get("ntupler/ElectronTree");
  
  TCanvas *c1 = new TCanvas("c1","",10,10,600,600);
  gStyle->SetOptStat(0);

  TH2F *hsig = new TH2F("hsig","",50,0,50, 75, -5, 10);
  hsig->GetXaxis()->SetTitle("rho");
  hsig->GetYaxis()->SetTitle("ISO_{neu.had}+ISO_{pho}, GeV");
  hsig->GetYaxis()->SetTitleOffset(1.2);
  treeSignal->Draw("isoNeutralHadrons+isoPhotons:rho>>hsig","pt>20 && abs(etaSC)<0.8 && isTrue==1", "colz", maxEvents);
  TLine *line1 = new TLine(0, cutVal, 50, cutVal+50*0.0973);
  line1->SetLineWidth(2);
  line1->SetLineColor(kRed+1);
  if(drawLine) line1->Draw();
  c1->Update();

  // Signal corrected
  TCanvas *c2 = new TCanvas("c2","",10,10,600,600);
  TH2F *hsigCorr = (TH2F*)hsig->Clone("hsigCorr");
  hsig->Reset();
  hsigCorr->GetYaxis()->SetTitle("ISO_{neu.had}+ISO_{pho}-#rhoA_{eff}, GeV");
  treeSignal->Draw("(isoNeutralHadrons+isoPhotons-0.0973*rho):rho>>hsigCorr",
		   "pt>20 && abs(etaSC)<0.8 && isTrue==1", "colz", maxEvents);
  TLine *line2 = new TLine(0, cutVal, 50, cutVal);
  line2->SetLineWidth(2);
  line2->SetLineColor(kRed+1);
  if(drawLine) line2->Draw();
  c2->Update();

  // Background uncorrected
  TFile *fbg = new TFile("/afs/cern.ch/user/r/rkamalie/workspace/public/TT_Spring15_Asympt50ns_24june2015.root");
  TTree *treeBg     = (TTree*)fbg->Get("ntupler/ElectronTree");
  
  TCanvas *c3 = new TCanvas("c3","",10,10,600,600);
  gStyle->SetOptStat(0);

  TH2F *hbg = new TH2F("hbg","",50,0,50, 75, -5, 10);
  hbg->GetXaxis()->SetTitle("rho");
  hbg->GetYaxis()->SetTitle("ISO_{neu.had}+ISO_{pho}, GeV");
  hbg->GetYaxis()->SetTitleOffset(1.2);
  treeBg->Draw("isoNeutralHadrons+isoPhotons:rho>>hbg",
	       "pt>20 && abs(etaSC)<0.8 && (isTrue==0 || isTrue==3)", "colz", maxEvents);
  TLine *line3 = new TLine(0, cutVal, 50, cutVal+50*0.0973);
  line3->SetLineWidth(2);
  line3->SetLineColor(kRed+1);
  if(drawLine) line3->Draw();
  c3->Update();

  // Background corrected
  TCanvas *c4 = new TCanvas("c4","",10,10,600,600);
  TH2F *hbgCorr = (TH2F*)hbg->Clone("hbgCorr");
  hbg->Reset();
  hbgCorr->GetYaxis()->SetTitle("ISO_{neu.had}+ISO_{pho}-#rhoA_{eff}, GeV");
  treeBg->Draw("(isoNeutralHadrons+isoPhotons-0.0973*rho):rho>>hbgCorr",
	       "pt>20 && abs(etaSC)<0.8 && (isTrue==0 || isTrue==3)", "colz", maxEvents);
  TLine *line4 = new TLine(0, cutVal, 50, cutVal);
  line4->SetLineWidth(2);
  line4->SetLineColor(kRed+1);
  if(drawLine) line4->Draw();
  c4->Update();

}
