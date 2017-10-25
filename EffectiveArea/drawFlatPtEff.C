{

  const TString finName1 = "/eos/user/i/ikrav/ElectronID/92X/DoubleEleFlat1to300_92X_v1.root";
  const TString finName2 = "/eos/user/i/ikrav/ElectronID/92X/DoubleEleFlat300to6500_92X_v1.root";
  const TString treeName = "ntupler/ElectronTree";

  TFile *fin1 = new TFile(finName1);
  TTree *tree1 = (TTree*)fin1->Get(treeName);

  TFile *fin2 = new TFile(finName2);
  TTree *tree2 = (TTree*)fin2->Get(treeName);

  //TString cutsGeneral = "abs(etaSC)<1.4442 && pt>10";
  TString cutsGeneral = "abs(etaSC)>1.566 && abs(etaSC)<2.5 && pt>10 && hOverE>0";

  // Flat cut
  //TString cutHOE = "hOverE < 0.100";
  // Scaled cut
  // rho 95%-based, E 95%-based
  //TString cutHOE = "hOverE < 7.01/eSC + 0.201*rho/eSC + 0.0155";
  // rho 95%-based, E 90%-based
  //TString cutHOE = "hOverE < 5.18/eSC + 0.201*rho/eSC + 0.0275";
  // rho 95%-based, E mean-based
  TString cutHOE = "hOverE < 2.65/eSC + 0.201*rho/eSC + 0.0465";

  // Random
  //TString cutHOE = "hOverE < 2.00/eSC + 0.0322*rho/eSC + 0.03";

  TString cutsWithHOE = cutsGeneral;
  cutsWithHOE += " && ";
  cutsWithHOE += cutHOE;

  TString cutsDen1 = TString::Format("1.0*(%s)",cutsGeneral.Data());
  TString cutsNum1 = TString::Format("1.0*(%s)",cutsWithHOE.Data());

  // Add weight for relative normalization of two flat samples
  TString cutsDen2 = TString::Format("21.70*(%s)",cutsGeneral.Data());
  TString cutsNum2 = TString::Format("21.70*(%s)",cutsWithHOE.Data());

  // Variable binning
  //  1 - 300 GeV, bins every 2 GeV
  // 300 - 2000 GeV, bins every 20 GeV 
  const int nbins = 235;
  double xlimits[nbins+1];
  xlimits[0] = 0;
  for(int ibin=1; ibin <= nbins; ibin++){
    if( ibin <= 150 )
      xlimits[ibin] = ibin * 2;
    else
      xlimits[ibin] = 300 + (ibin - 150)*20;
  }

  TH1F *hDen = new TH1F("hDen", "", nbins, xlimits);
  TH1F *hNum = new TH1F("hNum", "", nbins, xlimits);

  gStyle->SetOptStat(0);
  tree1->Draw("pt>>hDen" ,cutsDen1, "");
  tree2->Draw("pt>>+hDen",cutsDen2, "");

  tree1->Draw("pt>>hNum" ,cutsNum1, "");
  tree2->Draw("pt>>+hNum",cutsNum2, "");

  TH1F *hEff = (TH1F*)hNum->Clone("hEff");
  hEff->Divide(hDen);
  
  hEff->Draw();

  printf("Eff is %f\n", hNum->GetSumOfWeights()/hDen->GetSumOfWeights());


  TFile *fout = new TFile("effPt.root", "recreate");
  hEff->Write("eff");
  fout->Close();

}
