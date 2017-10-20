{

  const TString finName1 = "/eos/user/i/ikrav/ElectronID/92X/DoubleEleFlat1to300_92X_v1.root";
  const TString finName2 = "/eos/user/i/ikrav/ElectronID/92X/DoubleEleFlat300to6500_92X_v1.root";
  const TString treeName = "ntupler/ElectronTree";

  TFile *fin1 = new TFile(finName1);
  TTree *tree1 = (TTree*)fin1->Get(treeName);

  TFile *fin2 = new TFile(finName2);
  TTree *tree2 = (TTree*)fin2->Get(treeName);

  TString cutsGeneral = "abs(etaSC)<1.4442 && pt>10";
  //TString cuts = "abs(etaSC)>1.566 && abs(etaSC)<2.5 && pt>10 && hOverE>0";

  // Flat cut
  //TString cutHOE = "hOverE < 0.100";
  // Scaled cut
  //TString cutHOE = "hOverE < 3.91/eSC + 0.0322*rho/eSC + 0.003";
  TString cutHOE = "hOverE < 2.00/eSC + 0.0322*rho/eSC + 0.03";

  TString cuts = cutsGeneral;
  cuts += " && ";
  cuts += cutHOE;

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
  tree1->Draw("pt>>hDen",cutsGeneral, "");
  tree2->Draw("pt>>+hDen",cutsGeneral, "");

  tree1->Draw("pt>>hNum",cuts, "");
  tree2->Draw("pt>>+hNum",cuts, "");

  TH1F *hEff = (TH1F*)hNum->Clone("hEff");
  hEff->Divide(hDen);
  
  hEff->Draw();

  printf("Eff is %f\n", hNum->GetSumOfWeights()/hDen->GetSumOfWeights());

}
