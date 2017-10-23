{

  // Signal sample
  //const TString finName = "~/DYJetsToLL_cutID_tuning_92X_v1.root";

  // Background sample
  const TString finName1 = "~/TTJets_cutID_92X_v1.root";

  const TString treeName = "ntupler/ElectronTree";

  TFile *fin1 = new TFile(finName1);
  TTree *tree1 = (TTree*)fin1->Get(treeName);


  // Barrel cuts
  //TString cutsGeneral = "isTrue==1 && abs(etaSC)<1.4442 && pt>10";
  TString cutsGeneral = "(isTrue==0 || isTrue==3) && abs(etaSC)<1.4442 && pt>10";

  // Endcap cuts
  //TString cuts = "isTrue==1 && abs(etaSC)>1.566 && abs(etaSC)<2.5 && pt>10";
  //TString cuts = "((=isTrue==0 || isTrue==3) && abs(etaSC)>1.566 && abs(etaSC)<2.5 && pt>10";

  // Flat cut
  //TString cutHOE = "hOverE < 0.100";
  // Scaled cut
  // rho 90%-based, E 95%-based
  TString cutHOE = "hOverE < 3.91/eSC + 0.0368*rho/eSC + 0.002";
  // rho 90%-based, E 90%-based
  //TString cutHOE = "hOverE < 3.09/eSC + 0.0368*rho/eSC + 0.0165";
  // rho 90%-based, E mean-based
  //TString cutHOE = "hOverE < 1.63/eSC + 0.0368*rho/eSC + 0.0455";

  // Random
  //TString cutHOE = "hOverE < 2.00/eSC + 0.0322*rho/eSC + 0.03";

  TString cutsWithHOE = cutsGeneral;
  cutsWithHOE += " && ";
  cutsWithHOE += cutHOE;

  TString cutsDen1 = TString::Format("1.0*(%s)",cutsGeneral.Data());
  TString cutsNum1 = TString::Format("1.0*(%s)",cutsWithHOE.Data());

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

  printf("Denominator cuts: \n %s\n", cutsDen1.Data());
  printf("Numirator cuts: \n %s\n", cutsNum1.Data());

  gStyle->SetOptStat(0);
  tree1->Draw("pt>>hDen" ,cutsDen1, "");

  printf("Bin 10 denominator: %f\n", hDen->GetBinContent(10));

  tree1->Draw("pt>>hNum" ,cutsNum1, "");

  printf("Bin 10 numerator: %f\n", hNum->GetBinContent(10));

  TH1F *hEff = (TH1F*)hNum->Clone("hEff");
  hEff->Divide(hDen);
  
  hEff->Draw();

  printf("Eff is %f\n", hNum->GetSumOfWeights()/hDen->GetSumOfWeights());

  TFile *fout = new TFile("effPt.root", "recreate");
  hEff->Write("eff");
  fout->Close();

}
