{

  const TString finName1 = "~/DoubleEleFlat1to300_92X_v1.root";
  const TString finName2 = "~/DoubleEleFlat300to6500_92X_v1.root";
  const TString treeName = "ntupler/ElectronTree";

  TFile *fin1 = new TFile(finName1);
  TTree *tree1 = (TTree*)fin1->Get(treeName);

  TFile *fin2 = new TFile(finName2);
  TTree *tree2 = (TTree*)fin2->Get(treeName);

  TString cuts = "(abs(etaSC)<1.4442 && pt>20 && pt<100 && hOverE>0)";
  //TString cuts = "(abs(etaSC)>1.566 && abs(etaSC)<2.5 && pt>20 && hOverE>0)";
  
  const float rhoMin = -0.5;
  const float rhoMax = 50.5;
  const int rhoBins = 51;
  
  const float hMin = 0.0;
  const float hMax = 1000;
  const float hBins = 1000;

  const TString hName = "hVsRho";
  TH2F *hVsRho = new TH2F(hName, "",rhoBins, rhoMin, rhoMax, 
			  hBins, hMin, hMax);

  TString command = "hOverE*eSC:rho>>";
  command += hName;
  tree1->Draw(command, cuts, "colz");

  command = "hOverE*eSC:rho>>+";
  command += hName;
  // Empirical weight to have right relative normalization
  // between the 1-300 and 300-6500 flat samples
  cuts = "21.70*" + cuts;
  tree2->Draw(command, cuts, "colz");

  TFile *fout = new TFile("hVsRho.root","recreate");
  fout->cd();
  hVsRho->Write();
  fout->Close();

}
