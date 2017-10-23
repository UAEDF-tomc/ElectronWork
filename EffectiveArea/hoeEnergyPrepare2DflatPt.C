{

  const TString finName1 = "~/DoubleEleFlat1to300_92X_v1.root";
  const TString finName2 = "~/DoubleEleFlat300to6500_92X_v1.root";
  const TString treeName = "ntupler/ElectronTree";

  TFile *fin1 = new TFile(finName1);
  TTree *tree1 = (TTree*)fin1->Get(treeName);

  TFile *fin2 = new TFile(finName2);
  TTree *tree2 = (TTree*)fin2->Get(treeName);

  TString cuts = "(abs(etaSC)<1.4442 && pt>10 && hOverE>0)";
  //TString cuts = "(abs(etaSC)>1.566 && abs(etaSC)<2.5 && pt>10 && hOverE>0)";
  
  const float eMin  = 0;
  const float eMax  = 1000;
  const int   eBins = 1000;
  
  const float hoeMin = 0.0;
  const float hoeMax = 0.5;
  const float hoeBins = 1000;

  const TString hName = "hoeVsE";
  TH2F *hoeVsE = new TH2F(hName, "",eBins, eMin, eMax, 
			    hoeBins, hoeMin, hoeMax);

  TString command = "hOverE:eSC>>";
  command += hName;
  tree1->Draw(command, cuts, "colz");

  command = "hOverE:eSC>>+";
  command += hName;
  // Empirical weight to have right relative normalization
  // between the 1-300 and 300-6500 flat samples
  cuts = "21.70*" + cuts;
  tree2->Draw(command, cuts, "colz");


  TFile *fout = new TFile("hoeVsE.root","recreate");
  fout->cd();
  hoeVsE->Write();
  fout->Close();

}
