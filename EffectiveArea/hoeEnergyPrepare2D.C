{

  //const TString finName = "/eos/user/i/ikrav/ElectronID/92X/DYJetsToLL_cutID_tuning_92X_v1.root";
  const TString finName = "~/DYJetsToLL_cutID_tuning_92X_v1.root";
  const TString treeName = "ntupler/ElectronTree";

  TFile *fin = new TFile(finName);
  TTree *tree = (TTree*)fin->Get(treeName);

  //TString cuts = "isTrue==1 && abs(etaSC)<1.4442 && pt>10 && hOverE>0";
  TString cuts = "isTrue==1 && abs(etaSC)>1.566 && abs(etaSC)<2.5 && pt>10 && hOverE>0";
  
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
  
  tree->Draw(command, cuts, "colz");

  TFile *fout = new TFile("hoeVsE.root","recreate");
  fout->cd();
  hoeVsE->Write();
  fout->Close();

}
