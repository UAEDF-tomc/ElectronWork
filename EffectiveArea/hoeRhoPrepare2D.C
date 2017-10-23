{

  //const TString finName = "/eos/user/i/ikrav/ElectronID/92X/DYJetsToLL_cutID_tuning_92X_v1.root";
  const TString finName = "~/DYJetsToLL_cutID_tuning_92X_v1.root";
  const TString treeName = "ntupler/ElectronTree";

  TFile *fin = new TFile(finName);
  TTree *tree = (TTree*)fin->Get(treeName);

  //TString cuts = "isTrue==1 && abs(etaSC)<1.4442 && pt>20 && hOverE>0";
  TString cuts = "isTrue==1 && abs(etaSC)>1.566 && abs(etaSC)<2.5 && pt>20 && hOverE>0";
  
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
  
  tree->Draw(command, cuts, "colz");

  TFile *fout = new TFile("hVsRho.root","recreate");
  fout->cd();
  hVsRho->Write();
  fout->Close();

}
