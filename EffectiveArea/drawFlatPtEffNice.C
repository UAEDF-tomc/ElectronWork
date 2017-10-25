{

  const bool isBarrel = false;
  const int nFiles = 4;
  TString fileNamesSig[nFiles] = {
    "effPt_flat_EE_plain_sig.root",
    "effPt_flat_EE_95percent_sig.root",
    "effPt_flat_EE_90percent_sig.root",
    "effPt_flat_EE_mean_sig.root"
  };
  
  const int colors[nFiles] = {kBlack, kRed, kBlue, kMagenta};

  TString region = "Barrel";
  if( !isBarrel )
    region = "Endcap";

  //
  // ========== Draw signal curves ==================
  //

  TFile *finSig[nFiles];
  TH1F *histSig[nFiles];
  
  for(int i=0; i<nFiles; i++){
    finSig[i] = new TFile(fileNamesSig[i]);
    histSig[i] = (TH1F*)finSig[i]->Get("eff");
    // histSig[i]->SetMarkerStyle(20);
    // histSig[i]->SetMarkerSize(0.7);
    // histSig[i]->SetMarkerColor(colors[i]);
    histSig[i]->SetLineWidth(2);
    histSig[i]->SetLineColor(colors[i]);
  }

  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,500);
  c1->cd();
  gStyle->SetOptStat(0);

  histSig[0]->GetXaxis()->SetTitle("p_{T} [GeV]");
  histSig[0]->GetXaxis()->SetTitleOffset(1.2);
  histSig[0]->GetYaxis()->SetTitle("efficiency");
  histSig[0]->GetYaxis()->SetTitleOffset(1.2);
  // histSig[0]->GetXaxis()->SetRangeUser(0,300);
  // histSig[0]->GetYaxis()->SetRangeUser(0.8,1);

  histSig[0]->Draw("E2");
  histSig[0]->Draw("L same hist");
  for(int i=1; i<nFiles; i++){
    histSig[i]->Draw("E2 same");
    histSig[i]->Draw("L same hist");
  }

  TLegend *leg1 = new TLegend(0.2, 0.2, 0.8, 0.4);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(histSig[0], "flat H/E cut", "l");
  leg1->AddEntry(histSig[1], "95%-based scaled H/E cut", "l");
  leg1->AddEntry(histSig[2], "90%-based scaled H/E cut", "l");
  leg1->AddEntry(histSig[3], "mean-based scaled H/E cut", "l");
  leg1->Draw();

  TLatex *lat1 = new TLatex(0.2, 0.6, region);
  lat1->SetNDC(kTRUE);
  lat1->Draw();


}
