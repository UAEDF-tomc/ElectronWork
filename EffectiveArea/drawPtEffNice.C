{

  const bool isBarrel = true;
  const int nFiles = 4;
  TString fileNamesSig[nFiles] = {
    "effPt_DY_EB_plain_sig.root",
    "effPt_DY_EB_95percent_sig.root",
    "effPt_DY_EB_90percent_sig.root",
    "effPt_DY_EB_mean_sig.root"
  };
  TString fileNamesBg[nFiles] = {
    "effPt_DY_EB_plain_bg.root",
    "effPt_DY_EB_95percent_bg.root",
    "effPt_DY_EB_90percent_bg.root",
    "effPt_DY_EB_mean_bg.root"
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
    histSig[i]->SetMarkerStyle(20);
    histSig[i]->SetMarkerSize(0.7);
    histSig[i]->SetMarkerColor(colors[i]);
    histSig[i]->SetLineWidth(2);
    histSig[i]->SetLineColor(colors[i]);
  }

  TCanvas *c1 = new TCanvas("c1","c1",10,10,500,500);
  c1->cd();
  gStyle->SetOptStat(0);

  histSig[0]->GetXaxis()->SetTitle("p_{T} [GeV]");
  histSig[0]->GetXaxis()->SetTitleOffset(1.2);
  histSig[0]->GetYaxis()->SetTitle("efficiency");
  histSig[0]->GetYaxis()->SetTitleOffset(1.2);
  histSig[0]->GetXaxis()->SetRangeUser(0,300);
  histSig[0]->GetYaxis()->SetRangeUser(0.8,1);

  histSig[0]->Draw("L");
  for(int i=1; i<nFiles; i++){
    histSig[i]->Draw("L,same");
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

  //
  // ========== Draw background curves ==================
  //

  TFile *finBg[nFiles];
  TH1F *histBg[nFiles];
  
  for(int i=0; i<nFiles; i++){
    finBg[i] = new TFile(fileNamesBg[i]);
    histBg[i] = (TH1F*)finBg[i]->Get("eff");
    histBg[i]->SetMarkerStyle(20);
    histBg[i]->SetMarkerSize(0.7);
    histBg[i]->SetMarkerColor(colors[i]);
    histBg[i]->SetLineWidth(2);
    histBg[i]->SetLineColor(colors[i]);
  }

  TCanvas *c2 = new TCanvas("c2","c2",10,10,500,500);
  c2->cd();
  gStyle->SetOptStat(0);

  histBg[0]->GetXaxis()->SetTitle("p_{T} [GeV]");
  histBg[0]->GetXaxis()->SetTitleOffset(1.2);
  histBg[0]->GetYaxis()->SetTitle("efficiency");
  histBg[0]->GetYaxis()->SetTitleOffset(1.2);
  histBg[0]->GetXaxis()->SetRangeUser(0,300);
  histBg[0]->GetYaxis()->SetRangeUser(0.,1);

  histBg[0]->Draw("L");
  for(int i=1; i<nFiles; i++){
    histBg[i]->Draw("L,same");
  }

  TLegend *leg2 = new TLegend(0.2, 0.2, 0.8, 0.4);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(histBg[0], "flat H/E cut", "l");
  leg2->AddEntry(histBg[1], "95%-based scaled H/E cut", "l");
  leg2->AddEntry(histBg[2], "90%-based scaled H/E cut", "l");
  leg2->AddEntry(histBg[3], "mean-based scaled H/E cut", "l");
  leg2->Draw();

  TLatex *lat2 = new TLatex(0.2, 0.6, region);
  lat2->SetNDC(kTRUE);
  lat2->Draw();


}
