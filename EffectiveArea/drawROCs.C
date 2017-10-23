{

  TFile *fin1 = new TFile("roc_plain_EB_DY_pt20to1000.root");
  TGraph *gr1 = (TGraph*)fin1->Get("rocGraph");

  TFile *fin2 = new TFile("roc_scaled_mean_EB_DY_pt20to1000.root");
  TGraph *gr2 = (TGraph*)fin2->Get("rocGraph");

  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);
  c1->cd();
  gStyle->SetOptStat(0);

  TH2F *dummy = new TH2F("dummy","",100, 0.8, 1.0, 100, 0.0, 1.0);
  dummy->GetXaxis()->SetTitle("signal efficiency");
  dummy->GetYaxis()->SetTitle("background rejection");
  dummy->GetXaxis()->SetTitleOffset(1.2);
  dummy->GetYaxis()->SetTitleOffset(1.2);

  dummy->Draw();

  gr1->SetLineWidth(2);
  gr1->SetLineColor(kBlue);

  gr2->SetLineWidth(2);
  gr2->SetLineColor(kRed);

  gr1->Draw("L,same");
  gr2->Draw("L,same");

  TLegend *leg = new TLegend(0.2, 0.2, 0.8 ,0.4);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gr1, "traditional H/E cut", "l");
  leg->AddEntry(gr2, "H/E cut with scaling", "l");
  leg->Draw();

  TLatex *lat = new TLatex(0.2, 0.6, "Barrel");
  lat->SetNDC(kTRUE);
  lat->Draw();

}

