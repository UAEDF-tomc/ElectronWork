{

  TFile *fin1 = new TFile("roc_EE_flat_20_1000.root");
  TGraph *gr1 = (TGraph*)fin1->Get("rocGraph");

  TFile *fin2 = new TFile("roc_EE_scaled_20_1000.root");
  TGraph *gr2 = (TGraph*)fin2->Get("rocGraph");

  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);
  c1->cd();
  gStyle->SetOptStat(0);

  TH2F *dummy = new TH2F("dummy","",100, 0.8, 1.0, 100, 0.0, 1.0);
  dummy->GetXaxis()->SetTitle("signal efficiency");
  dummy->GetYaxis()->SetTitle("background rejection");
  dummy->Draw();

  gr1->SetLineWidth(2);
  gr1->SetLineColor(kBlue);

  gr2->SetLineWidth(2);
  gr2->SetLineColor(kRed);

  gr1->Draw("L,same");
  gr2->Draw("L,same");

  TLegend *leg = new TLegend(0.2, 0.2, 0.4 ,0.4);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gr1, "traditional H/E cut", "l");
  leg->AddEntry(gr2, "H/E cut with scaling", "l");
  leg->Draw();

  TLatex *lat = new TLatex(0.2, 0.6, "Endcap");
  lat->SetNDC(kTRUE);
  lat->Draw();

}

