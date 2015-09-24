// .x hptComparison.C
// Plots rB for different hpt cuts. Change file names within file as necessary

{

  char name[1000];
  sprintf(name,"/Users/zach/Research/pythia/ptHatTemplate/FFOutput/Sep24_0p3_1par_FIT.root");
  TFile *f31 = new TFile(name,"READ");
   sprintf(name,"/Users/zach/Research/pythia/ptHatTemplate/FFOutput/Sep24_0p5_1par_FIT.root");
  TFile *f51 = new TFile(name,"READ");
  sprintf(name,"/Users/zach/Research/pythia/ptHatTemplate/FFOutput/Sep24_0p3_2par_FIT.root");
  TFile *f32 = new TFile(name,"READ");
  sprintf(name,"/Users/zach/Research/pythia/ptHatTemplate/FFOutput/Sep24_0p5_2par_FIT.root");
  TFile *f52 = new TFile(name,"READ");

  // get the shared by all plots
  TGraphErrors* FONLL    = (TGraphErrors*)f31->Get("FONLL");
  TGraphErrors* FONLLmax = (TGraphErrors*)f31->Get("FONLLmax");
  TGraphErrors* FONLLmin = (TGraphErrors*)f31->Get("FONLLmin");
  TGraphErrors* prevData = (TGraphErrors*)f31->Get("PreviousData");
  // get the 0.3 hpt cut 2 param fit
  TGraphErrors* HT03    = (TGraphErrors*)f32->Get("HT0");
  TGraphErrors* HT23    = (TGraphErrors*)f32->Get("HT2");
  // get the 0.5 hpt cut 2 param fit
  TGraphErrors* HT05    = (TGraphErrors*)f52->Get("HT0");
  TGraphErrors* HT25    = (TGraphErrors*)f52->Get("HT2");

  TCanvas* c = new TCanvas("c","Bottom Contribution",150,0,1150,1000);
  c->cd();

  // Draw things that don't need to be changed
  prevData->SetTitle("Bottom Contribution to NPE");
  prevData->GetXaxis()->SetRangeUser(0,10);
  prevData->Draw("AP");
  FONLL->Draw("same");
  FONLLmax->Draw("same");
  FONLLmin->Draw("same");

  // Set Marker Styles
  HT03->SetMarkerStyle(22); // Up Triangle Closed
  HT23->SetMarkerStyle(21); // Square Closed
  HT05->SetMarkerStyle(23); // Down Triangle Closed
  HT25->SetMarkerStyle(29); // Star Closed
  HT03->SetMarkerSize(1);
  HT23->SetMarkerSize(1);
  HT05->SetMarkerSize(1);
  HT25->SetMarkerSize(1.5);
  HT03->SetMarkerColor(4); // Colors human readable in line style 
  HT05->SetMarkerColor(6);  
  HT23->SetMarkerColor(2);  
  HT25->SetMarkerColor(4); 

  // Set Line Styles
  HT03->SetLineColor(4);           // Green
  HT05->SetLineColor(6);           // Blue
  HT23->SetLineColor(2);           // Red
  HT25->SetLineColor(4);           // Magenta

  // Draw (use "same P" to remove connecting line)
  HT03->Draw("same P");
  HT23->Draw("same P");
  HT05->Draw("same P");
  HT25->Draw("same P");

  TLegend* leg2 = new TLegend(0.15,0.68,0.4,0.85);
  leg2->AddEntry(HT03,"High Tower 0 Trigs, hPt > 0.3","pe");
  leg2->AddEntry(HT23,"High Tower 2 Trigs, hPt > 0.3","pe");
  leg2->AddEntry(HT05,"High Tower 0 Trigs, hPt > 0.5","pe");
  leg2->AddEntry(HT25,"High Tower 2 Trigs, hPt > 0.5","pe");
  //leg2->AddEntry(grC,"Combined Trigs","pe");
  leg2->AddEntry(prevData,"Run 5/6 Analysis (Stat Uncertainty)","pe");
  // leg2->AddEntry(grPr,"Run 5/6 Refit (new Template)","pe");
  //leg2->AddEntry(grPPr,"Run 5/6 Refit (prev Template)","pe");
  leg2->AddEntry(FONLL,"FONLL (Uncertainty: Scale Only)","l");
  leg2->Draw("same");
}
