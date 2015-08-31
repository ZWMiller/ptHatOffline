// Offline Plots - Z. Miller July 24, 2015
//
// .L fractionFit.C
// fractionFit() 
// takes output of offline.C (pythia version, and data version) as inputs.
// Copy current best data and templates as "current%s.root" {B,C,Data}

void fractionFit()
{
   char name[1000];
  sprintf(name,"/Users/zach/Research/pythia/ptHatTemplate/outputs/currentB.root");
  TFile *fB = new TFile(name,"READ");
  sprintf(name,"/Users/zach/Research/pythia/ptHatTemplate/outputs/currentC.root");
  TFile *fC = new TFile(name,"READ");
   sprintf(name,"/Users/zach/Research/rootFiles/run12NPEhPhi/currentData.root");
  TFile *fD = new TFile(name,"READ");
  if (fB->IsOpen()==kFALSE || fC->IsOpen()==kFALSE)
    { std::cout << "!!!!!! Either B,C, or Data File not found !!!!!!" << std::endl
		<< "Looking for currentB.root, currentC.root, and currentData.root" << std::endl;
      exit(1); }
  
  // Set constants and projection bins
  const Int_t numPtBins = 10;
  Float_t lowpt[14] ={2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.5,10.,14.0};
  Float_t highpt[14]={3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.5,10.,14.,200.};
  Float_t hptCut=0.5;
  Double_t p00[numPtBins],p01[numPtBins],p20[numPtBins],p21[numPtBins];
  Double_t e00[numPtBins],e01[numPtBins],e20[numPtBins],e21[numPtBins];
  Double_t Rb0[numPtBins],Rb2[numPtBins],pT[numPtBins];

  // Make Canvases
  TCanvas* deltaPhi  = new TCanvas("deltaPhi","Pythia Delta Phi",150,0,1150,1000);
  TCanvas* fitResult0 = new TCanvas("fitResult0","RB Extraction HT0",150,0,1150,1000);
  TCanvas* fitResult2 = new TCanvas("fitResult2","RB Extraction HT2",150,0,1150,1000);
  deltaPhi->Divide(4,3);
  fitResult0->Divide(4,3);
  fitResult2->Divide(4,3);

  // Make histos
  TH1D* projB[numPtBins];
  TH1D* projC[numPtBins];
  TH1D* projData0[numPtBins];
  TH1D* projData2[numPtBins];
  
  // Get and Draw histos
  TPaveText* lbl[numPtBins];
  char textLabel[100];
  Int_t plotbin;

  for(Int_t ptbin=0; ptbin<numPtBins; ptbin++)
    {
      plotbin = ptbin;
      // Init necessary plotting tools
      lbl[ptbin] = new TPaveText(.25,.76,.5,.82,Form("NB NDC%i",ptbin));
      sprintf(textLabel,"%.1f < P_{T,e} < %.1f",lowpt[ptbin],highpt[ptbin]);
      lbl[ptbin]->AddText(textLabel);
      lbl[ptbin]->SetFillColor(kWhite);

      projB[ptbin] = (TH1D*)fB->Get(Form("delPhi_%i",ptbin));
      projC[ptbin] = (TH1D*)fC->Get(Form("delPhi_%i",ptbin));
      projData0[ptbin]= (TH1D*)fD->Get(Form("NPEhDelPhi_0_%i",ptbin));
      projData2[ptbin]= (TH1D*)fD->Get(Form("NPEhDelPhi_2_%i",ptbin));
      Int_t RB = 2;
      projB[ptbin]->Rebin(RB);
      projC[ptbin]->Rebin(RB);
      
      // Draw Templates on own plots
      deltaPhi->cd(plotbin+1);
      projData0[ptbin]->SetLineColor(kBlue);
      projData2[ptbin]->SetLineColor(kGreen+3);
      projB[ptbin]->SetLineColor(kRed);
      projC[ptbin]->SetLineColor(kBlack);
      projC[ptbin]->GetYaxis()->SetRangeUser(0.,1.);
      projC[ptbin]    -> Draw();
      projB[ptbin]    -> Draw("same");
      projData0[ptbin]-> Draw("same");
      projData2[ptbin]-> Draw("same");
      lbl[ptbin]      -> Draw("same");

      TLegend* leg = new TLegend(0.5,0.73,0.85,0.85);
      leg->AddEntry(projB[ptbin],"b#bar{b}->NPE","lpe");
      leg->AddEntry(projC[ptbin],"c#bar{c}->NPE","lpe");
      leg->AddEntry(projData0[ptbin],"HT0","lpe");
      leg->AddEntry(projData2[ptbin],"HT2","lpe");
      leg->Draw();

      // Do the actual fit
      TObjArray *mc = new TObjArray(2);   // MC histograms are put in this array
      mc->Add(projC[ptbin]);
      mc->Add(projB[ptbin]);

      fitResult0->cd(ptbin+1);
      TFractionFitter* fit = new TFractionFitter(projData0[ptbin], mc,"V"); // initialise
      fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
      fit->SetRangeX(45,56);                    // use only the first 15 bins in the fit
      Int_t status = fit->Fit();               // perform the fit
      std::cout << "fit status: " << status << std::endl;
      if (status == 0) {                       // check on fit status
	TH1F* result = (TH1F*) fit->GetPlot();
	projData0[ptbin]->GetXaxis()->SetRangeUser(-3.5,3.5);
	projData0[ptbin]->GetYaxis()->SetRangeUser(0,0.5);
	projData0[ptbin]->Draw("Ep");
	result->SetLineColor(kBlack);
	result->SetFillColor(kWhite);
	result->Draw("same");
	lbl[ptbin]->Draw("same");
	fit->GetResult(0,p00[ptbin],e00[ptbin]);
	fit->GetResult(1,p01[ptbin],e01[ptbin]);
      }

      
      fitResult2->cd(ptbin+1);
      TFractionFitter* fit2 = new TFractionFitter(projData2[ptbin], mc); // initialise
      fit2->Constrain(0,0.0,1.0);              // constrain fraction 0
      fit2->Constrain(1,0.0,1.0);              // constrain fraction 1 to be between 0 and 1
      fit2->SetRangeX(45,56);                  // use only the first 15 bins in the fit
      Int_t status2 = fit2->Fit();             // perform the fit
      std::cout << "fit status: " << status2 << std::endl;
      if (status2 == 0) {                       // check on fit status
	TH1F* result2 = (TH1F*) fit2->GetPlot();
	projData2[ptbin]->GetXaxis()->SetRangeUser(-3.5,3.5);
	projData2[ptbin]->GetYaxis()->SetRangeUser(0,0.5);
	projData2[ptbin]->Draw("Ep");
	result2->SetLineColor(kBlack);
	result2->SetFillColor(kWhite);
	result2->Draw("same");
	lbl[ptbin]->Draw("same");
	fit2->GetResult(0,p20[ptbin],e20[ptbin]);
	fit2->GetResult(1,p21[ptbin],e21[ptbin]);
      }

      pT[ptbin] = (lowpt[ptbin]+highpt[ptbin])/2.;
      Rb0[ptbin] = p01[ptbin]/(p01[ptbin]+p00[ptbin]);
      Rb2[ptbin] = p21[ptbin]/(p21[ptbin]+p20[ptbin]);

    }
  TCanvas* c1 = new TCanvas("c1","Bottom Contribution",150,0,1150,1000);
  TGraphErrors *gr0  = new TGraphErrors(numPtBins,pT,Rb0);
  TGraphErrors *gr2  = new TGraphErrors(numPtBins,pT,Rb2);
  c1->cd(1);
  gr0->SetMarkerStyle(20);
  gr0->SetMarkerSize(1);
  gr0->SetLineColor(kBlue);
  gr0->SetMarkerColor(kBlue);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(1);
  gr2->SetLineColor(kGreen+3);
  gr2->SetMarkerColor(kGreen+3);
  gr0->Draw();
  gr2->Draw("same");
  
  
}
