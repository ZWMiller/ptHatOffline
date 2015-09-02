// Fraction Fit - Z. Miller Sep 1, 2015
//
// .L fractionFit.C
// fractionFit() 
// takes output of offline.C (pythia version, and data version) as inputs.
// Copy current best data and templates as "current%s.root" {B,C,Data}

void fractionFit()
{
  gStyle->SetOptFit(1111);
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
  const Int_t numPtBins = 12;
  Float_t lowpt[14] ={2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.5,10.,14.0};
  Float_t highpt[14]={3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.5,10.,14.,200.};
  Float_t hptCut=0.5;
  Double_t p00[numPtBins],p01[numPtBins],p20[numPtBins],p21[numPtBins];
  Double_t e00[numPtBins],e01[numPtBins],e20[numPtBins],e21[numPtBins];
  Double_t pC0[numPtBins],pC1[numPtBins],eC0[numPtBins],eC1[numPtBins];
  Double_t Rb0[numPtBins],Rb2[numPtBins],RbC[numPtBins],pT[numPtBins];
  Double_t eb0[numPtBins],eb2[numPtBins],ebC[numPtBins],dx[numPtBins];
  Int_t rangeLow = 23; //45 for no extra rebin
  Int_t rangeHigh = 28; //56
  

  // Make Canvases
  TCanvas* deltaPhi  = new TCanvas("deltaPhi","Pythia Delta Phi",150,0,1150,1000);
  TCanvas* fitResult0 = new TCanvas("fitResult0","RB Extraction HT0",150,0,1150,1000);
  TCanvas* fitResult2 = new TCanvas("fitResult2","RB Extraction HT2",150,0,1150,1000);
  TCanvas* fitResultC = new TCanvas("fitResultC","RB Extraction Combined Trigs",150,0,1150,1000);
  deltaPhi->Divide(4,3);
  fitResult0->Divide(4,3);
  fitResult2->Divide(4,3);
  fitResultC->Divide(4,3);

  // Make histos
  TH1D* projB[numPtBins];
  TH1D* projC[numPtBins];
  TH1D* projData0[numPtBins];
  TH1D* projData2[numPtBins];
  TH1D* combData[numPtBins];
  TH1F* histoNorms;
  
  // Get and Draw histos
  TPaveText* lbl[numPtBins];
  TPaveText* stat[3][numPtBins];
  char statLabel[100];
  char textLabel[100];
  Int_t plotbin;
  Float_t norm0,norm2;

  // Get ptbin independent hists
  histoNorms = (TH1F*)fD->Get("histoNorms");

  for(Int_t ptbin=0; ptbin<numPtBins; ptbin++)
    {
      norm0 = histoNorms->GetBinContent(histoNorms->GetBin(1,ptbin+1));
      norm2 = histoNorms->GetBinContent(histoNorms->GetBin(3,ptbin+1));
      plotbin = ptbin;
      // Init necessary plotting tools
      lbl[ptbin] = new TPaveText(.25,.8,.5,.88,Form("NB NDC%i",ptbin));
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
      Int_t RB2 = 2;
      projB[ptbin]->Rebin(RB);
      projC[ptbin]->Rebin(RB);
      projData0[ptbin]->Rebin(RB);
      projData2[ptbin]->Rebin(RB); 
      
      // Draw Templates on own plots
      deltaPhi->cd(plotbin+1);
      projData0[ptbin]->SetLineColor(kBlue);
      projData2[ptbin]->SetLineColor(kGreen+3);
      projB[ptbin]->SetLineColor(kRed);
      projC[ptbin]->SetLineColor(kBlack);
      // projC[ptbin]->GetYaxis()->SetRangeUser(0.,1.);
      projC[ptbin]->GetXaxis()->SetRangeUser(-3.5,3.5);
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

      combData[ptbin] = (TH1D*) projData0[ptbin]->Clone();
      combData[ptbin]->Add(projData2[ptbin]);

      // Do the actual fit
      TObjArray *mc = new TObjArray(2);   // MC histograms are put in this array
      mc->Add(projC[ptbin]);
      mc->Add(projB[ptbin]);

      fitResult0->cd(ptbin+1);
      TFractionFitter* fit = new TFractionFitter(projData0[ptbin], mc,"V"); // initialise
      fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
      fit->SetRangeX(rangeLow,rangeHigh);      // use only the first 15 bins in the fit
      Int_t status = fit->Fit();               // perform the fit
      std::cout << "fit status: " << status << std::endl;
      if (status == 0) {                       // check on fit status
	TH1F* result = (TH1F*) fit->GetPlot();
	projData0[ptbin]->GetXaxis()->SetRangeUser(-3.5,3.5);
	//projData0[ptbin]->GetYaxis()->SetRangeUser(0,1);
	projData0[ptbin]->Draw("Ep");
	result->SetLineColor(kBlack);
	result->SetFillColor(kWhite);
	result->Draw("sames");
	lbl[ptbin]->Draw("same");
	fit->GetResult(0,p00[ptbin],e00[ptbin]);
	fit->GetResult(1,p01[ptbin],e01[ptbin]);
	Double_t chi2 = fit->GetChisquare();
	Int_t ndf = fit->GetNDF();
	stat[0][ptbin] = new TPaveText(.25,.6,.75,.8,Form("NB NDC%i",ptbin));
	sprintf(statLabel,"Chi2: %f",chi2);
	stat[0][ptbin]->InsertText(statLabel);
	sprintf(statLabel,"Rc: %f; Rb: %f",p00[ptbin],p01[ptbin]);
	stat[0][ptbin]->InsertText(statLabel);
	sprintf(statLabel,"eC: %f; eB: %f",e00[ptbin],e01[ptbin]);
	stat[0][ptbin]->InsertText(statLabel);
	stat[0][ptbin]->SetFillColor(kWhite);
	stat[0][ptbin]->Draw("same");
      }

      
      fitResult2->cd(ptbin+1);
      TFractionFitter* fit2 = new TFractionFitter(projData2[ptbin], mc); // initialise
      fit2->Constrain(0,0.0,1.0);              // constrain fraction 0
      fit2->Constrain(1,0.0,1.0);              // constrain fraction 1 to be between 0 and 1
      fit2->SetRangeX(rangeLow,rangeHigh);     // use only the first 15 bins in the fit
      Int_t status2 = fit2->Fit();             // perform the fit
      std::cout << "fit status: " << status2 << std::endl;
      if (status2 == 0) {                       // check on fit status
	TH1F* result2 = (TH1F*) fit2->GetPlot();
	projData2[ptbin]->GetXaxis()->SetRangeUser(-3.5,3.5);
	//projData2[ptbin]->GetYaxis()->SetRangeUser(0,1);
	projData2[ptbin]->Draw("Ep");
	result2->SetLineColor(kBlack);
	result2->SetFillColor(kWhite);
	result2->Draw("sames");
	lbl[ptbin]->Draw("same");
	fit2->GetResult(0,p20[ptbin],e20[ptbin]);
	fit2->GetResult(1,p21[ptbin],e21[ptbin]);
	Double_t chi2 = fit2->GetChisquare();
	Int_t ndf = fit2->GetNDF();
	stat[1][ptbin] = new TPaveText(.25,.6,.75,.8,Form("NB NDC%i",ptbin));
	sprintf(statLabel,"Chi2: %f",chi2);
	stat[1][ptbin]->InsertText(statLabel);
	sprintf(statLabel,"Rc: %f; Rb: %f",p20[ptbin],p21[ptbin]);
	stat[1][ptbin]->InsertText(statLabel);
	sprintf(statLabel,"eC: %f; eB: %f",e20[ptbin],e21[ptbin]);
	stat[1][ptbin]->InsertText(statLabel);
	stat[1][ptbin]->SetFillColor(kWhite);
	stat[1][ptbin]->Draw("same");
      }

      fitResultC->cd(ptbin+1);
      TFractionFitter* fitC = new TFractionFitter(combData[ptbin], mc); // initialise
      fitC->Constrain(0,0.0,1.0);              // constrain fraction 0
      fitC->Constrain(1,0.0,1.0);              // constrain fraction 1 to be between 0 and 1
      fitC->SetRangeX(rangeLow,rangeHigh);     // use only the first 15 bins in the fit
      Int_t statusC = fitC->Fit();             // perform the fit
      std::cout << "fit status: " << statusC << std::endl;
      if (statusC == 0) {                       // check on fit status
	TH1F* resultC = (TH1F*) fitC->GetPlot();
	combData[ptbin]->GetXaxis()->SetRangeUser(-3.5,3.5);
	//combData[ptbin]->GetYaxis()->SetRangeUser(0,1);
	combData[ptbin]->SetTitle("");
	combData[ptbin]->Draw("Ep");
	resultC->SetLineColor(kBlack);
	resultC->SetFillColor(kWhite);
	resultC->Draw("sames");
	lbl[ptbin]->Draw("same");
	fitC->GetResult(0,pC0[ptbin],eC0[ptbin]);
	fitC->GetResult(1,pC1[ptbin],eC1[ptbin]);
	Double_t chi2 = fitC->GetChisquare();
	Int_t ndf = fitC->GetNDF();
	stat[2][ptbin] = new TPaveText(.25,.6,.75,.8,Form("NB NDC%i",ptbin));
	sprintf(statLabel,"Chi2: %f",chi2);
	stat[2][ptbin]->InsertText(statLabel);
	sprintf(statLabel,"Rc: %f; Rb: %f",pC0[ptbin],pC1[ptbin]);
	stat[2][ptbin]->InsertText(statLabel);
	sprintf(statLabel,"ec: %f; eb: %f",eC0[ptbin],eC1[ptbin]);
	stat[2][ptbin]->InsertText(statLabel);
	stat[2][ptbin]->SetFillColor(kWhite);
	stat[2][ptbin]->Draw("same");
      }

      pT[ptbin] = (lowpt[ptbin]+highpt[ptbin])/2.;
      dx[ptbin] = 0;
      Rb0[ptbin] = p01[ptbin]/(p01[ptbin]+p00[ptbin]);
      Rb2[ptbin] = p21[ptbin]/(p21[ptbin]+p20[ptbin]);
      RbC[ptbin] = pC1[ptbin]/(pC1[ptbin]+pC0[ptbin]);

      eb0[ptbin] = sqrt(pow((p00[ptbin]/pow((p00[ptbin]+p01[ptbin]),2)*e00[ptbin]),2)+pow((p01[ptbin]/pow((p00[ptbin]+p01[ptbin]),2)*e01[ptbin]),2));
      eb2[ptbin] = sqrt(pow((p20[ptbin]/pow((p20[ptbin]+p21[ptbin]),2)*e20[ptbin]),2)+pow((p21[ptbin]/pow((p20[ptbin]+p21[ptbin]),2)*e21[ptbin]),2));
      ebC[ptbin] = sqrt(pow((pC0[ptbin]/pow((pC0[ptbin]+pC1[ptbin]),2)*eC0[ptbin]),2)+pow((pC1[ptbin]/pow((pC0[ptbin]+pC1[ptbin]),2)*eC1[ptbin]),2));
      
    }
  TCanvas* c1 = new TCanvas("c1","Bottom Contribution",150,0,1150,1000);
  TGraphErrors *gr0  = new TGraphErrors(numPtBins,pT,Rb0,dx,eb0);
  TGraphErrors *gr2  = new TGraphErrors(numPtBins,pT,Rb2,dx,eb2);
  TGraphErrors *grC  = new TGraphErrors(numPtBins,pT,RbC,dx,ebC);
  c1->cd(1);

  gr0->SetTitle("Bottom Contribution");
  gr0->GetXaxis()->SetTitle("p_{T,e}");
  gr0->GetYaxis()->SetTitle("#frac{r_{B}}{(r_{B}+r_{C})}");
  gr0->SetMarkerStyle(20);
  gr0->SetMarkerSize(1);
  gr0->SetLineColor(kBlue);
  gr0->SetMarkerColor(kBlue);
  gr0->GetXaxis()->SetRangeUser(2,10);
  gr0->GetYaxis()->SetRangeUser(0,1);
  gr2->SetMarkerStyle(22);
  gr2->SetMarkerSize(1);
  gr2->SetLineColor(kGreen+3);
  gr2->SetMarkerColor(kGreen+3);
  grC->SetMarkerStyle(21);
  grC->SetMarkerSize(1);
  grC->SetLineColor(kRed);
  grC->SetMarkerColor(kRed);
  
  gr0->Draw("AP");
  gr2->Draw("same P");
  // grC->Draw("same P");

  TLegend* leg2 = new TLegend(0.15,0.73,0.35,0.85);
  leg2->AddEntry(gr0,"High Tower 0 Trigs","pe");
  leg2->AddEntry(gr2,"High Tower 2 Trigs","pe");
  // leg2->AddEntry(grC,"Combined Trigs","pe");
  leg2->Draw("same");
  
  
}
