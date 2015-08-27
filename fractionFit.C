// Offline Plots - Z. Miller July 24, 2015
//
// .L fractionFit.C
// fractionFit() 
// takes output of offline.C (pythia version, and data version) as inputs.
// Copy current best data and templates as "current%s.root" {B,C,Data}

void fractionFit()
{
   char name[1000];
  sprintf(name,"/Users/zach/Research/pythia/npeTemplate/outputs/currentB.root");
  TFile *fB = new TFile(name,"READ");
  sprintf(name,"/Users/zach/Research/pythia/npeTemplate/outputs/currentC.root");
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
      lbl[ptbin] = new TPaveText(.2,.76,.5,.82,Form("NB NDC%i",ptbin));
      sprintf(textLabel,"%.1f < P_{T,e} < %.1f",lowpt[ptbin],highpt[ptbin]);
      lbl[ptbin]->AddText(textLabel);
      lbl[ptbin]->SetFillColor(kWhite);

      projB[ptbin] = (TH1D*)fB->Get(Form("projDelPhi_%i",ptbin));
      projC[ptbin] = (TH1D*)fC->Get(Form("projDelPhi_%i",ptbin));
      projData0[ptbin]= (TH1D*)fD->Get(Form("NPEhDelPhi_0_%i",ptbin));
      projData2[ptbin]= (TH1D*)fD->Get(Form("NPEhDelPhi_2_%i",ptbin));
      Int_t RB = 8;
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
      fit->SetRangeX(46,56);                    // use only the first 15 bins in the fit
      Int_t status = fit->Fit();               // perform the fit
      std::cout << "fit status: " << status << std::endl;
      if (status == 0) {                       // check on fit status
	TH1F* result = (TH1F*) fit->GetPlot();
	projData0[ptbin]->Draw("Ep");
	result->Draw("same");
	}

      
      fitResult2->cd(ptbin+1);
      TFractionFitter* fit2 = new TFractionFitter(projData2[ptbin], mc); // initialise
      fit2->Constrain(0,0.0,1.0);              // constrain fraction 0
      fit2->Constrain(1,0.0,1.0);              // constrain fraction 1 to be between 0 and 1
      fit2->SetRangeX(23,29);                  // use only the first 15 bins in the fit
      Int_t status2 = fit2->Fit();             // perform the fit
      std::cout << "fit status: " << status2 << std::endl;
      if (status2 == 0) {                       // check on fit status
	TH1F* result2 = (TH1F*) fit2->GetPlot();
	projData2[ptbin]->Draw("Ep");
	result2->Draw("same");
      }
    }
}
