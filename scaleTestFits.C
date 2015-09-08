// Fraction Fit - Z. Miller Sep 1, 2015
//
// .L scaleTestFits.C
// scaleTestFits.C() 
// takes output of offline.C (pythia version, and data version) as inputs.
// Copy current best data and templates as "current%s.root" {B,C,Data}
// Loops through various scales for MC, plot fit error and result for each scale.

#include "anaConst.h"

Bool_t checkMakePDF();
char FileName[100];

void scaleTestFits(Float_t maxScale=10000, Float_t increment=1000)
{
  
  gStyle->SetOptFit(1111);

   Bool_t makePDF = checkMakePDF();
  
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
  
  // Set constants and projection bins (from header file anaConst, analysis constants)
  const Int_t numPtBins = anaConst::nPtBins;
  Float_t lowpt[numPtBins],highpt[numPtBins];
  for(Int_t c=0; c< numPtBins; c++){
    lowpt[c] = anaConst::lpt[c];
    highpt[c] = anaConst::hpt[c];
  }
  Float_t hptCut=anaConst::hptCut;
  Int_t numScale = maxScale/increment + 1;
  Double_t p00[numPtBins],p01[numPtBins],p20[numPtBins],p21[numPtBins];
  Double_t e00[numPtBins],e01[numPtBins],e20[numPtBins],e21[numPtBins];
  Double_t pC0[numPtBins],pC1[numPtBins],eC0[numPtBins],eC1[numPtBins];
  Double_t Rb0[numPtBins][numScale],Rb2[numPtBins][numScale],RbC[numPtBins][numScale],x[numPtBins][numScale];
  Double_t eb0[numPtBins][numScale],eb2[numPtBins][numScale],ebC[numPtBins][numScale],dx[numPtBins][numScale];
  Double_t ptOFF1[numPtBins],ptOFF2[numPtBins];
  Int_t rangeLow  = 86;  //22-29 for near-side only
  Int_t rangeHigh = 116; //18-33 for -pi,pi
  Int_t plotCount0 = 0, plotCount2 = 0, plotCount = 0;
  
  // Make Canvases
  TCanvas* errors  = new TCanvas("errors","Error vs Scale Factor",150,0,1150,1000);
  TCanvas* results = new TCanvas("results","Result vs Scale Factor",150,0,1150,1000);
  errors  -> Divide(3,2);
  results -> Divide(3,2);

  // Make histos
  TH1D* projB[numPtBins];
  TH1D* projC[numPtBins];
  TH1D* projData0[numPtBins];
  TH1D* projData2[numPtBins];
  TH1D* combData[numPtBins];
  TH1D* plotD0[numPtBins];
  TH1D* plotD2[numPtBins];
  TH1D* plotC[numPtBins];
  TH1D* plotB[numPtBins];
  TH1F* histoNorms;
  TH1F* bPtNorms;
  TH1F* cPtNorms;
  
  // Get and Draw histos
  TPaveText* lbl[numPtBins];
  TPaveText* stat[3][numPtBins];
  char statLabel[100];
  char textLabel[100];
  Int_t plotbin;
  Float_t norm0,norm2,normB,normC;

  // Get ptbin independent hists
  histoNorms = (TH1F*)fD->Get("histoNorms");
  bPtNorms   = (TH1F*)fB->Get("ptNorm");
  cPtNorms   = (TH1F*)fC->Get("ptNorm");

  for(Int_t ptbin=0; ptbin<6; ptbin++)
    {
      plotCount=0;
      for(Float_t SCALE=100; SCALE < maxScale; SCALE += increment)
	{
     	  norm0 = histoNorms->GetBinContent(histoNorms->GetBin(1,ptbin+1));
	  norm2 = histoNorms->GetBinContent(histoNorms->GetBin(3,ptbin+1));
	  normB = bPtNorms->GetBinContent(bPtNorms->GetBin(ptbin+1));
	  normC = cPtNorms->GetBinContent(cPtNorms->GetBin(ptbin+1));
	  
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

	  // Do any rebinning
	  Int_t RB = 1;
	  projB[ptbin]->Rebin(RB);
	  projC[ptbin]->Rebin(RB);
	  projData0[ptbin]->Rebin(RB);
	  projData2[ptbin]->Rebin(RB);
	  
	  // Clone to make plots without effecting fits
	  plotD0[ptbin] = (TH1D*) projData0[ptbin]->Clone();
	  plotD2[ptbin] = (TH1D*) projData2[ptbin]->Clone();
	  plotB[ptbin]  = (TH1D*) projB[ptbin]->Clone();
	  plotC[ptbin]  = (TH1D*) projC[ptbin]->Clone();
	  
	  // Set features that are the same in plots
	  projData0[ptbin]->SetLineColor(kBlue);
	  projData2[ptbin]->SetLineColor(kGreen+3);
	  projB[ptbin]->SetLineColor(kRed);
	  projC[ptbin]->SetLineColor(kBlack);
	  projC[ptbin]->GetXaxis()->SetRangeUser(-3.5,3.5);
	  plotD0[ptbin]->SetLineColor(kBlue);
	  plotD2[ptbin]->SetLineColor(kGreen+3);
	  plotB[ptbin]->SetLineColor(kRed);
	  plotC[ptbin]->SetLineColor(kBlack);
	  plotC[ptbin]->GetXaxis()->SetRangeUser(-3.5,3.5);
	  
	  // Normalize
	  Float_t nScale = SCALE;
	  projB[ptbin] ->Scale(nScale/normB);
	  projC[ptbin] ->Scale(nScale/normC);
	  plotD0[ptbin]->Scale(1/norm0);
	  plotD2[ptbin]->Scale(1/norm2);
	  plotB[ptbin] ->Scale(1/normB);
	  plotC[ptbin] ->Scale(1/normC);
	  
	  combData[ptbin] = (TH1D*) projData0[ptbin]->Clone();
	  combData[ptbin] -> Add(projData2[ptbin]);

	  // Do the actual fit
	  TObjArray *mc = new TObjArray(2);   // MC histograms are put in this array
	  mc->Add(projC[ptbin]);
	  mc->Add(projB[ptbin]);
	  
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
	    resultC->SetLineColor(kBlack);
	    resultC->SetFillColor(kWhite);
	    fitC->GetResult(0,pC0[ptbin],eC0[ptbin]);
	    fitC->GetResult(1,pC1[ptbin],eC1[ptbin]);
	  }
	  
	  x[ptbin][plotCount] = SCALE;
	  dx[ptbin][plotCount] = 0;
	  RbC[ptbin][plotCount] = pC1[ptbin]/(pC1[ptbin]+pC0[ptbin]);
	  ebC[ptbin][plotCount] = eC1[ptbin];
	  plotCount++;
	}
    
      TGraph *grCe    = new TGraph(plotCount-1,x[ptbin],ebC[ptbin]);
      TGraph *grCr    = new TGraph(plotCount-1,x[ptbin],RbC[ptbin]);

      errors->cd(ptbin+1);
      grCe->SetLineColor(kBlack);
      grCe->SetMarkerStyle(20);
      grCe->SetMarkerColor(kBlack);
      grCe->GetYaxis()->SetRangeUser(0,1);
      grCe->GetYaxis()->SetTitle("Fitting Error");
      grCe->GetXaxis()->SetTitle("MC Scale Factor");
      grCe->SetTitle("Error vs Scale Factor");
      grCe->Draw("AP");
      lbl[ptbin]->Draw("same");
      
      results->cd(ptbin+1);
      grCr->SetLineColor(kBlack);
      grCr->SetMarkerStyle(20);
      grCr->SetMarkerColor(kBlack);
      grCr->GetYaxis()->SetRangeUser(0,1);
      grCr->GetYaxis()->SetTitle("Fit Result");
      grCr->GetXaxis()->SetTitle("MC Scale Factor");
      grCr->SetTitle("Fit vs Scale Factor");
      grCr->Draw("AP");
      lbl[ptbin]->Draw("same");
    }
  
  // Make PDF with output canvases
  if(makePDF)
    {
      //Set front page
      TCanvas* fp = new TCanvas("fp","Front Page",100,0,1000,900);
      fp->cd();
      TBox *bLabel = new TBox(0.01, 0.88, 0.99, 0.99);
      bLabel->SetFillColor(38);
      bLabel->Draw();
      TLatex tl;
      tl.SetNDC();
      tl.SetTextColor(kWhite);
      tl.SetTextSize(0.033);
      char tlName[100];
      char tlName2[100];
      
      TString titlename = FileName;
      int found = titlename.Last('/');
      if(found >= 0){
	titlename.Replace(0, found+1, "");
      } 
      sprintf(tlName, "RUN 12 NPE-h   #Delta#phi Correlations");
      tl.SetTextSize(0.05);
      tl.SetTextColor(kWhite);
      tl.DrawLatex(0.05, 0.92,tlName);
      
      TBox *bFoot = new TBox(0.01, 0.01, 0.99, 0.12);
      bFoot->SetFillColor(38);
      bFoot->Draw();
      tl.SetTextColor(kWhite);
      tl.SetTextSize(0.05);
      tl.DrawLatex(0.05, 0.05, (new TDatime())->AsString());
      tl.SetTextColor(kBlack);
      tl.SetTextSize(0.03);
      tl.DrawLatex(0.1, 0.14, titlename);
      sprintf(tlName,"TEST");
      tl.DrawLatex(0.1, 0.8,tlName);
      
      // Place canvases in order
      TCanvas* temp = new TCanvas();
      sprintf(name, "FFOutput/scaleTest/%s.pdf[", FileName);
      temp->Print(name);
      sprintf(name, "FFOutput/scaleTest/%s.pdf", FileName);

      temp = errors; 
      temp->Print(name);
      temp = results;
      temp->Print(name);
            
      sprintf(name, "FFOutput/scaleTest/%s.pdf]", FileName);
      temp->Print(name);
    }
}

Bool_t checkMakePDF(){

  // Set option for pdf creation
  Int_t number = 2; Bool_t fmakePDF = kTRUE;
  while(number > 1 || number < 0){
    std::cout << "Make PDF? [default: 1]: ";
    std::string input;
    std::getline( std::cin, input );
    if ( !input.empty() ){
      std::istringstream stream( input );
      stream >> number;
      if(number == 0)
	fmakePDF = kFALSE;
      if(number == 1)
	fmakePDF = kTRUE;
    }
    else
      number = 1; 
  }
  if(fmakePDF) // need a file name if making pdf
    {
      cout << "Need FileName (no ext.): ";
      std::string input2;
      std::getline( std::cin, input2 );
      if ( !input2.empty() ){
	std::istringstream stream2( input2 );
	string s = stream2.str();
	sprintf(FileName,"%s",s.c_str());
      }
      else
	{
	  sprintf(FileName, "scaleTest");
	}
    }

  return fmakePDF;
}
