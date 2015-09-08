// Offline Plots - Z. Miller July 24, 2015
//
// .L offline.C
// offline("FILENAME") # Without .root Extension
// takes output of offline.C (pythia version) as inputs.
// Requires input from both a C and B data run with same name

// Plots templates for "currentC" and "currentB" files,
// after offline has properly weighted them.

#include "anaConst.h"

void plotTemplates()
{
   char name[1000];
   sprintf(name,"/Users/zach/Research/pythia/ptHatTemplate/outputs/currentB.root");
   TFile *fB = new TFile(name,"READ");
   sprintf(name,"/Users/zach/Research/pythia/ptHatTemplate/outputs/currentC.root");
   TFile *fC = new TFile(name,"READ");
   if (fB->IsOpen()==kFALSE || fC->IsOpen()==kFALSE)
     { std::cout << "!!!!!! Either B or C File not found !!!!!!" << std::endl
		 << "Looking for currentC.root and currentB.root." << std::endl;
	 exit(1); }
   
   // Set constants and projection bins
   const Int_t numPtBins = anaConst::nPtBins;
   Float_t lowpt[numPtBins],highpt[numPtBins];
   for(Int_t c=0; c< numPtBins; c++){
     lowpt[c] = anaConst::lpt[c];
     highpt[c] = anaConst::hpt[c];
   }
  Float_t hptCut=anaConst::hptCut;
   
   // Make Canvases
   TCanvas* deltaPhi = new TCanvas("deltaPhi","Pythia Delta Phi",150,0,1150,1000);
   deltaPhi->Divide(4,3);
   TCanvas* singlePlot = new TCanvas("deltaPhiSP","Pythia Delta Phi",150,0,1150,1000);
   
   // Make histos
   TH1D* projB[numPtBins];
   TH1D* projC[numPtBins];
   TH1F* bPtNorms;
   TH1F* cPtNorms;
   Float_t norm0,norm2,normB,normC;

   // Get ptbin independent hists
   bPtNorms   = (TH1F*)fB->Get("ptNorm");
   cPtNorms   = (TH1F*)fC->Get("ptNorm");


  // Get and Draw histos
  TPaveText* lbl[numPtBins];
  char textLabel[100];

  for(Int_t ptbin=0; ptbin<numPtBins; ptbin++)
    {
      // Init necessary plotting tools
      lbl[ptbin] = new TPaveText(.2,.76,.5,.82,Form("NB NDC%i",ptbin));
      sprintf(textLabel,"%.1f < P_{T,e} < %.1f",lowpt[ptbin],highpt[ptbin]);
      lbl[ptbin]->AddText(textLabel);
      lbl[ptbin]->SetFillColor(kWhite);

      projB[ptbin] = (TH1D*)fB->Get(Form("delPhi_%i",ptbin));
      projC[ptbin] = (TH1D*)fC->Get(Form("delPhi_%i",ptbin));

      // Get Normalizations
      normB = bPtNorms->GetBinContent(bPtNorms->GetBin(ptbin+1));
      normC = cPtNorms->GetBinContent(cPtNorms->GetBin(ptbin+1));
      //projB[ptbin] ->Scale(1./normB);
      //projC[ptbin] ->Scale(1./normC);
      
      deltaPhi->cd(ptbin+1);
      projB[ptbin]->SetLineColor(kRed);
      projC[ptbin]->SetLineColor(kBlack);
      //projC[ptbin]->GetYaxis()->SetRangeUser(0.,1.5);
      projB[ptbin]->Draw();
      projC[ptbin]->Draw("same");
      lbl[ptbin]  ->Draw("same");

      TLegend* leg = new TLegend(0.5,0.73,0.85,0.85);
      leg->AddEntry(projB[ptbin],"b#bar{b}->NPE","lpe");
      leg->AddEntry(projC[ptbin],"c#bar{c}->NPE","lpe");
      leg->Draw();

      if(ptbin == 1)
	{
	  singlePlot->cd();
	  projC[ptbin]->Draw();
	  projB[ptbin]->Draw("same");
	  lbl[ptbin]  ->Draw("same");
	  leg->Draw("same");
	}
    }

}
