// Offline Plots - Z. Miller Sep 1, 2015 
// Updated to take ptHatBin Pythia Output
//
// .L offline.C
// offline("FILENAME", mode) # Without .root Extension
// modes: 0 = no ID about type of input, 1 = c/cbar, 2 = b/bbar

#include "anaConst.h"

// Declare functions
void checkBatchMode();
Bool_t checkMakePDF();
Bool_t checkMakeRoot();


void offline(const char* FileName="test", Int_t mode = 0)
{
  if (strcmp(FileName, "") == 0 || mode == 0 || mode > 2)
    {
      cout << "Error in input of offline('fileName',mode):" << endl
	   << "mode 1: c/cbar; mode 2: b/bbar." << endl
	   << "Need File Name: ''pythia_tree_Aug##_#''" << endl;
      abort();
    }
  
  // Set Style parameters for this macro
  //gStyle->SetOptTitle(1); // Show Title (off by default for cleanliness)
  gErrorIgnoreLevel = kError; // Set Verbosity Level (kPrint shows all)

   // Set Output options
  Int_t number;
  checkBatchMode();
  Bool_t makePDF = checkMakePDF();
  Bool_t makeROOT= checkMakeRoot();

  // Use mode input to decide whether C or B templates to work on
  char type[10] = "X";
  if(mode == 1)
    sprintf(type, "C");
  if(mode == 2)
    sprintf(type, "B");

  // Open output file
  char fname[100];
  TFile* file;
  if(makeROOT){
    sprintf(fname,"/Users/zach/Research/pythia/ptHatTemplate/%s_%s_processed.root",FileName,type);
    file = new TFile(fname,"RECREATE");
    if (file->IsOpen()==kFALSE)
      {
	std::cout << "!!! Outfile Not Opened !!!" << std::endl;
	makeROOT = kFALSE;
      }
  }

  // Initialize Histos for Summing and other global vars
  const Int_t numPtHatBins = 8;
  const Int_t numPtBins = anaConst::nPtBins;
  Float_t lowpt[numPtBins],highpt[numPtBins];
  for(Int_t c=0; c< numPtBins; c++){
    lowpt[c] = anaConst::lpt[c];
    highpt[c] = anaConst::hpt[c];
  }
  Float_t hptCut=anaConst::hptCut;
  Float_t hptMax=25; // Set max above range to allow overflow

  TH1F* ptHat     = new TH1F("pThat", "" ,1500, 0, 150);
  TH1F* ptHatCorr = new TH1F("pThatCorrected", "" ,1500, 0, 150);
  TH3F* mh3delPhi;
  TH2F* mh2npePt;
  TH1F* hStats;
  TH2F* mh2ptHatPt;
  TH1D* projpthatall;
  char hist[100];
  TH1F* delPhi[numPtBins];
  TH1F* NpeY[numPtBins];
  TH1F* ptNorm;
  TH1D* projDelPhi[numPtBins];
  TH1D* projNpeY[numPtBins];
  TH1D* projptHat[numPtBins];
  for(Int_t ptbin=0; ptbin<numPtBins; ptbin++) // initialize all before the actual sorting
    { delPhi[ptbin]= new TH1F(Form("delPhi_%i",ptbin), "Delta Phi" ,200, -10, 10);
      delPhi[ptbin]->Sumw2();
      NpeY[ptbin] = new TH1F(Form("NpeY_%i",ptbin),"NpeY",60,-3,3);
    }
  ptNorm = new TH1F("ptNorm", "pT Norm" ,200, 0, 20);
  ptNorm ->Sumw2();
      
  Float_t totalNorm[numPtBins]={0.};
  Double_t wt=0.;
   
  Int_t pthatlow[numPtHatBins] = {0,1,2,4,8,16,32,64};
  Int_t pthathigh[numPtHatBins]= {1,2,4,8,16,32,64,128};

  // Make Canvases
  TCanvas* deltaPhi = new TCanvas("deltaPhi","Pythia Delta Phi",150,0,1150,1000);
  deltaPhi -> Divide(4,3);
  TCanvas* ptHatC = new TCanvas("ptHatC","ptHat Stitching Comparison",150,0,1150,1000);
  ptHatC   -> Divide(1,2);

  TPaveText* lbl[numPtBins];
  char textLabel[100];
  char name[1000];
 
  // Loop over all ptHat bins
  for(Int_t pthBin=0; pthBin < numPtHatBins; pthBin++)
    {
      
      // Open ROOT File (example: output/pythia_tree_Aug31_1_C2_4.root)
      sprintf(name,"/Users/zach/Research/pythia/ptHatTemplate/%s_%s%i_%i.root",FileName,type,pthatlow[pthBin],pthathigh[pthBin]); 
      TFile *f = new TFile(name,"READ");
      if (f->IsOpen()==kFALSE)
	{ std::cout << "!!! File Not Found !!!" << std::endl;
	  exit(1); }
      else
	{ cout << name << " is open!" << endl;}
            
      char histName[100];
      // Get Histos from run output
     
      sprintf(hist, "histo3D%s0", type);
      mh3delPhi    = (TH3F*)f->Get(hist);
      sprintf(hist, "histos2D%s1", type);
      mh2npePt     = (TH2F*)f->Get(hist);
      sprintf(hist, "histos2D%s10", type);
      mh2ptHatPt   = (TH2F*)f->Get(hist);
      sprintf(hist, "hStatistics");
      hStats       = (TH1F*)f->Get(hist);
          
      // Calculate Weight factors
      wt = 1e9*1e-3*(hStats->GetBinContent(1)/hStats->GetBinContent(2)); // Taken from Zhenyu's method. The 1e# factors are luminosity(?) corrections?
      projpthatall = mh2ptHatPt->ProjectionY("test",0,-1);
      ptHat -> Add(projpthatall);
      ptHatCorr -> Add(projpthatall,wt);
   
      // Analyze each ptH bin individually, adding to the overall hists
      for(Int_t ptbin=0; ptbin<numPtBins; ptbin++)
	{
	  // DEBUGcout << "pthbin: " << pthBin << " ptbin: " << ptbin << endl;
	  projDelPhi[ptbin] = mh3delPhi->ProjectionZ(Form("projDelPhi_%i",ptbin),mh3delPhi->GetXaxis()->FindBin(lowpt[ptbin]),mh3delPhi->GetXaxis()->FindBin(highpt[ptbin])-1,mh3delPhi->GetYaxis()->FindBin(hptCut),mh3delPhi->GetYaxis()->FindBin(hptMax));
	  projNpeY[ptbin]   = mh2npePt->ProjectionY(Form("projNpeY_%i",ptbin),mh2npePt->GetXaxis()->FindBin(lowpt[ptbin]),mh2npePt->GetXaxis()->FindBin(highpt[ptbin])-1);
	  projptHat[ptbin]  = mh2ptHatPt->ProjectionY(Form("projPtHat_%i",ptbin),mh2ptHatPt->GetXaxis()->FindBin(lowpt[ptbin]),mh2ptHatPt->GetXaxis()->FindBin(highpt[ptbin])-1);
	
	  delPhi[ptbin] -> Add(projDelPhi[ptbin],wt);
	  NpeY[ptbin] -> Add(projNpeY[ptbin],wt);
	  
	  // Calculate scaling Factor
	  Double_t Norm = NpeY[ptbin]->Integral();
	  ptNorm->SetBinContent(ptNorm->GetBin(ptbin+1),Norm);
	  totalNorm[ptbin] += Norm;
	 
	}
    }

  // For making plots

  ptHatC->cd(1);
  gPad-> SetLogy();
  ptHat->GetXaxis()->SetTitle("pT-Hat (GeV/c)");
  ptHat->SetTitle("Raw pT Hat");
  ptHat->Draw();
  ptHatC->cd(2);
  gPad-> SetLogy();
  ptHatCorr->GetXaxis()->SetTitle("pT-Hat (GeV/c)");
  ptHatCorr->SetTitle("Weighted pT Hat");
  ptHatCorr->Draw();
  
  for(Int_t ptbin=0; ptbin<numPtBins; ptbin++)
    {
      // Init necessary plotting tools
      lbl[ptbin] = new TPaveText(.2,.8,.5,.85,Form("NB NDC%i",ptbin));
      sprintf(textLabel,"%.1f < P_{T,e} < %.1f",lowpt[ptbin],highpt[ptbin]);
      lbl[ptbin]->AddText(textLabel);
      lbl[ptbin]->SetFillColor(kWhite);

      deltaPhi->cd(ptbin+1);
      delPhi[ptbin]->GetXaxis()->SetTitle("#Delta#phi_{eh}");
      delPhi[ptbin]->Sumw2();
      //cout << totalNorm[ptbin] << endl;
      //delPhi[ptbin]->Scale(wt);
      delPhi[ptbin]->GetYaxis()->SetTitle("1/N_{NPE} #upoint dN/d(#Delta)#phi");
      delPhi[ptbin]->GetXaxis()->SetRangeUser(-3.5,3.5);
      if(ptbin == 0)
	{
	  if(mode == 1)
	    delPhi[ptbin]->SetTitle("Pythia NPE-had #Delta#phi - c/#bar{c}");
	  if(mode == 2)
	    delPhi[ptbin]->SetTitle("Pythia NPE-had #Delta#phi - b/#bar{b}");
	}
      else
	delPhi[ptbin]->SetTitle("");
      if(ptbin < 13){
	delPhi[ptbin]->Draw("E");
	lbl[ptbin]->Draw("same");
      }
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
      sprintf(tlName, "RUN 12 NPE-h   #Delta#phi Pythia Templates");
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
      sprintf(name, "%s.pdf[", FileName);
      temp->Print(name);
      sprintf(name, "%s.pdf", FileName);
      temp = fp; // print front page
      temp->Print(name);
      temp = ptHatC;
      temp->Print(name);
      temp = deltaPhi;
      temp->Print(name);
      sprintf(name, "%s.pdf]", FileName);
      temp->Print(name);
    }

  if(makeROOT)
    {
      file->Write();
      file->Close();
    }
}

void checkBatchMode(){

 // sets batch mode, so don't draw canvas
  Int_t number = 2;
  while(number > 1 || number < 0){
    std::cout << "Batch Mode? [default: 1]: ";
    std::string input;
    std::getline( std::cin, input );
    if ( !input.empty() ) {
      std::istringstream stream( input );
      stream >> number;
      if(number == 0)
	gROOT->SetBatch(kFALSE);
      if(number == 1)
	gROOT->SetBatch(kTRUE);
    }
    else
      {
	number = 1;
	gROOT->SetBatch(kTRUE);
      }
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
  return fmakePDF;
}

Bool_t checkMakeRoot(){

  // Set option for .root creation
  Int_t number = 2; Bool_t fmakeROOT = kTRUE;
  while(number > 1 || number < 0){
    std::cout << "Make .root? [default: 1]: ";
    std::string input;
    std::getline( std::cin, input );
    if ( !input.empty() ){
      std::istringstream stream( input );
      stream >> number;
      if(number == 0)
	fmakeROOT = kFALSE;
      if(number == 1){
	fmakeROOT = kTRUE;
      }
    }
    else
      number = 1; 
  }
  return fmakeROOT;
}

