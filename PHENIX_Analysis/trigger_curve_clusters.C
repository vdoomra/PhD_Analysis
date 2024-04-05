#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TCanvas.h>
#include "/direct/phenix+u/vdoomra/install/include/DileptonAnalysis/MyEvent.h"

//At the vertex reconstruction module level I only store good quality electron tracks with a minimum pT of 400 MeV/c.

using namespace std;
using namespace DileptonAnalysis;
const int centbin = 1;
const double pi = TMath::ACos(-1);

const int nsect = 8;
const char* sectNames[nsect] = {"E4", "E3", "E2", "E1", "W4", "W3", "W2", "W1"};

TH1D* ecore_spectra_MB[nsect] = {NULL};
TH1D* ecore_spectra_ERT[nsect] = {NULL};

void trigger_curve_clusters(const char* inFileMB, const char* inFileERT, const char* outFile)
{

  for(int isect = 0; isect < nsect; isect++){

    ecore_spectra_MB[isect] = new TH1D(Form("ecore_spectra_mb_%s", sectNames[isect]),"Ecore Distribution",20,0,6);
    ecore_spectra_MB[isect]->Sumw2();

    ecore_spectra_ERT[isect] = new TH1D(Form("ecore_spectra_ert_%s", sectNames[isect]),"Ecore Disribution",20,0,6);
    ecore_spectra_ERT[isect]->Sumw2();
  }

  TFile* inputMB = new TFile(inFileMB,"READ");
  if(!(inputMB))
  {
    cout << "no input MB file" << endl;
    exit(1);
  }

  TFile* inputERT = new TFile(inFileERT,"READ");
  if(!(inputERT))
  {
    cout << "no input ERT file" << endl;
    exit(1);
  }

  TTree* T_MB = (TTree*)inputMB->Get("T");
  TBranch* br_MB = T_MB->GetBranch("MyEvent");
  MyEvent* event_MB = 0;
  br_MB->SetAddress(&event_MB);

  TTree* T_ERT = (TTree*)inputERT->Get("T");
  TBranch* br_ERT = T_ERT->GetBranch("MyEvent");
  MyEvent* event_ERT = 0;
  br_ERT->SetAddress(&event_ERT);

  cout << "Trees read!" << endl; 

  TFile* output = new TFile(outFile,"RECREATE");
 
  int nevt_MB = T_MB->GetEntries();
  int nevt_ERT = T_ERT->GetEntries();

  cout << "=====Minimum Bias=========" << endl;

  for (int ievent_MB = 0; ievent_MB < nevt_MB; ievent_MB++)
  {
    if (ievent_MB%50000==0) cout << "Event: " << ievent_MB << " / " << nevt_MB << endl;

    event_MB->ClearEvent();
    br_MB->GetEntry(ievent_MB);

    if(fabs(event_MB->GetPreciseZ()) > 10) continue;
    int ncluster_MB = event_MB->GetNcluster();

    for (int iclust_MB = 0; iclust_MB < ncluster_MB; iclust_MB++)
    {
      MyCluster myclust_MB = event_MB->GetClusterEntry(iclust_MB);

      double ecore_MB = myclust_MB.GetEcore();
      if(ecore_MB < 0.3) continue;
      int sect_MB = myclust_MB.GetSect();
      if(sect_MB < 0) continue;
      
      ecore_spectra_MB[sect_MB]->Fill(ecore_MB);

    }
  }

  cout << "==========ERT============" << endl;
  
  for (int ievent_ERT = 0; ievent_ERT < nevt_ERT; ievent_ERT++)
  {
    if (ievent_ERT%50000==0) cout << "Event: " << ievent_ERT << " / " << nevt_ERT << endl;

    event_ERT->ClearEvent();
    br_ERT->GetEntry(ievent_ERT);

    if(fabs(event_ERT->GetPreciseZ()) > 10) continue;
    int ncluster_ERT = event_ERT->GetNcluster();

    for (int iclust_ERT = 0; iclust_ERT < ncluster_ERT; iclust_ERT++)
    {
      MyCluster myclust_ERT = event_ERT->GetEntry(iclust_ERT);

      double ecore_ERT = myclust_ERT.GetEcore();
      if(ecore_ERT < 0.3) continue;
      int sect_ERT = myclust_ERT.GetSect();
      int firedERT = myclust_ERT.GetisERT();
      if(sect_ERT<0)continue;
      if(!firedERT) continue;

      ecore_spectra_ERT[sect_ERT]->Fill(myclust_ERT.GetEcore());

    }
  }


  TCanvas* c1 = new TCanvas();
  c1->Divide(4,2);

  for(int isect = 0; isect<nsect; isect++){  

    int bin_lo = ecore_spectra_ERT[isect]->FindBin(3.5);
    int bin_hi = ecore_spectra_ERT[isect]->FindBin(6.0);

    if(isect == 6 || isect == 7) {  

      bin_lo = ecore_spectra_ERT[isect]->FindBin(4.0);
      bin_hi = ecore_spectra_ERT[isect]->FindBin(6.0);

    }

    double scale = ecore_spectra_MB[isect]->Integral(bin_lo, bin_hi)/ecore_spectra_ERT[isect]->Integral(bin_lo, bin_hi);
    ecore_spectra_ERT[isect]->Scale(scale);

    c1->cd(isect+1);
    gPad->SetLogy();
    ecore_spectra_MB[isect]->Draw();
    ecore_spectra_ERT[isect]->Draw("same");

  }

  TCanvas* c2 = new TCanvas();
  c2->Divide(4,2);

  for(int isect = 0; isect<nsect; isect++){  

    TH1D* hist = (TH1D*)ecore_spectra_ERT[isect]->Clone();
    hist->SetName(Form("trigger_turn_on_sector_%s",sectNames[isect]));

    hist->Divide(ecore_spectra_MB[isect]);

    c2->cd(isect+1);
    hist->Draw();

    hist->Reset("ICESM");

  }

  

  output->cd();
  c1->Write();
  c2->Write();
  output->Close();

} 
