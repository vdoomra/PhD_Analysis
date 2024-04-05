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
#include "/direct/phenix+u/vdoomra/install/include/DileptonAnalysis/MyEvent.h"

using namespace std;
using namespace DileptonAnalysis;

const double pi = TMath::ACos(-1);
const double step_size = 0.01;
const int total_nlayers = 8;
const int nlayers = 4;
const int ncharge = 2;
const int narms = 2;

TH1D* hist = NULL;

const float radii[total_nlayers] = { 2.64, 5.22, 10.41, 11.80, 12.90, 15.55, 16.74, 17.91};

double get_mean_phi_data(int layer, int charge, int arm)
{

  double dilep_mean_phi[nlayers][ncharge][narms] = {-999};

  dilep_mean_phi[0][0][0] = -0.000724;
  dilep_mean_phi[0][0][1] = -0.001465;
  dilep_mean_phi[0][1][0] =  0.000609;
  dilep_mean_phi[0][1][1] =  0.001143;

  dilep_mean_phi[1][0][0] =  -0.000603;
  dilep_mean_phi[1][0][1] =  -0.001290;
  dilep_mean_phi[1][1][0] =   0.000552;
  dilep_mean_phi[1][1][1] =   0.001375;

  dilep_mean_phi[2][0][0] =   -0.000699;
  dilep_mean_phi[2][0][1] =   -0.001283;
  dilep_mean_phi[2][1][0] =    0.000876;
  dilep_mean_phi[2][1][1] =    0.001479;

  dilep_mean_phi[3][0][0] =   -0.000876;
  dilep_mean_phi[3][0][1] =   -0.001094;
  dilep_mean_phi[3][1][0] =   0.000933;
  dilep_mean_phi[3][1][1] =   0.001498;

  return dilep_mean_phi[layer][charge][arm];

}

double get_sigma_phi_data(int layer, int arm, double pt)
{
  double dilep_par0_sigma_phi[nlayers][narms] = {-999};
  double dilep_par1_sigma_phi[nlayers][narms] = {-999};
  double dilep_par2_sigma_phi[nlayers][narms] = {-999};

  // Parameter 0 for eveluating sigma phi

  dilep_par0_sigma_phi[0][0] =  0.00568;
  dilep_par0_sigma_phi[0][1] =  0.00669;

  dilep_par0_sigma_phi[1][0] =  0.00429;
  dilep_par0_sigma_phi[1][1] =  0.00466;

  dilep_par0_sigma_phi[2][0] =  0.00355;
  dilep_par0_sigma_phi[2][1] =  0.00412;

  dilep_par0_sigma_phi[3][0] =  0.00340;
  dilep_par0_sigma_phi[3][1] =  0.00365;

  // Parameter 1 for eveluating sigma phi

  dilep_par1_sigma_phi[0][0] =  0.02474;
  dilep_par1_sigma_phi[0][1] =  0.01965;

  dilep_par1_sigma_phi[1][0] =  0.02431;
  dilep_par1_sigma_phi[1][1] =  0.01777;

  dilep_par1_sigma_phi[2][0] =  0.02238;
  dilep_par1_sigma_phi[2][1] =  0.01578;

  dilep_par1_sigma_phi[3][0] =  0.02234;
  dilep_par1_sigma_phi[3][1] =  0.01350;

  // Parameter 2 for eveluating sigma phi

  dilep_par2_sigma_phi[0][0] =  -2.59428;
  dilep_par2_sigma_phi[0][1] =  -2.32567;

  dilep_par2_sigma_phi[1][0] =  -2.46073;
  dilep_par2_sigma_phi[1][1] =  -1.91219;

  dilep_par2_sigma_phi[2][0] =  -2.30089;
  dilep_par2_sigma_phi[2][1] =  -1.74849;

  dilep_par2_sigma_phi[3][0] =  -2.57653;
  dilep_par2_sigma_phi[3][1] =  -1.69776;

  return dilep_par0_sigma_phi[layer][arm] + dilep_par1_sigma_phi[layer][arm]*exp(dilep_par2_sigma_phi[layer][arm] * pt);


}



void tree_QA(const char* inFile, const char* outFile){

  TFile* input = new TFile(inFile,"READ");
  if(!(input))
  {
    cout << "no input file" << endl;
    exit(1);
  }

  TTree* T = (TTree*)input->Get("T");
  TBranch* br = T->GetBranch("MyEvent");
  MyEvent* event = 0;
  br->SetAddress(&event);

  cout << "Trees read!" << endl; 

  TFile* output = new TFile(outFile,"RECREATE");
 
  int nevt = T->GetEntries();

  hist = new TH1D("hist","Debug Histogram", 12, -0.5, 11.5);

  TFile* fmag = new TFile("field_map.root","READ");
  TH2F* field_map_bz = (TH2F*)fmag->Get("hist_bz");
  TH2F* field_map_br = (TH2F*)fmag->Get("hist_br");

  for (int ievent = 0; ievent < nevt; ievent++)
  {
    if (ievent%50000==0) cout << "Event: " << ievent << " / " << nevt << endl;

    event->ClearEvent();
    br->GetEntry(ievent);

    int ntrack = event->GetNtrack();
    int nvtxhits = event->GetNVTXhit();

    hist->Fill(0);

    int nelectrons = 0;
    int nhadrons = 0;
    int nhadrons_onehit = 0;
    int nhadrons_twohit = 0;
    int nelectrons_onehit = 0;
    int nelectrons_twohit = 0;

    for(int itrk = 0; itrk < ntrack; itrk++){

      bool electron = true;
      bool hadron = true;

      MyTrack mytrk = event->GetEntry(itrk);
      double ecore = mytrk.GetEcore();
      double px = mytrk.GetPx();
      double py = mytrk.GetPy();
      double pz = mytrk.GetPz();
      double p = sqrt(px*px + py*py + pz*pz);

      if ( mytrk.GetN0() < 1 ) electron = false; 
      if ( fabs(mytrk.GetEmcdphi()) >= 0.05 ) electron = false;
      if ( fabs(mytrk.GetEmcdz()) >= 20 ) electron = false;
      if ( ecore/p <= 0.5 ) electron = false;
      if ( mytrk.GetProb() <= 0.01 ) electron = false;

		  // Hadron Identification Cuts

      if ( mytrk.GetN0() >= 0 ) hadron = false; 
      if ( ecore > 0.3) hadron = false; 

      if(electron) nelectrons++;
      if(hadron) nhadrons++;

      double phi0_trk_proj = mytrk.GetPhi0Prime();
      double the0_trk_proj = mytrk.GetThe0Prime();
	    pz = mytrk.GetPtPrime()*(TMath::Cos(mytrk.GetThe0Prime()))/(TMath::Sin(mytrk.GetThe0Prime()));

      double rp = sqrt(  pow(event->GetPreciseX(), 2) + pow(event->GetPreciseY(), 2) );
      double zp = event->GetPreciseZ();

      double dilep_phi_projection[total_nlayers] = {-999};
      double dilep_the_projection[total_nlayers] = {-999};
      double dilep_hit_index[total_nlayers] = {-999};
      double dilep_hit_counter[total_nlayers] = {-999};

      int charge_index = -999;

      for(int p=1; p<1800; p++){

		    for(int l=0; l < total_nlayers; l++){ 

          if( fabs(rp - radii[l]) < step_size && dilep_phi_projection[l] == -999 ){ dilep_phi_projection[l] = phi0_trk_proj; dilep_the_projection[l] = the0_trk_proj; } 

        }

        int rbin = field_map_bz->GetXaxis()->FindBin(rp);
        int zbin = field_map_bz->GetYaxis()->FindBin(zp);

        double bz = field_map_bz->GetBinContent(rbin, zbin)/10000; // Conversion from Gaus to Tesla

        double delta_phi0 = (mytrk.GetChargePrime()*0.3*step_size*bz)/(2*mytrk.GetPtPrime()*100);
        phi0_trk_proj += delta_phi0;

        double bradial = field_map_br->GetBinContent(rbin, zbin)/10000; // Conversion from Gaus to Tesla

        //Bend in the z direction does not depend upon the charge.

        double delta_the0 = 0.3*bradial*(step_size*TMath::Tan(pi/2 - the0_trk_proj))/(2*pz*100);
        if( mytrk.GetThe0Prime() > pi/2) the0_trk_proj -= delta_the0;
        else the0_trk_proj += delta_the0;

        zp += step_size*TMath::Tan(pi/2 - the0_trk_proj);
        rp += step_size;


      }

      if( mytrk.GetChargePrime() > 0 ) charge_index = mytrk.GetChargePrime();
      else charge_index = mytrk.GetChargePrime() + 1;
	
      int DCArm = mytrk.GetArm();
      TVector3 hitpoint;

      for(int ilayer = 0; ilayer < total_nlayers; ilayer++){

		    double min = -999;  int hit_flag_counter = 0;
                    
        for(int ihit = 0; ihit < nvtxhits ; ihit++){

          hitpoint.SetXYZ(0, 0, 0);
          MyVTXHit svxhit = event->GetVTXHitEntry(ihit);
          int layer = svxhit.GetLayer();
          double xhit = svxhit.GetXHit();
          double yhit = svxhit.GetYHit();
          double zhit = svxhit.GetZHit();
          double rr = sqrt(xhit*xhit + yhit*yhit);
          bool found = false;

          if( layer < 0 || layer > 3 ) continue;
          if( (ilayer == 0 || ilayer == 1) && layer != ilayer ) continue;
          if( !(rr > 10.20 && rr < 10.70) && ilayer == 2) continue;
          if( !(rr > 11.60 && rr < 12.00) && ilayer == 3) continue;
          if( !(rr > 12.70 && rr < 13.40) && ilayer == 4) continue;
          if( !(rr > 15.20 && rr < 15.80) && ilayer == 5) continue;
          if( !(rr > 16.50 && rr < 17.00) && ilayer == 6) continue;
          if( !(rr > 17.60 && rr < 18.20) && ilayer == 7) continue;
                        
          hitpoint.SetXYZ(xhit-event->GetPreciseX(), yhit-event->GetPreciseY(), zhit-event->GetPreciseZ());
          double phi_hit = hitpoint.Phi();
          if(phi_hit < -pi/2) phi_hit += 2*pi;
          double theta_hit = hitpoint.Theta();

          // Now defining our initial search windows

          double sigma_phi_value = -999;   double mean_phi_value = -999;

          sigma_phi_value = get_sigma_phi_data(layer, DCArm, mytrk.GetPtPrime());
          mean_phi_value = get_mean_phi_data(layer, charge_index, DCArm);

		      if((dilep_phi_projection[ilayer] - phi_hit) > (mean_phi_value - 10*sigma_phi_value) &&  (dilep_phi_projection[ilayer] - phi_hit) < (mean_phi_value + 10*sigma_phi_value)){ found = true; }

          if(found){ 
                
            hit_flag_counter++;
			      double diff = sqrt( pow( (dilep_phi_projection[ilayer] - phi_hit)/sigma_phi_value, 2) );

            if(hit_flag_counter == 1) { min = diff; dilep_hit_index[ilayer] = ihit; dilep_hit_counter[ilayer] = hit_flag_counter; }
            else{  if(diff < min){  min = diff; dilep_hit_index[ilayer] = ihit; dilep_hit_counter[ilayer] = hit_flag_counter; }   }  

          }
        }
      }

      int counter = 0;

      for(int i = 0; i < total_nlayers; i++){

        if(dilep_hit_index[i] == -999) continue;
        else counter++;

      }

      if(hadron && counter > 0) nhadrons_onehit++;
      if(hadron && counter > 1) nhadrons_twohit++;
      if(electron && counter > 0) nelectrons_onehit++;
      if(electron && counter > 1) nelectrons_twohit++;


    }

    if(nelectrons == 0) hist->Fill(1);
    if(nelectrons == 1) hist->Fill(2);
    if(nelectrons > 1) hist->Fill(3);

    if(nhadrons == 0) hist->Fill(4);
    if(nhadrons == 1) hist->Fill(5);
    if(nhadrons > 1) hist->Fill(6);

    if(nhadrons_onehit > 0) hist->Fill(7);
    if(nhadrons_twohit > 0) hist->Fill(8);

    if(nelectrons_onehit > 0) hist->Fill(9);
    if(nelectrons_twohit > 0) hist->Fill(10);

    if(nelectrons > 1 && nhadrons > 0) hist->Fill(11);

  }

  output->cd();
  hist->Write();
  output->Close();


}
