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
#include <TROOT.h>
#include "/phenix/plhf/vdoomra/dc_dead_map.h"

using namespace std;
using namespace DileptonAnalysis;
const double pi = TMath::ACos(-1);

const int nlayers = 4;
const int ncharge = 2;
const int narms = 2;
const int total_nlayers = 8;
const double step_size = 0.01;
const int rungroups = 10;
const double DCENTERCUT = 0.05; // RICH ghost cut 

int hit_index[total_nlayers] = {-999};

int nhits_trk = -999;
double theta_offset = -999;
int charge_index = -999;

double dilep_phi_projection[total_nlayers] = {-999};
double dilep_the_projection[total_nlayers] = {-999};
double dilep_hit_index[total_nlayers] = {-999};
double dilep_hit_counter[total_nlayers] = {-999};

double dilep_par0_theta[narms] = {-999}; 
double dilep_par1_theta[narms] = {-999}; 
double dilep_par2_theta[narms] = {-999};

const int dilep_run_lo[rungroups] = {421707, 423185, 425166, 427015, 429028, 430018, 430916, 431433, 432643, 435338}; 
const int dilep_run_hi[rungroups] = {423184, 425165, 427014, 429027, 430017, 430915, 431432, 432642, 435337, 436647}; 

const float dilep_par0_theta_east[rungroups] = { -0.021069, -0.014197, -0.011367, -0.007547, -0.013311, -0.021897, -0.035560, -0.031717, -0.013401, -0.032010};
const float dilep_par1_theta_east[rungroups] = { -0.009336, -0.008532, -0.009411, -0.009225, -0.007298, -0.007153, -0.008014, -0.008963, -0.008799, -0.009921};
const float dilep_par2_theta_east[rungroups] = { 0.017248, 0.010177, 0.007039, 0.003201, 0.008813, 0.016531, 0.030555, 0.027331, 0.009183, 0.028588};

const float dilep_par0_theta_west[rungroups] = { 0.039102, 0.012262, -0.002676, 0.016520, 0.000141, 0.032950, 0.031889, 0.012795, 0.017084, 0.008524};
const float dilep_par1_theta_west[rungroups] = { -0.008559, -0.009898, -0.009688, -0.010004, -0.008944, -0.008177, -0.007762, -0.009938, -0.011261, -0.012498};
const float dilep_par2_theta_west[rungroups] = { -0.031952, -0.005881, 0.009267, -0.009857, 0.006383, -0.025689, -0.025197, -0.006265, -0.010365, 0.000195};

TH1D* hist_diff_bbcz = NULL;
TH1D* hist_diff_vtxz = NULL;
TH1D* hist = NULL;
TH1D* hist_stag = NULL;
TH1D* debug = NULL;

const char* layer_label[nlayers] = {"B0", "B1"};
const char* arm_label[narms] = {"east", "west"};

TH2D* phi_z[narms][nlayers] = {NULL};

const int nbins_phi = 550;
const int nbins_z = 375;

int VTXDeadMap[narms][nlayers][nbins_phi][nbins_z] = {-999};

const float radii[total_nlayers] = { 2.64, 5.22, 10.41, 11.80, 12.90, 15.55, 16.74, 17.91};

bool single_cut(MyTrack mytrk)
{

  int arm = mytrk.GetArm();
  double alpha = mytrk.GetAlpha();
  double phi = mytrk.GetPhiDC();
  double zed = mytrk.GetZDC();
  double board = -999;
  int emcid = mytrk.GetEmcId();
  if(emcid<0) return false;
  if(arm == 0) board = (3.72402-phi+0.008047*cos(phi+0.87851))/0.01963496;
  else if(arm == 1) board = (0.573231+phi-0.0046*cos(phi+0.05721))/0.01963496;
  if(board < -99) return false; 
  bool dead_area = DCdeadArea(alpha, phi, board, zed); //Apply the DC Dead Map for Run15
  if(dead_area) return false;
  
  if (!(mytrk.GetTrkQuality() == 31 || mytrk.GetTrkQuality() == 51 || mytrk.GetTrkQuality() == 63) ) return false;
  if ( mytrk.GetN0() < 1 ) return false; 
  if ( fabs(mytrk.GetEmcdphi()) >= 0.05 ) return false;
  if ( fabs(mytrk.GetEmcdz()) >= 20 ) return false;
  if ( mytrk.GetProb() <= 0.01 ) return false;
  if ( mytrk.GetPtPrime() < 0.4 ) return false;
  if ( mytrk.GetDisp() >= 5) return false;
  if ( mytrk.GetChi2Npe0() >= 10) return false;

  return true;
}

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

double get_mean_theta_data(int layer, int charge, int arm)
{

  double dilep_mean_theta[nlayers][ncharge][narms] = {-999};

  dilep_mean_theta[0][0][0] = 0.001742;
  dilep_mean_theta[0][0][1] = 0.002142;
  dilep_mean_theta[0][1][0] = 0.002339;
  dilep_mean_theta[0][1][1] = 0.003140;

  dilep_mean_theta[1][0][0] =  0.001251;
  dilep_mean_theta[1][0][1] =  0.000626;
  dilep_mean_theta[1][1][0] =  0.001417;
  dilep_mean_theta[1][1][1] =  0.00140;

  dilep_mean_theta[2][0][0] =  0.000958;
  dilep_mean_theta[2][0][1] =  0.000369;
  dilep_mean_theta[2][1][0] =  0.001469;
  dilep_mean_theta[2][1][1] =  0.00117;

  dilep_mean_theta[3][0][0] =  0.000873;
  dilep_mean_theta[3][0][1] =  6.657e-5;
  dilep_mean_theta[3][1][0] =  0.001244;
  dilep_mean_theta[3][1][1] =  0.000805;

  return dilep_mean_theta[layer][charge][arm];

}

double get_sigma_theta_data(int layer, int arm, double pt)
{
  double dilep_par0_sigma_theta[nlayers][narms] = {-999};
  double dilep_par1_sigma_theta[nlayers][narms] = {-999};
  double dilep_par2_sigma_theta[nlayers][narms] = {-999};

  // Parameter 0 for eveluating sigma theta

  dilep_par0_sigma_theta[0][0] =  0.02952;
  dilep_par0_sigma_theta[0][1] =  0.02953;

  dilep_par0_sigma_theta[1][0] =  0.01571;
  dilep_par0_sigma_theta[1][1] =  0.01461;

  dilep_par0_sigma_theta[2][0] =  0.00887;
  dilep_par0_sigma_theta[2][1] =  0.00876;

  dilep_par0_sigma_theta[3][0] =  0.006560;
  dilep_par0_sigma_theta[3][1] =  0.005887;

  // Parameter 1 for eveluating sigma theta

  dilep_par1_sigma_theta[0][0] =  0.02705;
  dilep_par1_sigma_theta[0][1] =  0.05433;

  dilep_par1_sigma_theta[1][0] =  0.04767;
  dilep_par1_sigma_theta[1][1] =  0.03149;

  dilep_par1_sigma_theta[2][0] =  0.02942;
  dilep_par1_sigma_theta[2][1] =  0.03333;

  dilep_par1_sigma_theta[3][0] =  0.02545;
  dilep_par1_sigma_theta[3][1] =  0.02172;

  // Parameter 2 for eveluating sigma theta

  dilep_par2_sigma_theta[0][0] =  -8.3222;
  dilep_par2_sigma_theta[0][1] =  -4.4107;

  dilep_par2_sigma_theta[1][0] =  -4.6296;
  dilep_par2_sigma_theta[1][1] =  -2.9823;

  dilep_par2_sigma_theta[2][0] =  -2.8237;
  dilep_par2_sigma_theta[2][1] =  -3.0109;

  dilep_par2_sigma_theta[3][0] =  -2.7438;
  dilep_par2_sigma_theta[3][1] =  -2.2456;

  return dilep_par0_sigma_theta[layer][arm] + dilep_par1_sigma_theta[layer][arm]*exp(dilep_par2_sigma_theta[layer][arm] * pt);

}

void vertex_reconstruction(const char* inFile, const char* outFile, const char* inFileText, const int apply_deadmap, const int mode){


  TFile* input = new TFile(inFile,"READ");
  if(!(input))
  {
    cout << "no input file" << endl;
    exit(1);
  }

  TTree* T = (TTree*)input->Get("T");
  TBranch* br = T->GetBranch("MyEvent");
  MyEvent *event = 0;
  br->SetAddress(&event);

  cout << "Trees read!" << endl;  
  int nevt = T->GetEntries();

  hist_diff_bbcz = new TH1D("hist_diff_bbcz","BBCZ - Calculated Vertex", 600, -30, 30);
  hist_diff_vtxz = new TH1D("hist_diff_vtxz","VTXZ - Calculated Vertex", 600, -30, 30);
  hist = new TH1D("hist", "Vertex Distribution", 120, -30, 30);
  hist_stag = new TH1D("hist_stag", "Vertex Distribution", 121, -30.25, 30.25);
  debug = new TH1D("Debug", "Debug Histogram", 5, -0.5, 4.5);

  if(apply_deadmap){

    int a = -999; int l = -999; 
    int bin_phi = -999; int bin_z = -999;

    std::ifstream deadmap(inFileText);
    while(deadmap >> a >> l >> bin_phi >> bin_z){
      VTXDeadMap[a][l][bin_phi][bin_z] = 1;
    }

    for(int i=0; i<narms; i++){
      for(int j=0; j<nlayers; j++){

        if(j==0){
      
          if(i==0) phi_z[i][j] = new TH2D(Form("phi_z_%s_%s", arm_label[i], layer_label[j]), Form("phi_z_%s_%s", arm_label[i], layer_label[j]), 275, pi/2, 3*pi/2, 375, -21.25, 21.25);
          else phi_z[i][j] = new TH2D(Form("phi_z_%s_%s", arm_label[i], layer_label[j]), Form("phi_z_%s_%s", arm_label[i], layer_label[j]), 275, -pi/2, pi/2, 375, -21.25, 21.25);
          phi_z[i][j]->Sumw2();

        }

        if(j==1){
      
          if(i==0) phi_z[i][j] = new TH2D(Form("phi_z_%s_%s", arm_label[i], layer_label[j]), Form("phi_z_%s_%s", arm_label[i], layer_label[j]), 550, pi/2, 3*pi/2, 375, -21.25, 21.25);
          else phi_z[i][j] = new TH2D(Form("phi_z_%s_%s", arm_label[i], layer_label[j]), Form("phi_z_%s_%s", arm_label[i], layer_label[j]), 550, -pi/2, pi/2, 375, -21.25, 21.25);
          phi_z[i][j]->Sumw2();

        }

      }
    }
  }

  TFile* fmag = new TFile("field_map.root","READ");
  TH2F* field_map_bz = (TH2F*)fmag->Get("hist_bz");
  TH2F* field_map_br = (TH2F*)fmag->Get("hist_br");

  TFile* fout = new TFile(outFile, "RECREATE");
  TTree* newT = new TTree("T", "New Tree with Updated Vertex");
  MyEvent evt;
  newT->Branch("MyEvent", &evt);

  for (int ievent = 0; ievent < nevt; ievent++){

    int ghost[40], real[40];
    int all[40];
    int nghost = 0, nreal = 0, nall = 0;

    hist->Reset("ICESM");
    hist_stag->Reset("ICESM");

    if (ievent%10000==0) cout << "Event: " << ievent << " / " << nevt << endl;

    evt.ClearEvent();
    event->ClearEvent();
    br->GetEntry(ievent);

    int ntrack = event->GetNtrack();
    int nsvxhits = event->GetNVTXhit();
    int nclust = event->GetNcluster();
    int nmctrack = event->GetNMcTrack();
    double vtxx = event->GetPreciseX();
    double vtxy = event->GetPreciseY();
    double vtxz = event->GetPreciseZ();
    double calc_zvertex = -999;
    if(fabs(event->GetVtxZ()) > 30) continue;

    if(mode > 0){

      vtxx = gRandom->Gaus(event->GetPreciseX(), 0.0200);
      vtxy = gRandom->Gaus(event->GetPreciseY(), 0.0200);
      vtxz = gRandom->Gaus(event->GetPreciseZ(), 0.0751);

    }

    debug->Fill(0);

    if(mode == 0){

      for (int iprimary = 0; iprimary < ntrack; iprimary++){

        MyTrack mytrk = event->GetEntry(iprimary);
        if (!(mytrk.GetTrkQuality() == 31 || mytrk.GetTrkQuality() == 51 || mytrk.GetTrkQuality() == 63) ) continue;
	if ( mytrk.GetPtPrime() < 0.3 ) continue;

        for(int jj=0; jj<total_nlayers; jj++){ dilep_phi_projection[jj] = -999; }

        double phi0_trk_proj = mytrk.GetPhi0Prime();
        double the0_trk_proj = mytrk.GetThe0();
        double pz = mytrk.GetPtPrime()*(TMath::Cos(mytrk.GetThe0()))/(TMath::Sin(mytrk.GetThe0()));
        charge_index = -999;

        double rp = sqrt(  pow(vtxx, 2) + pow(vtxy, 2) );
        double zp = event->GetVtxZ();

        for(int p=1; p<1800; p++){

		      for(int l = 0; l < total_nlayers; l++){ 

            if( fabs(rp - radii[l]) < step_size && dilep_phi_projection[l] == -999 ){ dilep_phi_projection[l] = phi0_trk_proj; } 

          }

          int rbin = field_map_bz->GetXaxis()->FindBin(rp);
          int zbin = field_map_bz->GetYaxis()->FindBin(zp);

          double bz = field_map_bz->GetBinContent(rbin, zbin)/10000; // Conversion from Gaus to Tesla

          double delta_phi0 = (mytrk.GetChargePrime()*0.3*step_size*bz)/(2*mytrk.GetPtPrime()*100);
          phi0_trk_proj += delta_phi0;

          double bradial = field_map_br->GetBinContent(rbin, zbin)/10000; // Conversion from Gaus to Tesla

          //Bend in the z direction does not depend upon the charge.

          double delta_the0 = 0.3*bradial*(step_size*TMath::Tan(pi/2 - the0_trk_proj))/(2*pz*100);
          if( mytrk.GetThe0() > pi/2) the0_trk_proj -= delta_the0;
          else the0_trk_proj += delta_the0;

          zp += step_size*TMath::Tan(pi/2 - the0_trk_proj);
          rp += step_size;

        }

        if( mytrk.GetChargePrime() > 0 ) charge_index = mytrk.GetChargePrime();
        else charge_index = mytrk.GetChargePrime() + 1;

        int DCArm = mytrk.GetArm();
        TVector3 hitpoint;

        for(int ilayer = 0; ilayer < total_nlayers; ilayer++){

          for(int ihit = 0; ihit < nsvxhits ; ihit++){

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
                        
            hitpoint.SetXYZ(xhit-event->GetPreciseX(), yhit-event->GetPreciseY(), zhit-event->GetVtxZ());
            double phi_hit = hitpoint.Phi();
            if(phi_hit < -pi/2) phi_hit += 2*pi;

            double sigma_phi_value = -999;   double mean_phi_value = -999;

            sigma_phi_value = get_sigma_phi_data(layer, DCArm, mytrk.GetPtPrime());
            mean_phi_value = get_mean_phi_data(layer, charge_index, DCArm);

     if( (dilep_phi_projection[ilayer] - phi_hit) > (mean_phi_value - 7*sigma_phi_value) &&  (dilep_phi_projection[ilayer] - phi_hit) < (mean_phi_value + 7*sigma_phi_value)) found = true; 

            if(found && layer == 0){ 

              double xhit = svxhit.GetXHit();
              double yhit = svxhit.GetYHit();
              double zhit = svxhit.GetZHit();
              double r = sqrt(xhit*xhit + yhit*yhit);

              double zvertex = zhit - r/TMath::Tan(mytrk.GetThe0());
              hist->Fill(zvertex);
              hist_stag->Fill(zvertex);

            }
          }

        }

      } // End first track Loop

      if(hist->GetEntries() == 0){ debug->Fill(1); continue; }

      double max = hist->GetMaximum();
      calc_zvertex = hist->GetBinCenter(hist->GetMaximumBin());
      int nbins = hist->GetNbinsX();

      int counter = 0;

      for(int ibin = 1; ibin < nbins+1; ibin++){

        double content = hist->GetBinContent(ibin);
        if(content == 0) continue;
        if(content == max) counter++;

      }

      if(counter > 1){

        debug->Fill(2);
        counter = 0;
        nbins = hist_stag->GetNbinsX();
        max = hist_stag->GetMaximum();

        for(int ibin = 1; ibin < nbins+1; ibin++){

          double content = hist_stag->GetBinContent(ibin);
          if(content == 0) continue;
          if(content == max) counter++;
        }

        if(counter > 1){ debug->Fill(3); continue; }
        else calc_zvertex = hist_stag->GetBinCenter(hist_stag->GetMaximumBin());

      }

      debug->Fill(4);
      hist_diff_bbcz->Fill(event->GetVtxZ() - calc_zvertex);
      hist_diff_vtxz->Fill(event->GetPreciseZ() - calc_zvertex);

    }

    evt.SetPreciseX(vtxx);
    evt.SetPreciseY(vtxy);
    if(mode == 0) evt.SetPreciseZ(calc_zvertex);
    else evt.SetPreciseZ(vtxz);
    evt.SetVtxZ(event->GetVtxZ());
    evt.SetRunNumber(event->GetRunNumber());
    evt.SetEvtNo(event->GetEvtNo());
    

    // Get The VTX Hits Information Here

    for(int ihit = 0; ihit < nsvxhits ; ihit++){

      MyVTXHit svxhit = event->GetVTXHitEntry(ihit);
      MyVTXHit newHit;

      int layer = svxhit.GetLayer();

      if( (layer == 0 || layer ==1) && apply_deadmap ){

        int side = -999;
        double zhit = svxhit.GetZHit();
        double phi_hit = TMath::ATan2(svxhit.GetYHit(), svxhit.GetXHit());

        if(phi_hit < -pi/2) phi_hit += 2*pi;

        if(fabs(phi_hit) < pi/2) side = 1; //west
        else side = 0; //east

        int phi_bin = phi_z[side][layer]->GetXaxis()->FindBin(phi_hit);
        int z_bin = phi_z[side][layer]->GetYaxis()->FindBin(zhit);

        if(VTXDeadMap[side][layer][phi_bin][z_bin] == 1) continue;

      }

      newHit.SetClustId( svxhit.GetClustId() );
      newHit.SetLayer( svxhit.GetLayer() );
      newHit.SetLadder( svxhit.GetLadder() );
      newHit.SetSensor(svxhit.GetSensor() );
      newHit.SetXHit( svxhit.GetXHit() );
      newHit.SetYHit( svxhit.GetYHit() );
      newHit.SetZHit( svxhit.GetZHit() );
      evt.AddVTXHit( newHit );

    }

    // Now add the clusters to the new tree

    for(int iclust = 0; iclust < nclust ; iclust++){

      MyCluster myclust = event->GetClusterEntry(iclust);
      int cluster_id = myclust.GetID();
      bool asso_trk = false;
      MyCluster newcluster;

      for(int itrk = 0; itrk < ntrack; itrk++){

        MyTrack mytrk = event->GetEntry(itrk);
        int emcid = mytrk.GetEmcId();
        if(emcid == cluster_id) asso_trk = true;
      }

      if(asso_trk) continue;

      newcluster.SetID( myclust.GetID() );
      newcluster.SetX( myclust.GetX() );
      newcluster.SetY( myclust.GetY() );
      newcluster.SetZ( myclust.GetZ() );
      newcluster.SetEcore( myclust.GetEcore() );
      newcluster.SetProb( myclust.GetProb() );
      newcluster.SetChi2( myclust.GetChi2() );
      newcluster.SetTof( myclust.GetTof() );
      newcluster.SetTower( myclust.GetSect(), myclust.GetIY(), myclust.GetIZ() );
      evt.AddCluster( newcluster );

    }

    if(mode > 0){

      for(int imctrk = 0; imctrk < nmctrack; imctrk++){

	McTrack mctrk = event->GetMcEntry(imctrk);
        McTrack newMcTrack;

	newMcTrack.SetMcIndex( mctrk.GetMcIndex() );
	newMcTrack.SetParticleID( mctrk.GetParticleID() );
	newMcTrack.SetParentID( mctrk.GetParentID()  );
	newMcTrack.SetPrimaryID( mctrk.GetPrimaryID() );
	newMcTrack.SetGen( mctrk.GetGen() );

	newMcTrack.SetPx( mctrk.GetPx() );
	newMcTrack.SetPy( mctrk.GetPy() );
	newMcTrack.SetPz( mctrk.GetPz() );

	newMcTrack.SetParentPx( mctrk.GetParentPx() );
	newMcTrack.SetParentPy( mctrk.GetParentPy() );
	newMcTrack.SetParentPz( mctrk.GetParentPz() );

	newMcTrack.SetPrimaryPx( mctrk.GetPrimaryPx() );
	newMcTrack.SetPrimaryPy( mctrk.GetPrimaryPy() );
	newMcTrack.SetPrimaryPz( mctrk.GetPrimaryPz() );

	newMcTrack.SetVertexx( mctrk.GetVertexx() );
	newMcTrack.SetVertexy( mctrk.GetVertexy() );
	newMcTrack.SetVertexz( mctrk.GetVertexz() );

	newMcTrack.SetParentVertexx( mctrk.GetParentVertexx() );
	newMcTrack.SetParentVertexy( mctrk.GetParentVertexy() );
	newMcTrack.SetParentVertexz( mctrk.GetParentVertexz() );

	newMcTrack.SetPrimaryVertexx( mctrk.GetPrimaryVertexx() );
	newMcTrack.SetPrimaryVertexy( mctrk.GetPrimaryVertexy() );
	newMcTrack.SetPrimaryVertexz( mctrk.GetPrimaryVertexz() );

	newMcTrack.SetPhiDC( mctrk.GetPhiDC() );
	newMcTrack.SetPhi0( mctrk.GetPhi0() );
	newMcTrack.SetThe0( mctrk.GetThe0() );
	newMcTrack.SetZDC( mctrk.GetZDC() );
	newMcTrack.SetAlpha( mctrk.GetAlpha() );

	evt.AddMcTrack( newMcTrack );

      }
    }


    if(mode == 0){

      for(int kk = 0; kk < narms; kk++){ dilep_par0_theta[kk] = -999; dilep_par1_theta[kk] = -999; dilep_par2_theta[kk] = -999; }

      for(int l=0; l<rungroups; l++){

        if(evt.GetRunNumber() >= dilep_run_lo[l] && evt.GetRunNumber() <= dilep_run_hi[l]){

          dilep_par0_theta[0] = dilep_par0_theta_east[l]; dilep_par1_theta[0] = dilep_par1_theta_east[l]; dilep_par2_theta[0] = dilep_par2_theta_east[l];
          dilep_par0_theta[1] = dilep_par0_theta_west[l]; dilep_par1_theta[1] = dilep_par1_theta_west[l]; dilep_par2_theta[1] = dilep_par2_theta_west[l];

        }

      }


    }

    
    for (int k1 = 0; k1 < ntrack; ++k1)
      {
	  all[nall] = k1;
	  ++nall;
	  int ghost_flag = 0;

	  for (int k2 = 0; k2 < ntrack; ++k2)
	    {
	      if ( k1==k2 ) continue;
	      if ( fabs((event->GetEntry(k1)).GetCrkphi() - (event->GetEntry(k2)).GetCrkphi() ) < DCENTERCUT ) { ghost_flag = 1; break; }
	    }
	  if ( ghost_flag ) { ghost[nghost] = k1; ++nghost; }
	  else { real[nreal] = k1; ++nreal; }
      }

    for(int itrk = 0; itrk < nreal; itrk++){

      MyTrack mytrk = event->GetEntry(real[itrk]);
      if(!single_cut(mytrk)) continue;

      MyTrack newtrk;
    
      nhits_trk = 0;
      charge_index = -999;
      theta_offset = -999;

      for(int jj=0; jj<total_nlayers; jj++){

        dilep_hit_index[jj] = -999;
        dilep_hit_counter[jj] = -999;
        dilep_phi_projection[jj] = -999;
        dilep_the_projection[jj] = -999;

      }

      if(mode == 0){

        double new_the0 = mytrk.GetThe0() - ((evt.GetVtxZ() - evt.GetPreciseZ())/220)*TMath::Sin( mytrk.GetThe0()); // First order correction in theta0
        if( newtrk.GetArm() == 0) theta_offset = dilep_par0_theta[0]*TMath::Sin(new_the0) + dilep_par1_theta[0]*TMath::Cos(new_the0) + dilep_par2_theta[0];
        else theta_offset = dilep_par0_theta[1]*TMath::Sin(new_the0) + dilep_par1_theta[1]*TMath::Cos(new_the0) + dilep_par2_theta[1];
        newtrk.SetThe0Prime( new_the0 - theta_offset);

      }

      else newtrk.SetThe0Prime( mytrk.GetThe0());

      double ecore = mytrk.GetEcore();
      double px = mytrk.GetPtPrime()*TMath::Cos(mytrk.GetPhi0Prime());
      double py = mytrk.GetPtPrime()*TMath::Sin(mytrk.GetPhi0Prime());
      double pz = mytrk.GetPtPrime()*(TMath::Cos(newtrk.GetThe0Prime()))/(TMath::Sin(newtrk.GetThe0Prime()));
      double p = sqrt(px*px + py*py + pz*pz);

      if(ecore/p < 0.5) continue;

      newtrk.SetPhiDC( mytrk.GetPhiDC() );
      newtrk.SetPhi0( mytrk.GetPhi0() );
      newtrk.SetPhi0Prime( mytrk.GetPhi0Prime() );
      newtrk.SetPt( mytrk.GetPt() );
      newtrk.SetPtPrime( mytrk.GetPtPrime() );
      newtrk.SetThe0( mytrk.GetThe0() );
      newtrk.SetZDC( mytrk.GetZDC() );
      newtrk.SetAlpha( mytrk.GetAlpha() );
      newtrk.SetAlphaPrime( mytrk.GetAlphaPrime() );
      newtrk.SetQ( mytrk.GetCharge() );
      newtrk.SetQPrime( mytrk.GetChargePrime() );
      newtrk.SetEmcId( mytrk.GetEmcId() );
      newtrk.SetEmcTOF( mytrk.GetEmcTOF() );
      newtrk.SetEcore( mytrk.GetEcore() );
      newtrk.SetN0( mytrk.GetN0() );
      newtrk.SetDep( mytrk.GetDep() );
      newtrk.SetProb( mytrk.GetProb() );
      newtrk.SetEmcdz( mytrk.GetEmcdz() );
      newtrk.SetEmcdphi( mytrk.GetEmcdphi() );
      newtrk.SetCrkphi( mytrk.GetCrkphi() );
      newtrk.SetCrkz( mytrk.GetCrkz() );
      newtrk.SetArm( mytrk.GetArm() );
      newtrk.SetDisp( mytrk.GetDisp() );
      newtrk.SetChi2Npe0( mytrk.GetChi2Npe0() );
      newtrk.SetisERT( mytrk.GetisERT() );
      newtrk.SetMcId( mytrk.GetMcId() );
      newtrk.SetSect( mytrk.GetSect() );

      double phi0_trk_proj = newtrk.GetPhi0Prime();
      double the0_trk_proj = newtrk.GetThe0Prime();

      double rp = sqrt(  pow(evt.GetPreciseX(), 2) + pow(evt.GetPreciseY(), 2) );
      double zp = evt.GetPreciseZ();

      for(int p=1; p<1800; p++){

		    for(int l = 0; l < total_nlayers; l++){ 

          if( fabs(rp - radii[l]) < step_size && dilep_phi_projection[l] == -999 ){ dilep_phi_projection[l] = phi0_trk_proj; dilep_the_projection[l] = the0_trk_proj; } 

        }

        int rbin = field_map_bz->GetXaxis()->FindBin(rp);
        int zbin = field_map_bz->GetYaxis()->FindBin(zp);

        double bz = field_map_bz->GetBinContent(rbin, zbin)/10000; // Conversion from Gaus to Tesla

        double delta_phi0 = (newtrk.GetChargePrime()*0.3*step_size*bz)/(2*newtrk.GetPtPrime()*100);
        phi0_trk_proj += delta_phi0;

        double bradial = field_map_br->GetBinContent(rbin, zbin)/10000; // Conversion from Gaus to Tesla

        //Bend in the z direction does not depend upon the charge.

        double delta_the0 = 0.3*bradial*(step_size*TMath::Tan(pi/2 - the0_trk_proj))/(2*pz*100);
        if( newtrk.GetThe0Prime() > pi/2) the0_trk_proj -= delta_the0;
        else the0_trk_proj += delta_the0;

        zp += step_size*TMath::Tan(pi/2 - the0_trk_proj);
        rp += step_size;


      }

      

      if( newtrk.GetChargePrime() > 0 ) charge_index = newtrk.GetChargePrime();
      else charge_index = newtrk.GetChargePrime() + 1;

      int DCArm = newtrk.GetArm();
      TVector3 hitpoint;
      double chi2 = -999;

      for(int ilayer = 0; ilayer < total_nlayers; ilayer++){

		    double min = -999;  int hit_flag_counter = 0;

        for(int ihit = 0; ihit < nsvxhits ; ihit++){

          hitpoint.SetXYZ(0, 0, 0);
          MyVTXHit svxhit = evt.GetVTXHitEntry(ihit);
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
                        
          hitpoint.SetXYZ(xhit-evt.GetPreciseX(), yhit-evt.GetPreciseY(), zhit-evt.GetPreciseZ());
          double phi_hit = hitpoint.Phi();
          if(phi_hit < -pi/2) phi_hit += 2*pi;
          double theta_hit = hitpoint.Theta();

          // Now defining our initial search windows

          double sigma_phi_value = -999;   double mean_phi_value = -999;
          double sigma_theta_value = -999;   double mean_theta_value = -999;

          sigma_phi_value = get_sigma_phi_data(layer, DCArm, newtrk.GetPtPrime());
          mean_phi_value = get_mean_phi_data(layer, charge_index, DCArm);

          sigma_theta_value = get_sigma_theta_data(layer, DCArm, newtrk.GetPtPrime());
          mean_theta_value = get_mean_theta_data(layer, charge_index, DCArm);

	        if( (dilep_phi_projection[ilayer] - phi_hit) > (mean_phi_value - 7*sigma_phi_value) &&  (dilep_phi_projection[ilayer] - phi_hit) < (mean_phi_value + 7*sigma_phi_value) && (dilep_the_projection[ilayer] - theta_hit) > (mean_theta_value - 5*sigma_theta_value) &&  (dilep_the_projection[ilayer] - theta_hit) < (mean_theta_value + 5*sigma_theta_value)){ found = true; }

          if(found){ 
	    
            hit_flag_counter++;
	    double diff = sqrt( pow( (dilep_phi_projection[ilayer] - phi_hit)/sigma_phi_value, 2) + pow( (dilep_the_projection[ilayer] - theta_hit)/sigma_theta_value, 2) );
            if(hit_flag_counter == 1) { min = diff; dilep_hit_index[ilayer] = ihit; dilep_hit_counter[ilayer] = hit_flag_counter; }
            else{  if(diff < min){  min = diff; dilep_hit_index[ilayer] = ihit; dilep_hit_counter[ilayer] = hit_flag_counter; }   }

          }
        }

        if(dilep_hit_index[ilayer] != -999){  nhits_trk++; if(chi2 == -999){ chi2 = min; } else chi2 += min; }

      }

      if(nhits_trk > 0) newtrk.SetChi2(chi2/nhits_trk); 

      newtrk.SetNHits(nhits_trk);
      newtrk.SetHitIndexL1(dilep_hit_index[0]);
      newtrk.SetHitIndexL2(dilep_hit_index[1]);
      newtrk.SetHitIndexL3(dilep_hit_index[2]);
      newtrk.SetHitIndexL4(dilep_hit_index[3]);
      newtrk.SetHitIndexL5(dilep_hit_index[4]);
      newtrk.SetHitIndexL6(dilep_hit_index[5]);
      newtrk.SetHitIndexL7(dilep_hit_index[6]);
      newtrk.SetHitIndexL8(dilep_hit_index[7]);
      newtrk.SetHitCounterL1(dilep_hit_counter[0]);
      newtrk.SetHitCounterL2(dilep_hit_counter[1]);
      newtrk.SetHitCounterL3(dilep_hit_counter[2]);
      newtrk.SetHitCounterL4(dilep_hit_counter[3]);
      newtrk.SetHitCounterL5(dilep_hit_counter[4]);
      newtrk.SetHitCounterL6(dilep_hit_counter[5]);
      newtrk.SetHitCounterL7(dilep_hit_counter[6]);
      newtrk.SetHitCounterL8(dilep_hit_counter[7]);
      evt.AddTrack(newtrk);
    
    }

    newT->Fill();

  }// End Event A loop

  fout->cd();
  hist_diff_bbcz->Write();
  hist_diff_vtxz->Write();
  newT->Write();
  debug->Write();
  fout->Close();

} 
