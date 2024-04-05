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
#include "/direct/phenix+u/vdoomra/install/include/DileptonAnalysis/MyEvent.h"

using namespace std;
using namespace DileptonAnalysis;
const int centbin = 1;
const double pi = TMath::ACos(-1);

const double DCENTERCUT = 0.05; // RICH ghost cut
const double PC1_DPHI_CUT = 0.02, PC1_DZ_CUT = 0.5; // PC1 ghost cut 

const int nlayers = 4;
const int ncharge = 2;
const int nsides = 2;
const int total_nlayers = 8;

static const double Me = 0.000510998918;
static const double Me2 = Me*Me;

TH2D* h2d_ee_FG_temp[centbin] = {NULL};
TH2D* h2d_ee_FG[centbin] = {NULL};
TH1D* h1d_ee_FG[centbin] = {NULL};

TH2D* h2d_ee_FG_like_temp[centbin] = {NULL};
TH2D* h2d_ee_FG_like[centbin] = {NULL};
TH1D* h1d_ee_FG_like[centbin] = {NULL};

TH2D* h2d_pt_spectra_mass_temp[centbin] = {NULL};
TH2D* h2d_pt_spectra_mass[centbin] = {NULL};
TH2D* h2d_delta_phi_mass_temp[centbin] = {NULL};
TH2D* h2d_delta_phi_mass[centbin] = {NULL};

const int ntype = 6;
const char* description[ntype] = {"pair_dca_from_single_tracks_abs_diff", "pair_dca_from_single_tracks_diff", "pair_dca_from_single_tracks_abs_average", "pair_dca_from_single_tracks_average", "pair_dca_using_quadrature", "pair_dca_using_quadrature_add"};

TH3D* hist_temp[ntype][ncharge] = {NULL};
TH3D* hist[ntype][ncharge] = {NULL};

const int nbins_mass = 38;
const double mass_range[nbins_mass+1] = { 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75, 0.80, 0.87, 0.95, 1.05, 1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.00, 3.05, 3.10, 3.15, 3.25, 3.35, 3.50, 3.65, 3.75, 3.85, 4.00, 4.5};


const int nbins_dca = 100;
const double dca_range[nbins_dca+1] = {-2000,-1960,-1920,-1880,-1840,-1800,-1760,-1720,-1680,-1640,-1600,-1560,-1520,-1480,-1440,-1400,-1360,-1320,-1280,-1240,-1200,-1160,-1120,-1080,-1040,-1000,-960,-920,-880,-840,-800,-760,-720,-680,-640,-600,-560,-520,-480,-440,-400,-360,-320,-280,-240,-200,-160,-120,-80,-40,0,40,80,120,160,200,240,280,320,360,400,440,480,520,560,600,640,680,720,760,800,840,880,920,960,1000,1040,1080,1120,1160,1200,1240,1280,1320,1360,1400,1440,1480,1520,1560,1600,1640,1680,1720,1760,1800,1840,1880,1920,1960,2000};

const int nbins_pt = 100;
const double pt_range[nbins_pt+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.0};

const int nbins_mass_prof = 10;
const double mass_range_prof[nbins_mass_prof+1] = { 1.35, 1.65, 1.95, 2.25, 2.55, 2.85, 3.05, 3.25, 3.65, 3.85, 4.2 };

TProfile* hist_prof[ntype][ncharge] = {NULL};

const char* charge_type[ncharge] = {"like", "unlike"};

// Windows for finding the photon cluster //

const double mean_par0 = -0.000623;
const double mean_par1 = 0.03898;

const double sigma_par0 = 0.00159;
const double sigma_par1 = 0.00829;

// ===================================== //

double veto_window_mean(int layer, int charge, int delphi_side, double pt){

  double veto_window_mean_par0[nlayers-1][ncharge][nsides] = {-999};
  double veto_window_mean_par1[nlayers-1][ncharge][nsides] = {-999};
  double veto_window_mean_par2[nlayers-1][ncharge][nsides] = {-999};

  // Parameter 0

  veto_window_mean_par0[0][0][0] = -0.001876;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE
  veto_window_mean_par0[0][0][1] =  0.003013;
  veto_window_mean_par0[0][1][0] = -0.003025;
  veto_window_mean_par0[0][1][1] =  0.001911;

  veto_window_mean_par0[1][0][0] = -0.004279;
  veto_window_mean_par0[1][0][1] =  0.006654;
  veto_window_mean_par0[1][1][0] = -0.006633;
  veto_window_mean_par0[1][1][1] =  0.004403;

  veto_window_mean_par0[2][0][0] = -0.006293;
  veto_window_mean_par0[2][0][1] =  0.009369;
  veto_window_mean_par0[2][1][0] = -0.01025;
  veto_window_mean_par0[2][1][1] =  0.00720;
    

  // Parameter 1 

  veto_window_mean_par1[0][0][0] = -0.20232;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE
  veto_window_mean_par1[0][0][1] = 0.06833;
  veto_window_mean_par1[0][1][0] = -0.10529;
  veto_window_mean_par1[0][1][1] = 0.03370;

  veto_window_mean_par1[1][0][0] = -0.0252;
  veto_window_mean_par1[1][0][1] = 0.0338;
  veto_window_mean_par1[1][1][0] = -0.03311;
  veto_window_mean_par1[1][1][1] = 0.02312;

  veto_window_mean_par1[2][0][0] = -0.02067;
  veto_window_mean_par1[2][0][1] = 0.05404;
  veto_window_mean_par1[2][1][0] = -0.05873;
  veto_window_mean_par1[2][1][1] = 0.02255;

  // Parameter 2

  veto_window_mean_par2[0][0][0] = -9.1789;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE
  veto_window_mean_par2[0][0][1] = -8.58896;
  veto_window_mean_par2[0][1][0] = -9.4364;
  veto_window_mean_par2[0][1][1] = -6.8837;

  veto_window_mean_par2[1][0][0] = -2.1499;
  veto_window_mean_par2[1][0][1] = -1.9659;
  veto_window_mean_par2[1][1][0] = -1.8537;
  veto_window_mean_par2[1][1][1] = -1.7946;

  veto_window_mean_par2[2][0][0] = -1.5888;
  veto_window_mean_par2[2][0][1] = -1.7085;
  veto_window_mean_par2[2][1][0] = -1.8974;
  veto_window_mean_par2[2][1][1] = -1.9511;

  return veto_window_mean_par0[layer][charge][delphi_side] + veto_window_mean_par1[layer][charge][delphi_side]*exp(veto_window_mean_par2[layer][charge][delphi_side] * pt);
}

double veto_window_sigma(int layer, int charge, int delphi_side, double pt){

  double veto_window_sigma_par0[nlayers-1][ncharge][nsides] = {-999};
  double veto_window_sigma_par1[nlayers-1][ncharge][nsides] = {-999};
  double veto_window_sigma_par2[nlayers-1][ncharge][nsides] = {-999};

  // Parameter 0

  veto_window_sigma_par0[0][0][0] = 0.0002732;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE
  veto_window_sigma_par0[0][0][1] = 0.0008276;
  veto_window_sigma_par0[0][1][0] = 0.0008218;
  veto_window_sigma_par0[0][1][1] = 0.0002933;

  veto_window_sigma_par0[1][0][0] = 0.00041;
  veto_window_sigma_par0[1][0][1] = 0.001677;
  veto_window_sigma_par0[1][1][0] = 0.001626;
  veto_window_sigma_par0[1][1][1] = 0.000314;

  veto_window_sigma_par0[2][0][0] = 0.001299;
  veto_window_sigma_par0[2][0][1] = 0.001982;
  veto_window_sigma_par0[2][1][0] = 0.002024;
  veto_window_sigma_par0[2][1][1] = 0.000925;
    

  // Parameter 1 

  veto_window_sigma_par1[0][0][0] = 0.008219;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE
  veto_window_sigma_par1[0][0][1] = 0.11697;
  veto_window_sigma_par1[0][1][0] = 0.05595;
  veto_window_sigma_par1[0][1][1] = 0.01959;

  veto_window_sigma_par1[1][0][0] = 0.00903;
  veto_window_sigma_par1[1][0][1] = 0.005199;
  veto_window_sigma_par1[1][1][0] = 0.004186;
  veto_window_sigma_par1[1][1][1] = 0.01498;

  veto_window_sigma_par1[2][0][0] = 0.01592;
  veto_window_sigma_par1[2][0][1] = 0.006797;
  veto_window_sigma_par1[2][1][0] = 0.007108;
  veto_window_sigma_par1[2][1][1] = 0.009907;

  // Parameter 2

  veto_window_sigma_par2[0][0][0] = -6.0461;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE
  veto_window_sigma_par2[0][0][1] = -11.8997;
  veto_window_sigma_par2[0][1][0] = -10.3209;
  veto_window_sigma_par2[0][1][1] = -7.5985;

  veto_window_sigma_par2[1][0][0] = -2.3919;
  veto_window_sigma_par2[1][0][1] = -2.1199;
  veto_window_sigma_par2[1][1][0] = -1.6409;
  veto_window_sigma_par2[1][1][1] = -2.8119;

  veto_window_sigma_par2[2][0][0] = -3.5943;
  veto_window_sigma_par2[2][0][1] = -1.9312;
  veto_window_sigma_par2[2][1][0] = -1.9647;
  veto_window_sigma_par2[2][1][1] = -2.5570;

  return veto_window_sigma_par0[layer][charge][delphi_side] + veto_window_sigma_par1[layer][charge][delphi_side]*exp(veto_window_sigma_par2[layer][charge][delphi_side] * pt);
}


bool single_cut(MyTrack* mytrk)
{

  float px = mytrk->GetPtPrime()*TMath::Cos(mytrk->GetPhi0Prime());
  float py = mytrk->GetPtPrime()*TMath::Sin(mytrk->GetPhi0Prime());
  float pz = mytrk->GetPtPrime()*(TMath::Cos(mytrk->GetThe0Prime()))/(TMath::Sin(mytrk->GetThe0Prime()));
  float p = sqrt(px*px + py*py + pz*pz);
  float pt = sqrt( px*px + py*py );

  if(mytrk->GetN0() < 3) return false;
  if((mytrk->GetEcore()/p) < 0.8) return false;
  if(pt < 0.4) return false;

  int indexL1 = mytrk->GetHitIndexL1();
  int indexL2 = mytrk->GetHitIndexL2();
  int indexL3 = mytrk->GetHitIndexL3();
  int indexL4 = mytrk->GetHitIndexL4();
  int indexL5 = mytrk->GetHitIndexL5();
  int indexL6 = mytrk->GetHitIndexL6();
  int indexL7 = mytrk->GetHitIndexL7();
  int indexL8 = mytrk->GetHitIndexL8();

  if( !((indexL1 != -999) && (indexL2 != -999 || indexL3 != -999 || indexL4 != -999 || indexL5 != -999 || indexL6 != -999 || indexL7 != -999 || indexL8 != -999 )) ) return false;

  return true;
}

void DCA_data(const char* inFile, const char* outFile, const int isert, const int add_photon, const int additional_rejection)
{
  for(int i=0; i < ntype; i++){
    for(int j=0; j<ncharge; j++){
      hist_temp[i][j] = new TH3D(Form("hist_temp_%s_%s", description[i], charge_type[j]),"DCA Distribution", nbins_mass, mass_range, nbins_pt, pt_range, nbins_dca, dca_range);
      hist_temp[i][j]->Sumw2();
      hist_prof[i][j] = new TProfile(Form("hist_prof_%s_%s", description[i], charge_type[j]),"DCA Distribution Profile Histogram", nbins_mass_prof, mass_range_prof);
      hist_prof[i][j]->Sumw2();
    }
  }

  for(int icent = 0; icent < centbin; icent++){

    h2d_pt_spectra_mass_temp[icent] = new TH2D(Form("h2d_pt_spectra_mass_temp_cent%d",icent), "Pt spectra vs Pair Mass", nbins_mass, mass_range, 100, 0, 10);
    h2d_pt_spectra_mass_temp[icent]->Sumw2();

    h2d_delta_phi_mass_temp[icent] = new TH2D(Form("h2d_delta_phi_mass_temp_cent%d",icent), "Delta Phi vs Pair Mass", nbins_mass, mass_range, 100, 0, pi);
    h2d_delta_phi_mass_temp[icent]->Sumw2();
   
    h2d_ee_FG_temp[icent] = new TH2D(Form("h2d_ee_FG_temp_cent%d", icent), "ee mass distribution FG", nbins_mass, mass_range ,100, 0, 10);
    h2d_ee_FG_temp[icent]->Sumw2();

    h2d_ee_FG_like_temp[icent] = new TH2D(Form("h2d_ee_FG_like_temp_cent%d", icent), "ee mass distribution Like FG", nbins_mass, mass_range ,100, 0, 10);
    h2d_ee_FG_like_temp[icent]->Sumw2();
  }


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

  for (int ievent = 0; ievent < nevt; ievent++)
  {
    if (ievent%100000==0) cout << "Event: " << ievent << " / " << nevt << endl;

    event->ClearEvent();
    br->GetEntry(ievent);
    
    int ntrack = event->GetNtrack();
    int nvtxhits = event->GetNVTXhit();
    int nclust = event->GetNcluster();
    float vtxz = event->GetPreciseZ();
    float vtxx = event->GetPreciseX();
    float vtxy = event->GetPreciseY();
    if( fabs(vtxz) > 10 ) continue;

    int ghost[40], real[40];
    int all[40];
    int nghost = 0, nreal = 0, nall = 0;
   
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

    if(nreal<1) continue;

    int nprimary = 0;
    int primary_id[10] = {-999};

    if(additional_rejection){

    for(int ireal=0; ireal<nreal; ireal++){

      MyTrack mytrk = event->GetEntry(real[ireal]);
      float px = mytrk.GetPtPrime()*TMath::Cos(mytrk.GetPhi0Prime());
      float py = mytrk.GetPtPrime()*TMath::Sin(mytrk.GetPhi0Prime());
      float pz = mytrk.GetPtPrime()*(TMath::Cos(mytrk.GetThe0Prime()))/(TMath::Sin(mytrk.GetThe0Prime()));
      float pt = sqrt( px*px + py*py );
      
      
      int charge = mytrk.GetChargePrime();

      int charge_index = -999;
      if(charge<0) charge_index = 0;
      else charge_index = 1;

      int hit_index[total_nlayers] = {-999};

      hit_index[0] = mytrk.GetHitIndexL1();
      hit_index[1] = mytrk.GetHitIndexL2();
      hit_index[2] = mytrk.GetHitIndexL3();
      hit_index[3] = mytrk.GetHitIndexL4();
      hit_index[4] = mytrk.GetHitIndexL5();
      hit_index[5] = mytrk.GetHitIndexL6();
      hit_index[6] = mytrk.GetHitIndexL7();
      hit_index[7] = mytrk.GetHitIndexL8();

      bool veto_L1 = false; bool veto_L2 = false; bool veto_L3 = false;

      for(int ilayer=1; ilayer<total_nlayers; ilayer++){

        if(hit_index[ilayer] == -999) continue;

        int layer_index = -999;
        if(ilayer == 1) layer_index = 0;
        else if(ilayer == 2 || ilayer == 3 || ilayer == 4) layer_index = 1;
        else layer_index = 2;

        MyVTXHit vtxhit = event->GetVTXHitEntry(hit_index[ilayer]);
        double zhit = vtxhit.GetZHit();
        double phi_hit = TMath::ATan2(vtxhit.GetYHit(), vtxhit.GetXHit());
        if(phi_hit < -pi/2) phi_hit += 2*pi;

        for(int ihit = 0; ihit < nvtxhits; ihit++){

          if(ihit == hit_index[ilayer]) continue;
          int dphi_index = -999;
          MyVTXHit vtx_near_hit = event->GetVTXHitEntry(ihit);
          int layer = vtx_near_hit.GetLayer();
          if( (ilayer==1) && ilayer != layer) continue;
          if( (ilayer==2 || ilayer==3 || ilayer==4) && layer !=2) continue;
          if( (ilayer==5 || ilayer==6 || ilayer==7) && layer !=3) continue;
          double phi_near_hit = TMath::ATan2(vtx_near_hit.GetYHit(), vtx_near_hit.GetXHit());
          if(phi_near_hit < -pi/2) phi_near_hit += 2*pi;
          double z_near_hit = vtx_near_hit.GetZHit();

          double dphi = phi_near_hit - phi_hit;
          double dz = z_near_hit - zhit;

          if(layer==1 && fabs(dz) > 0.2) continue; 
          else if((layer==2 || layer==3) && fabs(dz) > 0.4 ) continue;

          if(dphi<0) dphi_index = 0;
          else dphi_index = 1;

          double mean = veto_window_mean(layer_index, charge_index, dphi_index, pt);
          double sigma = veto_window_sigma(layer_index, charge_index, dphi_index, pt);

          if( fabs(dphi) > 0.001 && fabs(dphi) < fabs(mean)+15*sigma && layer==1) veto_L1 = true;  
          if( fabs(dphi) > 0.001 && fabs(dphi) < fabs(mean)+15*sigma && layer==2) veto_L2 = true; 
          if( fabs(dphi) > 0.001 && fabs(dphi) < fabs(mean)+15*sigma && layer==3) veto_L3 = true;

        }
      }

      if( (veto_L1 || veto_L2 || veto_L3) ) continue;
      primary_id[nprimary] = real[ireal];
      nprimary++;

    }

    if(nprimary<1) continue;

    }

    else{

      nprimary = nreal;
      for(int ii=0; ii<nreal; ii++) { primary_id[ii] = real[ii];  }

    }


    for (int iprimary_A1 = 0; iprimary_A1 < nprimary; iprimary_A1++)
    {
      MyTrack mytrk_A1 = event->GetEntry(primary_id[iprimary_A1]);
      if(!(single_cut(&mytrk_A1))) continue;

      double pt_A1 = mytrk_A1.GetPtPrime();
      double px_A1 = mytrk_A1.GetPtPrime()*TMath::Cos(mytrk_A1.GetPhi0Prime());
      double py_A1 = mytrk_A1.GetPtPrime()*TMath::Sin(mytrk_A1.GetPhi0Prime());
      double pz_A1 = mytrk_A1.GetPtPrime()*(TMath::Cos(mytrk_A1.GetThe0Prime()))/(TMath::Sin(mytrk_A1.GetThe0Prime()));
      double phi0_A1 = mytrk_A1.GetPhi0Prime();
      double phi_A1 = mytrk_A1.GetPhiDC();
      double zed_A1 = mytrk_A1.GetZDC();

      int charge_A1 = mytrk_A1.GetChargePrime();
      int emcid_A1 = mytrk_A1.GetEmcId();
      int arm_A1 = mytrk_A1.GetArm();
      int isERT_A1 = mytrk_A1.GetisERT();
      double ecore_A1 = sqrt(pt_A1*pt_A1 + pz_A1*pz_A1 + Me2);
      bool isert_fired_A1 = false;

      double energy_threshold_A1 = -999;

      if(mytrk_A1.GetArm() == 1){  // West Arm

	      if (mytrk_A1.GetSect() == 0) { energy_threshold_A1 = 1.652; }
	      else if (mytrk_A1.GetSect() == 1) { energy_threshold_A1 = 1.511; }
	      else if (mytrk_A1.GetSect() == 2) { energy_threshold_A1 = 1.541; }
	      else { energy_threshold_A1 = 1.497; }

      }


      else {  // East Arm                                                                                                                               

	      if (mytrk_A1.GetSect() == 0) { energy_threshold_A1 = 2.239; }
	      else if (mytrk_A1.GetSect() == 1) { energy_threshold_A1 = 2.561; }
	      else if(mytrk_A1.GetSect() == 2) { energy_threshold_A1 = 1.676; }
	      else { energy_threshold_A1 = 1.465; }

      }

      isert_fired_A1 = isERT_A1 && ecore_A1 > energy_threshold_A1;

      TLorentzVector trk_A1;
      trk_A1.SetPxPyPzE(px_A1, py_A1, pz_A1, ecore_A1);

      TVector3 trk3_A1;
      trk3_A1.SetXYZ(px_A1, py_A1, pz_A1);

      if(add_photon){

	    double highest_energy_cluster_A1 = -999;
	    double ecore_highest_energy_cluster_A1 = -999;

      for (int iclust_A1 = 0; iclust_A1 < nclust; iclust_A1++)
      {
        MyCluster myclust_A1 = event->GetClusterEntry(iclust_A1);
        int id_clust_A1 = myclust_A1.GetID();
        if(emcid_A1 == id_clust_A1) continue;
        if(myclust_A1.GetChi2() > 3) continue;
        int arm_clust_A1 = myclust_A1.GetArm();
        double ecore_clust_A1 = myclust_A1.GetEcore();
        if(ecore_clust_A1 < 0.3) continue;
        if(arm_A1 != arm_clust_A1) continue;

        double xclust_A1 = myclust_A1.GetX() - event->GetPreciseX();
        double yclust_A1 = myclust_A1.GetY() - event->GetPreciseY();
        double zclust_A1 = myclust_A1.GetZ() - event->GetPreciseZ(); // subtract the zvertex

        int sect_clust_A1 = myclust_A1.GetSect();

        if (sect_clust_A1<0) continue;

        TLorentzVector emcPhoton_A1;
        emcPhoton_A1.SetX( ecore_clust_A1*xclust_A1/sqrt(xclust_A1*xclust_A1 + yclust_A1*yclust_A1 + zclust_A1*zclust_A1) );
        emcPhoton_A1.SetY( ecore_clust_A1*yclust_A1/sqrt(xclust_A1*xclust_A1 + yclust_A1*yclust_A1 + zclust_A1*zclust_A1) );
        emcPhoton_A1.SetZ( ecore_clust_A1*zclust_A1/sqrt(xclust_A1*xclust_A1 + yclust_A1*yclust_A1 + zclust_A1*zclust_A1) );
        emcPhoton_A1.SetE( ecore_clust_A1 );

        double phi0_clust_A1 = emcPhoton_A1.Phi();
        if(phi0_clust_A1 < -pi/2) phi0_clust_A1 += 2*pi;

	      double window_lo = -999;
	      double window_hi = -999;

	      if(charge_A1 == -1){ window_lo = 0.0; window_hi = (mean_par0 + mean_par1/pt_A1) + 2*(sigma_par0 + sigma_par1/pt_A1);  }
	      else { window_hi = 0.0; window_lo = -1*((mean_par0 + mean_par1/pt_A1) + 2*(sigma_par0 + sigma_par1/pt_A1));   }

        if( (phi0_A1 - phi0_clust_A1) > window_lo && (phi0_A1 - phi0_clust_A1) < window_hi){ 

          if(highest_energy_cluster_A1 == -999){ highest_energy_cluster_A1 = iclust_A1; ecore_highest_energy_cluster_A1 = ecore_clust_A1; }
          if(highest_energy_cluster_A1 != -999 && ecore_clust_A1 > ecore_highest_energy_cluster_A1){ highest_energy_cluster_A1 = iclust_A1; ecore_highest_energy_cluster_A1 = ecore_clust_A1; }

        }

      }

      if(highest_energy_cluster_A1 != -999){

        MyCluster myclust_found = event->GetClusterEntry(highest_energy_cluster_A1);
        double ecore_clust_found = myclust_found.GetEcore();
        double xclust_found = myclust_found.GetX() - event->GetPreciseX();
        double yclust_found = myclust_found.GetY() - event->GetPreciseY();
        double zclust_found = myclust_found.GetZ() - event->GetPreciseZ(); // subtract the zvertex

        TLorentzVector emcPhoton_found;
        emcPhoton_found.SetX( ecore_clust_found*xclust_found/sqrt(xclust_found*xclust_found + yclust_found*yclust_found + zclust_found*zclust_found) );
        emcPhoton_found.SetY( ecore_clust_found*yclust_found/sqrt(xclust_found*xclust_found + yclust_found*yclust_found + zclust_found*zclust_found) );
        emcPhoton_found.SetZ( ecore_clust_found*zclust_found/sqrt(xclust_found*xclust_found + yclust_found*yclust_found + zclust_found*zclust_found) );
        emcPhoton_found.SetE( ecore_clust_found );

        trk_A1 += emcPhoton_found;

      }

      }

      double R_A1 = trk_A1.Pt()/(0.003*0.90); // Radius in cm

      int hit_index_A1[total_nlayers] = {-999};

      hit_index_A1[0] = mytrk_A1.GetHitIndexL1();
      hit_index_A1[1] = mytrk_A1.GetHitIndexL2();
      hit_index_A1[2] = mytrk_A1.GetHitIndexL3();
      hit_index_A1[3] = mytrk_A1.GetHitIndexL4();
      hit_index_A1[4] = mytrk_A1.GetHitIndexL5();
      hit_index_A1[5] = mytrk_A1.GetHitIndexL6();
      hit_index_A1[6] = mytrk_A1.GetHitIndexL7();
      hit_index_A1[7] = mytrk_A1.GetHitIndexL8();

      double x1_A1 = -999; double y1_A1 = -999;
      double x2_A1 = -999; double y2_A1 = -999;

      for(int ilayer=0; ilayer<total_nlayers; ilayer++){

        if(hit_index_A1[ilayer] == -999) continue;

        MyVTXHit vtxhit = event->GetVTXHitEntry(hit_index_A1[ilayer]);
        if(x1_A1 == -999){ x1_A1 = vtxhit.GetXHit();  y1_A1 = vtxhit.GetYHit(); continue; }
        if(x2_A1 == -999){ x2_A1 = vtxhit.GetXHit();  y2_A1 = vtxhit.GetYHit(); }

      }

      double a_A1 = ((x1_A1-x2_A1) -(y1_A1*y1_A1 - y2_A1*y2_A1)/(x1_A1-x2_A1))*0.5;
      double b_A1 = (y1_A1-y2_A1)/(x1_A1-x2_A1);

      double yc1_A1 = (y1_A1-a_A1*b_A1 + sqrt(R_A1*R_A1*(1+b_A1*b_A1) - (a_A1 + b_A1*y1_A1)*(a_A1 + b_A1*y1_A1)))/(1+b_A1*b_A1);
      double yc2_A1 = (y1_A1-a_A1*b_A1 - sqrt(R_A1*R_A1*(1+b_A1*b_A1) - (a_A1 + b_A1*y1_A1)*(a_A1 + b_A1*y1_A1)))/(1+b_A1*b_A1);

      double xc1_A1 = x1_A1 - a_A1 - b_A1*yc1_A1;
      double xc2_A1 = x1_A1 - a_A1 - b_A1*yc2_A1;

      double xc_A1 = -999; double yc_A1 = -999;

      if(mytrk_A1.GetPhi0Prime() > -pi/2 && mytrk_A1.GetPhi0Prime() < pi/2 && mytrk_A1.GetChargePrime() == -1){ 
	
	      if(yc1_A1<0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
	      else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      if(mytrk_A1.GetPhi0Prime() > -pi/2 && mytrk_A1.GetPhi0Prime() < pi/2 && mytrk_A1.GetChargePrime() == 1){

        if(yc1_A1>0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
	      else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      if(mytrk_A1.GetPhi0Prime() > pi/2 && mytrk_A1.GetPhi0Prime() < 3*pi/2 && mytrk_A1.GetChargePrime() == -1){

        if(yc1_A1>0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
	      else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      if(mytrk_A1.GetPhi0Prime() > pi/2 && mytrk_A1.GetPhi0Prime() < 3*pi/2 && mytrk_A1.GetChargePrime() == 1){

        if(yc1_A1<0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
        else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      double L_A1 = sqrt((xc_A1-vtxx)*(xc_A1-vtxx) + (yc_A1-vtxy)*(yc_A1-vtxy));
      double dca_A1 = (L_A1-R_A1)*10000; // In Micro Meters

      for (int iprimary_A2 = iprimary_A1 +1; iprimary_A2 < nprimary; iprimary_A2++)
      {
        MyTrack mytrk_A2 = event->GetEntry(primary_id[iprimary_A2]);
        if(!(single_cut(&mytrk_A2))) continue;

        double pt_A2 = mytrk_A2.GetPtPrime();
        double px_A2 = mytrk_A2.GetPtPrime()*TMath::Cos(mytrk_A2.GetPhi0Prime());
        double py_A2 = mytrk_A2.GetPtPrime()*TMath::Sin(mytrk_A2.GetPhi0Prime());
        double pz_A2 = mytrk_A2.GetPtPrime()*(TMath::Cos(mytrk_A2.GetThe0Prime()))/(TMath::Sin(mytrk_A2.GetThe0Prime()));

        int charge_A2 = mytrk_A2.GetChargePrime();
        int emcid_A2 = mytrk_A2.GetEmcId();
        int arm_A2 = mytrk_A2.GetArm();
        int isERT_A2 = mytrk_A2.GetisERT();
	      double ecore_A2 = sqrt(pt_A2*pt_A2 + pz_A2*pz_A2 + Me2);
        double phi0_A2 = mytrk_A2.GetPhi0Prime();
	      double phi_A2 = mytrk_A2.GetPhiDC();
	      double zed_A2 = mytrk_A2.GetZDC();

        bool isert_fired_A2 = false;

	      double energy_threshold_A2 = -999;

        if(mytrk_A2.GetArm() == 1){  // West Arm

	        if (mytrk_A2.GetSect() == 0) { energy_threshold_A2 = 1.652; }
	        else if (mytrk_A2.GetSect() == 1) { energy_threshold_A2 = 1.511; }
	        else if (mytrk_A2.GetSect() == 2) { energy_threshold_A2 = 1.541; }
	        else { energy_threshold_A2 = 1.497; }

        }


        else {  // East Arm                                                                                                                               

	        if (mytrk_A2.GetSect() == 0) { energy_threshold_A2 = 2.239; }
	        else if (mytrk_A2.GetSect() == 1) { energy_threshold_A2 = 2.561; }
	        else if(mytrk_A2.GetSect() == 2) { energy_threshold_A2 = 1.676; }
	        else { energy_threshold_A2 = 1.465; }

        }

        isert_fired_A2 = isERT_A2 && ecore_A2 > energy_threshold_A2;

        if (isert && !(isert_fired_A1 || isert_fired_A2)) continue;
	      if (fabs(phi_A1-phi_A2)<PC1_DPHI_CUT && fabs(zed_A1-zed_A2)<PC1_DZ_CUT) continue;

        TLorentzVector trk_A2;
        trk_A2.SetPxPyPzE(px_A2, py_A2, pz_A2, ecore_A2);

	      TVector3 trk3_A2;
	      trk3_A2.SetXYZ(px_A2, py_A2, pz_A2);

	      if(add_photon){

          double highest_energy_cluster_A2 = -999;
          double ecore_highest_energy_cluster_A2 = -999;

          for (int iclust_A2 = 0; iclust_A2 < nclust; iclust_A2++)
          {
            MyCluster myclust_A2 = event->GetClusterEntry(iclust_A2);
            int id_clust_A2 = myclust_A2.GetID();
            if(emcid_A2 == id_clust_A2) continue;
            if(myclust_A2.GetChi2() > 3) continue;
            int arm_clust_A2 = myclust_A2.GetArm();
            double ecore_clust_A2 = myclust_A2.GetEcore();
	          if(ecore_clust_A2 < 0.3) continue;
            if(arm_A2 != arm_clust_A2) continue;

            double xclust_A2 = myclust_A2.GetX() - event->GetPreciseX();
            double yclust_A2 = myclust_A2.GetY() - event->GetPreciseY();
            double zclust_A2 = myclust_A2.GetZ() - event->GetPreciseZ(); // subtract the zvertex

            int sect_clust_A2 = myclust_A2.GetSect();

            if (sect_clust_A2<0) continue;

            TLorentzVector emcPhoton_A2;
            emcPhoton_A2.SetX( ecore_clust_A2*xclust_A2/sqrt(xclust_A2*xclust_A2 + yclust_A2*yclust_A2 + zclust_A2*zclust_A2) );
            emcPhoton_A2.SetY( ecore_clust_A2*yclust_A2/sqrt(xclust_A2*xclust_A2 + yclust_A2*yclust_A2 + zclust_A2*zclust_A2) );
            emcPhoton_A2.SetZ( ecore_clust_A2*zclust_A2/sqrt(xclust_A2*xclust_A2 + yclust_A2*yclust_A2 + zclust_A2*zclust_A2) );
            emcPhoton_A2.SetE( ecore_clust_A2 );

            double phi0_clust_A2 = emcPhoton_A2.Phi();
            if(phi0_clust_A2 < -pi/2) phi0_clust_A2 += 2*pi;

	          double window_lo = -999;
	          double window_hi = -999;

	          if(charge_A2 == -1){ window_lo = 0.0; window_hi = (mean_par0 + mean_par1/pt_A2) + 2*(sigma_par0 + sigma_par1/pt_A2);  }
	          else { window_hi = 0.0; window_lo = -1*((mean_par0 + mean_par1/pt_A2) + 2*(sigma_par0 + sigma_par1/pt_A2));   }

            if( (phi0_A2 - phi0_clust_A2) > window_lo && (phi0_A2 - phi0_clust_A2) < window_hi){ 

              if(highest_energy_cluster_A2 == -999){ highest_energy_cluster_A2 = iclust_A2; ecore_highest_energy_cluster_A2 = ecore_clust_A2; }
              if(highest_energy_cluster_A2 != -999 && ecore_clust_A2 > ecore_highest_energy_cluster_A2){ highest_energy_cluster_A2 = iclust_A2; ecore_highest_energy_cluster_A2 = ecore_clust_A2; }

            }

          }

          if(highest_energy_cluster_A2 != -999){

            MyCluster myclust_found = event->GetClusterEntry(highest_energy_cluster_A2);
            double ecore_clust_found = myclust_found.GetEcore();
            double xclust_found = myclust_found.GetX() - event->GetPreciseX();
            double yclust_found = myclust_found.GetY() - event->GetPreciseY();
            double zclust_found = myclust_found.GetZ() - event->GetPreciseZ(); // subtract the zvertex

            TLorentzVector emcPhoton_found;
            emcPhoton_found.SetX( ecore_clust_found*xclust_found/sqrt(xclust_found*xclust_found + yclust_found*yclust_found + zclust_found*zclust_found) );
            emcPhoton_found.SetY( ecore_clust_found*yclust_found/sqrt(xclust_found*xclust_found + yclust_found*yclust_found + zclust_found*zclust_found) );
            emcPhoton_found.SetZ( ecore_clust_found*zclust_found/sqrt(xclust_found*xclust_found + yclust_found*yclust_found + zclust_found*zclust_found) );
            emcPhoton_found.SetE( ecore_clust_found );

            trk_A2 += emcPhoton_found;

          }
	      }

        double R_A2 = trk_A2.Pt()/(0.003*0.90); // Radius in cm

        int hit_index_A2[total_nlayers] = {-999};

        hit_index_A2[0] = mytrk_A2.GetHitIndexL1();
        hit_index_A2[1] = mytrk_A2.GetHitIndexL2();
        hit_index_A2[2] = mytrk_A2.GetHitIndexL3();
        hit_index_A2[3] = mytrk_A2.GetHitIndexL4();
        hit_index_A2[4] = mytrk_A2.GetHitIndexL5();
        hit_index_A2[5] = mytrk_A2.GetHitIndexL6();
        hit_index_A2[6] = mytrk_A2.GetHitIndexL7();
        hit_index_A2[7] = mytrk_A2.GetHitIndexL8();

        double x1_A2 = -999; double y1_A2 = -999;
        double x2_A2 = -999; double y2_A2 = -999;

        for(int ilayer=0; ilayer<total_nlayers; ilayer++){

          if(hit_index_A2[ilayer] == -999) continue;

          MyVTXHit vtxhit = event->GetVTXHitEntry(hit_index_A2[ilayer]);
          if(x1_A2 == -999){ x1_A2 = vtxhit.GetXHit();  y1_A2 = vtxhit.GetYHit(); continue; }
          if(x2_A2 == -999){ x2_A2 = vtxhit.GetXHit();  y2_A2 = vtxhit.GetYHit(); }

        }

        double a_A2 = ((x1_A2-x2_A2) -(y1_A2*y1_A2 - y2_A2*y2_A2)/(x1_A2-x2_A2))*0.5;
        double b_A2 = (y1_A2-y2_A2)/(x1_A2-x2_A2);

        double yc1_A2 = (y1_A2-a_A2*b_A2 + sqrt(R_A2*R_A2*(1+b_A2*b_A2) - (a_A2 + b_A2*y1_A2)*(a_A2 + b_A2*y1_A2)))/(1+b_A2*b_A2);
        double yc2_A2 = (y1_A2-a_A2*b_A2 - sqrt(R_A2*R_A2*(1+b_A2*b_A2) - (a_A2 + b_A2*y1_A2)*(a_A2 + b_A2*y1_A2)))/(1+b_A2*b_A2);

        double xc1_A2 = x1_A2 - a_A2 - b_A2*yc1_A2;
        double xc2_A2 = x1_A2 - a_A2 - b_A2*yc2_A2;

        double xc_A2 = -999; double yc_A2 = -999;

        if(mytrk_A2.GetPhi0Prime() > -pi/2 && mytrk_A2.GetPhi0Prime() < pi/2 && mytrk_A2.GetChargePrime() == -1){ 
	
	        if(yc1_A2<0){ yc_A2 = yc1_A2; xc_A2 = xc1_A2; }
	        else{ yc_A2 = yc2_A2; xc_A2 = xc2_A2; }

        }

        if(mytrk_A2.GetPhi0Prime() > -pi/2 && mytrk_A2.GetPhi0Prime() < pi/2 && mytrk_A2.GetChargePrime() == 1){

          if(yc1_A2>0){ yc_A2 = yc1_A2; xc_A2 = xc1_A2; }
	        else{ yc_A2 = yc2_A2; xc_A2 = xc2_A2; }

        }

        if(mytrk_A2.GetPhi0Prime() > pi/2 && mytrk_A2.GetPhi0Prime() < 3*pi/2 && mytrk_A2.GetChargePrime() == -1){

          if(yc1_A2>0){ yc_A2 = yc1_A2; xc_A2 = xc1_A2; }
	        else{ yc_A2 = yc2_A2; xc_A2 = xc2_A2; }

        }

        if(mytrk_A2.GetPhi0Prime() > pi/2 && mytrk_A2.GetPhi0Prime() < 3*pi/2 && mytrk_A2.GetChargePrime() == 1){

          if(yc1_A2<0){ yc_A2 = yc1_A2; xc_A2 = xc1_A2; }
          else{ yc_A2 = yc2_A2; xc_A2 = xc2_A2; }

        }

        double L_A2 = sqrt((xc_A2-vtxx)*(xc_A2-vtxx) + (yc_A2-vtxy)*(yc_A2-vtxy));
        double dca_A2 = (L_A2-R_A2)*10000; // In Micro Meters

        TLorentzVector pair = trk_A1 + trk_A2;
        double pair_mass = pair.M();
        double pair_pt = pair.Pt();

	      double dca1 = sqrt( fabs(dca_A1*dca_A1 - dca_A2*dca_A2) );
	      double dca2 = sqrt( dca_A1*dca_A1 + dca_A2*dca_A2 );

        double angle = trk3_A1.Angle(trk3_A2);

        if(dca1 > 1500) continue;
        if(pair_pt > 5.0) continue;

	      if(charge_A1 == charge_A2){

          h2d_ee_FG_like_temp[0]->Fill(pair_mass, pair_pt); 

          hist_temp[0][0]->Fill(pair_mass, pair_pt, fabs(dca_A1)-fabs(dca_A2)); 
          hist_prof[0][0]->Fill(pair_mass, fabs(dca_A1)-fabs(dca_A2));
          hist_temp[1][0]->Fill(pair_mass, pair_pt, dca_A1-dca_A2); 
          hist_prof[1][0]->Fill(pair_mass, dca_A1-dca_A2);
          hist_temp[2][0]->Fill(pair_mass, pair_pt, (fabs(dca_A1)+fabs(dca_A2))/2); 
          hist_prof[2][0]->Fill(pair_mass, (fabs(dca_A1)+fabs(dca_A2))/2);
          hist_temp[3][0]->Fill(pair_mass, pair_pt, (dca_A1+dca_A2)/2); 
          hist_prof[3][0]->Fill(pair_mass, (dca_A1+dca_A2)/2);
          hist_temp[4][0]->Fill(pair_mass, pair_pt, dca1);
          hist_prof[4][0]->Fill(pair_mass, dca1);	
          hist_temp[5][0]->Fill(pair_mass, pair_pt, dca2);
          hist_prof[5][0]->Fill(pair_mass, dca2);

	      }


	      else{ 

          h2d_ee_FG_temp[0]->Fill(pair_mass, pair_pt); 
          h2d_delta_phi_mass_temp[0]->Fill(pair_mass, angle); 
          h2d_pt_spectra_mass_temp[0]->Fill(pair_mass, pt_A1);
          h2d_pt_spectra_mass_temp[0]->Fill(pair_mass, pt_A2);

	        hist_temp[0][1]->Fill(pair_mass, pair_pt, fabs(dca_A1)-fabs(dca_A2));
	        hist_prof[0][1]->Fill(pair_mass, fabs(dca_A1)-fabs(dca_A2));
	        hist_temp[1][1]->Fill(pair_mass, pair_pt, dca_A1-dca_A2);
	        hist_prof[1][1]->Fill(pair_mass, dca_A1-dca_A2);
	        hist_temp[2][1]->Fill(pair_mass, pair_pt, (fabs(dca_A1)+fabs(dca_A2))/2);
	        hist_prof[2][1]->Fill(pair_mass, (fabs(dca_A1)+fabs(dca_A2))/2);
	        hist_temp[3][1]->Fill(pair_mass, pair_pt, (dca_A1+dca_A2)/2);
	        hist_prof[3][1]->Fill(pair_mass, (dca_A1+dca_A2)/2);
	        hist_temp[4][1]->Fill(pair_mass, pair_pt, dca1);
	        hist_prof[4][1]->Fill(pair_mass, dca1);
	        hist_temp[5][1]->Fill(pair_mass, pair_pt, dca2);
	        hist_prof[5][1]->Fill(pair_mass, dca2);

	      }

      } // End Track A2 Loop
    
    } // End Track A1 Loop
  } // End Event A loop

  for(int icent=0; icent<centbin; icent++){

      h2d_ee_FG[icent] = (TH2D*)h2d_ee_FG_temp[icent]->Clone();
      h2d_ee_FG[icent]->Reset("ICESM");
      h2d_ee_FG[icent]->SetName(Form("h2d_ee_FG_cent%d",icent));

      int nbinsX = h2d_ee_FG_temp[icent]->GetNbinsX();
      int nbinsY = h2d_ee_FG_temp[icent]->GetNbinsY();

      for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

        double binwidthX = h2d_ee_FG_temp[icent]->GetXaxis()->GetBinWidth(ibinX);

        for(int ibinY=1; ibinY < nbinsY+1; ibinY++){
 
          double binwidthY = h2d_ee_FG_temp[icent]->GetYaxis()->GetBinWidth(ibinY);
          double binwidth = binwidthX*binwidthY;
          double content = h2d_ee_FG_temp[icent]->GetBinContent(ibinX, ibinY);
          double err = h2d_ee_FG_temp[icent]->GetBinError(ibinX,ibinY);
          h2d_ee_FG[icent]->SetBinContent(ibinX, ibinY, content/binwidth);
          h2d_ee_FG[icent]->SetBinError(ibinX, ibinY, err/binwidth);

        }
      }
    }

    for(int icent=0; icent<centbin; icent++){

      h2d_ee_FG_like[icent] = (TH2D*)h2d_ee_FG_like_temp[icent]->Clone();
      h2d_ee_FG_like[icent]->Reset("ICESM");
      h2d_ee_FG_like[icent]->SetName(Form("h2d_ee_FG_like_cent%d",icent));

      int nbinsX = h2d_ee_FG_like_temp[icent]->GetNbinsX();
      int nbinsY = h2d_ee_FG_like_temp[icent]->GetNbinsY();

      for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

        double binwidthX = h2d_ee_FG_like_temp[icent]->GetXaxis()->GetBinWidth(ibinX);

        for(int ibinY=1; ibinY < nbinsY+1; ibinY++){

          double binwidthY = h2d_ee_FG_like_temp[icent]->GetYaxis()->GetBinWidth(ibinY);
          double binwidth = binwidthX*binwidthY;
          double content = h2d_ee_FG_like_temp[icent]->GetBinContent(ibinX, ibinY);
          double err = h2d_ee_FG_like_temp[icent]->GetBinError(ibinX,ibinY);
          h2d_ee_FG_like[icent]->SetBinContent(ibinX, ibinY, content/binwidth);
          h2d_ee_FG_like[icent]->SetBinError(ibinX, ibinY, err/binwidth);

        }
      }
    }

    for(int icent=0; icent<centbin; icent++){

      h2d_pt_spectra_mass[icent] = (TH2D*)h2d_pt_spectra_mass_temp[icent]->Clone();
      h2d_pt_spectra_mass[icent]->Reset("ICESM");
      h2d_pt_spectra_mass[icent]->SetName(Form("h2d_pt_spectra_mass_cent%d",icent));

      int nbinsX = h2d_pt_spectra_mass_temp[icent]->GetNbinsX();
      int nbinsY = h2d_pt_spectra_mass_temp[icent]->GetNbinsY();

      for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

        double binwidthX = h2d_pt_spectra_mass_temp[icent]->GetXaxis()->GetBinWidth(ibinX);

        for(int ibinY=1; ibinY < nbinsY+1; ibinY++){

          double binwidthY = h2d_pt_spectra_mass_temp[icent]->GetYaxis()->GetBinWidth(ibinY);
          double binwidth = binwidthX*binwidthY;
          double content = h2d_pt_spectra_mass_temp[icent]->GetBinContent(ibinX, ibinY);
          double err = h2d_pt_spectra_mass_temp[icent]->GetBinError(ibinX,ibinY);
          h2d_pt_spectra_mass[icent]->SetBinContent(ibinX, ibinY, content/binwidth);
          h2d_pt_spectra_mass[icent]->SetBinError(ibinX, ibinY, err/binwidth);

        }
      }
    }

     for(int icent=0; icent<centbin; icent++){

      h2d_delta_phi_mass[icent] = (TH2D*)h2d_delta_phi_mass_temp[icent]->Clone();
      h2d_delta_phi_mass[icent]->Reset("ICESM");
      h2d_delta_phi_mass[icent]->SetName(Form("h2d_delta_phi_mass_cent%d",icent));

      int nbinsX = h2d_delta_phi_mass_temp[icent]->GetNbinsX();
      int nbinsY = h2d_delta_phi_mass_temp[icent]->GetNbinsY();

      for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

        double binwidthX = h2d_delta_phi_mass_temp[icent]->GetXaxis()->GetBinWidth(ibinX);

        for(int ibinY=1; ibinY < nbinsY+1; ibinY++){

          double binwidthY = h2d_delta_phi_mass_temp[icent]->GetYaxis()->GetBinWidth(ibinY);
          double binwidth = binwidthX*binwidthY;
          double content = h2d_delta_phi_mass_temp[icent]->GetBinContent(ibinX, ibinY);
          double err = h2d_delta_phi_mass_temp[icent]->GetBinError(ibinX,ibinY);
          h2d_delta_phi_mass[icent]->SetBinContent(ibinX, ibinY, content/binwidth);
          h2d_delta_phi_mass[icent]->SetBinError(ibinX, ibinY, err/binwidth);

        }
      }
    }

    for(int icent=0; icent<centbin; icent++){

      h2d_ee_FG[icent]->GetYaxis()->SetRangeUser(0.,10.);
      h1d_ee_FG[icent] = (TH1D*)h2d_ee_FG[icent]->ProjectionX();
      h1d_ee_FG[icent]->SetName(Form("h1d_ee_FG_cent%d",icent));

      h2d_ee_FG_like[icent]->GetYaxis()->SetRangeUser(0.,10.);
      h1d_ee_FG_like[icent] = (TH1D*)h2d_ee_FG_like[icent]->ProjectionX();
      h1d_ee_FG_like[icent]->SetName(Form("h1d_ee_FG_like_cent%d",icent));


    }


    for(int i=0; i < ntype; i++){
      for(int j=0; j<ncharge; j++){

        hist[i][j] = (TH3D*)hist_temp[i][j]->Clone();
        hist[i][j]->Reset("ICESM");
        hist[i][j]->SetName(Form("hist_%s_%s", description[i], charge_type[j]));

        int nbinsX = hist_temp[i][j]->GetNbinsX();
        int nbinsY = hist_temp[i][j]->GetNbinsY();
        int nbinsZ = hist_temp[i][j]->GetNbinsZ();

        for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

          double binwidthX = hist_temp[i][j]->GetXaxis()->GetBinWidth(ibinX);

          for(int ibinY=1; ibinY < nbinsY+1; ibinY++){

              double binwidthY = hist_temp[i][j]->GetYaxis()->GetBinWidth(ibinY);

            for(int ibinZ=1; ibinZ < nbinsZ+1; ibinZ++){

              double binwidthZ = hist_temp[i][j]->GetZaxis()->GetBinWidth(ibinZ);
              double binwidth = binwidthX*binwidthY*binwidthZ;
              double content = hist_temp[i][j]->GetBinContent(ibinX, ibinY, ibinZ);
              double err = hist_temp[i][j]->GetBinError(ibinX,ibinY, ibinZ);
              hist[i][j]->SetBinContent(ibinX, ibinY, ibinZ, content/binwidth);
              hist[i][j]->SetBinError(ibinX, ibinY, ibinZ, err/binwidth);

            }
          }
        }
      }
    }



  output->cd();
  for(int i=0; i<ntype; i++){
    for(int j=0; j<ncharge; j++){
      hist[i][j]->Write();
      hist_prof[i][j]->Write();
    }
  }

  for(int icent =0; icent < centbin; icent++){
    
    h2d_ee_FG[icent]->Write();
    h1d_ee_FG[icent]->Write();

    h2d_ee_FG_like[icent]->Write();
    h1d_ee_FG_like[icent]->Write();

    h2d_delta_phi_mass[icent]->Write();
    h2d_pt_spectra_mass[icent]->Write();
  }


  output->Close();

} 
