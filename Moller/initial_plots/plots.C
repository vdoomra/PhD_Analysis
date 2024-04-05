// Cretaed for the MOLLER Optics study by Vassu Doomra Date: April 17, 2022


#include "TH2D.h"

const int nFiles = 1;
const int nPlots = 29;
const char* filename[nFiles] = { "../slim-output/slim_moller_optics1DS_p5.root"};

const char* passNames[nFiles] = {"pass5"};
const int color[nFiles] = {kRed};
const char* plotName[nPlots] = { "main_xy","sieve_xy","gem1_xy","gem2_xy","gem3_xy","gem4_xy","gem1_r","gem1_r'","gem1_phi","gem1_phi'",
				                          "theta_tg_vs_r_gem1","theta_tg_vs_phi_gem1","theta_tg_vs_r'_gem1","theta_tg_vs_phi'_gem1" ,
                                  "tg_phi_vs_r_gem1","tg_phi_vs_phi_gem1","tg_phi_vs_r'_gem1","tg_phi_vs_phi'_gem1",
                                  "tg_p_vs_r_gem1","tg_p_vs_phi_gem1","tg_p_vs_r'_gem1","tg_p_vs_phi'_gem1" , "sieve_r",
				                          "sieve_phi", "gem1_k_gem1_r", "main_r_gem_r", "profile_sieve_gem_r", "sieve_ph_gem_ph", "gem_k_r"};

TCanvas* c[nPlots];
TCanvas* c1;

TH1D* gr[nFiles] = {NULL};    // r distribution at gem1 plane
TH1D* grp[nFiles] = {NULL};   // r' distribution at gem1 plane
TH1D* gphi[nFiles] = {NULL};  // phi distribution at gem1 plane
TH1D* gphip[nFiles] = {NULL}; // phi' distribution at gem1 plane
TH1D* siever[nFiles] = {NULL};    // r distribution at sieve plane
TH1D* sievephi[nFiles] = {NULL};   // phi distribution at sieve plane

TH2D* gxy1[nFiles] = {NULL};  // xy distribution at gem1 plane
TH2D* gxy2[nFiles] = {NULL};  // xy distribution at gem2 plane
TH2D* gxy3[nFiles] = {NULL};  // xy distribution at gem3 plane
TH2D* gxy4[nFiles] = {NULL};  // xy distribution at gem3 plane

TH2D* Etheta[nFiles] = {NULL};  // E vs theta_tg
TH2D* Etheta_weighted[nFiles] = {NULL};  // E vs theta_tg
TH2D* sievexy[nFiles] = {NULL}; // xy distribution at the sieve plane
TH2D* mainxy[nFiles] = {NULL};  // xy distribution at the main detector plane

TH2D* thetar[nFiles] = {NULL};  // theta_tg vs r at gem1
TH2D* thetarp[nFiles] = {NULL};  // theta_tg vs r' at gem1
TH2D* thetaphi[nFiles] = {NULL}; // theta_tg vs phi at gem1
TH2D* thetaphip[nFiles] = {NULL}; // theta_tg vs phi' at gem1

TH2D* phir[nFiles] = {NULL};  // phi_tg vs r at gem1
TH2D* phirp[nFiles] = {NULL};  // phi_tg vs r' at gem1
TH2D* phiphi[nFiles] = {NULL}; // phi_tg vs phi at gem1
TH2D* phiphip[nFiles] = {NULL}; // phi_tg vs phi' at gem1

TH2D* pr[nFiles] = {NULL};  // phi_tg vs r at gem1
TH2D* prp[nFiles] = {NULL};  // phi_tg vs r' at gem1
TH2D* pphi[nFiles] = {NULL}; // phi_tg vs phi at gem1
TH2D* pphip[nFiles] = {NULL}; // phi_tg vs phi' at gem1

TH2D* gemrphi[nFiles] = {NULL};
TH2D* gemkr[nFiles] = {NULL};

TH2D* sieve_gem_r[nFiles] = {NULL};
TProfile* prof_sieve_gem_r[nFiles] = {NULL};
TH2D* sieve_gem_ph[nFiles] = {NULL};


const double pi = acos(-1);

void plots(){

  for(int k=0; k<nPlots; k++){

    c[k] = new TCanvas(Form("%s",plotName[k])); c[k]->Divide(3,2);}

     c1 = new TCanvas("E_vs_theta"); c1->Divide(1,2);

  for(int i=0; i<nFiles; i++){

    TFile* fin = new TFile(Form("%s",filename[i]),"READ");

    TTree* T = (TTree*)fin->Get("newT");

    int nEntries = T->GetEntries();
 
    double main_x, main_y, main_r, sieve_x, sieve_y, sieve_r, sieve_ph;
    double gem1_x, gem1_y, gem2_x, gem2_y, gem3_x, gem3_y, gem4_x, gem4_y, gem1_k;
    double gem1_r, gem1_rp, gem1_phi, gem1_phip, gem1_px, gem1_py, gem1_pz;
    double tg_th, tg_ph; double tg_p; double rate;

    T->SetBranchAddress("main_x",&main_x);
    T->SetBranchAddress("main_r",&main_r);
    T->SetBranchAddress("main_y",&main_y);
    T->SetBranchAddress("gem1_r",&gem1_r);
    T->SetBranchAddress("gem1_ph",&gem1_phi);
    T->SetBranchAddress("gem1_x",&gem1_x);
    T->SetBranchAddress("gem1_k",&gem1_k);
    T->SetBranchAddress("gem1_y",&gem1_y);
    T->SetBranchAddress("gem2_x",&gem2_x);
    T->SetBranchAddress("gem2_y",&gem2_y);
    T->SetBranchAddress("gem3_x",&gem3_x);
    T->SetBranchAddress("gem3_y",&gem3_y);
    T->SetBranchAddress("gem4_x",&gem4_x);
    T->SetBranchAddress("gem4_y",&gem4_y);
    T->SetBranchAddress("gem1_px",&gem1_px);
    T->SetBranchAddress("gem1_py",&gem1_py);
    T->SetBranchAddress("gem1_pz",&gem1_pz);
    T->SetBranchAddress("sieve_x",&sieve_x);
    T->SetBranchAddress("sieve_y",&sieve_y);
    T->SetBranchAddress("tg_th",&tg_th);
    T->SetBranchAddress("tg_ph",&tg_ph);
    T->SetBranchAddress("tg_p",&tg_p);
    T->SetBranchAddress("rate",&rate);
    T->SetBranchAddress("sieve_r",&sieve_r);
    T->SetBranchAddress("sieve_ph",&sieve_ph);
   

    gr[i] = new TH1D(Form("gem1_r_%s",passNames[i]),"gem1_r_distribution;r[mm];Counts",100,500,1300);
    grp[i] = new TH1D(Form("gem1_r'_%s",passNames[i]),"gem1_r_prime_distribution;r';Counts",100,0,0.1);
    gphi[i] = new TH1D(Form("gem1_phi_%s",passNames[i]),"gem1_phi_distribution;#phi(rad);Counts",200,0,2*pi);
    gphip[i] = new TH1D(Form("gem1_phi'_%s",passNames[i]),"gem1_phi_prime_distribution;#phi';Counts",1000,-0.1,0.1);
    Etheta[i] = new TH2D(Form("Etheta_%s",passNames[i]),"E' vs #theta;#theta_{tg};E(MeV)",100,0,1.8,100,0,12);
    Etheta_weighted[i] = new TH2D(Form("Etheta_weighted_%s",passNames[i]),"E' vs #theta weighted;#theta_{tg};E(MeV)",100,0,1.8,100,0,12);
    sievexy[i] = new TH2D(Form("sieve_xy_%s",passNames[i]),"sieve_xy_distribution;x[mm];y[mm]",500,-100,100,500,-100,100);
    siever[i] = new TH1D(Form("sieve_r_%s",passNames[i]),"sieve_r_distribution;r[mm];counts",200,0,100);
    sievephi[i] = new TH1D(Form("sieve_phi_%s",passNames[i]),"sieve_phi_distribution;#phi(rad);counts",200,0,2*pi);
    mainxy[i] = new TH2D(Form("main_xy_%s",passNames[i]),"main_xy_distribution;x[mm];y[mm]",1000,-1500,1500,1000,-1500,1500);
    gxy1[i] = new TH2D(Form("gem1_xy_%s",passNames[i]),"gem1_xy_distribution;x[mm];y[mm]",1000,-1500,1500,1000,-1500,1500);
    gxy2[i] = new TH2D(Form("gem2_xy_%s",passNames[i]),"gem2_xy_distribution;x[mm];y[mm]",1000,-1500,1500,1000,-1500,1500);
    gxy3[i] = new TH2D(Form("gem3_xy_%s",passNames[i]),"gem3_xy_distribution;x[mm];y[mm]",1000,-1500,1500,1000,-1500,1500);
    gxy4[i] = new TH2D(Form("gem4_xy_%s",passNames[i]),"gem4_xy_distribution;x[mm];y[mm]",1000,-1500,1500,1000,-1500,1500);
    thetar[i] = new TH2D(Form("theta_r_gem1_%s",passNames[i]),"theta_r_gem1_distribution;r[mm];#theta_{tg}",100, 500, 1300, 50, 0, 1.8);
    thetarp[i] = new TH2D(Form("theta_r'_gem1_%s",passNames[i]),"theta_r_prime_gem1_distribution;r';#theta_{tg}",100, 0, 0.1, 50, 0, 1.8);
    thetaphi[i] = new TH2D(Form("theta_phi_sieve_%s",passNames[i]),"theta_phi_sieve_distribution;#phi;#theta_{tg}",200, 0, 2*pi, 50, 0, 1.8);
    thetaphip[i] = new TH2D(Form("theta_phi'_gem1_%s",passNames[i]),"theta_phi_prime_gem1_distribution;#phi';#theta_{tg}",1000, -0.1, 0.1, 50, 0, 1.8);
    phir[i] = new TH2D(Form("phi_r_gem1_%s",passNames[i]),"phi_r_gem1_distribution;#phi_{tg};r[mm]", 200, 0, 2*pi,100, 500, 1300);
    phirp[i] = new TH2D(Form("phi_r'_gem1_%s",passNames[i]),"phi_r_prime_gem1_distribution;#phi_{tg};r'", 200,0, 2*pi,100, 0, 0.1);
    phiphi[i] = new TH2D(Form("phi_phi_gem1_%s",passNames[i]),"phi_phi_gem1_distribution;#phi_{tg};#phi",200, 0, 2*pi, 100, -pi, pi);
    phiphip[i] = new TH2D(Form("phi_phi'_gem1_%s",passNames[i]),"phi_phi_prime_gem1_distribution;#phi_{tg};#phi'",200, 0, 2*pi,1000, -0.01, 0.01);

    pr[i] = new TH2D(Form("tg_p_r_gem1_%s",passNames[i]),"tg_p_r_gem1_distribution;r[mm];E'(GeV)",1000, 500, 1300, 100, 0, 12);
    prp[i] = new TH2D(Form("tg_p_r'_gem1_%s",passNames[i]),"tg_p_r_prime_gem1_distribution;r';E'(GeV)",1000, 0, 0.1, 100, 0, 12);
    pphi[i] = new TH2D(Form("tg_p_phi_gem1_%s",passNames[i]),"tg_p_phi_gem1_distribution;#phi;E'(GeV)", 200, 0, 2*pi, 100, 0, 12);
    pphip[i] = new TH2D(Form("tg_p_phi'_gem1_%s",passNames[i]),"tg_p_phi_prime_gem1_distribution;#phi';E'(GeV)",1000, -0.006, 0.006, 100, 0, 12);

    gemrphi[i] = new TH2D(Form("gem_r_phi_%s",passNames[i]),"gem_r_phi_distribution;r(mm);phi(rad)",800, 500, 1300, 100, 0, 7);
    gemkr[i] = new TH2D(Form("gem_r_k_%s",passNames[i]),"gem_r_k_distribution;r(mm);k(GeV)",800, 500, 1300, 100, 0,  12);
    sieve_gem_r[i] = new TH2D(Form("sieve_main_r_%s",passNames[i]),"sieve_main_r_distribution;r at sieve(mm);r at main(mm)",200, 0, 100, 800, 500, 1300);
    prof_sieve_gem_r[i] = new TProfile(Form("prof_sieve_gem_r_%s",passNames[i]),"sieve_gem_r_distribution;r at sieve(mm);r at gem(mm)",200, 0, 100, 500, 1300);
    sieve_gem_ph[i] = new TH2D(Form("sieve_gem_ph_%s",passNames[i]),"sieve_gem_ph_distribution;ph at sieve(rad);ph at gem(ph)",200, 0, 2*pi, 200, 0, 2*pi);

    for(int j=0; j<nEntries; j++){  

      T->GetEntry(j);
      //rate = 1; //Beam Gnerator

      if(gem1_r < 500) continue;
   
      double rprime = (gem1_x*gem1_px + gem1_y*gem1_py)/(gem1_r*gem1_pz);
      double phip = (-gem1_y*gem1_px + gem1_x*gem1_py)/(gem1_r*gem1_pz); 
      if(sieve_ph<0) sieve_ph = sieve_ph + 2*pi;
      if(gem1_phi < 0) gem1_phi += 2*pi;


      gr[i]->Fill(gem1_r,rate);  grp[i]->Fill(rprime,rate); gphi[i]->Fill(gem1_phi,rate); gphip[i]->Fill(phip,rate); 
      Etheta[i]->Fill(tg_th*(180/pi),tg_p/1000);  Etheta_weighted[i]->Fill(tg_th*(180/pi),tg_p/1000,rate); 
      sievexy[i]->Fill(sieve_x,sieve_y,tg_th*(180/pi));   mainxy[i]->Fill(main_x,main_y,rate);
      gxy1[i]->Fill(gem1_x,gem1_y,rate);  gxy2[i]->Fill(gem2_x,gem2_y,rate);  gxy3[i]->Fill(gem3_x,gem3_y,rate);  gxy4[i]->Fill(gem4_x,gem4_y,rate);
      thetar[i]->Fill(gem1_r,tg_th*(180/pi),rate); thetaphi[i]->Fill(sieve_ph,tg_th*(180/pi),rate);
      thetarp[i]->Fill(rprime,tg_th*(180/pi),rate); thetaphip[i]->Fill(phip,tg_th*(180/pi),rate);
      phir[i]->Fill(tg_ph,gem1_r,rate); phiphi[i]->Fill(tg_ph,gem1_phi,rate);
      phirp[i]->Fill(tg_ph,rprime,rate); phiphip[i]->Fill(tg_ph,phip,rate);

      pr[i]->Fill(gem1_r,tg_p/1000,rate); pphi[i]->Fill(gem1_phi,tg_p/1000,rate);
      prp[i]->Fill(rprime,tg_p/1000,rate); pphip[i]->Fill(phip,tg_p/1000,rate);

      siever[i]->Fill(sieve_r,rate); sievephi[i]->Fill(sieve_ph,rate);

      gemrphi[i]->Fill(gem1_r, gem1_phi, rate);
      gemkr[i]->Fill(gem1_r, sqrt(gem1_px*gem1_px + gem1_py*gem1_py + gem1_pz*gem1_pz)/1000, rate);

      sieve_gem_r[i]->Fill(sieve_r, main_r, rate);
      prof_sieve_gem_r[i]->Fill(sieve_r, gem1_r, rate);

      sieve_gem_ph[i]->Fill(sieve_ph, gem1_phi, rate);


    }

      c[0]->cd(i+1); mainxy[i]->Draw("colz"); 
      TLatex* t1 = new TLatex();
      t1->SetTextAlign(11);
      t1->SetTextSize(0.04);
      t1->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));
         

      c[1]->cd(i+1); sievexy[i]->Draw("colz"); 
      TLatex* t2 = new TLatex();
      t2->SetTextAlign(11);
      t2->SetTextSize(0.04);
      t2->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[2]->cd(i+1); gxy1[i]->Draw("colz"); 
      TLatex* t3 = new TLatex();
      t3->SetTextAlign(11);
      t3->SetTextSize(0.04);
      t3->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i])); 

      c[3]->cd(i+1); gxy2[i]->Draw("colz"); 
      TLatex* t4 = new TLatex();
      t4->SetTextAlign(11);
      t4->SetTextSize(0.04);
      t4->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[4]->cd(i+1); gxy3[i]->Draw("colz"); 
      TLatex* t5 = new TLatex();
      t5->SetTextAlign(11);
      t5->SetTextSize(0.04);
      t5->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i])); 
  
      c[5]->cd(i+1); gxy4[i]->Draw("colz"); 
      TLatex* t6 = new TLatex();
      t6->SetTextAlign(11);
      t6->SetTextSize(0.04);
      t6->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[6]->cd(i+1); gr[i]->Draw("hist"); 
      TLatex* t7 = new TLatex();
      t7->SetTextAlign(11);
      t7->SetTextSize(0.04);
      t7->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i])); 

      c[7]->cd(i+1); grp[i]->Draw("hist"); 
      TLatex* t8 = new TLatex();
      t8->SetTextAlign(11);
      t8->SetTextSize(0.04);
      t8->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i])); 

      c[8]->cd(i+1); gphi[i]->Draw("hist"); 
      TLatex* t9 = new TLatex();
      t9->SetTextAlign(11);
      t9->SetTextSize(0.04);
      t9->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i])); 
  
      c[9]->cd(i+1); gphip[i]->Draw("hist"); 
      TLatex* t10 = new TLatex();
      t10->SetTextAlign(11);
      t10->SetTextSize(0.04);
      t10->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i])); 

      c[10]->cd(i+1); thetar[i]->Draw("colz"); 
      TLatex* t11 = new TLatex();
      t11->SetTextAlign(11);
      t11->SetTextSize(0.04);
      t11->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[11]->cd(i+1); thetaphi[i]->Draw("colz"); 
      TLatex* t12 = new TLatex();
      t12->SetTextAlign(11);
      t12->SetTextSize(0.04);
      t12->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));  

      c[12]->cd(i+1); thetarp[i]->Draw("colz"); 
      TLatex* t13 = new TLatex();
      t13->SetTextAlign(11);
      t13->SetTextSize(0.04);
      t13->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i])); 
  
      c[13]->cd(i+1); thetaphip[i]->Draw("colz"); 
      TLatex* t14 = new TLatex();
      t14->SetTextAlign(11);
      t14->SetTextSize(0.04);
      t14->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[14]->cd(i+1); phir[i]->Draw("colz"); 
      TLatex* t15 = new TLatex();
      t15->SetTextAlign(11);
      t15->SetTextSize(0.04);
      t15->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[15]->cd(i+1); phirp[i]->Draw("colz"); 
      TLatex* t16 = new TLatex();
      t16->SetTextAlign(11);
      t16->SetTextSize(0.04);
      t16->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[16]->cd(i+1); phiphi[i]->Draw("colz"); 
      TLatex* t17 = new TLatex();
      t17->SetTextAlign(11);
      t17->SetTextSize(0.04);
      t17->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));


      c[17]->cd(i+1); phiphip[i]->Draw("colz"); 
      TLatex* t18 = new TLatex();
      t18->SetTextAlign(11);
      t18->SetTextSize(0.04);
      t18->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[18]->cd(i+1); pr[i]->Draw("colz"); 
      TLatex* t19 = new TLatex();
      t19->SetTextAlign(11);
      t19->SetTextSize(0.04);
      t19->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[19]->cd(i+1); prp[i]->Draw("colz"); 
      TLatex* t20 = new TLatex();
      t20->SetTextAlign(11);
      t20->SetTextSize(0.04);
      t20->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[20]->cd(i+1); pphi[i]->Draw("colz"); 
      TLatex* t21 = new TLatex();
      t21->SetTextAlign(11);
      t21->SetTextSize(0.04);
      t21->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[21]->cd(i+1); pphip[i]->Draw("colz"); 
      TLatex* t22 = new TLatex();
      t22->SetTextAlign(11);
      t22->SetTextSize(0.04);
      t22->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[22]->cd(i+1); siever[i]->Draw("hist"); 
      TLatex* t23 = new TLatex();
      t23->SetTextAlign(11);
      t23->SetTextSize(0.04);
      t23->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[23]->cd(i+1); sievephi[i]->Draw("hist"); 
      TLatex* t24 = new TLatex();
      t24->SetTextAlign(11);
      t24->SetTextSize(0.04);
      t24->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[24]->cd(i+1); gemrphi[i]->Draw("colz"); 
      TLatex* t25 = new TLatex();
      t25->SetTextAlign(11);
      t25->SetTextSize(0.04);
      t25->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[25]->cd(i+1); sieve_gem_r[i]->Draw("colz"); 
      TLatex* t26 = new TLatex();
      t26->SetTextAlign(11);
      t26->SetTextSize(0.04);
      t26->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[26]->cd(i+1); prof_sieve_gem_r[i]->Draw("colz"); 
      TLatex* t27 = new TLatex();
      t27->SetTextAlign(11);
      t27->SetTextSize(0.04);
      t27->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[27]->cd(i+1); sieve_gem_ph[i]->Draw("colz"); 
      TLatex* t28 = new TLatex();
      t28->SetTextAlign(11);
      t28->SetTextSize(0.04);
      t28->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c[28]->cd(i+1); gemkr[i]->Draw("colz"); 
      TLatex* t29 = new TLatex();
      t29->SetTextAlign(11);
      t29->SetTextSize(0.04);
      t29->DrawLatexNDC(0.19,0.85,Form("%s",passNames[i]));

      c1->cd(1); if(i==0) Etheta[i]->Draw("p");
      else Etheta[i]->Draw("same p");
      Etheta[i]->SetMarkerColor(color[i]);


      c1->cd(2); if(i==0) Etheta_weighted[i]->Draw("p");
      else Etheta_weighted[i]->Draw("same p");
      Etheta_weighted[i]->SetMarkerColor(color[i]); 

  }

  TFile* fout = new TFile("initial_plots_moller_optics1DS_p5.root","RECREATE");
  fout->cd();

  for(int i=0; i<nFiles; i++){

    Etheta[i]->Write(); }

  for(int i=0; i<nPlots; i++){

    c[i]->Write();
 
    }

  fout->Close();

  }

     
