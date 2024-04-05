#include <iostream>
#include <fstream>
#include <TF1.h>
#include <TF2.h>

const int npar = 6;
const double Me2 = 0.511*0.511;

// Theta Parameters
const double par[npar] = { -0.01688, 1.897e-4, -1.7333, -8.478e-8, -1.549e-4, 10.453 };

void plotting_residuals(){

    TH1D* hist = new TH1D("hist","Residual Distribution",1000,-0.5,0.5);
    TH1D* hist_true = new TH1D("hist_true","Residual Distribution",1000,-0.5,0.5);
    TH2D* h2d_residual_vs_sieve_r = new TH2D("h2d_residual_vs_sieve_r","Residuals vs sieve r;Residuals[rad];sieve_{r}[mm]", 200, -0.01, 0.01, 500, 30, 100);
    TH2D* h2d_residual_vs_gem_ph = new TH2D("h2d_residual_vs_gem_ph","Residuals vs gem phi;Residuals[rad];gem_{#phi}[rad]", 200, -0.01, 0.01, 500, -0.5, 0.5);
    
    TF2* func = new TF2("func","[0]+ [1]*x + [2]*y + [3]*x*x + [4]*x*y + [5]*y*y", 0, 1100, 0, 0.1);

    func->SetParameters(par); 

    std::ifstream inFile("csv_output/opticsDS_p2_non_radiative/C12_opticsDS_p2_merged.csv");
    int index;
    double tg_th, tg_ph, tg_vz, tg_p, gem1_r, gem1_rp, gem1_ph, gem1_php, gem1_ph_local, sieve_r, sieve_ph, rate;

    while(!inFile.eof()){

        inFile >> index >> tg_th >> tg_ph >> tg_vz >> tg_p >> gem1_r >> gem1_rp >> gem1_ph >>  gem1_php >> gem1_ph_local >> sieve_r >> sieve_ph >> rate ;

        hist->Fill(func->Eval(gem1_r, gem1_rp));
        hist_true->Fill(tg_th);

        h2d_residual_vs_gem_ph->Fill(tg_th-func->Eval(gem1_r, gem1_rp), gem1_ph_local);
        h2d_residual_vs_sieve_r->Fill(tg_th-func->Eval(gem1_r, gem1_rp), sieve_r);
        


    }

    TF1* g1 = new TF1("g1", "gaus", -0.0025, 0.0025);
    g1->SetParameter(1, 0.0);
    g1->SetParameter(2, 0.0005);
    hist->Fit(g1);

    TF1* g2 = new TF1("g2", "gaus", -0.0025, 0.0025);
    g2->SetParameter(1, 0.0);
    g2->SetParameter(2, 0.002);
    hist->Fit(g2);

    double gaus_par[6];
    g1->GetParameters(&gaus_par[0]);
    g2->GetParameters(&gaus_par[3]);

    TF1* fit_func = new TF1("fit_func", "gaus(0) + gaus(3)", -0.0025, 0.0025);
    fit_func->SetParameters(gaus_par);
    hist->Fit(fit_func, "R");
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    hist->Draw();
    hist_true->Draw("same");

    TCanvas* c2 = new TCanvas();
    c2->cd();
    h2d_residual_vs_sieve_r->Draw("colz");

    TCanvas* c3 = new TCanvas();
    c3->cd();
    h2d_residual_vs_gem_ph->Draw("colz");

    TFile* fout = new TFile("residuals_opticsDS_p1_non_radiative.root","RECREATE");
    fout->cd();
    hist->Write();
    h2d_residual_vs_sieve_r->Write();
    h2d_residual_vs_gem_ph->Write();
    fout->Close();




}