#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TH1.h>

//===============PPG037==================//

const int npts_ppg037 = 12;

double x_ppg037[npts_ppg037] = {0.45, 0.55, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.25, 2.75, 3.5, 4.5};
double y_ppg037[npts_ppg037] = {1.453e-02,6.142e-03,2.986e-03,8.341e-04,3.441e-04,1.251e-04,1.119e-04,3.875e-05,1.371e-05,2.581e-06,8.401e-07,1.895e-07};

double x_lohi_ppg037[npts_ppg037] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03};

double stat_lo_ppg037[npts_ppg037] = {3.322e-03,1.923e-03,7.183e-04,3.340e-04,1.659e-04,6.671e-05,2.952e-05,7.555e-06,2.380e-06,8.910e-07,2.714e-07,1.053e-07};
double stat_hi_ppg037[npts_ppg037] = {3.322e-03,1.923e-03,7.183e-04,3.340e-04,1.659e-04,6.671e-05,2.952e-05,7.555e-06,2.380e-06, 1.065e-06, 3.345e-07, 1.632e-07};

double sys_lo_ppg037[npts_ppg037] = {1.140e-02,5.006e-03,1.829e-03,5.339e-04,1.895e-04,7.211e-05,3.735e-05,1.230e-05,3.428e-06,7.133e-07,1.779e-07,3.536e-08};
double sys_hi_ppg037[npts_ppg037] = {1.140e-02,5.006e-03,1.829e-03,5.339e-04,1.895e-04,7.211e-05,3.735e-05,1.230e-05,3.428e-06,7.133e-07,1.779e-07,3.536e-08};

//===============PPG065==================//

const int npts_ppg065 = 28;

double x_ppg065[npts_ppg065] = {3.50e-01,4.50e-01,5.50e-01,6.50e-01,7.50e-01,8.50e-01,9.50e-01,1.10e+00,1.30e+00,1.50e+00,1.70e+00,1.90e+00,2.10e+00,2.30e+00,2.50e+00,
2.70e+00,2.90e+00,3.10e+00,3.30e+00,3.50e+00,3.70e+00,3.90e+00,4.25e+00,4.75e+00,5.50e+00,6.50e+00,7.50e+00,8.50e+00};
double y_ppg065[npts_ppg065] = {1.36e-02,6.05e-03,3.36e-03,1.81e-03,1.30e-03,7.50e-04,5.61e-04,3.18e-04,1.26e-04,6.58e-05,3.97e-05,1.99e-05,1.14e-05,6.83e-06,3.98e-06,
2.44e-06,1.63e-06,1.05e-06,7.21e-07,5.04e-07,3.45e-07,2.37e-07,1.33e-07,6.13e-08,2.06e-08,6.62e-09,1.70e-09,9.72e-10};

double x_lohi_ppg065[npts_ppg065] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03};

double stat_lo_ppg065[npts_ppg065] = {4.45e-03,1.70e-03,7.56e-04,3.98e-04,2.28e-04,1.35e-04,9.14e-05,3.77e-05,2.10e-05,1.25e-05,2.07e-06,1.18e-06,7.44e-07,4.96e-07,3.63e-07,6.49e-08,4.77e-08,3.53e-08,2.73e-08,2.15e-08,1.70e-08,1.36e-08,5.95e-09,3.73e-09,1.79e-09,9.26e-10,5.73e-10,5.16e-10};
double stat_hi_ppg065[npts_ppg065] = {4.45e-03,1.70e-03,7.56e-04,3.98e-04,2.28e-04,1.35e-04,9.14e-05,3.77e-05,2.10e-05,1.25e-05,2.07e-06,1.18e-06,7.44e-07,4.96e-07,3.63e-07,6.49e-08,4.77e-08,3.53e-08,2.73e-08,2.15e-08,1.70e-08,1.36e-08,5.95e-09,3.73e-09,1.79e-09,9.26e-10,5.73e-10,5.16e-10};

double sys_lo_ppg065[npts_ppg065] = {5.95e-03,2.39e-03,1.04e-03,5.12e-04,2.59e-04,1.35e-04,7.85e-05,3.64e-05,1.50e-05,6.56e-06,3.32e-06,1.65e-06,8.96e-07,5.09e-07,3.05e-07,
2.80e-07,1.77e-07,1.10e-07,7.31e-08,4.89e-08,3.30e-08,2.25e-08,1.22e-08,5.45e-09,1.85e-09,5.72e-10,1.60e-10,8.07e-11};
double sys_hi_ppg065[npts_ppg065] = {5.95e-03,2.39e-03,1.04e-03,5.12e-04,2.59e-04,1.35e-04,7.85e-05,3.64e-05,1.50e-05,6.56e-06,3.32e-06,1.65e-06,8.96e-07,5.09e-07,3.05e-07,
2.80e-07,1.77e-07,1.10e-07,7.31e-08,4.89e-08,3.30e-08,2.25e-08,1.22e-08,5.45e-09,1.85e-09,5.72e-10,1.60e-10,8.07e-11};

//===============PPG094==================//

const int npts_ppg094 = 20;

double x_ppg094[npts_ppg094] = {1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 4.7, 4.9, 5.5, 6.5, 7.5, 8.5};
double y_ppg094[npts_ppg094] = {1.9e-05, 1.05e-05, 6.29e-06, 3.96e-06, 2.36e-06, 1.52e-06, 1.06e-06, 6.7e-07, 4.79e-07, 3.09e-07, 2.1e-07,1.54e-07, 1.1e-07, 7.5e-08, 5.8e-08, 5.1e-08,
1.73e-08, 5.58e-09, 1.99e-09, 1.36e-09};

double x_lohi_ppg094[npts_ppg094] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

double error_lo_ppg094[npts_ppg094] = {2.7e-07, 9.9e-07, 4.7e-07, 3e-07, 1.2e-07, 1.1e-07, 1.9e-07, 5e-08, 2.4e-08, 2.3e-08, 1.6e-08, 1.2e-08, 8e-09, 7.3e-09, 5.7e-09, 3.8e-09, 2.1e-09, 5.5e-10, 3.353e-10, 2.202e-10};
double error_hi_ppg094[npts_ppg094] = {2.1e-07, 8e-07, 6.9e-07, 4.3e-07, 3.3e-07, 1.3e-07, 9e-08, 7.2e-08, 5.2e-08, 3.4e-08, 2.3e-08, 1.2e-08, 9e-09, 8.2e-09, 6.3e-09, 4.1e-09, 1.9e-09, 6.1e-10, 2.802e-10, 2.46e-10};

//===============PPG141==================//

const int npts_ppg141 = 18;

double x_ppg141[npts_ppg141] = {0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625,4.875};
double y_ppg141[npts_ppg141] = {0.00212	,0.000793,0.000278,0.000109,4.77E-05,2.34E-05,1.15E-05,6.05E-06,3.28E-06,1.82E-06,1.08E-06,6.2E-07,	
4.07E-07,2.42E-07,1.59E-07,1.07E-07,8.02E-08,5.38E-08};

double x_lohi_ppg141[npts_ppg141] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03};

double stat_lo_ppg141[npts_ppg141] = {4E-05	,1E-05,3E-06,2E-06,1E-06,1E-06,4E-07,2E-07,2E-07,1E-07,1E-07,4E-08,3E-08,2E-08,2E-08,1E-08,
1.1E-08,7E-09};
double stat_hi_ppg141[npts_ppg141] = {4E-05	,1E-05,3E-06,2E-06,1E-06,1E-06,4E-07,2E-07,2E-07,1E-07,1E-07,4E-08,3E-08,2E-08,2E-08,1E-08,
1.1E-08,7E-09};

double sys_lo_ppg141[npts_ppg141] = {0.0005,0.00011,4E-05,1E-05,6E-06,3E-06,1E-06,7E-07,4E-07,2E-07,1E-07,7E-08,5E-08,3E-08,2E-08,1E-08,
9E-09,6E-09};
double sys_hi_ppg141[npts_ppg141] = {0.0005,0.00011,4E-05,1E-05,6E-06,3E-06,1E-06,7E-07,4E-07,2E-07,1E-07,7E-08,5E-08,3E-08,2E-08,1E-08,
9E-09,6E-09};

void heavy_flavor_pt_spectra_PHENIX(){

    TGraphAsymmErrors* ppg037_stat = new TGraphAsymmErrors(npts_ppg037, x_ppg037, y_ppg037, 0, 0, stat_lo_ppg037, stat_hi_ppg037);
    ppg037_stat->SetMarkerStyle(20);
    ppg037_stat->SetMarkerSize(0.9);
    ppg037_stat->SetMarkerColor(kRed);
    ppg037_stat->SetLineColor(kRed);
    TGraphAsymmErrors* ppg037_sys = new TGraphAsymmErrors(npts_ppg037, x_ppg037, y_ppg037, x_lohi_ppg037, x_lohi_ppg037, sys_lo_ppg037, sys_hi_ppg037);
    ppg037_sys->SetFillColor(kRed-3);
    ppg037_sys->SetFillStyle(0);
    ppg037_sys->SetLineColor(kRed);

    TGraphAsymmErrors* ppg065_stat = new TGraphAsymmErrors(npts_ppg065, x_ppg065, y_ppg065, 0, 0, stat_lo_ppg065, stat_hi_ppg065);
    ppg065_stat->SetMarkerStyle(20);
    ppg065_stat->SetMarkerSize(0.9);
    ppg065_stat->SetMarkerColor(kBlue);
    ppg065_stat->SetLineColor(kBlue);
    TGraphAsymmErrors* ppg065_sys = new TGraphAsymmErrors(npts_ppg065, x_ppg065, y_ppg065, x_lohi_ppg065, x_lohi_ppg065, sys_lo_ppg065, sys_hi_ppg065);
    ppg065_sys->SetFillColor(kBlue-3);
    ppg065_sys->SetFillStyle(0);
    ppg065_sys->SetLineColor(kBlue);

    TGraph* ppg094 = new TGraph(npts_ppg094, x_ppg094, y_ppg094);
    ppg094->SetMarkerStyle(20);
    ppg094->SetMarkerSize(0.9);
    ppg094->SetMarkerColor(kGreen-3);

    TGraphAsymmErrors* ppg094_err = new TGraphAsymmErrors(npts_ppg094, x_ppg094, y_ppg094, x_lohi_ppg094, x_lohi_ppg094, error_lo_ppg094, error_hi_ppg094);
    ppg094_err->SetFillColor(kGreen-3);
    ppg094_err->SetFillStyle(0);
    ppg094_err->SetLineColor(kGreen);

    TGraphAsymmErrors* ppg141_stat = new TGraphAsymmErrors(npts_ppg141, x_ppg141, y_ppg141, 0, 0, stat_lo_ppg141, stat_hi_ppg141);
    ppg141_stat->SetMarkerStyle(20);
    ppg141_stat->SetMarkerSize(0.9);
    ppg141_stat->SetMarkerColor(kBlack);
    ppg141_stat->SetLineColor(kBlack);
    TGraphAsymmErrors* ppg141_sys = new TGraphAsymmErrors(npts_ppg141, x_ppg141, y_ppg141, x_lohi_ppg141, x_lohi_ppg141, sys_lo_ppg141, sys_hi_ppg141);
    ppg141_sys->SetFillColor(kBlack-3);
    ppg141_sys->SetFillStyle(0);
    ppg141_sys->SetLineColor(kBlack);


    TFile* f1 = new TFile("pt_spectra_ccbar_detroit.root","READ");
    TH1D* pythia8_detroit_ccbar = (TH1D*)f1->Get("pt_spectra");
    
    TFile* f2 = new TFile("pt_spectra_bbbar_detroit.root","READ");
    TH1D* pythia8_detroit_bbbar = (TH1D*)f2->Get("pt_spectra");

    TH1D* pythia8_detroit_hf = (TH1D*)pythia8_detroit_ccbar->Clone();
    pythia8_detroit_hf->SetName("pythia8_detroit_hf");
    pythia8_detroit_hf->Add(pythia8_detroit_bbbar);
    pythia8_detroit_hf->Rebin(2);
    pythia8_detroit_hf->Scale(1./(2*5000.));
    pythia8_detroit_hf->SetMarkerStyle(20);
    pythia8_detroit_hf->SetMarkerColor(kCyan);

    TFile* f3 = new TFile("pt_spectra_ccbar_monash.root","READ");
    TH1D* pythia8_monash_ccbar = (TH1D*)f3->Get("pt_spectra");
    
    TFile* f4 = new TFile("pt_spectra_bbbar_monash.root","READ");
    TH1D* pythia8_monash_bbbar = (TH1D*)f4->Get("pt_spectra");

    TH1D* pythia8_monash_hf = (TH1D*)pythia8_monash_ccbar->Clone();
    pythia8_monash_hf->SetName("pythia8_monash_hf");
    pythia8_monash_hf->Add(pythia8_monash_bbbar);
    pythia8_monash_hf->Rebin(2);
    pythia8_monash_hf->Scale(1./(2*5000.));
    pythia8_monash_hf->SetMarkerStyle(20);
    pythia8_monash_hf->SetMarkerColor(kPink-4);

    TFile* f5 = new TFile("pt_spectra_ccbar_pythia6.root","READ");
    TH1D* pythia6_ccbar = (TH1D*)f5->Get("pt_spectra");
    
    TFile* f6 = new TFile("pt_spectra_bbbar_pythia6.root","READ");
    TH1D* pythia6_bbbar = (TH1D*)f6->Get("pt_spectra");

    TH1D* pythia6_hf = (TH1D*)pythia6_ccbar->Clone();
    pythia6_hf->SetName("pythia6_hf");
    pythia6_hf->Add(pythia6_bbbar);
    pythia6_hf->Rebin(2);
    pythia6_hf->Scale(1./(2*5000.));
    pythia6_hf->SetMarkerStyle(20);
    pythia6_hf->SetMarkerColor(kRed-4);


    TCanvas* c1 = new TCanvas();
    c1->cd();

    ppg065_stat->Draw("AP");
    ppg065_sys->Draw("e2same");
    ppg037_stat->Draw("p same");
    ppg037_sys->Draw("e2same");
    ppg094->Draw("p same");
    ppg094_err->Draw("e2same");
    ppg141_stat->Draw("p same");
    ppg141_sys->Draw("e2same");
    pythia8_detroit_hf->Draw("p same");
    pythia8_monash_hf->Draw("p same");
    pythia6_hf->Draw("p same");

    TMultiGraph* mg = new TMultiGraph();
    
    mg->Add(ppg065_stat,"AP");
    mg->Add(ppg065_sys, "e2same");
    mg->Add(ppg037_stat, "p same");
    mg->Add(ppg037_sys, "e2same");
    mg->Add(ppg094, "p same");
    mg->Add(ppg094_err, "e2same");
    mg->Add(ppg141_stat, "p same");
    mg->Add(ppg141_sys, "e2same");
    mg->SetMaximum(0.1);
    mg->SetMinimum(1e-10);
  

    TF1* func = new TF1("func", "[0]/( (exp(-[1]*x - [2]*x*x) + x/[3])^[4] )", 0.2, 10);
    func->SetParameters(1e-2, -1.42849, 1.17713e-1, 3.26575e-1, 3.25844);
    mg->Fit(func,"R+");


    TCanvas* c2 = new TCanvas();
    c2->cd();
    mg->Draw("a");
    pythia8_detroit_hf->Draw("p same");
    pythia8_monash_hf->Draw("p same");
    pythia6_hf->Draw("p same");
   
    TH1D* h1 = new TH1D("h1", "Ratio: Data vs Pythia8 Detroit", 50, 0, 10);
    TH1D* h2 = new TH1D("h2", "Ratio: Data vs Pythia8 Monash", 50, 0, 10);
    TH1D* h3 = new TH1D("h3", "Ratio: Data vs Pythia6", 50, 0, 10);

    int nbins = pythia8_detroit_hf->GetNbinsX();

    for(int ibin = 1; ibin < nbins+1; ibin++){

        if(pythia8_detroit_hf->GetBinContent(ibin) == 0) continue;
        double x = pythia8_detroit_hf->GetBinCenter(ibin);
        double y = func->Eval(x)/pythia8_detroit_hf->GetBinContent(ibin);

        h1->SetBinContent(h1->FindBin(x), y);

    }

    for(int ibin = 1; ibin < nbins+1; ibin++){

        if(pythia8_monash_hf->GetBinContent(ibin) == 0) continue;
        double x = pythia8_monash_hf->GetBinCenter(ibin);
        double y = func->Eval(x)/pythia8_monash_hf->GetBinContent(ibin);

        h2->SetBinContent(h2->FindBin(x), y);

    }

    for(int ibin = 1; ibin < nbins+1; ibin++){

        if(pythia6_hf->GetBinContent(ibin) == 0) continue;
        double x = pythia6_hf->GetBinCenter(ibin);
        double y = func->Eval(x)/pythia6_hf->GetBinContent(ibin);

        h3->SetBinContent(h3->FindBin(x), y);

    }

    TCanvas* c3 = new TCanvas();
    c3->cd();
    h1->SetLineColor(kRed);
    h1->Draw();
    h2->SetLineColor(kBlue);
    h2->Draw("same");
    h3->SetLineColor(kBlack);
    h3->Draw("same");

    TF1* func_detroit = new TF1("func_detroit","pol0",7.2,10.0);
    TF1* func_monash = new TF1("func_monash","pol0",7.2,10.0);
    TF1* func_pythia6 = new TF1("func_pythia6","pol0",7.2,10.0);

    h1->Fit(func_detroit,"R+");
    h2->Fit(func_monash, "R+");
    h3->Fit(func_pythia6, "R+");
    
    TH1D* pythia8_detroit_hf_scaled = (TH1D*)pythia8_detroit_hf->Clone();
    pythia8_detroit_hf_scaled->SetName("pythia8_detroit_hf_scaled");

    TH1D* pythia8_monash_hf_scaled = (TH1D*)pythia8_monash_hf->Clone();
    pythia8_monash_hf_scaled->SetName("pythia8_monash_hf_scaled");

    TH1D* pythia6_hf_scaled = (TH1D*)pythia6_hf->Clone();
    pythia6_hf_scaled->SetName("pythia6_hf_scaled");

    pythia8_detroit_hf_scaled->Scale(func_detroit->GetParameter(0));
    pythia8_monash_hf_scaled->Scale(func_monash->GetParameter(0));
    pythia6_hf_scaled->Scale(func_pythia6->GetParameter(0));

    TCanvas* c4 = new TCanvas();
    c4->cd();
    mg->Draw("a");
    pythia8_detroit_hf_scaled->Draw("p same");
    pythia8_monash_hf_scaled->Draw("p same");
    pythia6_hf_scaled->Draw("p same");

}
