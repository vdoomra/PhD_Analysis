#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>

const int nbins = 3;
const double mass_lo[nbins] = {1.15, 2.05, 2.95};
const double mass_hi[nbins] = {2.05, 2.95, 3.4};
TH1D* h1d_dca_dist[nbins] = {NULL};

void plot(){

    TFile* fin = new TFile("DCA_ccbar_no_rejection_no_gencut_eff.root", "READ");
    TH3D* h3d_spectra = (TH3D*)fin->Get("hist_pair_dca_using_quadrature_unlike");
    TH2D* h2d_spectra = (TH2D*)h3d_spectra->Project3D("yx"); // Mass Along X
    h2d_spectra->SetName("h2d_spectra");

    TCanvas* c1 = new TCanvas();
    c1->cd();
    h2d_spectra->Draw("colz");

    h2d_spectra->GetYaxis()->SetRangeUser(0, 10);
    TH1D* h1d_temp_spectra = (TH1D*)h2d_spectra->ProjectionX();
    h1d_temp_spectra->SetName("h1d_temp_spectra");

    TH1D* h1d_spectra = (TH1D*)h1d_temp_spectra->Clone();
    h1d_spectra->Reset("ICESM");
    h1d_spectra->SetName("h1d_spectra");

    int nbinsX = h1d_temp_spectra->GetNbinsX();

    for(int ibin=1; ibin < nbinsX+1; ibin++){

        double binwidth = h1d_temp_spectra->GetBinWidth(ibin);
        double content = h1d_temp_spectra->GetBinContent(ibin);
        double err = h1d_temp_spectra->GetBinError(ibin);
        h1d_spectra->SetBinContent(ibin, content/binwidth);
        h1d_spectra->SetBinError(ibin, err/binwidth);

    }
    TCanvas* c2 = new TCanvas();
    c2->cd();
    h1d_spectra->Draw();

    TCanvas* c3 = new TCanvas();
    c3->cd();

    TH2D* h2d_dca_mass = (TH2D*)h3d_spectra->Project3D("zx"); // Mass Along X
    h2d_dca_mass->Draw("colz");

    TCanvas* c4 = new TCanvas();
    c4->SetLogy();
    c4->cd();

    for(int ibin = 0; ibin < nbins; ibin++){

        h2d_dca_mass->GetXaxis()->SetRangeUser(mass_lo[ibin], mass_hi[ibin]);
        h1d_dca_dist[ibin] = (TH1D*)h2d_dca_mass->ProjectionY();
        h1d_dca_dist[ibin]->SetName(Form("h1d_dac_dist_%d",ibin));

        if(ibin==0) h1d_dca_dist[ibin]->Draw();
        else h1d_dca_dist[ibin]->Draw("same");
    }

    TFile* fout = new TFile("DCA_dist_ccbar_eff.root","RECREATE");
    fout->cd();
    h2d_spectra->Write();
    h1d_spectra->Write();
    for(int ibin = 0; ibin < nbins; ibin++){ 
        h1d_dca_dist[ibin]->Rebin(2); 
        h1d_dca_dist[ibin]->Scale(1./h1d_dca_dist[ibin]->Integral()); 
        h1d_dca_dist[ibin]->Write(); 
    }
    fout->Close();

    for(int i = 0; i < nbins; i++){

        std::ofstream datafile;
        datafile.open(Form("outFile_%f_%f.txt",mass_lo[i], mass_hi[i]),std::ofstream::app);
        int nbinsX = h1d_dca_dist[i]->GetNbinsX();
        for(int j = 1; j < nbinsX+1; j++){

            if(h1d_dca_dist[i]->GetBinCenter(j) < 0) continue;
            datafile << h1d_dca_dist[i]->GetBinCenter(j) << "\t" << h1d_dca_dist[i]->GetBinContent(j) << "\t" << h1d_dca_dist[i]->GetBinError(j) << endl;

        }

    }


}