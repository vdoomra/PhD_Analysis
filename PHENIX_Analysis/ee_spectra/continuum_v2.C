#include <TFile.h>
#include <TH1.h>

const int nbins = 11;
const double mass_range[nbins+1] = {1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 3.05};

void continuum_v2(){

    TFile* fdata_unlike = new TFile("ee_spectra_data_FG_unlike.root","READ");
    TH1D* h1d_unlike_data = (TH1D*)fdata_unlike->Get("h1d_ee_FG_cent0");
    h1d_unlike_data->SetName("h1d_unlike_data");

    TFile* fdata_like = new TFile("ee_spectra_data_FG.root","READ");
    TH1D* h1d_like_data = (TH1D*)fdata_like->Get("h1d_ee_FG_like_cent0");
    h1d_like_data->SetName("h1d_like_data");


    TFile* fdata_nocuts = new TFile("ee_spectra_data_FG_nocuts.root","READ");
    TH1D* h1d_like_data_nocuts = (TH1D*)fdata_nocuts->Get("h1d_ee_FG_like_cent0");
    h1d_like_data_nocuts->SetName("h1d_like_data_nocuts");

    double scale1 = h1d_like_data->Integral(h1d_like_data->FindBin(1.0), h1d_like_data->FindBin(3.0))/h1d_like_data_nocuts->Integral(h1d_like_data_nocuts->FindBin(1.0), h1d_like_data_nocuts->FindBin(3.0));

    h1d_like_data_nocuts->Scale(scale1);

    TFile* fin_BG = new TFile("ee_spectra_data_BG_nocuts.root","READ");
    TH1D* h1d_unlike_BG = (TH1D*)fin_BG->Get("h1d_ee_BG_cent0");
    TH1D* h1d_like_BG = (TH1D*)fin_BG->Get("h1d_ee_BG_like_cent0");

    TH1D* h1d_correction = (TH1D*)h1d_unlike_BG->Clone();
    h1d_correction->SetName("h1d_correction");
    h1d_correction->Divide(h1d_like_BG);

    h1d_like_data_nocuts->Multiply(h1d_correction);
    h1d_like_data->Multiply(h1d_correction);

    TFile* fjpsi = new TFile("ee_spectra_jpsi_sim_correct_inflection_point.root","READ");
    TH1D* h1d_unlike_jpsi = (TH1D*)fjpsi->Get("h1d_ee_FG_cent0");
    h1d_unlike_jpsi->SetName("h1d_unlike_jpsi");

    double scale = h1d_unlike_data->GetBinContent(h1d_unlike_data->FindBin(3.10))/h1d_unlike_jpsi->GetBinContent(h1d_unlike_jpsi->FindBin(3.10));

    h1d_unlike_jpsi->Scale(scale);

    TCanvas* c1 = new TCanvas();
    c1->cd();
    h1d_unlike_data->SetMarkerColor(kRed);
    h1d_unlike_data->Draw();
    h1d_like_data_nocuts->SetMarkerColor(kGreen-7);
    h1d_like_data_nocuts->Draw("same");
    h1d_like_data->SetMarkerColor(kRed);
    h1d_like_data->Draw("same");
    h1d_unlike_jpsi->SetMarkerColor(kBlue);
    h1d_unlike_jpsi->Draw("same");

    TH1D* h1d_signal_temp = (TH1D*)h1d_unlike_data->Clone();
    h1d_signal_temp->Reset("ICESM");
    h1d_signal_temp->SetName("h1d_signal_temp");

    int nbinsX = h1d_signal_temp->GetNbinsX();

    for(int ibin=1; ibin < nbinsX+1; ibin++){

        double bin_content = h1d_unlike_data->GetBinContent(ibin) - h1d_unlike_jpsi->GetBinContent(ibin) - h1d_like_data_nocuts->GetBinContent(ibin);
        double bin_error = sqrt( pow(h1d_unlike_data->GetBinError(ibin), 2) + pow(h1d_unlike_jpsi->GetBinError(ibin), 2)  + pow(h1d_like_data_nocuts->GetBinError(ibin), 2));

        h1d_signal_temp->SetBinContent(ibin, bin_content);
        h1d_signal_temp->SetBinError(ibin, bin_error);


    }


    TH1D* h1d_signal = new TH1D("h1d_signal","h1d_signal", nbins, mass_range);
    h1d_signal->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_signal_temp->FindBin(mass_range[ibin]);
        int bin_hi = h1d_signal_temp->FindBin(mass_range[ibin+1]);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i < bin_hi; i++){

            bin_content += h1d_signal_temp->GetBinContent(i)*h1d_signal_temp->GetBinWidth(i);
            bin_error += pow(h1d_signal_temp->GetBinError(i)*h1d_signal_temp->GetBinWidth(i), 2);

        }

        h1d_signal->SetBinContent(ibin+1, bin_content/h1d_signal->GetBinWidth(ibin+1));
        h1d_signal->SetBinError(ibin+1, sqrt(bin_error)/h1d_signal->GetBinWidth(ibin+1));


    }


    TFile* fin_ccbar_sim = new TFile("ee_spectra_ccbar_sim_correct_inflection_point.root","READ");
    TH1D* h1d_ccbar_sim_temp = (TH1D*)fin_ccbar_sim->Get("h1d_ee_FG_cent0");
    h1d_ccbar_sim_temp->SetName("h1d_ccbar_sim_temp");

    TH1D* h1d_ccbar_sim = new TH1D("h1d_ccbar_sim","h1d_ccbar_sim", nbins, mass_range);
    h1d_ccbar_sim->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_ccbar_sim_temp->FindBin(mass_range[ibin]);
        int bin_hi = h1d_ccbar_sim_temp->FindBin(mass_range[ibin+1]);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i < bin_hi; i++){

            bin_content += h1d_ccbar_sim_temp->GetBinContent(i)*h1d_ccbar_sim_temp->GetBinWidth(i);
            bin_error += pow(h1d_ccbar_sim_temp->GetBinError(i)*h1d_ccbar_sim_temp->GetBinWidth(i), 2);

        }

        h1d_ccbar_sim->SetBinContent(ibin+1, bin_content/h1d_ccbar_sim->GetBinWidth(ibin+1));
        h1d_ccbar_sim->SetBinError(ibin+1, sqrt(bin_error)/h1d_ccbar_sim->GetBinWidth(ibin+1));


    }

    double scale_ccbar = h1d_signal->Integral(h1d_signal->FindBin(1.4), h1d_signal->FindBin(2.5))/h1d_ccbar_sim->Integral(h1d_ccbar_sim->FindBin(1.4), h1d_ccbar_sim->FindBin(2.5));

    h1d_ccbar_sim->Scale(scale_ccbar);

    TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetLogy();
    h1d_signal->Draw();
    h1d_ccbar_sim->Draw("same");
    

}