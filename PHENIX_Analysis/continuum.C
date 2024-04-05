#include <TFile.h>
#include <TH1.h>

const int nbins = 11;
const double mass_range[nbins+1] = {1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 3.05};

void continuum(){

    TH1D* h1d_like_data = new TH1D("h1d_like_data","Like Sign Distribution", 200, 1.10, 3.30);

    // Create the Like Sign Mass Histogram

    std::ifstream inFile("mlp_likesign_mass.txt");
    double x, y, err_y;
    while(!inFile.eof()){

        inFile >> x >> y >> err_y;

        int binX = h1d_like_data->FindBin(x);
        h1d_like_data->SetBinContent(binX, y);
        h1d_like_data->SetBinError(binX, err_y);

    }

    TFile* fdata = new TFile("unlike_sign_data_correct_energy_cut.root","READ");
    TH1D* h1d_unlike_data = (TH1D*)fdata->Get("h1d_ee_FG_cent0");

    TFile* fjpsi = new TFile("DCA_jpsi_correct_energy_cut.root","READ");
    TH1D* h1d_unlike_jpsi = (TH1D*)fjpsi->Get("h1d_ee_FG_cent0");
    h1d_unlike_jpsi->SetName("h1d_unlike_jpsi");

    double scale = h1d_unlike_data->GetBinContent(h1d_unlike_data->FindBin(3.10))/h1d_unlike_jpsi->GetBinContent(h1d_unlike_jpsi->FindBin(3.10));

    h1d_unlike_jpsi->Scale(scale);

    TH1D* h1d_unlike_data_minus_jpsi_temp = (TH1D*)h1d_unlike_data->Clone();
    h1d_unlike_data_minus_jpsi_temp->Reset("ICESM");
    h1d_unlike_data_minus_jpsi_temp->SetName("h1d_unlike_data_minus_jpsi_temp");

    int nbinsX = h1d_unlike_data_minus_jpsi_temp->GetNbinsX();

    for(int ibin=1; ibin < nbinsX+1; ibin++){

        double bin_content = h1d_unlike_data->GetBinContent(ibin) - h1d_unlike_jpsi->GetBinContent(ibin);
        double bin_error = sqrt( pow(h1d_unlike_data->GetBinError(ibin), 2) + pow(h1d_unlike_jpsi->GetBinError(ibin), 2) );

        h1d_unlike_data_minus_jpsi_temp->SetBinContent(ibin, bin_content);
        h1d_unlike_data_minus_jpsi_temp->SetBinError(ibin, bin_error);


    }

    TCanvas* c1 = new TCanvas();
    c1->cd();

    h1d_unlike_data_minus_jpsi_temp->SetMarkerColor(kBlue);
    h1d_unlike_data_minus_jpsi_temp->Draw("same");

    h1d_like_data->SetMarkerColor(kBlack);
    h1d_like_data->Draw("same");


    TH1D* h1d_unlike_data_minus_jpsi = new TH1D("h1d_unlike_data_minus_jpsi","h1d_unlike_data_minus_jpsi", nbins, mass_range);
    h1d_unlike_data_minus_jpsi->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_unlike_data_minus_jpsi_temp->FindBin(mass_range[ibin]);
        int bin_hi = h1d_unlike_data_minus_jpsi_temp->FindBin(mass_range[ibin+1]);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i < bin_hi; i++){

            bin_content += h1d_unlike_data_minus_jpsi_temp->GetBinContent(i)*h1d_unlike_data_minus_jpsi_temp->GetBinWidth(i);
            bin_error += pow(h1d_unlike_data_minus_jpsi_temp->GetBinError(i)*h1d_unlike_data_minus_jpsi_temp->GetBinWidth(i), 2);

        }

        h1d_unlike_data_minus_jpsi->SetBinContent(ibin+1, bin_content/h1d_unlike_data_minus_jpsi->GetBinWidth(ibin+1));
        h1d_unlike_data_minus_jpsi->SetBinError(ibin+1, sqrt(bin_error)/h1d_unlike_data_minus_jpsi->GetBinWidth(ibin+1));


    }

    

    TH1D* h1d_signal = new TH1D("h1d_signal","ccbar continuum distribution", nbins, mass_range);
    h1d_signal->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        double bin_content = h1d_unlike_data_minus_jpsi->GetBinContent(h1d_unlike_data_minus_jpsi->FindBin(mass_range[ibin])) - h1d_like_data->GetBinContent(h1d_like_data->FindBin(mass_range[ibin]));
        double bin_error = sqrt( pow(h1d_unlike_data_minus_jpsi->GetBinError(h1d_unlike_data_minus_jpsi->FindBin(mass_range[ibin])), 2) + pow(h1d_like_data->GetBinError(h1d_like_data->FindBin(mass_range[ibin])), 2) );

        h1d_signal->SetBinContent(ibin+1, bin_content);
        h1d_signal->SetBinError(ibin+1, bin_error);


    }

    TFile* fin_ccbar_sim_pthat0 = new TFile("DCA_ccbar_correct_energy_cut_pthat0.root","READ");
    TH1D* h1d_ccbar_sim_pthat0_temp = (TH1D*)fin_ccbar_sim_pthat0->Get("h1d_ee_FG_cent0");
    h1d_ccbar_sim_pthat0_temp->SetName("h1d_ccbar_sim_pthat0_temp");

    TH1D* h1d_ccbar_sim_pthat0 = new TH1D("h1d_ccbar_sim_pthat0","h1d_ccbar_sim_pthat0", nbins, mass_range);
    h1d_ccbar_sim_pthat0->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_ccbar_sim_pthat0_temp->FindBin(mass_range[ibin]);
        int bin_hi = h1d_ccbar_sim_pthat0_temp->FindBin(mass_range[ibin+1]);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i < bin_hi; i++){

            bin_content += h1d_ccbar_sim_pthat0_temp->GetBinContent(i)*h1d_ccbar_sim_pthat0_temp->GetBinWidth(i);
            bin_error += pow(h1d_ccbar_sim_pthat0_temp->GetBinError(i)*h1d_ccbar_sim_pthat0_temp->GetBinWidth(i), 2);

        }

        h1d_ccbar_sim_pthat0->SetBinContent(ibin+1, bin_content/h1d_ccbar_sim_pthat0->GetBinWidth(ibin+1));
        h1d_ccbar_sim_pthat0->SetBinError(ibin+1, sqrt(bin_error)/h1d_ccbar_sim_pthat0->GetBinWidth(ibin+1));


    }

    double scale_ccbar = h1d_signal->Integral(h1d_signal->FindBin(1.1), h1d_signal->FindBin(3.0))/h1d_ccbar_sim_pthat0->Integral(h1d_ccbar_sim_pthat0->FindBin(1.1), h1d_ccbar_sim_pthat0->FindBin(3.0));

    h1d_ccbar_sim_pthat0->Scale(scale_ccbar);

    TFile* fin_ccbar_sim_pthat2 = new TFile("DCA_ccbar_correct_energy_cut_pthat2.root","READ");
    TH1D* h1d_ccbar_sim_pthat2_temp = (TH1D*)fin_ccbar_sim_pthat2->Get("h1d_ee_FG_cent0");

    TH1D* h1d_ccbar_sim_pthat2 = new TH1D("h1d_ccbar_sim_pthat2","h1d_ccbar_sim_pthat2", nbins, mass_range);
    h1d_ccbar_sim_pthat2->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_ccbar_sim_pthat2_temp->FindBin(mass_range[ibin]);
        int bin_hi = h1d_ccbar_sim_pthat2_temp->FindBin(mass_range[ibin+1]);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i < bin_hi; i++){

            bin_content += h1d_ccbar_sim_pthat2_temp->GetBinContent(i)*h1d_ccbar_sim_pthat2_temp->GetBinWidth(i);
            bin_error += pow(h1d_ccbar_sim_pthat2_temp->GetBinError(i)*h1d_ccbar_sim_pthat2_temp->GetBinWidth(i), 2);

        }

        h1d_ccbar_sim_pthat2->SetBinContent(ibin+1, bin_content/h1d_ccbar_sim_pthat2->GetBinWidth(ibin+1));
        h1d_ccbar_sim_pthat2->SetBinError(ibin+1, sqrt(bin_error)/h1d_ccbar_sim_pthat2->GetBinWidth(ibin+1));


    }

    double scale_ccbar_pthat2 = h1d_signal->Integral(h1d_signal->FindBin(1.1), h1d_signal->FindBin(3.0))/h1d_ccbar_sim_pthat2->Integral(h1d_ccbar_sim_pthat2->FindBin(1.1), h1d_ccbar_sim_pthat2->FindBin(3.0));

    h1d_ccbar_sim_pthat2->Scale(scale_ccbar_pthat2);

    TFile* fin_ccbar_sim_pythia6 = new TFile("DCA_ccbar_pythia6.root","READ");
    TH1D* h1d_ccbar_sim_pythia6_temp = (TH1D*)fin_ccbar_sim_pythia6->Get("h1d_ee_FG_cent0");

    TH1D* h1d_ccbar_sim_pythia6 = new TH1D("h1d_ccbar_sim_pythia6","h1d_ccbar_sim_pythia6", nbins, mass_range);
    h1d_ccbar_sim_pythia6->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_ccbar_sim_pythia6_temp->FindBin(mass_range[ibin]);
        int bin_hi = h1d_ccbar_sim_pythia6_temp->FindBin(mass_range[ibin+1]);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i < bin_hi; i++){

            bin_content += h1d_ccbar_sim_pythia6_temp->GetBinContent(i)*h1d_ccbar_sim_pythia6_temp->GetBinWidth(i);
            bin_error += pow(h1d_ccbar_sim_pythia6_temp->GetBinError(i)*h1d_ccbar_sim_pythia6_temp->GetBinWidth(i), 2);

        }

        h1d_ccbar_sim_pythia6->SetBinContent(ibin+1, bin_content/h1d_ccbar_sim_pythia6->GetBinWidth(ibin+1));
        h1d_ccbar_sim_pythia6->SetBinError(ibin+1, sqrt(bin_error)/h1d_ccbar_sim_pythia6->GetBinWidth(ibin+1));


    }

    double scale_ccbar_pythia6 = h1d_signal->Integral(h1d_signal->FindBin(1.1), h1d_signal->FindBin(3.0))/h1d_ccbar_sim_pythia6->Integral(h1d_ccbar_sim_pythia6->FindBin(1.1), h1d_ccbar_sim_pythia6->FindBin(3.0));

    h1d_ccbar_sim_pythia6->Scale(scale_ccbar_pythia6);


    TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetLogy();
    h1d_signal->Draw();
    h1d_ccbar_sim_pthat0->Draw("same");
    h1d_ccbar_sim_pthat2->Draw("same");
    h1d_ccbar_sim_pythia6->Draw("same");

    TFile* fout = new TFile("continuum.root","RECREATE");
    fout->cd();
    h1d_unlike_data->Write();
    h1d_like_data->Write();
    h1d_unlike_jpsi->Write();
    h1d_signal->Write();
    h1d_ccbar_sim_pthat0->Write();
    h1d_ccbar_sim_pthat2->Write();
    h1d_ccbar_sim_pythia6->Write();
    fout->Close();
}