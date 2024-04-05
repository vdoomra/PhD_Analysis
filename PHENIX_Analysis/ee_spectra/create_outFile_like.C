#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>

const int nbins = 20;
const double mass_range[nbins+1] = {0.60, 0.75, 0.87, 0.95, 1.05, 1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.10, 3.25, 3.35, 3.50};

void create_outFile_like(){

    TCanvas* c1 = new TCanvas();
    c1->SetLogy();

    TFile* fin = new TFile("ee_spectra_pp_B0_and_B1orB2orB3_no_highpt_cutoff_correct_inflection_point.root", "READ");
    TH1D* h1d_spectra_like_temp = (TH1D*)fin->Get("h1d_ee_FG_like_cent0"); // Mass Along X

    TFile* fin_BG = new TFile("ee_spectra_pp_BG_B0_and_B1orB2orB3_no_highpt_cutoff_correct_inflection_point.root","READ");
    TH1D* h1d_unlike_BG = (TH1D*)fin_BG->Get("h1d_ee_BG_cent0");
    TH1D* h1d_like_BG = (TH1D*)fin_BG->Get("h1d_ee_BG_like_cent0");

    TH1D* h1d_correction = (TH1D*)h1d_unlike_BG->Clone();
    h1d_correction->SetName("h1d_correction");
    h1d_correction->Divide(h1d_like_BG);

    h1d_spectra_like_temp->Multiply(h1d_correction);

    TH1D* h1d_spectra_like = new TH1D("h1d_spectra_like","Like Sign Distribution",nbins, mass_range);
    h1d_spectra_like->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_spectra_like_temp->FindBin(mass_range[ibin]);
        int bin_hi = h1d_spectra_like_temp->FindBin(mass_range[ibin+1]);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i < bin_hi; i++){

            bin_content += h1d_spectra_like_temp->GetBinContent(i)*h1d_spectra_like_temp->GetBinWidth(i);
            bin_error += pow(h1d_spectra_like_temp->GetBinError(i)*h1d_spectra_like_temp->GetBinWidth(i), 2);

        }

        h1d_spectra_like->SetBinContent(ibin+1, bin_content/h1d_spectra_like->GetBinWidth(ibin+1));
        h1d_spectra_like->SetBinError(ibin+1, sqrt(bin_error)/h1d_spectra_like->GetBinWidth(ibin+1));


    }

    h1d_spectra_like->Draw();

    std::ofstream datafile;
    datafile.open("outFile_likesign.txt",std::ofstream::app);
    int nbinsX = h1d_spectra_like->GetNbinsX();
    for(int j = 1; j < nbinsX+1; j++){

        if(h1d_spectra_like->GetBinContent(j) == 0) continue;
        datafile << h1d_spectra_like->GetBinCenter(j) << "\t" << h1d_spectra_like->GetBinContent(j) << "\t" << h1d_spectra_like->GetBinError(j) << endl;

    }

}

1	238.186	168.423