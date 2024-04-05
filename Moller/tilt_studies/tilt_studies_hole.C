#include <TTree.h>
#include <TFile.h>

const int nsectors = 7;
const float pi = TMath::ACos(-1);
const int nholes = 3;

const char* hole_name[nholes] = {"hole1", "hole2", "hole3"};
const char* sector[nsectors] = {"sector1", "sector2", "sector3", "sector4", "sector5" ,"sector6", "sector7"};
double phi_lo[nsectors] = {0, 2*pi/7, 4*pi/7, 6*pi/7, 8*pi/7, 10*pi/7, 12*pi/7};
double phi_hi[nsectors] = {2*pi/7, 4*pi/7, 6*pi/7, 8*pi/7, 10*pi/7, 12*pi/7, 2*pi};

TH1D* deltaE[nsectors][nholes] = {NULL};
TH1D* delta_theta[nsectors][nholes] = {NULL};
TH2D* gemrphi_eloss[nsectors][nholes] = {NULL};
TH2D* gemrphi[nsectors][nholes] = {NULL};
TH2D* gemrphi_thetaloss[nsectors][nholes] = {NULL};

void tilt_studies_hole(const char* infile, const char* outfile){

    for(int i=6; i<nsectors; i++){ 
        for(int j=0; j<nholes; j++){
            deltaE[i][j] = new TH1D(Form("deltaE_%s_%s", sector[i], hole_name[j]), "Delta E Distribution", 70000, 0, 7000); 
            delta_theta[i][j] = new TH1D(Form("delta_theta_%s_%s", sector[i],hole_name[j]), "Delta theta Distribution", 2000, -100, 100);
            gemrphi[i][j] = new TH2D(Form("gemrphi_%s_%s", sector[i], hole_name[j]), "GEM r-phi Distribution", 1000, 600, 1100, 1000, 0, 7);
            gemrphi_eloss[i][j] = new TH2D(Form("gemrphi_eloss_%s_%s", sector[i], hole_name[j]), "GEM r-phi Distribution", 1000, 600, 1100, 1000, 0, 7);
            gemrphi_thetaloss[i][j] = new TH2D(Form("gemrphi_thetaloss_%s_%s", sector[i], hole_name[j]), "GEM r-phi Distribution", 1000, 600, 1100, 1000, 0, 7);
        }
    } 

    TFile* fin = new TFile(infile, "READ");
    TTree* T = (TTree*)fin->Get("newT");

    int nevts = T->GetEntries();

    double front_sieve_k, behind_sieve_k, sieve_r, rate, sieve_ph, gem_r, gem_ph, front_sieve_pz, behind_sieve_pz, front_sieve_px, behind_sieve_px, front_sieve_py, behind_sieve_py;

    T->SetBranchAddress("front_sieve_k", &front_sieve_k);
    T->SetBranchAddress("front_sieve_px", &front_sieve_px);
    T->SetBranchAddress("front_sieve_py", &front_sieve_py);
    T->SetBranchAddress("front_sieve_pz", &front_sieve_pz);
    T->SetBranchAddress("behind_sieve_k", &behind_sieve_k);
    T->SetBranchAddress("behind_sieve_px", &behind_sieve_px);
    T->SetBranchAddress("behind_sieve_py", &behind_sieve_py);
    T->SetBranchAddress("behind_sieve_pz", &behind_sieve_pz);
    T->SetBranchAddress("sieve_r", &sieve_r);
    T->SetBranchAddress("sieve_ph", &sieve_ph);
    T->SetBranchAddress("rate", &rate);
    T->SetBranchAddress("gem1_r", &gem_r);
    T->SetBranchAddress("gem1_ph", &gem_ph);

    for(int ievt = 0; ievt < nevts; ievt++){

        T->GetEntry(ievt);
        if(sieve_r < 26.5) continue;
        if(sieve_ph < 0) sieve_ph += 2*pi;
        if(gem_ph < 0) gem_ph += 2*pi;
        rate = rate/200;

        double front_theta = TMath::ACos(front_sieve_pz/sqrt(front_sieve_px*front_sieve_px + front_sieve_py*front_sieve_py + front_sieve_pz*front_sieve_pz));
        double behind_theta = TMath::ACos(behind_sieve_pz/sqrt(behind_sieve_px*behind_sieve_px + behind_sieve_py*behind_sieve_py + behind_sieve_pz*behind_sieve_pz));

        for(int j=6; j<nsectors; j++){ 

            if(sieve_ph > phi_lo[j] && sieve_ph < phi_hi[j]){

                // For sector 7

                if(sieve_r < 55){

                    delta_theta[j][0]->Fill((behind_theta-front_theta)*1000, rate);
                    gemrphi[j][0]->Fill(gem_r, gem_ph, rate); 
                    deltaE[j][0]->Fill((front_sieve_k - behind_sieve_k), rate);
                    if(fabs(front_sieve_k-behind_sieve_k) > 0.1) gemrphi_eloss[j][0]->Fill(gem_r, gem_ph, rate);
                    if(fabs(front_theta-behind_theta)*1000 > 0.1) gemrphi_thetaloss[j][0]->Fill(gem_r, gem_ph, rate);   

                }

                else if(sieve_r > 55 && sieve_r < 75){

                    delta_theta[j][1]->Fill((behind_theta-front_theta)*1000, rate);
                    gemrphi[j][1]->Fill(gem_r, gem_ph, rate); 
                    deltaE[j][1]->Fill((front_sieve_k - behind_sieve_k), rate);
                    if(fabs(front_sieve_k-behind_sieve_k) > 0.1) gemrphi_eloss[j][1]->Fill(gem_r, gem_ph, rate);
                    if(fabs(front_theta-behind_theta)*1000 > 0.1) gemrphi_thetaloss[j][1]->Fill(gem_r, gem_ph, rate);

                }

                else if(sieve_r > 75){

                    delta_theta[j][2]->Fill((behind_theta-front_theta)*1000, rate);
                    gemrphi[j][2]->Fill(gem_r, gem_ph, rate); 
                    deltaE[j][2]->Fill((front_sieve_k - behind_sieve_k), rate);
                    if(fabs(front_sieve_k-behind_sieve_k) > 0.1) gemrphi_eloss[j][2]->Fill(gem_r, gem_ph, rate);
                    if(fabs(front_theta-behind_theta)*1000 > 0.1) gemrphi_thetaloss[j][2]->Fill(gem_r, gem_ph, rate); 

                } 

                

            }


        }

    }

    TFile* fout = new TFile(outfile,"RECREATE");
    fout->cd();
    for(int k=6; k<nsectors; k++){ 
       for(int j=0; j<nholes; j++){
            deltaE[k][j]->Write(); 
            gemrphi[k][j]->Write();
            gemrphi_eloss[k][j]->Write();
            gemrphi_thetaloss[k][j]->Write();
            delta_theta[k][j]->Write();

        }
    }

    fout->Close();

}