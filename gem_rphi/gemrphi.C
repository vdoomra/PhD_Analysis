const int nFiles = 1;
const int nsectors = 7;
const char* filename[nFiles] = { "/volatile/halla/moller12gev/vdoomra/slim_new_sieve_beam_optics1DS_6.6/slim_new_sieve_beam_optics1DS_p3.root"};
const char* passNames[nFiles] = {"pass3"};
const char* sector[nsectors] = {"sector7", "sector6", "sector5", "sector4", "sector3" ,"sector2", "sector1"};
double phi_lo[nsectors] = {0, 1, 2, 2.8, 3.6, 4.5, 5.3};
double phi_hi[nsectors] = {1, 2, 2.8, 3.6, 4.5, 5.3, 6.28};

TH2D* rphi[nFiles] = {NULL};
TH2D* rrp[nFiles] = {NULL};
TH1D* radial[nFiles][nsectors] = {NULL};
TH2D* rrprime[nFiles][nsectors] = {NULL};

const double pi = acos(-1);

void gemrphi(){

  for(int i=0; i<nFiles; i++){ 

    TFile* fin = new TFile(Form("%s",filename[i]),"READ");

    TTree* T = (TTree*)fin->Get("newT");

    int nEntries = T->GetEntries();
 
  
    double gem1_r, gem1_phi, sieve_r, gem1_x, gem1_y, gem1_px, gem1_py, gem1_pz;
    double rate;
 
    T->SetBranchAddress("gem1_r",&gem1_r);
    T->SetBranchAddress("gem1_ph",&gem1_phi);
    T->SetBranchAddress("rate",&rate);
    T->SetBranchAddress("sieve_r",&sieve_r);
    T->SetBranchAddress("gem1_px",&gem1_px);
    T->SetBranchAddress("gem1_py",&gem1_py);
    T->SetBranchAddress("gem1_pz",&gem1_pz);
    T->SetBranchAddress("gem1_x",&gem1_x);
    T->SetBranchAddress("gem1_y",&gem1_y);
   
    rphi[i] = new TH2D(Form("gem_r_phi_%s",passNames[i]),"gem_r_phi_distribution;r(mm);phi(rad)",1000, 500, 1300, 1000, 0, 7);
    rrp[i] = new TH2D(Form("gem_r_rprime_%s",passNames[i]),"gem_r_rprime_distribution;r(mm);r'",1000, 500, 1300, 1000, 0, 0.1);

    for(int j=0; j<nsectors; j++){

      radial[i][j] = new TH1D(Form("r_%s_%s",passNames[i], sector[j]),"gem_r_distribution;r[mm];Counts",1000,500,1300); 
      rrprime[i][j] = new TH2D(Form("gem_r_rprime_%s_%s",passNames[i], sector[j]),"gem_r_r'_distribution;r[mm];r'",1000,500,1300, 1000, 0, 0.1);

  }
   

    for(int k=0; k<nEntries; k++){  

      T->GetEntry(k);  int index_sector = -999;

      rate = 1; //Beam Gnerator

      if(sieve_r < 26.5) continue;
      if(gem1_phi < 0) gem1_phi += 2*pi;
      double rprime = (gem1_x*gem1_px + gem1_y*gem1_py)/(gem1_r*gem1_pz);

      rphi[i]->Fill(gem1_r, gem1_phi, rate); rrp[i]->Fill(gem1_r, rprime, rate);
    
      for(int l=0; l<nsectors; l++){  if(gem1_phi > phi_lo[l] && gem1_phi < phi_hi[l]) index_sector = l; }

      radial[i][index_sector]->Fill(gem1_r, rate); 
      rrprime[i][index_sector]->Fill(gem1_r, rprime, rate);

     
    }

  }

  TFile* fout = new TFile("gemrphi_new_sieve_beam_optics1DS_p3.root","RECREATE");
  fout->cd();
  
  for(int i=0; i<nFiles; i++){  

    rphi[i]->Write(); rrp[i]->Write();

   for(int j=0; j<nsectors; j++){

     radial[i][j]->Write();  rrprime[i][j]->Write();

   }
 }

  fout->Close();

}
