const int nsides = 2;
const int nsectors = 2;
const char* side[nsides] = {"left", "right"};
double phi_lo[nsectors] = { 0.8, 3.6};
double phi_hi[nsectors] = { 1.8, 4.4};
double gem_rcut[nsectors] = { 950, 950};
const char* sector[nsectors] = {"sector6", "sector3"};

TH2D* th_p[nsectors][nsides];
TH1D* rsieve[nsectors][nsides];
TH1D* zvertex[nsectors][nsides];

const double pi = acos(-1);

void gemrphi_Ana(){

  TFile* fin = new TFile("slim_new_sieve_beam_p1.root","READ");
  TTree* T = (TTree*)fin->Get("newT");

  for(int i=0; i<nsectors; i++){

    for(int j=0; j<nsides; j++){
 
      rsieve[i][j] = new TH1D(Form("rsieve_%s_%s",sector[i], side[j]),"Radial location at the sieve plane", 400, 0, 200);
      zvertex[i][j] = new TH1D(Form("vz_%s_%s",sector[i], side[j]),"Scaterred Electron z vertex distribution", 1000, -10000, 0);
      th_p[i][j] = new TH2D(Form("th_p_%s_%s",sector[i], side[j]),"Scaterred Electron energy vs #theta distribution", 100,0,1.8,1000, 0, 3);

  }
 }

  int nEntries = T->GetEntries();

  cout << "Number of Events in the tree = " << nEntries << endl;

  double sieve_px, sieve_py, sieve_pz, sieve_r, rate, gem1_r, gem1_ph, gem1_vx, gem1_vy, gem1_vz;

  T->SetBranchAddress("sieve_px",&sieve_px);
  T->SetBranchAddress("sieve_py",&sieve_py);
  T->SetBranchAddress("sieve_pz",&sieve_pz);
  T->SetBranchAddress("sieve_r",&sieve_r);
  T->SetBranchAddress("gem1_r",&gem1_r);
  T->SetBranchAddress("rate",&rate);
  T->SetBranchAddress("gem1_ph",&gem1_ph);
  T->SetBranchAddress("gem1_vx",&gem1_vx);
  T->SetBranchAddress("gem1_vy",&gem1_vy);
  T->SetBranchAddress("gem1_vz",&gem1_vz);

  for(int i=0; i<nEntries; i++){

    T->GetEntry(i);  rate = 1; int index_side = -999; int index_sector = -999;

     if(sieve_r < 26.5) continue;
     if(gem1_ph < 0) gem1_ph += 2*pi;

     for( int k=0; k<nsectors; k++){ 

       if(gem1_ph > phi_lo[k] && gem1_ph < phi_hi[k]) index_sector = k;

       if(gem1_r < gem_rcut[k]) index_side = 0;

       else if(gem1_r > gem_rcut[k]) index_side = 1;  }

     if(index_sector == -999) continue;

    double p = sqrt(sieve_px*sieve_px + sieve_py*sieve_py + sieve_pz*sieve_pz);
    double angle = acos(sieve_pz/p)*(180/pi);

    rsieve[index_sector][index_side]->Fill(sieve_r);
    zvertex[index_sector][index_side]->Fill(gem1_vz);
    th_p[index_sector][index_side]->Fill(angle, p/1000);

  }

 
  TFile* fout = new TFile("gemrphi_Ana_new_sieve.root","RECREATE");
  fout->cd();
  for(int i=0; i<nsectors; i++){

    for(int j=0; j<nsides; j++){

     th_p[i][j]->Write();
     rsieve[i][j]->Write();
     zvertex[i][j]->Write(); } }

  fout->Close();

}
