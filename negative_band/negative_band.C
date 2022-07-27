TFile* fout;
const int nPlanes = 9;
const int nsectors = 7;
const double pi = acos(-1);
const char* PlaneNames[nPlanes] = {"sieve_xy", "dscoil3_xy", "dscoil6_xy", "dscoil9_xy", "dscoil12_xy", "dscoil13_xy", "collar1Ent_xy", "collar1Exit_xy", "main_xy" };
const char* sector[nsectors] = {"sector7", "sector6", "sector5", "sector4", "sector3" ,"sector2", "sector1"};
double phi_lo[nsectors] = {0, 2*pi/7, 4*pi/7, 6*pi/7, 8*pi/7, 10*pi/7, 12*pi/7};
double phi_hi[nsectors] = {2*pi/7, 4*pi/7, 6*pi/7, 8*pi/7, 10*pi/7, 12*pi/7, 2*pi};
double radial_cut[nPlanes] = { 26.5, 170, 200, 240, 320, 320, 500, 500, 500 };

TH2D* xy[nPlanes];
TH1D* r[nPlanes][nsectors];
std::map<int, int> dt{{270,1}, {51,2}, {54,3}, {57,4}, {5696,5}, {5697,6}, {70,7}, {71,8}, {28,9}};

long nTotEv(0);
string fin;

void initHisto();
long processOne(string);
void process(int);
void writeOutput();

void negative_band(const string& finName = "./remollout.root", int test_run=0){
  fin = finName;
  initHisto();
  process(test_run);
  writeOutput();
}

void process(int test_run){

  if(fin==""){
    cout<<"\t did not find input file. Quitting!"<<endl;
    return 2;
  }

  int nFiles(0);
  if( fin.find(".root") < fin.size() ){
    cout<<"Processing single file:\n\t"<<fin<<endl;
    nTotEv+=processOne(fin);
    nFiles=1;
  }

   else{
    cout<<"Attempting to process list of output from\n\t"<<fin<<endl;
    ifstream ifile(fin.c_str());
    string data;
    while(ifile>>data){
      cout<<" processing: "<<data<<endl;
      nTotEv+=processOne(data);
      nFiles++;
    
      if(test_run==1 && nFiles==10) break;

    }
  }
  
  cout<<"\nFinished processing a total of "<<nTotEv<<endl;
}

void initHisto(){

  fout = new TFile("negative_band.root", "RECREATE");

  for(int i=0; i<nPlanes; i++){
   
    xy[i] = new TH2D(Form("xy_%s",PlaneNames[i]), Form("xy_%s",PlaneNames[i]), 2000, -1200, 1200, 2000, -1200, 1200);

    for(int j=0; j<nsectors; j++){

      r[i][j] = new TH1D(Form("r_%s_%s",PlaneNames[i], sector[j]), Form("r_%s_%s",PlaneNames[i], sector[j]), 1000, 0, 800);

  }

  }

} 

long processOne(string fnm){

  TFile *fin=TFile::Open(fnm.c_str(),"READ");
  if(!fin->IsOpen() || fin->IsZombie()){
    cout<<"Problem: can't find file: "<<fnm<<endl;
    fin->Close();
    delete fin;
    return 0;
  }

  TTree *t=(TTree*)fin->Get("T");
  if (t == 0) return 0;
  Double_t rate=0;
  std::vector<remollGenericDetectorHit_t> *hit=0;
  t->SetBranchAddress("rate", &rate);
  t->SetBranchAddress("hit", &hit);
  
  long nEntries = t->GetEntries();

  for (Long64_t event = 0; event < nEntries; event++) {
   
    t->GetEntry(event); 
   
    if(event%10000==0) cout << "At tree entry = " << event << endl;

    if(std::isnan(rate) || std::isinf(rate)) continue;

    for(int j=0; j<hit->size(); j++){

      int det = hit->at(j).det;
      int pid = hit->at(j).pid;
      int trid = hit->at(j).trid;
      int mtrid = hit->at(j).mtrid;
      double pz = hit->at(j).pz;
      double px = hit->at(j).px;
      double py = hit->at(j).py;
      double rr = hit->at(j).r; 
      double ph = hit->at(j).ph;
      if(ph<0) ph = ph + 2*pi;
      if(rr < 30) continue;
      if( pz < 0 ) continue;
      int index = dt[det] - 1;
      int sec_index = -999;
      if(index == -1) continue; 
      if(!(pid==11 && trid==1 && mtrid==0)) continue;

      for(int k=0; k<nsectors; k++){ if(ph > phi_lo[k] && ph < phi_hi[k]) sec_index = k; }

      if( sec_index == -999) continue;

      xy[index]->Fill(hit->at(j).x, hit->at(j).y, rate);
      if(rr > radial_cut[index]) r[index][sec_index]->Fill(rr, rate);
    
    }

  }

  fin->Close();
  delete fin;
  return nEntries;
}
    

void writeOutput(){
  fout->cd();

  for(int i=0; i<nPlanes; i++){ 

     xy[i]->Write(); 

    for(int j=0; j<nsectors; j++){

      r[i][j]->Write();

      } 

  }

  fout->Close();

}
