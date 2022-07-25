
const int nsectors = 2;
const char* secName[nsectors] = {"sector1", "sector7"};
TFile *fout;
TH1D* r[nsectors];
const double pi = acos(-1);
long nTotEv(0);
string fin;

void initHisto();
long processOne(string);
void process(int);
void writeOutput();

void negative_band_Analysis(const string& finName = "./remollout.root", int test_run=0){
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
  string foutNm = "negative_band_analysis.root";

  fout = new TFile(foutNm.c_str(),"RECREATE");
 
  for(int i=0; i<nsectors; i++){
    r[i] = new TH1D(Form("r_%s",secName[i]), Form("r_%s",secName[i]), 100, 0, 100); }
 
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

    int hits_main = 0;
 
    for(int j=0;j<hit->size();j++){

      if(std::isnan(hit->at(j).x) || std::isinf(hit->at(j).x)) continue;
      if(std::isnan(hit->at(j).y) || std::isinf(hit->at(j).y)) continue;

      double ph = hit->at(j).ph;
      if(ph< 0) ph = ph +  2*pi; 
      double rr = hit->at(j).r;
      int det = hit->at(j).det;
      int pid = hit->at(j).pid;
      int trid = hit->at(j).trid;
      int mtrid = hit->at(j).mtrid;
    

      if(!(det==70 && pid==11 && trid==1 && mtrid==0 && rr > 500)) continue;

      
      if( (ph > 0 && ph < 2*pi/7) || (ph > 12*pi/7 && ph < 2*pi) ) hits_main++;
     
    }

    if(hits_main>0){

          for(int k=0;k<hit->size();k++){

               double ph = hit->at(k).ph;
               if(ph<0) ph = ph + 2*pi;
               double rr = hit->at(k).r;
               int det = hit->at(k).det;
               int pid = hit->at(k).pid;
               int trid = hit->at(k).trid;
               int mtrid = hit->at(k).mtrid;

               if(!(det == 270 && pid==11 && trid==1 && mtrid==0 && rr > 26.5)) continue; 
  
               int sec_index = -999;

               if(ph > 0 && ph < 2*pi/7) sec_index = 1;
               if(ph > 12*pi/7 && ph < 2*pi) sec_index = 0;
   
	       if( sec_index == -999 ) continue;
	       
               r[sec_index]->Fill(rr);               

	  }

    }

  }

  fin->Close();
  delete fin;
  return nEntries;
}
    

void writeOutput(){
  fout->cd();

  for(int i=0; i<nsectors; i++){ r[i]->Write(); }

  fout->Close();
}

  

