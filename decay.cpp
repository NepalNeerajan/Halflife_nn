#include <bitset>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TCanvas.h>

//using namespace std;


void ReadTree(int gEx=200, int gEy=200, Float_t Zmin=25.0, Float_t Zmax=40.0, Float_t AoQmin= 2.0, Float_t AoQmax=4.0, int Nbin=100, Double_t t0= 1000.0, Double_t t1=25000.0, int Nmax= -1, Long64_t dt_min= 200.0 /* usec */){


struct merge_data_struct{
    Long64_t ts;
    Float_t zet;
    Float_t aoq;
    int xm;
    int ym;
    int zm;
    int EX;
    int EY;
    int ID;
  };


  struct luckyaida{
   ULong64_t T;
   ULong64_t Tfast;
    double E;
    double EX;
    double EY;
    double x;
    double y;
    double z;
    int nx;
    int ny;
    int nz;
    unsigned char ID;
  };

  double implant[6][128][128];
  double decay_old[6][128][128];

  double at_no[6][128][128];
  double mass_no[6][128][128];

  for(int i=0;i<6;i++){
    for(int j=0;j<128;j++){
      for(int k=0;k<128;k++){
	implant[i][j][k]=-1;
	decay_old[i][j][k]=-1;

	at_no[i][j][k]=0;
	mass_no[i][j][k]=0;
      }
    }
  }
  TH2I * HPID;
  HPID = new TH2I("HPID","PID, aoq, zet",100,2.55,2.8,100,27,32);

  TH1I * Hdecay[6];
  TH1I * Hdecay_old[6];

  Hdecay[0]= new TH1I("Hdecay0","DSSD1;t",Nbin,t0,t1);
  Hdecay[1]= new TH1I("Hdecay1","DSSD2;t",Nbin,t0,t1);
  Hdecay[2]= new TH1I("Hdecay2","DSSD3;t",Nbin,t0,t1);
  Hdecay[3]= new TH1I("Hdecay3","DSSD4;t",Nbin,t0,t1);
  Hdecay[4]= new TH1I("Hdecay4","DSSD5;t",Nbin,t0,t1);
  Hdecay[5]= new TH1I("Hdecay5","DSSD6;t",Nbin,t0,t1);

  Hdecay_old[0]= new TH1I("Hdecay_old0","DSSD1;t (prev)",Nbin,t0,t1);
  Hdecay_old[1]= new TH1I("Hdecay_old1","DSSD2;t (prev)",Nbin,t0,t1);
  Hdecay_old[2]= new TH1I("Hdecay_old2","DSSD3;t (prev)",Nbin,t0,t1); 
  Hdecay_old[3]= new TH1I("Hdecay_old3","DSSD4;t (prev)",Nbin,t0,t1);
  Hdecay_old[4]= new TH1I("Hdecay_old4","DSSD5;t (prev)",Nbin,t0,t1);
  Hdecay_old[5]= new TH1I("Hdecay_old5","DSSD6;t (prev)",Nbin,t0,t1); 

  TH1I * HEx[6];
  HEx[0]= new TH1I("HE0x","Ex DSSD1 !E;ADC ch",512,0,16384);
  HEx[1]= new TH1I("HE1x","Ex DSSD2 !2;ADC ch",512,0,16384);
  HEx[2]= new TH1I("HE2x","Ex DSSD3 !E;ADC ch",512,0,16384);
  HEx[3]= new TH1I("HE3x","Ex DSSD4 !E;ADC ch",512,0,16384);
  HEx[4]= new TH1I("HE4x","Ex DSSD5 !2;ADC ch",512,0,16384);
  HEx[5]= new TH1I("HE5x","Ex DSSD6 !E;ADC ch",512,0,16384);
  TH1I * HEy[6];
  HEy[0]= new TH1I("HE0y","Ey DSSD1 !E;ADC ch",512,0,16384);
  HEy[1]= new TH1I("HE1y","Ey DSSD2 !2;ADC ch",512,0,16384);
  HEy[2]= new TH1I("HE2y","Ey DSSD3 !E;ADC ch",512,0,16384);
  HEy[3]= new TH1I("HE3y","Ey DSSD4 !E;ADC ch",512,0,16384);
  HEy[4]= new TH1I("HE4y","Ey DSSD5 !2;ADC ch",512,0,16384);
  HEy[5]= new TH1I("HE5y","Ey DSSD6 !E;ADC ch",512,0,16384);
  int Npixel=1;

  string runName;
  string runName1;
  int Nruns;

  TFile* input_file; 
  TTree * treeD;

  TFile* input_file2; 
  TTree * treeI;
  
  //tree from Converter step will data in a struct_entry_midas structure
  luckyaida data_d;
  merge_data_struct data_i;

  runName=  "/home/neerajan/HalfLife/161110_1940_aida27_3to18.root";
  input_file = new TFile(runName.data(),"read");
  
  //  if (input_file != 0){
  cout << "Opened file " << runName << endl;
  //read tree from file
  treeD= (TTree*) input_file->Get("aida");
  // treeD->Print();
  treeD->SetBranchAddress("aida",&data_d);

  runName1 = "/home/neerajan/HalfLife/testmerger2.root";
  input_file2 = new TFile(runName1.data(),"read");
  //read tree from file
  treeI= (TTree*) input_file2->Get("AB_implant");
  //treeI->Print();
  treeI->SetBranchAddress("implant",&data_i,0);

  //number of data points in TTree
  long long n_entries;
  long long N_d, N_i;
  n_entries = treeD->GetEntries();
    if(Nmax>0 && Nmax<n_entries) N_d = Nmax;
   else N_d = n_entries; 
 cout<<"N_d "<<N_d<<endl;

  n_entries = treeI->GetEntries();
   if(Nmax>0 && Nmax<n_entries) N_i = Nmax;
   else N_i = n_entries; 
 cout<<"N_i "<<N_i<<endl;
  
  Long64_t dt;
  int my_x, my_y;
  
  //  long long N_i=0;
  // long long N_d=0;
    
  Long64_t t_next_i=0;
  Long64_t t_next_d=0;
  
  bool b_implant= false;
  bool b_decay= false;
  
  bool b_end_implant= false;
  
  long long i_i=0;
  long long i_d=0;

  int cnt=0;

  bool b_break= false;
  //now looping only over first 99 events, or up to n_entries, whatever is less
  //   for(long long i=0; i<NloopD; i++){


 TFile *fOut;

  string foutname="out_decay.root"; 
  fOut= new TFile(foutname.data(),"recreate");

  int ddd=2000000;
  for(int i = 0; i<ddd; i++){
  //for(;;){

    b_implant= false;
    b_decay= false;
    
    // if(b_break) break;

    treeD->GetEntry(i_d);
    if(!b_end_implant) treeI->GetEntry(i_i);
    
      if(data_i.ts > data_d.T || b_end_implant){
    // if(b_end_implant){
      //if good decay
       if((data_d.ID)==5){
	//E cut
	if( data_d.EX > gEx && data_d.EY > gEy){	 
          b_implant= false;
	  b_decay= true;
	  //cout<<"b_decay"<<b_decay<<endl;
	}
      }
      i_d++;
      if(i_d==N_d){
	cout << "\n**** REACHED END OF decay FILE ***"<<endl;
	b_break=true;
	break;
      }
    }
    else {
      //if good implant
      if(data_i.ID==1){
	//check PID!!!!
	if( data_i.zet > Zmin && data_i.zet < Zmax){
	  if(data_i.aoq > AoQmin && data_i.aoq < AoQmax){
	    b_implant= true;
	    b_decay= false;
	    //cout<<"b_implant"<<b_implant<<endl;
	  }
	}
      }
      i_i++;
      if(i_i==N_i){
	cout << "\n**** REACHED END OF implant FILE ***"<<endl;
	b_end_implant = true;
      }
    }
    
    //if IMPLANT
    if(b_implant){ 
      implant[data_i.zm][data_i.xm][data_i.ym]= data_i.ts;
      HPID->Fill(data_i.aoq, data_i.zet);
      //  cout<<"ARRAY"<<decay_old[data_i.zm][data_i.xm][data_i.ym]<<endl;
      //??????????????/
      if((decay_old[data_i.zm][data_i.xm][data_i.ym])>0){
	  //  cout<<"TEST"<<data_i.zm<<endl;}
	dt= (data_i.ts - decay_old[data_i.zm][data_i.xm][data_i.ym])*4/1.e8; //in seconds
	//	cout<<"dt of implant: "<<dt<<endl;
	if(dt>dt_min/1.e6){
	  Hdecay_old[data_i.zm]->Fill(dt);
	  decay_old[data_i.zm][data_i.xm][data_i.ym]= -1;
	}
      } //??????
    }
    

    //IF DECAY
    else if(b_decay){
      if(implant[int(data_d.z)][int(data_d.x)][int(data_d.y)]>0){
	dt= (data_d.T-implant[int(data_d.z)][int(data_d.x)][int(data_d.y)])*4/1.e8; //in seconds
	//	cout<<"dt of decay: "<<dt<<endl;
	if(dt>dt_min/1.e6){  
	  Hdecay[int(data_d.z)]->Fill(dt);
	  int aa= int(data_d.z);	    
	  switch (aa){
	  case 0:
	  HEx[0]->Fill(data_d.EX);
	  HEy[0]->Fill(data_d.EY);
	  break;
	  case 1:
	  HEx[1]->Fill(data_d.EX);
	  HEy[1]->Fill(data_d.EY);
	  break;
	  case 2:
	  HEx[2]->Fill(data_d.EX);
	  HEy[2]->Fill(data_d.EY);
	  break;
	  case 3:
	  HEx[3]->Fill(data_d.EX);
	  HEy[3]->Fill(data_d.EY);
	  break;
	  case 4:
	  HEx[4]->Fill(data_d.EX);
	  HEy[4]->Fill(data_d.EY);
	  break;
	  case 5:
	  HEx[5]->Fill(data_d.EX);
	  HEy[5]->Fill(data_d.EY);
	  break;
	  }
	  //cout<<"data_d.EX "<<data_d.EX<<endl;
	  implant[int(data_d.z)][int(data_d.x)][int(data_d.y)]= -1;
	}
      }
      decay_old[int(data_d.z)][int(data_d.x)][int(data_d.y)]= data_d.T;
    }
 } //Done looing 
  cout << "\nDone looping through entries..."<<endl;
  //}
   
  TCanvas *cc;
  cc = new TCanvas("cc", "cc", 30,30,600,600); 
  HPID->Draw("cc");

  TCanvas *hdecay;
  hdecay = new TCanvas("hdecay","hdecay"); hdecay->Divide(4,3);
  hdecay->cd(1); Hdecay_old[0]->Draw("");
  hdecay->cd(2); Hdecay_old[1]->Draw(""); 
  hdecay->cd(3); Hdecay_old[2]->Draw("");
  hdecay->cd(4); Hdecay_old[3]->Draw("");
  hdecay->cd(5); Hdecay_old[4]->Draw("");
  hdecay->cd(6); Hdecay_old[5]->Draw("");
  
  hdecay->cd(7); Hdecay[0]->Draw("");
  hdecay->cd(8); Hdecay[1]->Draw(""); 
  hdecay->cd(9); Hdecay[2]->Draw("");
  hdecay->cd(10); Hdecay[3]->Draw("");
  hdecay->cd(11); Hdecay[4]->Draw("");
  hdecay->cd(12); Hdecay[5]->Draw("");

  TH1I * HDiff[6];
  HDiff[0]= new TH1I("HDiff_0","DSSD1 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiff[1]= new TH1I("HDiff_1","DSSD2 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiff[2]= new TH1I("HDiff_2","DSSD3 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiff[3]= new TH1I("HDiff_3","DSSD4 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiff[4]= new TH1I("HDiff_4","DSSD5 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiff[5]= new TH1I("HDiff_5","DSSD6 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  TCanvas *c0;
  c0 = new TCanvas("c0","c0"); c0->Divide(2,3);
  c0->cd(1); HDiff[0]->Draw();
  c0->cd(2); HDiff[1]->Draw();
  c0->cd(3); HDiff[2]->Draw();
  c0->cd(4); HDiff[3]->Draw();
  c0->cd(5); HDiff[4]->Draw();
  c0->cd(6); HDiff[5]->Draw();

  TCanvas *c3;
  c3 = new TCanvas("c3","c3",90,90, 600,600); c3->Divide(3,4);
  c3->cd(1); HEx[0]->Draw();
  // gPad->SetLogy(1);
  c3->cd(2); HEx[1]->Draw(); gPad->SetLogy(1);
  c3->cd(3); HEx[2]->Draw(); gPad->SetLogy(1);
  c3->cd(4); HEx[3]->Draw(); gPad->SetLogy(1);
  c3->cd(5); HEx[4]->Draw(); gPad->SetLogy(1);
  c3->cd(6); HEx[5]->Draw(); gPad->SetLogy(1);
  c3->cd(7); HEy[0]->Draw(); gPad->SetLogy(1);
  c3->cd(8); HEy[1]->Draw(); gPad->SetLogy(1);
  c3->cd(9); HEy[2]->Draw(); gPad->SetLogy(1);
  c3->cd(10); HEy[3]->Draw(); gPad->SetLogy(1);
  c3->cd(11); HEy[4]->Draw(); gPad->SetLogy(1);
  c3->cd(12); HEy[5]->Draw(); gPad->SetLogy(1);    
  
 

 
  for(int i=0;i<6;i++){
    //subtraction of background events estimated with backwards correlations (random correlations)
        HDiff[i]->Add(Hdecay[i],Hdecay_old[i],1,-1);
  }
  c0->Write();
  c3->Write();
  cc->Write();

  fOut->Close();
  //input_file.close();
}

