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
  double implant2[6][128][128];
  double implantE[6][128][128];

  double decay_old[6][128][128];
  double decay_old2[6][128][128];
  double decay_oldE[6][128][128];

  double decay_oldB[6][128][128];
  double decay_oldB2[6][128][128];
  double decay_oldBE[6][128][128];

  double at_no[6][128][128];
  double mass_no[6][128][128];

  for(int i=0;i<6;i++){
    for(int j=0;j<128;j++){
      for(int k=0;k<128;k++){
	implant[i][j][k]=-1;
	implant2[i][j][k]=-1;
	implantE[i][j][k]=-1;

	decay_old[i][j][k]=-1;
	decay_old2[i][j][k]=-1;
	decay_oldE[i][j][k]=-1;

	decay_oldB[i][j][k]=-1;
	decay_oldB2[i][j][k]=-1;
	decay_oldBE[i][j][k]=-1;

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

  TH1I * Hdecay2[6];
  TH1I * Hdecay_old2[6];

  Hdecay2[0]= new TH1I("Hdecay20","DSSD1;t",Nbin,t0,t1);
  Hdecay2[1]= new TH1I("Hdecay21","DSSD2;t",Nbin,t0,t1);
  Hdecay2[2]= new TH1I("Hdecay22","DSSD3;t",Nbin,t0,t1);
  Hdecay2[3]= new TH1I("Hdecay23","DSSD4;t",Nbin,t0,t1);
  Hdecay2[4]= new TH1I("Hdecay24","DSSD5;t",Nbin,t0,t1);
  Hdecay2[5]= new TH1I("Hdecay25","DSSD6;t",Nbin,t0,t1);

  Hdecay_old2[0]= new TH1I("Hdecay_old02","DSSD1;t (prev)",Nbin,t0,t1);
  Hdecay_old2[1]= new TH1I("Hdecay_old12","DSSD2;t (prev)",Nbin,t0,t1);
  Hdecay_old2[2]= new TH1I("Hdecay_old22","DSSD3;t (prev)",Nbin,t0,t1); 
  Hdecay_old2[3]= new TH1I("Hdecay_old32","DSSD4;t (prev)",Nbin,t0,t1);
  Hdecay_old2[4]= new TH1I("Hdecay_old42","DSSD5;t (prev)",Nbin,t0,t1);
  Hdecay_old2[5]= new TH1I("Hdecay_old52","DSSD6;t (prev)",Nbin,t0,t1); 

  TH1I * HdecayE[6];
  TH1I * Hdecay_oldE[6];

  HdecayE[0]= new TH1I("Hdecay0E","DSSD1;t",Nbin,t0,t1);
  HdecayE[1]= new TH1I("Hdecay1E","DSSD2;t",Nbin,t0,t1);
  HdecayE[2]= new TH1I("Hdecay2E","DSSD3;t",Nbin,t0,t1);
  HdecayE[3]= new TH1I("Hdecay3E","DSSD4;t",Nbin,t0,t1);
  HdecayE[4]= new TH1I("Hdecay4E","DSSD5;t",Nbin,t0,t1);
  HdecayE[5]= new TH1I("Hdecay5E","DSSD6;t",Nbin,t0,t1);

  Hdecay_oldE[0]= new TH1I("Hdecay_old0E","DSSD1;t (prev)",Nbin,t0,t1);
  Hdecay_oldE[1]= new TH1I("Hdecay_old1E","DSSD2;t (prev)",Nbin,t0,t1);
  Hdecay_oldE[2]= new TH1I("Hdecay_old2E","DSSD3;t (prev)",Nbin,t0,t1); 
  Hdecay_oldE[3]= new TH1I("Hdecay_old3E","DSSD4;t (prev)",Nbin,t0,t1);
  Hdecay_oldE[4]= new TH1I("Hdecay_old4E","DSSD5;t (prev)",Nbin,t0,t1);
  Hdecay_oldE[5]= new TH1I("Hdecay_old5E","DSSD6;t (prev)",Nbin,t0,t1); 

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
  //input_file->ls();
  //read tree from file
  treeD= (TTree*) input_file->Get("aida");
  // treeD->Print();
  treeD->SetBranchAddress("aida",&data_d);


 

 // runName=  "/Disk/ds-sopa-personal/aestrade/rootfiles/AIDA/May2015/aida_sort_"+inFile+".root";
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
  // cout << "\nNumber of entries in input tree= " << n_entries << ", in file " << inFile << endl;


  n_entries = treeI->GetEntries();
   if(Nmax>0 && Nmax<n_entries) N_i = Nmax;
   else N_i = n_entries; 
 cout<<"N_i "<<N_i<<endl;
    // cout << "\nNumber of entries in input tree D= " << n_entries << ", in file " << inFile2 << endl;
    
  
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
      // cout<<"CHECK2POINT"<<endl;
      decay_old[int(data_d.z)][int(data_d.x)][int(data_d.y)]= data_d.T;
    }
    
    
    // cout<<"CHECKPOINt3"<<endl;
    
 }
   //cout<<"Checkpoint4"<<endl;
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


  TH1I * HDiff2[6];
  HDiff2[0]= new TH1I("HDiff2_0","DSSD1 !2 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiff2[1]= new TH1I("HDiff2_1","DSSD2 !2 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiff2[2]= new TH1I("HDiff2_2","DSSD3 !2 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiff2[3]= new TH1I("HDiff2_3","DSSD4 !2 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiff2[4]= new TH1I("HDiff2_4","DSSD5 !2 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiff2[5]= new TH1I("HDiff2_5","DSSD6 !2 (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  TCanvas *c1;
  c1 = new TCanvas("c1","c1"); c1->Divide(2,3);
  c1->cd(1); HDiff2[0]->Draw();
  c1->cd(2); HDiff2[1]->Draw();
  c1->cd(3); HDiff2[2]->Draw();
  c1->cd(4); HDiff2[3]->Draw();
  c1->cd(5); HDiff2[4]->Draw();
  c1->cd(6); HDiff2[5]->Draw();

  TH1I * HDiffE[6];
  HDiffE[0]= new TH1I("HDiffE_0","DSSD1 !E (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiffE[1]= new TH1I("HDiffE_1","DSSD2 !E (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiffE[2]= new TH1I("HDiffE_2","DSSD3 !E (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiffE[3]= new TH1I("HDiffE_3","DSSD4 !E (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiffE[4]= new TH1I("HDiffE_4","DSSD5 !E (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  HDiffE[5]= new TH1I("HDiffE_5","DSSD6 !E (backgrnd. substr.);t (prev)",Nbin,t0,t1);
  TCanvas *c2;
  c2 = new TCanvas("c2","c2"); c2->Divide(2,3);
  c2->cd(1); HDiffE[0]->Draw();
  c2->cd(2); HDiffE[1]->Draw();
  c2->cd(3); HDiffE[2]->Draw();
  c2->cd(4); HDiffE[3]->Draw();
  c2->cd(5); HDiffE[4]->Draw();
  c2->cd(6); HDiffE[5]->Draw();




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

    //    c_d[i]= 1.*Hdecay[i]->Integral()/Hdecay_old[i]->Integral();
    //    HDiff[i]->Add(Hdecay[i],Hdecay_old[i],1,-1*c_d[i]);
        HDiff[i]->Add(Hdecay[i],Hdecay_old[i],1,-1);

    //    c_d2[i]= 1.*Hdecay2[i]->Integral()/Hdecay_old2[i]->Integral();
    //    HDiff2[i]->Add(Hdecay2[i],Hdecay_old2[i],1,-1*c_d2[i]);
      HDiff2[i]->Add(Hdecay2[i],Hdecay_old2[i],1,-1);

    //  c_de[i]= 1.*HdecayE[i]->Integral()/Hdecay_oldE[i]->Integral();
    //    HDiffE[i]->Add(HdecayE[i],Hdecay_oldE[i],1,-1*c_de[i]);
      HDiffE[i]->Add(HdecayE[i],Hdecay_oldE[i],1,-1);
  }
  c0->Write();
  c1->Write();
  c2->Write();
  c3->Write();
  cc->Write();

 
 

  fOut->Close();

  //input_file.close();
}

