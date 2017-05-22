#include <bitset>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom3.h>


using namespace std;

//void Plot(Long64_t time_window= 300, bool b_debug= false, int run0= 3002, Long64_t Nmax_a= -1, int Nmax_b= -1, Long64_t offset=0, int Nscale= 1, Long64_t ts_reset= 20000000000){
void Plot(bool b_debug= false, int run0= 3002, Long64_t Nmax_a= -1, int Nmax_b= -1, Long64_t offset=0, int Nscale= 1, Long64_t ts_reset= 20000000000){
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


  struct aida{
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

  TH1D * DeltaTS[3]; //for histogram of ta-tb or delta
  TCanvas * c1;  
  c1 = new TCanvas("c1","c1");
  DeltaTS[0]= new TH1D("DeltaTSZ0", "#Delta time-stamp !AIDA BigRIPS coincidence;t_{BigRIPS} - t_{AIDA} [10ns]", 10000, -10e5, 10e5);
  DeltaTS[1]= new TH1D("DeltaTSZ1", "#Delta time-stamp !AIDA behind;t_{BigRIPS} - t_{AIDA} [10ns]", 10000, -10e5, 10e5);
  DeltaTS[2]= new TH1D("DeltaTSZ2", "#Delta time-stamp !BigRIPS behind;t_{BigRIPS} - t_{AIDA} [10ns]", 10000, -10e5, 10e5);

  TFile* input_file1; 
  TFile* input_file2; 
  TFile* input_file3; 

  TTree * tree1;
  TTree * tree2;

  TTree * treeOut;

  aida event_data;

  merge_data_struct merged_data;


  TFile *fOut = new TFile("testmerger2.root","RECREATE");

  treeOut = new TTree("AB_implant","AB_implant");
  treeOut->Branch("implant",&merged_data,"ts/L:zet/F:aoq/F:xm/I:ym/I:zm/I:EX/I:EY/I:ID/I");




   const Int_t kMaxEventInfo = 1;
   const Int_t kMaxBigRIPSRaw = 10;



   Long64_t       ts;
   Long64_t       sts;
   Double_t        tof;
   Double_t        zet;
   Double_t        aoq;
   Double_t        f5x;
   Double_t        f11x;
   Double_t        f11y;
   Double_t        f11dt;
   Double_t        beta;
   
   // List of branches
   TBranch        *b_bigrips_ts;   //!
   TBranch        *b_bigrips_sts;   //!
   TBranch        *b_bigrips_tof;   //!
   TBranch        *b_bigrips_zet;   //!
   TBranch        *b_bigrips_aoq;   //!
   TBranch        *b_bigrips_f5x;   //!
   TBranch        *b_bigrips_f11x;   //!
   TBranch        *b_bigrips_f11y;   //!
   TBranch        *b_bigrips_f11dt;   //!
   TBranch        *b_bigrips_beta;   //!
   
  //AIDA FILE
  string runName;
  string path="/home/neerajan/HalfLife/";


  //PID FILE
  // string path2="/Disk/ds-sopa-personal/aestrade/rootfiles/AIDA/May2015/";
  char pid_file[1024];

  //number of data points in TTree
    Long64_t n_entries;
   Long64_t n_max;

   Long64_t delta_ts=0;
   Long64_t  ts_imp_prev=0;
   Long64_t ts_dec_prev=0;
   Long64_t ts_corr_prev=0;

   // Long64_t TW= time_window;
   // Long64_t minusTW= -1*time_window;
   Long64_t TW= 15000;
   Long64_t minusTW= -15000;

  int i_a=0;
  int i_b=0;
   Long64_t N_a;
   Long64_t N_b;

   Long64_t count=0;

   Long64_t t_a_scale=0; //time_externa * Nscale
   Long64_t t_b=0;

   Long64_t dt_prev;

   Long64_t dt_prev_a=9e12;
   Long64_t dt_prev_b=9e12;
   Long64_t t_prev_a=0;
   Long64_t t_prev_b=0;

  cout << "opening AIDA file...."<<endl;

  // runName= path+"aida_sort_"+runNum+".root";
  runName = path+"161110_1940_aida27_3to18.root";
  input_file1 = new TFile(runName.data(),"read");
  tree1= (TTree*) input_file1->Get("aida");
  //tree1->Print();
  tree1->SetBranchAddress("aida",&event_data,0);

  N_a = tree1->GetEntries();
  cout << "\nNumber of entries in AIDA tree= " << N_a << endl;
  if(Nmax_a>0 && N_a>Nmax_a) N_a= Nmax_a; //loop only through part of file?

  //  tree1->Show(24);

  //find first data point in AIDA with correlation scaler set
  bool b_corr_flag=false;
  for( Long64_t i=0;i<N_a;i++){

    tree1->GetEntry(i);


    if(event_data.T > 0){
      t_a_scale = event_data.T * Nscale;
      i_a=i;
      
      cout << "\nFound first AIDA entry with CORR SCALER set: i= "<< i_a << ", ts= "<< t_a_scale << endl;
      b_corr_flag= true;
      break;
    }
  }
  
  if(!b_corr_flag){
    cout << "\n!!! Did not find any entry in AIDA with CORR SCALER flag set: ts(last)= " <<  event_data.T<<" !!! "<<endl;			 
    cout << "Exittttttttttttttt..........."<<endl;
    exit(1);
  }


  cout << "\nopening BigRIPS file...."<<endl;

  //  sprintf(pid_file,"/Disk/ds-sopa-group/np/RIKEN/May2015/PID_data/go4_%d.root",run0);
  sprintf(pid_file,"/home/neerajan/HalfLife/161110_1938_bigrips0444.root",run0);
  input_file2 = new TFile(pid_file,"read");
  
  tree2= (TTree*) input_file2->Get("tree");
  if(!tree2){
    cout <<"\n\nNOOOOOOOOOOOO TREEEEEEEEEEEEEEEEEEEEEEEEE"<<endl;
    return;
  }

  tree2->SetMakeClass(1);
  tree2->SetBranchAddress("ts", &ts, &b_bigrips_ts);
  tree2->SetBranchAddress("sts", &sts, &b_bigrips_sts);
  tree2->SetBranchAddress("tof", &tof, &b_bigrips_tof);
  tree2->SetBranchAddress("zet", &zet, &b_bigrips_zet);
  tree2->SetBranchAddress("aoq", &aoq, &b_bigrips_aoq);
  tree2->SetBranchAddress("f5x", &f5x, &b_bigrips_f5x);
  tree2->SetBranchAddress("f11x", &f11x, &b_bigrips_f11x);
  tree2->SetBranchAddress("f11y", &f11y, &b_bigrips_f11y);
  tree2->SetBranchAddress("f11dt", &f11dt, &b_bigrips_f11dt);
  tree2->SetBranchAddress("beta", &beta, &b_bigrips_beta);

  N_b = tree2->GetEntries();
  cout << "\nNumber of entries in PID tree= " << N_b << endl; 
  if(Nmax_b>0 && N_b>Nmax_b) N_b= Nmax_b; //loop only through part of file?

  tree2->GetEntry(1);

  cout << " \n***ZET, AOQ: " << zet << "  " << aoq << endl;
  tree2->Show(1);


  bool b_new=false;
  bool b_coinc= false;

  /****** read one entry from each file *****/

  i_a++;
  if(i_a==N_a){
    cout << "\n**** REACHED END OF AIDA FILE ***"<<endl;
    //break;
  }
  else {
    tree1->GetEntry(i_a);
    
    b_new= true;
	
    t_prev_a =  event_data.T * Nscale;
    t_a_scale = event_data.T * Nscale;
	
    //FILL AAA
	
  }
  
  
	
  i_b++;
  if(i_b==N_b){
    cout << "\n**** REACHED END OF BigRIPS FILE ***"<<endl;
    //break;
  }
  else {
    tree2->GetEntry(i_b);
    /*  if( EventInfo_fBit[0]==2){ */
	    
    //	  if(ts>=t_prev_b){ //monotonically increasing
	      
    b_new= true;
    t_b=ts ;
    t_prev_b= ts; 


    cout << "initial values for bigrips: " << t_b <<endl;	    
  }
   

  /******************/



   for(;;){
  // for(int iii=0 ;iii<60000 ;iii++){
     // for(int iii=0 ;iii<1500000 ;iii++){


    delta_ts = t_b - t_a_scale + offset; 

    //If have new data to evaluate....
    if(b_new){


      b_new= false;
      
      //coincidence!
      if( delta_ts < TW && delta_ts > minusTW){
	DeltaTS[0]->Fill(delta_ts);
	
	//	cout << " option zero! coincidence... ta=" << t_a_scale << "  tb= "<< t_a_scale << " delta_ts=" << delta_ts << endl;

	b_coinc= true;
	
	//Fill TTree...
	merged_data.ts = t_a_scale;
	merged_data.zet=zet;
	merged_data.aoq=aoq;
	merged_data.xm= event_data.x;
	merged_data.ym= event_data.y;
	merged_data.zm= event_data.z;
	merged_data.EX= event_data.EX;
	merged_data.EY= event_data.EY;
	merged_data.ID= 1;
	treeOut->Fill();
	
	//get one of each
	i_a++;// cout  << "count i_a 1" << endl;
	if(i_a==N_a){
	  cout << "\n**** REACHED END OF AIDA FILE ***"<<endl;
	  break;
	}
	else {
	  tree1->GetEntry(i_a);
	  
	  if(event_data.ID == 4 && event_data.T > 0){
	    
	    if( event_data.T >= t_prev_a){ //if monotonically increasing
	      
	      b_new= true;
	      
	      t_prev_a = t_a_scale; // event_data.T;
	      t_a_scale = event_data.T * Nscale;
	      
	    }
	    else{
	      cout << "AIDA not monotonically increasing! i: "<<i_a<<":  ts_prev, t_this: "<<	t_prev_a <<"  "<< event_data.T<< endl;
	    }
	  }
	}
	
	i_b++;
	if(i_b==N_b){
	  cout << "\n**** REACHED END OF BigRIPS FILE ***"<<endl;
	  break;
	}
	else {
	  tree2->GetEntry(i_b);
	  /***/
	  //what would be EventInfo_fbit??
	  //if( EventInfo_fBit[0]==2){
	    
	  if(ts>=t_prev_b){ //monotonically increasing
	      
	    //if(b_new) b_new= true; //if both A+B are ok data points for next comparison
	    b_new= true; //if both A+B are ok data points for next comparison
	    t_b=ts;
	    t_prev_b=ts;
	  }
	  else{
	    cout << "BigRIPS not monotonically increasing! i: "<<i_a<<":  ts_prev, t_this: "<<	t_prev_b <<"  "<< ts<< endl;
	  }
	    //} /***/
	}
	
      }//if inside coincidence window
      else{ //one of them outside... 
	
	//if AIDA behind, get new AIDA entry for implant
	if( ts > t_a_scale){ //delta_ts>0
	  DeltaTS[1]->Fill(delta_ts);
	 
 //DeltaTS_Merge->Fill(delta_ts);
	  // cout << "option2:  ia, ib: " << i_a << "  " << i_b << endl;
	  //	  cout<<"check"<<endl;


	  //	  cout << " done with option 2 ...." << endl;
	}
	//if BigRIPS behind, get new BIG RIPS
	else{
	  DeltaTS[2]->Fill(delta_ts);
	 
	  //	  cout << "option3:  ia, ib: " << i_a << "  " << i_b << endl;
	  //DeltaTS_Merge->Fill(delta_ts);

	}
	
      }//if outside coincidence

    }
    else{ //if need new data

      //      b_new= false; //need new data.... somewhat implicit?

      //if AIDA behind, get new AIDA entry for implant
      if( t_b > t_a_scale){ //delta_ts>0

	i_a++; //cout << " count i_a 2 " << endl;
	if(i_a==N_a){
	  cout << "\n**** REACHED END OF AIDA FILE ***"<<endl;
	  break;
	}
	else {
	  tree1->GetEntry(i_a);


	  //	 

	  if(event_data.ID == 4 && event_data.T > 0){

	    if( event_data.T >= t_prev_a){ //if monotonically increasing

	      b_new= true;
	      
	      t_prev_a =  t_a_scale;
	      t_a_scale = event_data.T * Nscale;

	       //FILL AAA

	    }
	    else{
	      cout << "AIDA not monotonically increasing! i: "<<i_a<<":  ts_prev, t_this: "<<	t_prev_a <<"  "<< event_data.T<< endl;
	    }
	  }
	}
      }//if A behind
      //if BigRIPS behind, get new BIG RIPS
      else{
	
	//	cout << " big rips needs new data! " << endl;
	i_b++;
	if(i_b==N_b){
	  cout << "\n**** REACHED END OF BigRIPS FILE ***"<<endl;
	  break;
	}
	else {
	  tree2->GetEntry(i_b);
	  /*  if( EventInfo_fBit[0]==2){ */
	   
	 
	  if(ts>=t_prev_b){ //monotonically increasing
	      
	    //	    if( (iii%2)==0)  cout << " *** values BIGRIPS: " << ts << "  " << t_prev_a << endl;
	    b_new= true;
	    t_b=ts ;
	    t_prev_b= ts; 
	    
	    //FILL BBB
	    
	    
	  }
	  else{
	    cout << "BigRIPS not monotonically increasing! i: "<<i_a<<":  ts_prev, t_this: "<<	t_prev_b << " " << ts << endl;
	  }
	  /*****  }*/
	}
      }//if B behind
    }//if need new data

    //  Delta->SetLineColor(1); //histogram for ta-tb or delta_ts
   

    //cout<<" final values of timestamps = "<<t_b<<", taida= "<<t_a_scale<<endl; 
 }// FOR EVER AND EVER

    DeltaTS[0]->SetLineColor(kBlue);
    DeltaTS[0]->Draw();
    DeltaTS[1]->SetLineColor(kBlue);
    DeltaTS[1]->Draw("same");
    DeltaTS[2]->SetLineColor(kBlue);
    DeltaTS[2]->Draw("same");
  cout << " final values of timestamps... ts=" << t_b<<",  taida=" << t_a_scale << endl;



  // c1->Write();
  fOut->cd();
  c1->Write();
  treeOut->Write();
  fOut->Close();

  input_file1->Close();
  input_file2->Close();

  cout << "\n out written to file..."<<endl;

}



//  LocalWords:  hEE
