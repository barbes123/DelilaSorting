#include "TChain.h"
#include "TProof.h"
#include <sstream>
#include <string>
#include <iostream> 
using namespace std;


string data_dir = "/data/2022_w3/root_files/";


void Delila_selector_chain(UInt_t first_run=195,  UInt_t last_run=195, UInt_t vol0=1, UInt_t vol1=1){

 
 std::cout<<" Chain_selector.C is running "<<std::endl;
 std::cout<<" vol0 "<<vol0<<" vol1 "<<vol1<<std::endl;
  
 TChain *ch = new TChain("ELIADE_Tree","ELIADE_Tree");

 for(UInt_t run=first_run;run<=last_run;++run){     
        for (UInt_t vol=vol0;vol<=vol1;++vol){
         
        //TChain *ch = new TChain("ELIADE_Tree","ELIADE_Tree");
        string szRun, szVol;

        szRun = Form("%i",run);szVol = Form("%i",vol);

        std::stringstream ifile;
        ifile<<Form("%s/run%s_%s_ssgant1.root", data_dir.c_str(), szRun.c_str(),szVol.c_str());
        std::cout<<"File "<<ifile.str().c_str()<<std::endl;  
        ch->Add(Form("%s/run%s_%s_ssgant1.root", data_dir.c_str(), szRun.c_str(),szVol.c_str()));
        //double beta=0.0; 
        //std::ostringstream options;
        //options<<run<<","<<vol<<","<<beta<<","<<","<<0<<","<<"0";
        //ch->Process("~/DelilaSorting/DelilaSelector.C+",options.str().c_str());
         
        };  
       	double beta=0.0;  
        std::ostringstream options;
        options<<first_run<<","<<vol1<<","<<beta<<","<<","<<0<<","<<"0";
        ch->Process("~/DelilaSorting/DelilaSelector.C+",options.str().c_str());
  };
}

