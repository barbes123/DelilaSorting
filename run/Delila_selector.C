#include "TChain.h"
#include "TProof.h"
#include <sstream>
#include <string>
#include <iostream> 
using namespace std;

void Delila_selector(UInt_t first_run=195,  UInt_t last_run=195, UInt_t vol0=1, UInt_t vol1=1){

 
 std::cout<<" Chain_selector.C is running "<<std::endl;
 std::cout<<" vol0 "<<vol0<<" vol1 "<<vol1<<std::endl;
  
 for(UInt_t run=first_run;run<=last_run;++run){     
 	for (UInt_t vol=vol0;vol<=vol1;++vol){
         
        TChain *ch = new TChain("ELIADE_Tree","ELIADE_Tree");
        string szRun, szVol;

        szRun = Form("%i",run);szVol = Form("%i",vol);

        std::stringstream ifile;

        ifile<<Form("/rosphere/ROSPHERE/DELILA/root_files_sorted/run%s_%s_ssgant1.root",  szRun.c_str(),szVol.c_str());
        std::cout<<"File "<<ifile.str().c_str()<<std::endl;  
        ch->Add(Form("/rosphere/ROSPHERE/DELILA/root_files_sorted/run%s_%s_ssgant1.root",  szRun.c_str(),szVol.c_str()));
         
        std::ostringstream options;
        options<<run<<","<<vol<<","<<","<<","<<0<<","<<"0";
        ch->Process("~/DelilaSorting/DelilaSelector.C+",options.str().c_str());
         
	};   
  };
}

