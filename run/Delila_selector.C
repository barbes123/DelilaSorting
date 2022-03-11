#include "TChain.h"
#include "TProof.h"
#include <sstream>
#include <string>
#include <iostream> 
#include <TROOT.h>
using namespace std;


string data_dir = "/rosphere/2022_w10/root_files";

/*bool is_file_exist(std::stringstream fileName)
{
    std::ifstream in_file(fileName.c_str().str());
    return in_file.good();
}
*/

void Delila_selector(UInt_t first_run=195,  UInt_t last_run=195, UInt_t vol0=1, UInt_t vol1=1, UInt_t numberofevents=0){

 
 std::cout<<" Chain_selector.C is running "<<std::endl;
 std::cout<<" vol0 "<<vol0<<" vol1 "<<vol1<<std::endl;
 
//  double beta=0.0669; 

 for(UInt_t run=first_run;run<=last_run;++run){     
        for (UInt_t vol=vol0;vol<=vol1;++vol){
         
        TChain *ch = new TChain("ELIADE_Tree","ELIADE_Tree");
        string szRun, szVol;

        szRun = Form("%i",run);szVol = Form("%i",vol);

        std::stringstream ifile;
        ifile<<Form("%s/run%s_%s_ssgant1.root", data_dir.c_str(), szRun.c_str(),szVol.c_str());
//        std::cout<<"File "<<ifile.str().c_str()<<std::endl;  
        
        
	 if(gSystem->AccessPathName(ifile.str().c_str())){
	        std::cout << "File "<<ifile.str().c_str()<<" does not exist, skipping" << std::endl;
	        continue;
    	} else {
	        std::cout << "File "<<ifile.str().c_str()<<"  exists! " << std::endl;
	        ch->Add(Form("%s/run%s_%s_ssgant1.root", data_dir.c_str(), szRun.c_str(),szVol.c_str()));
	        std::ostringstream options;
//         options<<run<<","<<vol<<","<<beta<<","<<0<<","<<"0";
	        options<<run<<","<<vol<<","<<0<<","<< 0<<","<<"100";
	       std::cout<<"I will start DelilaSelectorElifant with the options: "<<options.str().c_str()<<std::endl;
        	if (numberofevents == 0){
        	        ch->Process("~/DelilaSorting/DelilaSelector.C+",options.str().c_str());
                }
                else {ch->Process("~/DelilaSorting/DelilaSelector.C+",options.str().c_str(),numberofevents);};    
		};
             };
      };
}

