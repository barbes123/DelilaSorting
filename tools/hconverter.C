//#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TString.h"
#include "TRegexp.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h> 

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <sstream>
#include "TFile.h"

#include "TROOT.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h> 
#include <TSystem.h>
#include <vector>
#include <utility>
#include <ctime>
#include "TF1.h"
#include "TMath.h"
#include <string>
#include <TObjArray.h>
#include <TObjString.h>

/*#ifndef CalibrationCoefficient_h
#define CalibrationCoefficient_h
class CalibrationCoefficient {
public:
int test;
};
#endif
*/
// using namespace std;
// TList *FileList;
// TFile *Target;
// 


//const TString frootname = "run_selector_1521_1530_beta_0.015";
//const TString frootname = "run_selector_1530_beta_0.015";
//const TString frootname = "83Ga";
//const TString foldername = "";  
//const TString foldername = "Counts/_beta_";  
//const TString foldername = "Coinc";  
TString frootname;
TString foldername;
//TString foldername = "_beta_";
TString fout;

Int_t maxbin;
UInt_t rebin_const=0;
 
bool tkt = true;

void ConvertTH1F (TH1* h);
void ConvertTH2F (TH2* m);


void hconverter(TString file="toconvertfile.root", TString hname="mEliade_raw", Double_t maxbin_=16384, UInt_t rebin=0, TString foldername1="", TString foldername2="",int fortkt=1)

{   
  if (fortkt==0) tkt=false;
  TObjArray *toks;
  foldername = foldername1;
  maxbin = maxbin_;
  UInt_t nrun=0;
  rebin_const = rebin;
  std::stringstream InputFile;
  
   
  if (nrun!=0){InputFile<<"run_"<<nrun;}
  else{   
    toks = file.Tokenize(".");
    TString name = ((TObjString* )toks->At(0))->GetString();
    InputFile<< name;
    std::cout<<"----- name -----"<<  InputFile.str().c_str()<<std::endl;     
  }
  fout = InputFile.str().c_str();
  frootname = InputFile.str().c_str();
  fout.Append("/");  
  fout.Append(foldername);     
  std::cout<<"Out Path "  <<fout<<std::endl;
  
  struct stat st;
  if(stat(frootname,&st) == 0){std::cout<<"yes"<<std::endl;}
  else mkdir(frootname, 0777); 
  
  struct stat st1;
  if(stat(fout,&st1) == 0){std::cout<<"yes"<<std::endl;}
  else mkdir(fout, 0777);
  
  if (foldername2 != "")
  {
    
  fout.Append("/");  
  fout.Append(foldername2);
  foldername.Append("/");  
  foldername.Append(foldername2);     
  
  struct stat st2;
  if(stat(fout,&st1) == 0){std::cout<<"yes"<<std::endl;}
  else mkdir(fout, 0777);    
  };
  std::cout<<"Out Path "  <<fout<<std::endl;
//   if (nrun!=0)
    InputFile<<".root";
  TString fname = InputFile.str().c_str();
  std::cout<<"fname "  <<fname<<std::endl;
  TFile *_file0 = TFile::Open(fname);
   _file0->cd(foldername);
//   _file0->cd();
// _file0->cd(foldername);
   TDirectory *current_sourcedir = gDirectory;
   
   TChain *globChain = 0;
   TIter nextkey( current_sourcedir->GetListOfKeys() );
   TKey *key, *oldkey=0;
   
//    key = (TKey*)nextkey();   
   while ( (key = (TKey*)nextkey())) {
      TObject *obj = key->ReadObj();     
      if (!obj){continue;};
      if ( obj->IsA()->InheritsFrom( TH1F::Class() ) ) {
	std::cout<<"TH1::";//<<endl;
	TH1 *h1 = (TH1*)obj;
	std::cout<<h1->GetName()<<std::endl;
 	ConvertTH1F(h1);
	delete h1;
      }
      else
        if ( obj->IsA()->InheritsFrom( TH1D::Class() ) ) {
	std::cout<<"TD1::";//<<std::endl;
	TH1 *h1 = (TH1D*)obj;
	std::cout<<h1->GetName()<<std::endl;
 	ConvertTH1F(h1);
	delete h1;
      }
      else 
      if ( obj->IsA()->InheritsFrom( TH2F::Class() ) ) {
	TH2 *h2 = (TH2*)obj;
	std::cout<<"TH2";//<<std::endl;
	std::cout<<h2->GetName()<<std::endl;
		
	if ((h2->GetName() == hname )||(hname == "all"))
	{
		ConvertTH2F(h2); 
		std::cout<<" TH2 HistName "<<h2->GetName()<<std::endl;
	};
	
//	ConvertTH2F(h2);
	delete h2;
      }
      else {std::cout<<" type of hist in unknown \n"; continue;};     
   };
}

void ConvertTH1F (TH1* h)
{
  TString str;
  TString oname;  
  str = h->GetName();
  //For CoMPASS to get rid of stupid symbols
  std::cout<<str<<std::endl;
  str.Replace(str.First("@"),1,"_");
  str.Replace(str.First("/"),1,"_");
  std::cout<<str<<std::endl;
  
  
  
  if (rebin_const!=0){h->Rebin(rebin_const);}
  oname=frootname;
  oname.Append("/");
  oname.Append(foldername);
  oname.Append("/");
  oname.Append(str);  
  if (tkt) {oname.Append(".spe");}
  else {oname.Append(".dat");};
  
  
  ofstream fout(oname);  
  std::cout<<oname<<" n bins "<<h->GetNbinsX()<<std::endl;
  Double_t scale = 1;
   
//   if ((oname.Contains("_beta_"))||((foldername.Contains("_n_nim_")))) scale = 1e-6;
        
  
  if (tkt){
    for (int i=1; i<=h->GetNbinsX(); i++){
     Double_t value = h->GetBinContent(i);
    if (h->GetBinCenter(i)<=maxbin)
    {
    	fout<<value<<std::endl;
    	//std::cout<<value<<std::endl;
    };
    }
   }    
  else{
    for (int i=1; i<=h->GetNbinsX(); i++){
       Double_t value = h->GetBinContent(i);
    if (h->GetBinCenter(i)<=maxbin)
    {
     fout<<h->GetBinCenter(i)*scale<<" "<<value<< " "<< sqrt(value) <<std::endl;
//     std::cout<<h->GetBinCenter(i)*scale<<" "<<value<<std::std::endl;
    };
   }; 
  };      
  fout.close();
  std::cout<<" Completed!! \n";
}

void ConvertTH2F (TH2* m)
{
  TString oname;  
  TString str = m->GetName();
//  if (rebin_const!=0){m->Rebin(rebin_const);}

  int n_hists=m->GetXaxis()->GetNbins();
  for (int j=0;j<=n_hists;j++){
      TH1D *_py=m->ProjectionY("_py",j,j);
      if (_py->GetEntries() == 0) continue;
      
  oname=frootname;
  oname.Append("/");
  oname.Append(foldername);
  oname.Append("/");
  oname.Append(str);
  oname.Append("_py_");    
  int domain = m->GetXaxis()->GetBinCenter(j);
  oname.Append(Form("%i",domain));  
  //if (tkt) {oname.Append(".spe");}
  if (tkt) {oname.Append(".spe");}
  else {oname.Append(".dat");};      
  
  ofstream fout(oname);  
  std::cout<<oname<<" n bins "<<_py->GetNbinsX()<<std::endl;
  Double_t scale = 1;
  
  
  if (tkt){
    for (int i=1; i<=_py->GetNbinsX(); i++){
     Double_t value = _py->GetBinContent(i);
    if (_py->GetBinCenter(i)<=maxbin){fout<<value<<std::endl;};
    };
   }    
    else{
    for (int i=1; i<=_py->GetNbinsX(); i++){
       Double_t value = _py->GetBinContent(i);
       if (_py->GetBinCenter(i)<=maxbin){
//       fout<<_py->GetBinCenter(i)*scale<<" "<<value<< " "<< sqrt(value) <<std::endl;
       fout<<_py->GetBinCenter(i)*scale<<" "<<value <<std::endl;
//     std::cout<<h->GetBinCenter(i)*scale<<" "<<value<<std::std::endl;
        };
       }
      };
      fout.close();
  };         
}
