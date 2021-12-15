#include <iostream>
#include <fstream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
void TimeAlignement(TH2 *matrix, float alignment_pos=0.)
{

  std::cout << "Welcome to the automatic TimeAlignement macro \n"
	    << "I will calculate by default the offset to put in the LookUpTable\n"
	    << "to put the coincidence peak in "<<alignment_pos <<" , except if you give a second parameter"
	    << std::endl;
  fstream outputFile;
  outputFile.open("TimeCalib.dat", ios_base::out);		    
  int channel_offset=0;
  TString mat_name = matrix->GetName();
  Int_t nb_bin_mat  = matrix->GetYaxis()->GetNbins();
  Float_t range_min = matrix->GetYaxis()->GetXmin();
  Float_t range_max = matrix->GetYaxis()->GetXmax();
 // std::cout<<"range_min "<<range_min<<" range_max "<<range_max<<"\n";
  TF1 *gaus = new TF1("gaus","gaus", -45e3, 45e3);
  TCanvas *c1 = new TCanvas("c1","c1");
  TH1D *proj_y =new TH1D("projection","proj",nb_bin_mat,range_min,range_max);
   for(Int_t j=0 ; j<=matrix->GetXaxis()->GetNbins() ; j++){
    matrix->ProjectionY(proj_y->GetName(),j,j);
    if (proj_y->GetEntries() == 0) continue;
    int max_bin = proj_y->GetMaximumBin();
    float max_value = proj_y->GetXaxis()->GetBinCenter(max_bin);
    c1->cd();
    proj_y->Fit(gaus,"MR");
    proj_y->Draw();
    gaus->SetParameter(1,max_value);
    std::cout<<"max_value "<<max_value<<"\n";
//     c1->WaitPrimitive();
    std::cout << "   Offset for coincID " << j+channel_offset << "  " << gaus->GetParameter(1) << std::endl;;
    outputFile << j+channel_offset-1 << "  " << gaus->GetParameter(1) <<"\n";;
    proj_y->Reset();
   }
    outputFile.close();
}
