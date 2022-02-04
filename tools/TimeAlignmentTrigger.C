#include <iostream>
#include <fstream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
void TimeAlignementTrigger(TH2 *matrix, float alignment_pos=0.)
{

  std::cout << "Welcome to the automatic TimeAlignement macro \n"
	    << "I will calculate by default the offset to put in the LookUpTable\n"
	    << "to put the coincidence peak in "<<alignment_pos <<" , except if you give a second parameter"
	    << std::endl;
  fstream outputFile;
  outputFile.open("TimeCalib.dat", ios_base::out);		    

  TString mat_name = matrix->GetName();
  Int_t nb_bin_mat  = matrix->GetYaxis()->GetNbins();
  Float_t range_min = matrix->GetYaxis()->GetXmin();
  Float_t range_max = matrix->GetYaxis()->GetXmax();
 // std::cout<<"range_min "<<range_min<<" range_max "<<range_max<<"\n";
  //TF1 *gaus = new TF1("gaus","gaus", -45e3, 45e3);//for mTimeCalib
 TF1 *gaus = new TF1("gaus","gaus", 80e3, 150e3);
  //    gaus->SetParLimits(1,100e3,300e3);
    //gaus->SetParLimits(0,30,100);
    //gaus->SetParLimits(2,4500,7500);
  TCanvas *c1 = new TCanvas("c1","c1");
  TH1D *proj_y =new TH1D("projection","proj",nb_bin_mat,range_min,range_max);

   for(int jj=0 ; jj<=matrix->GetXaxis()->GetNbins() ; jj++){
    matrix->ProjectionY(proj_y->GetName(),jj,jj);
    if (proj_y->GetEntries() == 0) continue;
    proj_y->GetXaxis()->SetRangeUser(50e3,250e3);
    int max_bin = proj_y->GetMaximumBin();
    float max_value = proj_y->GetXaxis()->GetBinCenter(max_bin);
    c1->cd();
//    gaus->SetParameter(1,alignment_pos);
    //gaus->SetParLimits(1,4.0e+03,5e+03);
//    proj_y->Fit(gaus,"MR");
    proj_y->Draw();

//    std::cout<<"max_value "<<max_value<<" max bin "<< proj_y-> GetBinCenter(max_bin-1) <<"\n";

//    std::cout << "   Offset for domain " << j+channel_offset << "  " << gaus->GetParameter(1) << std::endl;;
    std::cout << " domain " << jj-1<<" offset " << proj_y-> GetBinCenter(max_bin) - alignment_pos << " maxbin_pos " << proj_y-> GetBinCenter(max_bin) << std::endl;;
    //outputFile << j+channel_offset-1 << "  " << gaus->GetParameter(1) <<"\n";;
    //outputFile << j+channel_offset-1 << "  " << max_value<< " offset: "<< alignment_pos - max_value <<"\n";;
    outputFile << jj-1 << "  " << proj_y-> GetBinCenter(max_bin) - alignment_pos <<"\n";;
    
    //c1->WaitPrimitive();
    proj_y->Reset();
   }
    outputFile.close();
}
