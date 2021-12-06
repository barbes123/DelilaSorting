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
//  if(mat_name.Contains("Eu")) channel_offset=64;
  Int_t nb_bin_mat  = matrix->GetYaxis()->GetNbins();
  Float_t range_min = matrix->GetYaxis()->GetXmin();
  Float_t range_max = matrix->GetYaxis()->GetXmax();
  TF1 *gaus = new TF1("gaus","gaus");
  TCanvas *c1 = new TCanvas("c1","c1");
  TH1D *proj_y =new TH1D("projection","proj",nb_bin_mat,range_min,range_max);
   for(Int_t j=0 ; j<=matrix->GetXaxis()->GetNbins() ; j++){
    matrix->ProjectionY(proj_y->GetName(),j,j);
    //if(proj_y->Integral()<100)continue;
    if (proj_y->GetEntries() == 0) continue;
    int max_bin = proj_y->GetMaximumBin();
    float max_value = proj_y->GetBinCenter(max_bin);
    c1->cd();
    proj_y->GetXaxis()->SetRangeUser(max_value-30,max_value+30);
    proj_y->Fit(gaus,"MQ");
    proj_y->Draw();
  //  c1->WaitPrimitive();
    std::cout << "   Offset for coincID " << j+channel_offset << "  " << gaus->GetParameter(1) << std::endl;;
    outputFile << j+channel_offset << "  " << gaus->GetParameter(1) <<"\n";;
    proj_y->Reset();
   }
    outputFile.close();
}
