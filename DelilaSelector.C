#define DelilaSelector_cxx
// The class definition in data.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after hTimeDiffCoreCoreBegin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("DelilaSelector.C")
// root> T->Process("DelilaSelector.C","some options")
// root> T->Process("DelilaSelector.C+")
//


#include "DelilaSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TObjString.h>
#include <unordered_set>
#include <iomanip>      // std::setwsorted
using namespace std;


////////////////////////////////Please, modify if needed////////////////////////////////////////////
bool blGammaGamma = true;
bool blGammaGammaCS = true;
bool blCS = true;
bool blIsTrigger = true; //trigger signal is present
bool blTriggerMode = false;//to reset the queues after each trigger - to be implemented
bool blOutTree = false;
// bool blTimeEnergy = false;

////////////////////////////////Please, DO NOT modify ////////////////////////////////////////////
int addBackMode = 0; //0 - no addback; 1- addback;//not in use for ELIFANT
// bool blPulserTimeAllignement = false; // for normal use should be false

bool debug = false;

ULong64_t trigger_cnt = 0;

// const int NumberOfClovers = 2;
const int max_domain = 100;
const int nbr_of_ch = 200;
// ULong64_t lastTimeStampTrigger = 0;

// const int zero_channel = 101; //for time allignement
// ULong64_t lastTime_pulser = 0;
// ULong64_t lastTime_dom0 = 0;
// ULong64_t lastTimeStamp = 0;

double beta = 0;

std::stringstream OutputFile;


void DelilaSelector::Read_ELIADE_LookUpTable() {
  std::cout << "I am Reading ELIADE LookUpTable ... ";
//  std::stringstream CUTFile;
//  CUTFile << CUTG_Directory.Data() << "cut_EeEw_galileo.root";
//  TFile *file_EeEw = TFile::Open(CUTFile.str().c_str());

 char* pLUT_Path;
  pLUT_Path = getenv ("DELILA_LUT");
  if (pLUT_Path!=NULL)
    printf ("The LookUpTable path is: %s \n",pLUT_Path);


  std::stringstream LUTFile;
  LUTFile << pLUT_Path <<"/"<<"LUT_DELILA.dat";
  //LUTFile << LUT_Directory << "LUT_ELIADE.dat";
  std::ifstream lookuptable(LUTFile.str().c_str());

  if (!lookuptable.good()) {
    std::ostringstream os;
    os << "Could not open " << LUTFile.str().c_str()
       << " and I need it ;(\n";
    Abort(os.str().c_str());
  } else {
    while (lookuptable.good()) {
      std::string oneline;
      std::getline(lookuptable, oneline);
      if (!lookuptable.good()) continue;
      if (oneline[0] == '#') continue; // ignore lines stating with #
      if (oneline.empty())   continue; // ignore empty lines
      TDelilaDetector curDet;
      Float_t theta(-1.), phi(-1.);
      int upperThreshold = 1e6;
      std::istringstream is(oneline);
      if (debug) std::cout << is.str().c_str() << std::endl;
//       is >> curDet.ch >> curDet.dom >> curDet.theta >> curDet.phi >> curDet.TimeOffset >> curDet.upperThreshold;
      is >> curDet.ch >> curDet.dom >> curDet.detType >> curDet.serial >> curDet.TimeOffset >> curDet.theta >> curDet.phi >> curDet.upperThreshold >> curDet.cs_dom;
    //  std::cout<<" curDfalseet.ch  "<<curDet.ch <<" curDet.TimeOffset " <<curDet.TimeOffset<<std::endl;
      
      if (curDet.ch >= 0) {
          curDet.theta *= TMath::DegToRad();
          curDet.phi *= TMath::DegToRad();
	//theta *= TMath::DegToRad();
	//phi *= TMath::DegToRad();
//	TVector3 DetPos;
//	curDet.direction.SetMagThetaPhi(210, theta, phi);
	int pol_order = 0;
	//Now we sorted_run_354.roottry to get the EeEw selection with a simple line
 	float offset_gate(0.),slope_gate(1.);
    is >> offset_gate;
// 	is >> offset_gate >> slope_gate;
	//curDet.rejectionEeEw = new TF1(Form("Ge_%2i_EeEw",curDet.domain),
	//			       "pol1");
	//curDet.rejectionEeEw->FixParameter(0,offset_gate);
	//curDet.rejectionEeEw->FixParameter(1,slope_gate);
	is >> pol_order;
	curDet.pol_order = pol_order;
	if (debug) std::cout << "Cal order " << pol_order << "  ";
	std::vector<float> DetCal_sub0;
	for (int k = 0; k < pol_order; k++) {
	  float par = -FLT_MAX;
	  is >> par;
	  if (par > -FLT_MAX) {
	    if (debug) std::cout << par << "  ";
	    curDet.calibE.push_back(par);
	  }
	}
	LUT_DELILA[curDet.ch] = curDet;
      }
    }
  }
  lookuptable.close();
  std::cout << " done" << std::endl;
  //  std::exit(1);
}


void DelilaSelector::Read_TimeAlignment_LookUpTable() {
  std::cout << "I am Reading TimeAlignment LookUpTable ... ";
 
  char* pLUT_Path;
  pLUT_Path = getenv ("DELILA_LUT");
  if (pLUT_Path!=NULL)
    printf ("The LookUpTable path is: %s \n",pLUT_Path);


  std::stringstream LUTFile;
  LUTFile << pLUT_Path <<"/"<<"LUT_TA.dat";
  //LUTFile << LUT_Directory << "LUT_ELIADE.dat";
  const int nbr_of_ch = 200;
  std::ifstream lookuptable(LUTFile.str().c_str());

  if (!lookuptable.good()) {
    std::ostringstream os;
    os << "Could not open " << LUTFile.str().c_str()
       << " and I need it ;(\n";
    Abort(os.str().c_str());
  } else {
    while (lookuptable.good()) {
      std::string oneline;
      std::getline(lookuptable, oneline);
      if (!lookuptable.good()) continue;
      if (oneline[0] == '#') continue; // ignore lines stating with #
      if (oneline.empty())   continue; // ignore empty lines

      std::istringstream is(oneline);
 
      int coinc_id = 0; int time_corr = 0;
      is >> coinc_id >> time_corr;
      LUT_TA[coinc_id] = time_corr;

 
  }
  lookuptable.close();
  }
  std::cout << " done" << std::endl;
  //  std::exit(1);
}

void DelilaSelector::Read_TimeAlignment_Trigger() {
  std::cout << "I am Reading Read_TimeAlignment_Trigger LookUpTable ... ";
 
  char* pLUT_Path;
  pLUT_Path = getenv ("DELILA_LUT");
  if (pLUT_Path!=NULL)
    printf ("The LookUpTable path is: %s \n",pLUT_Path);


  std::stringstream LUTFile;
  LUTFile << pLUT_Path <<"/"<<"LUT_TRIGGER.dat";
  //LUTFile << LUT_Directory << "LUT_ELIADE.dat";
  const int nbr_of_ch = 200;
  std::ifstream lookuptable(LUTFile.str().c_str());

  if (!lookuptable.good()) {
    std::ostringstream os;
    os << "Could not open " << LUTFile.str().c_str()
       << " and I need it ;( but i can survive without time alignment \n";
    //Abort(os.str().c_str());
  } else {
    while (lookuptable.good()) {
      std::string oneline;
      std::getline(lookuptable, oneline);
      if (!lookuptable.good()) continue;
      if (oneline[0] == '#') continue; // ignore lines stating with #
      if (oneline.empty())   continue; // ignore empty lines

      std::istringstream is(oneline);
 
      int det_id = 0; int time_corr = 0;
      is >> det_id >> time_corr;
      LUT_TA_TRG[det_id] = time_corr;

 
  }
  lookuptable.close();
  }
  std::cout << " done" << std::endl;
  //  std::exit(1);
}

void DelilaSelector::Read_CoincGates() {
  std::cout << "I am Reading Coic Gates LookUpTable ... ";
 
  char* pLUT_Path;
  pLUT_Path = getenv ("DELILA_LUT");
  if (pLUT_Path!=NULL)
    printf ("The LookUpTable path is: %s \n",pLUT_Path);


  std::stringstream LUTFile;
  LUTFile << pLUT_Path <<"/"<<"LUT_COINC.dat";
  std::ifstream lookuptable(LUTFile.str().c_str());

  if (!lookuptable.good()) {
    std::ostringstream os;
    os << "Could not open " << LUTFile.str().c_str()
       << " and I need it ;( but will can continue with default values\n";
//     Abort(os.str().c_str());
  } else {
    while (lookuptable.good()) {
      std::string oneline;
      std::getline(lookuptable, oneline);
      if (!lookuptable.good()) continue;
      if (oneline[0] == '#') continue; // ignore lines stating with #
      if (oneline.empty())   continue; // ignore empty lines

      std::istringstream is(oneline);
      TString coinc_name;
      int coinc_id = 0; int gate = 0;
      is >> coinc_name>> coinc_id >> gate;
      
      std::cout<<coinc_name<<" coin_id " << coinc_id <<" gate "<<gate <<" ps \n";
      LUT_COINC[coinc_id] = gate; 
  }
  lookuptable.close();
  }
  std::cout << " done" << std::endl;
}



void DelilaSelector::Print_ELIADE_LookUpTable()
{
    std::cout<<"Print_ELIADE_LookUpTable \n";		
    std::map<unsigned int, TDelilaDetector > ::iterator it__ = LUT_DELILA.begin();
    for (; it__ != LUT_DELILA.end(); ++it__) {
     // is >> curDet.ch >> curDet.dom >> theta >> phi >> curDet.TimeOffset >> curDet.upperThreshold;
	std::cout<<" Ch "<<LUT_DELILA[it__->first].ch<<" Dom "<< LUT_DELILA[it__->first].dom<<" "<< LUT_DELILA[it__->first].theta<<" "<< LUT_DELILA[it__->first].phi <<" offset "<< LUT_DELILA[it__->first].TimeOffset<<" Thr "<< LUT_DELILA[it__->first].upperThreshold<<" serial "<<LUT_DELILA[it__->first].serial<<" theta" <<LUT_DELILA[it__->first].theta<<" phi "<<LUT_DELILA[it__->first].phi <<" cs_dom: "<<LUT_DELILA[it__->first].cs_dom<<" pol_order: " <<LUT_DELILA[it__->first].pol_order <<std::endl;
    }
};

void DelilaSelector::Print_TimeAlignment_LookUpTable()
{
    std::cout<<"Print_TimeAlignment_LookUpTable \n";		
    std::map<int, int > ::iterator it__ = LUT_TA.begin();
    for (; it__ != LUT_TA.end(); ++it__) {
     // is >> curDet.ch >> curDet.dom >> theta >> phi >> curDet.TimeOffset >> curDet.upperThreshold;
	std::cout<<" coinc_id "<<it__->first<<" time_corr "<< it__->second<<std::endl;
    }
};



float DelilaSelector::CalibDet(float val, int ch)
{

 Float_t Randomize = 0.;
 std::vector<float> cal_par;
 cal_par = LUT_DELILA[ch].calibE;  
 float cal_value = 0.;
 float flou = (float)rand () / RAND_MAX - 0.5; 
  for (UInt_t j = 0; j < cal_par.size(); j++) {
    cal_value += cal_par[j] * TMath::Power(val + flou, (Int_t) j);
  }
  return cal_value;

};

void DelilaSelector::Begin(TTree * tree)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

  
  TString option = GetOption();
  toks = option.Tokenize(",");
  TString RunID = ((TObjString*) toks->At(0))->GetString();
  beta = ((TObjString*) toks->At(2))->GetString().Atof();
  addBackMode = atoi(((TObjString*) toks->At(3))->GetString());
  std::cout << "addBackMode  " << addBackMode <<std::endl;
  
  std::cout<<"Beta is "<<beta<<" % \n";

//   TString VolID = ((TObjString*) toks->At(1))->GetString();
// 
//   std::stringstream OutputFile;
//   OutputFile << "selectror_run" << "_" << RunID <<"_"<<VolID<< ".root";
//   std::cout << "OUTFILE  run" << "_" << RunID<<"_"<<VolID<< ".root"<<std::endl;
  
   if(!tree) {std::cout<<" TTree NOT found "<<std::endl; return;};
  
  std::cout<<" Delila Sorting "<<std::endl;
  std::cout<<" TTree found "<<std::endl;
  fChain = tree;
  fChain->SetMakeClass(1);  
  fChain->SetBranchAddress("ChargeLong", 	&DelilaEvent.fEnergy, 	       &b_energ);
  fChain->SetBranchAddress("ChargeShort", 	&DelilaEvent.fEnergyShort, 	   &b_energ_short);
  fChain->SetBranchAddress("FineTS", 	    &DelilaEvent.fTimeStampFS, 	   &b_tstmp_fine);
  fChain->SetBranchAddress("TimeStamp",     &DelilaEvent.fTimeStamp,       &b_tstmp); 
  fChain->SetBranchAddress("Ch", 	        &DelilaEvent.fChannel,         &b_channel);
  fChain->SetBranchAddress("Mod", 	        &DelilaEvent.fMod,             &b_mod);
    
  auto index = new TTreeIndex(fChain,"TimeStamp", "Ch");
  fChain->SetTreeIndex(index);
  ULong64_t lastStamp = 0;
  ULong64_t timeStamp = 0;
  const auto nEvents = fChain->GetEntries();

  for(auto iEve = 0; iEve < nEvents; iEve++) {
    auto local = fChain->GetEntryNumberWithIndex(iEve);
    fChain->GetEntry(local);

    if(lastStamp > timeStamp){
      std::cout << "hit: " << lastStamp <<"\t"<< timeStamp << std::endl;
    }
    //std::cout << "hit: " << timeStamp  << std::endl;
    lastStamp = timeStamp;
  }

  
  Read_ELIADE_LookUpTable();
  Print_ELIADE_LookUpTable();
  
  Read_TimeAlignment_LookUpTable();
  Print_TimeAlignment_LookUpTable();
  
//   Read_TimeAlignment_Trigger();
//   Read_CoincGates();

  
   
}




void DelilaSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

 nevents = 0;
 nevents_reset=0;
 reset_counter = 0;

   hTimeSort = new TH1F("hTimeSort", "hTimeSort", 1e3, -1e5,1e5);
   fOutput->Add(hTimeSort);
    
   hChannelHit = new TH1F("hChannelHit", "hChannelHit",3216,0,3216);
   fOutput->Add(hChannelHit);
   
   hSegmentHit = new TH1F("hSegmentHit", "hSegmentHit",20,0,20);
   fOutput->Add(hSegmentHit);
   
   hDomainHit = new TH1F("hDomainHit", "hDomainHit",max_domain,0,max_domain);
   fOutput->Add(hDomainHit);
   
   hDetTypeHit = new TH1F("hDetTypeHit", "hDetTypeHit",20,0,20);
   fOutput->Add(hDetTypeHit);
   
   hMult_gg = new TH1F("hMult_gg", "hMult_gg",20,0,20);
   fOutput->Add(hMult_gg);
   
   hMult_gg_CS = new TH1F("hMult_gg_CS", "hMult_gg_CS",20,0,20);
   fOutput->Add(hMult_gg_CS);
   
/*        
   hDelila = new TH1F("hDelila", "hDelila", 4096, -0.5, 16383.5);
   hDelila->GetYaxis()->SetTitle("counts");
   hDelila->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelila);
   
   hDelilaDC = new TH1F("hDelilaDC", "hDelilaDC", 4096, -0.5, 16383.5);
   hDelilaDC->GetYaxis()->SetTitle("counts");
   hDelilaDC->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelilaDC);
   
   hDelilaCS = new TH1F("hDelilaCS", "hDelilaCS", 4096, -0.5, 16383.5);
   hDelilaCS->GetYaxis()->SetTitle("counts");
   hDelilaCS->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelilaCS);
   
   hDelilaCS_DC = new TH1F("hDelilaCS_DC", "hDelilaCS_DC", 4096, -0.5, 16383.5);
   hDelilaCS_DC->GetYaxis()->SetTitle("counts");
   hDelilaCS_DC->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelilaCS_DC);*/
   
   hTimeDiffPulser = new TH1F("hTimeDiffPulser", "hTimeDiffPulser", 1000, -99.5, 899.5);
   fOutput->Add(hTimeDiffPulser);

   mDelila_raw = new TH2F("mDelila_raw", "mDelila_raw", max_domain, -0.5, max_domain-0.5, 16384, -0.5, 16383.5);
   mDelila_raw->GetXaxis()->SetTitle("domain");
   mDelila_raw->GetYaxis()->SetTitle("ADC channels");   
   fOutput->Add(mDelila_raw);
   
   mDelila = new TH2F("mDelila", "mDelila", max_domain, 0, max_domain, 16384, -0.5, 16383.5);
   mDelila->GetXaxis()->SetTitle("domain");
   mDelila->GetYaxis()->SetTitle("keV");
   fOutput->Add(mDelila);
   
   mDelilaDC = new TH2F("mDelilaDC", "mDelilaDC", max_domain, 0, max_domain, 16384, -0.5, 16383.5);
   mDelilaDC->GetXaxis()->SetTitle("domain");
   mDelilaDC->GetYaxis()->SetTitle("keV");
   fOutput->Add(mDelilaDC);   
   
   mDelilaCS = new TH2F("mDelilaCS", "mDelilaCS", max_domain, 0, max_domain, 16384, -0.5, 16383.5);
   mDelilaCS->GetXaxis()->SetTitle("domain");
   mDelilaCS->GetYaxis()->SetTitle("keV");
   fOutput->Add(mDelilaCS);
   
   mDelilaCS_DC = new TH2F("mDelilaCS_DC", "mDelilaCS_DC", max_domain, 0, max_domain, 16384, -0.5, 16383.5);
   mDelilaCS_DC->GetXaxis()->SetTitle("domain");
   mDelilaCS_DC->GetYaxis()->SetTitle("keV");
   fOutput->Add(mDelilaCS_DC);
   
   ////////////////
   mDelila_long = new TH2F("mDelila_long", "mDelila_long", max_domain, 0, max_domain, 4096, -0.5,  65535.5);
   mDelila_long->GetXaxis()->SetTitle("domain");
   mDelila_long->GetYaxis()->SetTitle("keV");
   fOutput->Add(mDelila_long);
   
   mDelilaDC_long = new TH2F("mDelilaDC_long", "mDelilaDC_long", max_domain, 0, max_domain, 4096, -0.5,  65535.5);
   mDelilaDC_long->GetXaxis()->SetTitle("domain");
   mDelilaDC_long->GetYaxis()->SetTitle("keV");
   fOutput->Add(mDelilaDC_long);   
   
   mDelilaCS_long = new TH2F("mDelilaCS_long", "mDelilaCS_long", max_domain, 0, max_domain, 4096, -0.5,  65535.5);
   mDelilaCS_long->GetXaxis()->SetTitle("domain");
   mDelilaCS_long->GetYaxis()->SetTitle("keV");
   fOutput->Add(mDelilaCS_long);
   
   mDelilaCS_DC_long = new TH2F("mDelilaCS_DC_long", "mDelilaCS_DC_long", max_domain, 0, max_domain, 4096, -0.5,  65535.5);
   mDelilaCS_DC_long->GetXaxis()->SetTitle("domain");
   mDelilaCS_DC_long->GetYaxis()->SetTitle("keV");
   fOutput->Add(mDelilaCS_DC_long);
   
   ///////////////////////
   
   
   
   
   
   
   mThetaPhi = new TH2F("mThetaPhi", "mThetaPhi", 90,60,150,360,0,360);
   mThetaPhi->GetXaxis()->SetTitle("theta, degrees");
   mThetaPhi->GetYaxis()->SetTitle("phi, degrees");
   fOutput->Add(mThetaPhi);
      
   mGammaGamma = new TH2F("mGammaGamma", "mGammaGamma", 4096, -0.5, 16383.5, 4096, -0.5, 16383.5);
   mGammaGamma->GetXaxis()->SetTitle("keV");
   mGammaGamma->GetYaxis()->SetTitle("keV");
   fOutput->Add(mGammaGamma);
   
   mGammaGammaCS = new TH2F("mGammaGammaCS", "mGammaGammaCS", 4096, -0.5, 16383.5, 4096, -0.5, 16383.5);
   mGammaGammaCS->GetXaxis()->SetTitle("keV");
   mGammaGammaCS->GetYaxis()->SetTitle("keV");
   fOutput->Add(mGammaGammaCS);
   
   mGammaGammaDC = new TH2F("mGammaGammaDC", "mGammaGammaDC", 4096, -0.5, 16383.5, 4096, -0.5, 16383.5);
   mGammaGammaDC->GetXaxis()->SetTitle("keV");
   mGammaGammaDC->GetYaxis()->SetTitle("keV");
   fOutput->Add(mGammaGammaDC);
   
   mGammaGammaCS_DC = new TH2F("mGammaGammaCS_DC", "mGammaGammaCS_DC", 4096, -0.5, 16383.5, 4096, -0.5, 16383.5);
   mGammaGammaCS_DC->GetXaxis()->SetTitle("keV");
   mGammaGammaCS_DC->GetYaxis()->SetTitle("keV");
   fOutput->Add(mGammaGammaCS_DC);
   
   coinc_gates[11]=600000;//hpge-hpge 40 ns
   coinc_gates[15]=100000;//hpge-bgo 40 ns
//    coinc_gates[51]=10000;//bgo-hpge 40 ns
   
   coinc_gates[33]=10000;//labr-labr 40 ns
   coinc_gates[35]=10000;//labr-bgo 40 ns
//    coinc_gates[53]=10000;//hpge-bgo 40 ns
   
   coinc_gates[13]=1000000;//hpge-bgo 40 ns
//    coinc_gates[31]=10000;//hpge-bgo 40 ns
   
   gg_coinc_id[11]="mgg_hpge_hpge";
   gg_coinc_id[13]="mgg_labr_hpge";
   gg_coinc_id[33]="mgg_labr_labr";
   
   detector_name[1]="HPGe";
   detector_name[3]="LabBr";

   std::map<UInt_t,std::string>::iterator itna =  gg_coinc_id.begin();

  for(;itna!=gg_coinc_id.end();++itna){
      
  mGG[itna->first] = new TH2F(Form("%s",itna->second.c_str()), Form("%s",itna->second.c_str()), 4096, -0.5, 16383.5, 4096, -0.5, 16383.5);
  mGG[itna->first]->GetXaxis()->SetTitle("keV");
  mGG[itna->first]->GetYaxis()->SetTitle("keV");
  fOutput->Add(mGG[itna->first]);
  
  mGG_CS[itna->first] = new TH2F(Form("%s_CS",itna->second.c_str()), Form("%s_CS",itna->second.c_str()), 4096, -0.5, 16383.5, 4096, -0.5, 16383.5);
  mGG_CS[itna->first]->GetXaxis()->SetTitle("keV");
  mGG_CS[itna->first]->GetYaxis()->SetTitle("keV");
  fOutput->Add(mGG_CS[itna->first]);
    
  mGG_DC[itna->first] = new TH2F(Form("%s_DC",itna->second.c_str()), Form("%s_DC",itna->second.c_str()), 4096, -0.5, 16383.5, 4096, -0.5, 16383.5);
  mGG_DC[itna->first]->GetXaxis()->SetTitle("keV");
  mGG_DC[itna->first]->GetYaxis()->SetTitle("keV");
  fOutput->Add(mGG_DC[itna->first]);
   
  mGG_CS_DC[itna->first] = new TH2F(Form("%s_CS_DC",itna->second.c_str()), Form("%s_CS_DC",itna->second.c_str()), 4096, -0.5, 16383.5, 4096, -0.5, 16383.5);
  mGG_CS_DC[itna->first]->GetXaxis()->SetTitle("keV");
  mGG_CS_DC[itna->first]->GetYaxis()->SetTitle("keV");
  fOutput->Add(mGG_CS_DC[itna->first]);
  
  //////////////////////////////
  mGG_long[itna->first] = new TH2F(Form("%s_long",itna->second.c_str()), Form("%s_long",itna->second.c_str()), 4096, -0.5,  65535.5, 4096, -0.5,  65535.5);
  mGG_long[itna->first]->GetXaxis()->SetTitle("keV");
  mGG_long[itna->first]->GetYaxis()->SetTitle("keV");
  fOutput->Add(mGG_long[itna->first]);
  
  mGG_CS_long[itna->first] = new TH2F(Form("%s_CS_long",itna->second.c_str()), Form("%s_CS_long",itna->second.c_str()), 4096, -0.5,  65535.5, 4096, -0.5,  65535.5);
  mGG_CS_long[itna->first]->GetXaxis()->SetTitle("keV");
  mGG_CS_long[itna->first]->GetYaxis()->SetTitle("keV");
  fOutput->Add(mGG_CS_long[itna->first]);
    
  mGG_DC_long[itna->first] = new TH2F(Form("%s_DC_long",itna->second.c_str()), Form("%s_DC_long",itna->second.c_str()), 4096, -0.5,  65535.5, 4096, -0.5,  65535.5);
  mGG_DC_long[itna->first]->GetXaxis()->SetTitle("keV");
  mGG_DC_long[itna->first]->GetYaxis()->SetTitle("keV");
  fOutput->Add(mGG_DC_long[itna->first]);
   
  mGG_CS_DC_long[itna->first] = new TH2F(Form("%s_CS_DC_long",itna->second.c_str()), Form("%s_CS_DC_long",itna->second.c_str()), 4096, -0.5,  65535.5, 4096, -0.5,  65535.5);
  mGG_CS_DC_long[itna->first]->GetXaxis()->SetTitle("keV");
  mGG_CS_DC_long[itna->first]->GetYaxis()->SetTitle("keV");
  fOutput->Add(mGG_CS_DC_long[itna->first]);
  
  ///////////////////////////////
  
  
   if (itna->first == 13){
       mGG[itna->first]->GetXaxis()->SetTitle("LaBr, keV"); mGG[itna->first]->GetYaxis()->SetTitle("HPGe, keV");
       mGG_CS[itna->first]->GetXaxis()->SetTitle("LaBr, keV"); mGG_CS[itna->first]->GetYaxis()->SetTitle("HPGe, keV");
       mGG_DC[itna->first]->GetXaxis()->SetTitle("LaBr, keV"); mGG_DC[itna->first]->GetYaxis()->SetTitle("HPGe, keV");
       mGG_CS_DC[itna->first]->GetXaxis()->SetTitle("LaBr, keV"); mGG_CS_DC[itna->first]->GetYaxis()->SetTitle("HPGe, keV");

       mGG_long[itna->first]->GetXaxis()->SetTitle("LaBr, keV"); mGG_long[itna->first]->GetYaxis()->SetTitle("HPGe, keV");
       mGG_CS_long[itna->first]->GetXaxis()->SetTitle("LaBr, keV"); mGG_CS_long[itna->first]->GetYaxis()->SetTitle("HPGe, keV");
       mGG_DC_long[itna->first]->GetXaxis()->SetTitle("LaBr, keV"); mGG_DC_long[itna->first]->GetYaxis()->SetTitle("HPGe, keV");
       mGG_CS_DC_long[itna->first]->GetXaxis()->SetTitle("LaBr, keV"); mGG_CS_DC_long[itna->first]->GetYaxis()->SetTitle("HPGe, keV");              
   };
   
   if (itna->first == 11){
//     mGG_time_diff[itna->first] = new TH2F(Form("%s_time_diff",itna->second.c_str()), Form("%s_time_diff",itna->second.c_str()), 100, 0, 100, 10e3, -2e6, 2e6);//was tuned like that
       mGG_time_diff[itna->first] = new TH2F(Form("%s_time_diff",itna->second.c_str()), Form("%s_time_diff",itna->second.c_str()), 100, 0, 100, 4e4, -2e6, 2e6);
//     mGG_time_diff[itna->first]->GetXaxis()->SetTitle("domain"); mGG_DC[itna->first]->GetYaxis()->SetTitle("ps");
//     fOutput->Add(mGG_time_diff[itna->first]);   

   }else{
     mGG_time_diff[itna->first] = new TH2F(Form("%s_time_diff",itna->second.c_str()), Form("%s_time_diff",itna->second.c_str()), 100, 0, 100, 4e4, -2e6, 2e6);
   };
     mGG_time_diff[itna->first]->GetXaxis()->SetTitle("domain"); mGG_DC[itna->first]->GetYaxis()->SetTitle("ps");
     fOutput->Add(mGG_time_diff[itna->first]);
//    };
//   std::cout<<itna->first<<" "<< Form("%s",itna->second.c_str())<<" Initialized \n" ;
};
   
    
   
  std::map<UInt_t,std::string>::iterator itna1 =  detector_name.begin();

  for(;itna1!=detector_name.end();++itna1){
   hDelila[itna1->first] = new TH1F(Form("%s",itna1->second.c_str()), Form("%s",itna1->second.c_str()), 4096, -0.5, 16383.5);
   hDelila[itna1->first]->GetYaxis()->SetTitle("counts");
   hDelila[itna1->first]->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelila[itna1->first]);
   
   hDelilaCS[itna1->first] = new TH1F(Form("%s_CS",itna1->second.c_str()), Form("%s_CS",itna1->second.c_str()), 4096, -0.5, 16383.5);
   hDelilaCS[itna1->first]->GetYaxis()->SetTitle("counts");
   hDelilaCS[itna1->first]->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelilaCS[itna1->first]);
   
   hDelilaDC[itna1->first] = new TH1F(Form("%s_DC",itna1->second.c_str()), Form("%s_DC",itna1->second.c_str()), 4096, -0.5, 16383.5);
   hDelilaDC[itna1->first]->GetYaxis()->SetTitle("counts");
   hDelilaDC[itna1->first]->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelilaDC[itna1->first]);
   
   hDelilaCS_DC[itna1->first] = new TH1F(Form("%s_CS_DC",itna1->second.c_str()), Form("%s_CS_DC",itna1->second.c_str()), 4096, -0.5, 16383.5);
   hDelilaCS_DC[itna1->first]->GetYaxis()->SetTitle("counts");
   hDelilaCS_DC[itna1->first]->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelilaCS_DC[itna1->first]);
      
  };
  
  
  itna1 =  detector_name.begin();

  for(;itna1!=detector_name.end();++itna1){
   hDelila_long[itna1->first] = new TH1F(Form("%s_long",itna1->second.c_str()), Form("%s_long",itna1->second.c_str()), 4096, -0.5, 65535.5);
   hDelila_long[itna1->first]->GetYaxis()->SetTitle("counts");
   hDelila_long[itna1->first]->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelila_long[itna1->first]);
   
   hDelilaCS_long[itna1->first] = new TH1F(Form("%s_CS_long",itna1->second.c_str()), Form("%s_CS_long",itna1->second.c_str()), 4096, -0.5, 65535.5);
   hDelilaCS_long[itna1->first]->GetYaxis()->SetTitle("counts");
   hDelilaCS_long[itna1->first]->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelilaCS_long[itna1->first]);
   
   hDelilaDC_long[itna1->first] = new TH1F(Form("%s_DC_long",itna1->second.c_str()), Form("%s_DC_long",itna1->second.c_str()), 4096, -0.5, 65535.5);
   hDelilaDC_long[itna1->first]->GetYaxis()->SetTitle("counts");
   hDelilaDC_long[itna1->first]->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelilaDC_long[itna1->first]);
   
   hDelilaCS_DC_long[itna1->first] = new TH1F(Form("%s_CS_DC_long",itna1->second.c_str()), Form("%s_CS_DC_long",itna1->second.c_str()), 4096, -0.5, 65535.5);
   hDelilaCS_DC_long[itna1->first]->GetYaxis()->SetTitle("counts");
   hDelilaCS_DC_long[itna1->first]->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelilaCS_DC_long[itna1->first]);
      
  };
   
   
   
   mEnergyTimeDiff_trigger = new TH2F("mEnergyTimeDiff_trigger", "mEnergyTimeDiff_trigger", 16384, -0.5, 16383.5, 10e3, 0, 10e6);
   mEnergyTimeDiff_trigger->GetXaxis()->SetTitle("Energy, keV");
   mEnergyTimeDiff_trigger->GetYaxis()->SetTitle("Time diff, ps");
   fOutput->Add(mEnergyTimeDiff_trigger);
   
   mDomainTimeDiff_trigger = new TH2F("mDomainTimeDiff_trigger", "mDomainTimeDiff_trigger", max_domain, 0, max_domain, 10e3, 0, 10e6);
//    mDomainTimeDiff_trigger = new TH2F("mDomainTimeDiff_trigger", "mDomainTimeDiff_trigger", max_domain, 0, max_domain, 3e2, 0, 3e5);
   mDomainTimeDiff_trigger->GetXaxis()->SetTitle("domain");
   mDomainTimeDiff_trigger->GetYaxis()->SetTitle("Time diff, ps");
   fOutput->Add(mDomainTimeDiff_trigger);
   
   hTriggerTrigger = new TH1F("hTriggerTrigger", "hTriggerTrigger", 10e4, 0, 10e6);
   hTriggerTrigger->GetYaxis()->SetTitle("Counts");
   hTriggerTrigger->GetXaxis()->SetTitle("Time diff, ps");
   hTriggerTrigger->SetTitle("Time between two trigger signals");
   fOutput->Add(hTriggerTrigger);
   
   mTimeCalib = new TH2F("mTimeCalib", "mTimeCalib", 10000, 0, 10000, 2e3, -1e6, 1e6);
   mTimeCalib->GetXaxis()->SetTitle("coinc ID");
   mTimeCalib->GetYaxis()->SetTitle("ps");
   mTimeCalib->SetTitle("Sci time diff");
   fOutput->Add(mTimeCalib);
   
   mTimeCalibDomain0 = new TH2F("mTimeCalibDomain0", "mTimeCalibDomain0", 100, 0, 100, 4e4, 0, 2e12);
   mTimeCalibDomain0->GetXaxis()->SetTitle("coinc ID");
   mTimeCalibDomain0->GetYaxis()->SetTitle("ps");
   mTimeCalibDomain0->SetTitle("Sci time diff");
   fOutput->Add(mTimeCalibDomain0);
   
   mLaBr_BGO_time_diff = new TH2F("mLaBr_BGO_time_diff", "mLaBr_BGO_time_diff", 100, 0, 100, 4e4, -2e6, 2e6);
   mLaBr_BGO_time_diff->GetXaxis()->SetTitle("domin");
   mLaBr_BGO_time_diff->GetYaxis()->SetTitle("ps");
   mLaBr_BGO_time_diff->SetTitle("LaBr_BGO time diff");
   fOutput->Add(mLaBr_BGO_time_diff);
   
//    mLaBr_LabBr_time_diff = new TH2F("mLaBr_LabBr_time_diff", "mLaBr_LabBr_time_diff", 100, 0, 100, 4e3, -2e5, 2e5);
   mLaBr_LabBr_time_diff = new TH2F("mLaBr_LabBr_time_diff", "mLaBr_LabBr_time_diff", 100, 0, 100, 4e3, -4e4, 4e4);
   mLaBr_LabBr_time_diff->GetXaxis()->SetTitle("domain");
   mLaBr_LabBr_time_diff->GetYaxis()->SetTitle("ps");
   mLaBr_LabBr_time_diff->SetTitle("LaBr_LabBr time diff");
   fOutput->Add(mLaBr_LabBr_time_diff);

   mPulser0TimeDiff = new TH2F("mPulser0TimeDiff", "mPulser0TimeDiff", 100, 0.5, 100.5, 6e3, -3e6, 3e6);
   mPulser0TimeDiff->GetXaxis()->SetTitle("domain");
   mPulser0TimeDiff->GetYaxis()->SetTitle("ps");
   mPulser0TimeDiff->SetTitle("PulsePulser time diff");
   fOutput->Add(mPulser0TimeDiff);//time_diff relevant to the 1st channel (101), i.e. ch 101 is a trigger
   
   mTimeDiff_gg = new TH2F("mTimeDiff_gg", "mTimeDiff_gg", 200, 0.5, 200.5, 500, -99.5, 399.5);
   fOutput->Add(mTimeDiff_gg);
   
   mTimeDiff_gg_CS = new TH2F("mTimeDiff_gg_CS", "mTimeDiff_gg_CS", 200, 0.5, 200.5, 500, -99.5, 399.5);
   fOutput->Add(mTimeDiff_gg_CS);
       
   mPulserPulser = new TH2F("mPulserPulser", "mPulserPulser",4096, -0.5, 8191.5, 4096, -0.5, 8195.5);
   fOutput->Add(mPulserPulser);
    
     for (int j = 0;j<=50;j++){
          mEnergy_time_diff[j] = new TH2F(Form("mEnergy_time_diff_dom%i",j), Form("mEnergy_time_diff_dom%i",j), 4096, -0.5, 65535.5, 10e3, 0, 10e6);
          mEnergy_time_diff[j] ->GetXaxis()->SetTitle("16 keV / bin");
          mEnergy_time_diff[j] ->GetYaxis()->SetTitle("1000 ps / bin");
          mEnergy_time_diff[j]->SetTitle(" Time Diff Event.Time - Trigger.Time ");
           fOutput->Add(mEnergy_time_diff[j]); 
     };
   
  TString option = GetOption();
  toks = option.Tokenize(",");
  TString RunID = ((TObjString*) toks->At(0))->GetString();
  TString VolID = ((TObjString*) toks->At(1))->GetString();
  TString ServerID = ((TObjString*) toks->At(3))->GetString();
    
  //std::stringstream OutputFile;
  OutputFile.str(std::string());
  OutputFile << "selected_run" << "_" << RunID <<"_"<<VolID;
  if (atoi(ServerID) != 0) {OutputFile<<"_eliadeS"<<ServerID;};
  OutputFile << ".root";
//   std::cout <<"ServerID "<<ServerID<<" "<< OutputFile.str().c_str() <<std::endl;
   lastTime = 0;
   
  foutFile = new TFile (OutputFile.str().c_str(),"recreate"); 
   
   std::cout <<"ServerID "<<ServerID<<" "<< OutputFile.str().c_str() <<std::endl;
    start = std::clock();
    
    lastDelilaEvent.fTimeStamp = 0;
//     lastEliadeZeroEvent.fTimeStamp = 0;
    LastTriggerEvent.Time = 0;
    PulserEvent.Time = 0;
    DomainZeroEvent.Time = 0;

    
    outputQu.clear();
    
    /*std::map<UInt_t, TDelilaEvent> ::iterator it__ = LUT_ELIADE.begin();
    for (; it__ != LUT_ELIADE.end(); ++it__) {
    //      particle_id[it__->second] = it__->first;
	std::cout<<LUT_ELIADE[it__->first].ch<<" "<< LUT_ELIADE[it__->first].pol_order <<std::endl;
    }*/

}

Bool_t DelilaSelector::Process(Long64_t entry)
{

  
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  // fReader.SetLocalEntry(entry);
    GetEntry(entry);
    
//     std::cout<<"Warning: .HEREEEEE \n";
    
//   if(nb_entries == 0) {
    nb_entries = GetEntries();
    
    if (DelilaEvent.fTimeStamp == 0) {hTimeZero->Fill(DelilaEvent.fChannel);};
    
	int mod = DelilaEvent.fMod;
	int ch = DelilaEvent.fChannel;
    int daq_ch = (mod)*100+ch;
 	hChannelHit->Fill(daq_ch);
    
    //Check that daq_ch is defined in LUT
    std::map<unsigned int, TDelilaDetector >::iterator it = LUT_DELILA.find(daq_ch);
    if(it == LUT_DELILA.end()){return kTRUE;};
    
    DelilaEvent.channel = daq_ch;
    DelilaEvent.det_def = LUT_DELILA[daq_ch].detType;
    DelilaEvent.domain = LUT_DELILA[daq_ch].dom;   
    int domain = DelilaEvent.domain;
    
    mDelila_raw->Fill(domain,DelilaEvent.fEnergy);
    hDomainHit->Fill(domain);


//Check if the tree is time sorted
//     if (lastDelilaEvent.fTimeStamp > DelilaEvent.fTimeStamp){std::cout<<"Warning: .fTimeStamp TTree may be not sorted by time \n";return kTRUE;};
    
//      DelilaEvent.fTimeStamp = DelilaEvent.fTimeStamp + LUT_DELILA[daq_ch].TimeOffset;
//      DelilaEvent.Time = 1000*DelilaEvent.fTimeStamp + DelilaEvent.fTimeStampFS;
//         DelilaEvent.Time = DelilaEvent.fTimeStamp*1.0;
     DelilaEvent.Time = DelilaEvent.fTimeStampFS;
//      if (DelilaEvent.det_def == 1)  DelilaEvent.Time= DelilaEvent.Time / 4;
//      if (DelilaEvent.det_def == 9)  DelilaEvent.Time= DelilaEvent.Time / 4;
     
//      double_t time_diff_trigger = DelilaEvent.Time - LastTriggerEvent.Time;
     double_t time_diff_trigger = DelilaEvent.Time+LUT_TA_TRG[DelilaEvent.domain] - LastTriggerEvent.Time;//+LUT_TA_TRG[LastTriggerEvent.domain];
     
     if (DelilaEvent.det_def == 99){
        hTriggerTrigger->Fill(time_diff_trigger);
        LastTriggerEvent = DelilaEvent;
        trigger_cnt++;
        return kTRUE;
        
    };


    double time_diff_last = DelilaEvent.Time - lastDelilaEvent.Time;
//      if (lastDelilaEvent.Time > DelilaEvent.Time){std::cout<<"Warning: .Time  TTree may be not sorted by time \n";};
      if (time_diff_last<0){std::cout<<"Warning: .Time  TTree may be not sorted by time"<< time_diff_last<<" \n";};
     hTimeSort->Fill(DelilaEvent.Time - lastDelilaEvent.Time);
    
    
    //Check Pulser TimeAlignment
    if (DelilaEvent.det_def == 9){
        CheckPulserAllignement(90);
         return kTRUE;
    };
//    if (blPulserTimeAllignement)       return kTRUE;
    DelilaEvent.trg = trigger_cnt;
    
//     if (DelilaEvent.det_def == 9){return kTRUE;};
//     if ((DelilaEvent.fEnergy < LUT_DELILA[daq_ch].upperThreshold))  {/*std::cout<<DelilaEvent.fEnergy<< " "<< LUT_ELIADE[daq_ch].upperThreshold<<std::endl; */return kTRUE;}
    
    
    
	DelilaEvent.EnergyCal = CalibDet(DelilaEvent.fEnergy, daq_ch);
    
    double costheta = TMath::Cos(LUT_DELILA[daq_ch].theta);
    DelilaEvent.EnergyDC = DelilaEvent.EnergyCal*(1./sqrt(1 - beta*beta) * (1 - beta*costheta));
	DelilaEvent.cs_domain = LUT_DELILA[daq_ch].cs_dom;
    
    DelilaEvent.theta= LUT_DELILA[daq_ch].theta;
    DelilaEvent.phi= LUT_DELILA[daq_ch].phi;
        
    hDetTypeHit->Fill(DelilaEvent.det_def);
    mThetaPhi->Fill(DelilaEvent.theta, DelilaEvent.phi);
	
//     if (blTimeEnergy){mapTimeEnergy[domain]->Fill(DelilaEvent.Time, DelilaEvent.EnergyCal);}
  	
	
	mDelila->Fill(domain,DelilaEvent.EnergyCal);
    mDelila_long->Fill(domain,DelilaEvent.EnergyCal);
    mDelilaDC->Fill(domain,DelilaEvent.EnergyDC);
    mDelilaDC_long->Fill(domain,DelilaEvent.EnergyDC);
    
    if (((DelilaEvent.det_def == 3)||(DelilaEvent.det_def == 1))&&(blIsTrigger)){
//     if ((DelilaEvent.det_def == 3)&&(blIsTrigger)){
        mDomainTimeDiff_trigger->Fill(domain,time_diff_trigger);
//         mEnergyTimeDiff_trigger->Fill(DelilaEvent.EnergyCal,time_diff_trigger);
        mEnergy_time_diff[domain]->Fill(DelilaEvent.EnergyCal,time_diff_trigger);
    };
    
        
//        time_alignment();   
//     time_alignment_pulser(90);
       time_alignment_symmetric();   //use this one
//     time_alignment_dom0(0);

     //Doing Gamma-Gamma matrix     
     if (blGammaGamma)
      {
     if (DelilaEvent.det_def == 3) gamma_gamma();
        };

    if (blCS){
        cs();
        };
        
//     if ((DelilaEvent.det_def ==1 )||(DelilaEvent.det_def ==3 )) gamma_gamma_LaBr_HPGe();
//         if (DelilaEvent.det_def ==1) gamma_gamma_LaBr_HPGe();   
   
  if ((entry) % int(nb_entries / 100) == 0 || (entry) % 100000 == 0) {
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    double eventsPerSec = entry / duration;
    double eta = (nb_entries - entry) / eventsPerSec;
    std::cout << "                      \r" << entry << " / " << nb_entries
	      << " ====> " << round((float) entry / nb_entries * 100.)
	      << " % " << " (" << eventsPerSec << " ev/s, " << round(eta / 60)
	      << ":" << std::setw(2) << std::setfill('0')
	      << round(((int) eta) % 60) << std::setw(8) << " min ETA)";
    std::cout.flush();
  }
   
   lastDelilaEvent = DelilaEvent;
   nevents++;
   nevents_reset++;
   return kTRUE;
}

void DelilaSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void DelilaSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the clieint coreID, nt, it can be used to present
   // the results graphically or save the results to file.
    
  std::cout<<"I  will terminate soon... "<<std::endl;  

  TIter iter(fOutput);
      

  std::cout << std::endl << "Finished after " << round(duration / 60) << ":"
	    << std::setw(2) << round(((int) duration) % 60) << std::setw(8)
	    << ". Write spectra to file" << std::endl;


      TObject *obj; 
      
      foutFile->mkdir("GammaGamma","GammaGamma");
      foutFile->mkdir("long","long");
      foutFile->mkdir("Energy_time_diff","Energy_time_diff");


      
       while ((obj = iter())) {
           TString name = obj->GetName();
//             std::cout<<"name "<<name<<std::endl;

           if(name.Contains("mDelila_long")){
           foutFile->cd(Form("%s:/", OutputFile.str().c_str()));
           }else if(name.Contains("mEnergy_time_diff")){
               foutFile->cd(Form("%s:/Energy_time_diff", OutputFile.str().c_str()));      
           }else if (name.Contains("mgg_")){
             foutFile->cd(Form("%s:/GammaGamma", OutputFile.str().c_str()));      
//              std::cout<<Form("%s:/GammaGamma", OutputFile.str().c_str())<<std::endl;
            }else if (name.Contains("long")){
             foutFile->cd(Form("%s:/long", OutputFile.str().c_str()));      
            } else {
            foutFile->cd(Form("%s:/", OutputFile.str().c_str()));
            };  
            
           
            
           if ( obj->IsA()->InheritsFrom(TH1F::Class())){
            TH1 *h1 = (TH1*)obj;
                if (h1->GetEntries()>0) obj->Write();
            }else if (obj->IsA()->InheritsFrom(TH2F::Class())){
            TH2 *h2 = (TH2*)obj;
                if (h2->GetEntries()>0) obj->Write();
            }
        };
        
       outputTree->Write();
       foutFile->Close();
   

}

int DelilaSelector::GetCoincTimeCorrection(int dom1, int dom2)
{
 int coin_id = GetCoincID(dom1, dom2);
//  if (dom1 >= dom2 ){coin_id = dom1*100+dom2;}else {coin_id = dom2*100+dom1;};
    
//  int coin_id = dom1*100+dom2;
 int time_corr = 0;
 std::map<int, int >::iterator it = LUT_TA.find(coin_id);
 if(it != LUT_TA.end()){time_corr = LUT_TA[coin_id];};
 return time_corr;
}

void DelilaSelector::cs()
{
    int gate = 20000;
    double_t time_diff_cs = 0;
    int cs_dom = DelilaEvent.cs_domain;
    
    if (DelilaEvent.det_def == 1) gate = 100000;
    
//     int cc_id = GetCoinc_det_def(ggLaBr_HPGe_Qu.front().det_def,DelilaEvent.det_def);

//     if (DelilaEvent.det_def == 3) {gate = 20000; }
//     else if (DelilaEvent.det_def == 1) {gate = 200000; };
    
    
    //if (DelilaEvent.det_def == 3) {waitingQu_gamma[cs_dom].push_back(DelilaEvent); hDelila->Fill(DelilaEvent.EnergyCal);hDelilaDC->Fill(DelilaEvent.EnergyDC);}
    if ((DelilaEvent.det_def == 3)||(DelilaEvent.det_def == 1)) {
        waitingQu_gamma[cs_dom].push_back(DelilaEvent); 
        hDelila[DelilaEvent.det_def]->Fill(DelilaEvent.EnergyCal);
        hDelilaDC[DelilaEvent.det_def]->Fill(DelilaEvent.EnergyDC);
        hDelila_long[DelilaEvent.det_def]->Fill(DelilaEvent.EnergyCal);
        hDelilaDC_long[DelilaEvent.det_def]->Fill(DelilaEvent.EnergyDC);
    }
    if (DelilaEvent.det_def == 5) waitingQu_bgo[cs_dom].push_back(DelilaEvent);
        
        if ((!waitingQu_gamma[cs_dom].empty())&&(!waitingQu_bgo.empty()))
        {
            std::deque<TDelilaEvent>  ::iterator it_g__ = waitingQu_gamma[cs_dom].begin();
            std::deque<TDelilaEvent>  ::iterator it_bgo__ = waitingQu_bgo[cs_dom].begin();
            for (; it_g__ != waitingQu_gamma[cs_dom].end();++it_g__){
                for (; it_bgo__ != waitingQu_bgo[cs_dom].end();++it_bgo__){   
                time_diff_cs = it_g__->Time - it_bgo__->Time -  GetCoincTimeCorrection(it_g__->domain,it_bgo__->domain);
                mLaBr_BGO_time_diff->Fill(it_g__->domain, time_diff_cs);
                if (abs(time_diff_cs)<gate){
                    if (it_g__->CS == 1) continue;
                    it_g__->CS = 1;   
                    it_g__->bgo_time_diff = time_diff_cs;
//                      if (it_g__->det_def == 1) {std::cout<<time_diff_cs<<" time_diff_cs \n";};
                    //mLaBr_BGO_time_diff->Fill(it_g__->domain, time_diff_cs);
                    }
                }
            }
        }
        
         if (!waitingQu_bgo[cs_dom].empty())
         {
             std::deque<TDelilaEvent>  ::iterator it1_ = waitingQu_bgo[cs_dom].begin();
             for (; it1_ != waitingQu_bgo[cs_dom].end();)
             {
              if (abs(DelilaEvent.Time - it1_->Time)>gate) {it1_=waitingQu_bgo[cs_dom].erase(it1_);}
                  else ++it1_;
             }
         };
         
         if (!waitingQu_gamma[cs_dom].empty())
         {
             std::deque<TDelilaEvent>  ::iterator it2_ = waitingQu_gamma[cs_dom].begin();
             for (; it2_ != waitingQu_gamma[cs_dom].end();)
             {
              if (abs(DelilaEvent.Time - it2_->Time)>gate) 
              {
                 // output_pQu.push(*it2_);                  
                  if (it2_->CS ==0){
                      hDelilaCS[it2_->det_def]->Fill(it2_->EnergyCal);
                      hDelilaCS_long[it2_->det_def]->Fill(it2_->EnergyCal);
                      mDelilaCS->Fill(it2_->cs_domain, it2_->EnergyCal);
                      gamma_gamma_cs(*it2_);
                      hDelilaCS_DC[it2_->det_def]->Fill(it2_->EnergyDC);
                      hDelilaCS_DC_long[it2_->det_def]->Fill(it2_->EnergyDC);
                      mDelilaCS_DC->Fill(it2_->cs_domain, it2_->EnergyDC);
                };
                  
                 if (blOutTree) {DelilaEventCS = *it2_; outputTree->Fill();};
                  
                  it2_=waitingQu_gamma[cs_dom].erase(it2_);
              }
              else ++it2_;
             };
         };
};

void DelilaSelector::gamma_gamma()
{
 if (gammagammaQu.empty()){gammagammaQu.push_back(DelilaEvent);/*std::cout<<"Empty Coic \n";*/}
    else
         {
                           
//              double time_diff = DelilaEvent.Time - coincQu_TA.front().Time - GetCoincTimeCorrection(DelilaEvent.domain,coincQu_TA.front().domain);
             double_t time_diff_gg = DelilaEvent.Time - gammagammaQu.front().Time - GetCoincTimeCorrection(DelilaEvent.domain, gammagammaQu.front().domain);
             //mLaBr_LabBr_time_diff->Fill(DelilaEvent.domain, time_diff_gg);
//              mLaBr_LabBr_time_diff->Fill(gammagammaQu.front().domain, time_diff_gg);
             
             
             if (std::abs(time_diff_gg) < 10000)//2000 ps
             {
                 gammagammaQu.push_back(DelilaEvent);
             }
             else
             {
                hMult_gg->Fill(gammagammaQu.size());             
                std::deque<TDelilaEvent>  ::iterator it1__ = gammagammaQu.begin();
                std::deque<TDelilaEvent>  ::iterator it2__ = gammagammaQu.begin(); 
                for (; it1__ != gammagammaQu.end(); ++it1__){
                     it2__ = gammagammaQu.begin(); 
                  for (; it2__ != gammagammaQu.end(); ++it2__){
                      if (it1__ == it2__) continue;
                      mGammaGamma->Fill((*it1__).EnergyCal, (*it2__).EnergyCal);
                      mGammaGammaDC->Fill((*it1__).EnergyCal, (*it2__).EnergyCal);
                      mLaBr_LabBr_time_diff->Fill((*it1__).domain,(*it1__).Time - (*it2__).Time -  GetCoincTimeCorrection((*it1__).domain,(*it2__).domain));
                  };
                };
                gammagammaQu.clear();
                gammagammaQu.push_back(DelilaEvent);
             };
         };
         return;
};

void DelilaSelector::gamma_gamma_cs(TDelilaEvent &ev_)
{
    
 if (ev_.det_def != 3  ) return;
    
 if (gammagammaQu_CS.empty()){gammagammaQu_CS.push_back(ev_);/*std::cout<<"Empty Coic \n";*/}
    else
         {
             double_t time_diff_gg = ev_.Time - gammagammaQu_CS.front().Time - GetCoincTimeCorrection(ev_.domain, gammagammaQu_CS.front().domain);
             if (std::abs(time_diff_gg) < 2000) 
             {
                 gammagammaQu_CS.push_back(ev_);
             }
             else
             {
                hMult_gg->Fill(gammagammaQu_CS.size());             
                std::deque<TDelilaEvent>  ::iterator it1__ = gammagammaQu_CS.begin();
                std::deque<TDelilaEvent>  ::iterator it2__ = gammagammaQu_CS.begin(); 
                for (; it1__ != gammagammaQu_CS.end(); ++it1__){
                  for (; it2__ != gammagammaQu_CS.end(); ++it2__){
                      if (it1__ == it2__) continue;
                      mGammaGammaCS->Fill((*it1__).EnergyCal, (*it2__).EnergyCal);
                      mGammaGammaCS_DC->Fill((*it1__).EnergyCal, (*it2__).EnergyCal);
                  };
                };
                gammagammaQu_CS.clear();
                gammagammaQu_CS.push_back(ev_);
             };
         };
         return;
};
    
void DelilaSelector::time_alignment()
{
 ///To produce matrix for to Check Time Alignment
        if (coincQu_TA.empty()){coincQu_TA.push_back(DelilaEvent);}
         else
         {
//              double time_diff = DelilaEvent.Time - coincQu_TA.front().Time - GetCoincTimeCorrection(DelilaEvent.domain,coincQu_TA.front().domain);
              double time_diff = DelilaEvent.Time - coincQu_TA.front().Time + GetCoincTimeCorrection(DelilaEvent.domain,coincQu_TA.front().domain);

             if (std::abs(time_diff) < 400000) //ps was 4000000
              {
                 coincQu_TA.push_back(DelilaEvent);
              }
              else
              {
                 std::deque<TDelilaEvent>  ::iterator it1__ = coincQu_TA.begin();
                 std::deque<TDelilaEvent>  ::iterator it2__ = coincQu_TA.begin(); 
                 for (; it1__ != coincQu_TA.end(); ++it1__){
                   for (; it2__ != coincQu_TA.end(); ++it2__){
                       
                         double time_diff_temp = it1__->Time - it2__->Time - GetCoincTimeCorrection(it1__->domain,it2__->domain) - GetCoincTimeCorrection(DelilaEvent.domain,coincQu_TA.front().domain);
                         mTimeCalib->Fill(it1__->domain*100+it2__->domain, time_diff_temp); 
                   };
                 };
                 coincQu_TA.clear();
                 coincQu_TA.push_back(DelilaEvent);
              }
         }     
}

void DelilaSelector::time_alignment_symmetric()
{
 ///To produce matrix for to Check Time Alignment
        if (coincQu_TA.empty()){coincQu_TA.push_back(DelilaEvent);}
         else
         {
              double time_diff = DelilaEvent.Time - coincQu_TA.front().Time - GetCoincTimeCorrection(DelilaEvent.domain,coincQu_TA.front().domain);
//               double time_diff = DelilaEvent.Time - coincQu_TA.front().Time + GetCoincTimeCorrection(DelilaEvent.domain,coincQu_TA.front().domain);

             if (std::abs(time_diff) < 1000000) //ps was 4000000
              {
                 coincQu_TA.push_back(DelilaEvent);
              }
              else
              {
                 std::deque<TDelilaEvent>  ::iterator it1__ = coincQu_TA.begin();
                 std::deque<TDelilaEvent>  ::iterator it2__ = coincQu_TA.begin(); 
                 for (; it1__ != coincQu_TA.end(); ++it1__){
                   for (; it2__ != coincQu_TA.end(); ++it2__){
                         if (it1__ == it2__) continue;
                         double time_diff_temp = it1__->Time - it2__->Time - GetCoincTimeCorrection(it1__->domain,it2__->domain);
//                          if (time_diff_temp > 0) std::cout<<"!!! "<<time_diff_temp<<"\n";
                         mTimeCalib->Fill(GetCoincID((*it1__).domain, (*it2__).domain), time_diff_temp); 
                         //mTimeCalib->Fill(it1__->domain*100+it2__->domain, time_diff_temp); 
                   };
                 };
                 coincQu_TA.clear();
                 coincQu_TA.push_back(DelilaEvent);
              }
         }     
}

 void DelilaSelector::time_alignment_dom0(int zero_dom)
 {
    double time_diff_dom0;
    int cur_dom = DelilaEvent.domain;
    
 //    if (cur_dom != zero_dom){
        time_diff_dom0 =  DelilaEvent.Time - DomainZeroEvent.Time;// - GetCoincTimeCorrection(PulserEvent.domain, DelilaEvent.domain);;
        mTimeCalibDomain0->Fill(cur_dom, time_diff_dom0);
 //        std::cout<<" time_alignment_pulser "<<cur_dom <<"  DelilaEvent.domain "<<DelilaEvent.Time<<" DomainZeroEvent.Time "<<DomainZeroEvent.Time<<" dt ="<<time_diff_dom0<<" \n";
        if (cur_dom == zero_dom) DomainZeroEvent = DelilaEvent;
 //    }else {DomainZeroEvent = DelilaEvent;};//else {DomainZeroEvent = DelilaEvent; std::cout<<" time_alignment_pulser "<<cur_dom <<"  DelilaEvent.domain "<<DelilaEvent.Time<<" DomainZeroEvent.Time "<<DomainZeroEvent.Time<<" dt ="<<time_diff_dom0<<" \n";exit(1);};
 }

void DelilaSelector::gamma_gamma_LaBr_HPGe()
{
    //if ((waitingQu_gg_[1].empty()) && ((waitingQu_gg_[3].empty())) {waitingQu_gg_[DelilaEvent.domain].push_back();};
    if (ggLaBr_HPGe_Qu.empty()) {ggLaBr_HPGe_Qu.push_back(DelilaEvent);return;};
    
    //int cc_id = ggLaBr_HPGe_Qu.front().det_def*10+DelilaEvent.det_def;//std::cout<<cc_id<<" cc_id1 \n";
    int cc_id = GetCoinc_det_def(ggLaBr_HPGe_Qu.front().det_def,DelilaEvent.det_def);
    double_t time_diff_gg = DelilaEvent.Time - ggLaBr_HPGe_Qu.front().Time - GetCoincTimeCorrection(DelilaEvent.domain, ggLaBr_HPGe_Qu.front().domain);
//     mGG_time_diff[cc_id] ->Fill(DelilaEvent.domain,time_diff_gg);
//      std::cout<<" cc_id "<<cc_id<<std::endl;
     //if (std::abs(time_diff_gg) < 2000)
//      if (std::abs(time_diff_gg) < 20000)
    if (std::abs(time_diff_gg) < coinc_gates[cc_id]) 
              {
                  ggLaBr_HPGe_Qu.push_back(DelilaEvent);
              }
              else
              {
 //                 hMult_gg->Fill(gammagammaQu_CS.size());             
                  std::deque<TDelilaEvent>  ::iterator it1__ = ggLaBr_HPGe_Qu.begin();
                  std::deque<TDelilaEvent>  ::iterator it2__ = ggLaBr_HPGe_Qu.begin(); 
                  for (; it1__ != ggLaBr_HPGe_Qu.end(); ++it1__){
                   for (; it2__ != ggLaBr_HPGe_Qu.end(); ++it2__){
                        if (it1__ == it2__) continue;
                        cc_id=GetCoinc_det_def((*it1__).det_def, (*it2__).det_def);
//                         cc_id = (*it1__).det_def*10+(*it2__).det_def;
 //                       if ((*it1__).det_def == 3) continue;
 //                       if ((*it2__).det_def == 1) continue;
 //                       mgg_labr_hpge->Fill((*it1__).EnergyCal, (*it2__).EnergyCal);
 //                       mgg_labr_hpgeDC->Fill((*it1__).EnergyCal, (*it2__).EnergyCal);
                        
                        if (cc_id == 13){
                            std::deque<TDelilaEvent>  ::iterator it_hpge__;
                            std::deque<TDelilaEvent>  ::iterator it_labr__;
                            if ((*it1__).det_def == 1) {it_hpge__ = it1__; it_labr__ = it2__;}
                            else {it_hpge__ = it2__; it_labr__ = it1__;};
//                             mGG_time_diff[cc_id]->Fill((*it_hpge__).domain,(*it_hpge__).Time - (*it_labr__).Time -  GetCoincTimeCorrection((*it_hpge__).domain,(*it_labr__).domain));
//                             mGG[cc_id]->Fill((*it_hpge__).EnergyCal, (*it_labr__).EnergyCal);
//                             mGG_DC[cc_id]->Fill((*it_hpge__).EnergyDC, (*it_labr__).EnergyDC);
                            mGG_time_diff[cc_id]->Fill((*it_labr__).domain,(*it_labr__).Time - (*it_hpge__).Time -  GetCoincTimeCorrection((*it_labr__).domain,(*it_hpge__).domain));
                            mGG[cc_id]->Fill((*it_labr__).EnergyCal, (*it_hpge__).EnergyCal);
                            mGG_DC[cc_id]->Fill((*it_labr__).EnergyDC, (*it_hpge__).EnergyDC);
                            mGG_long[cc_id]->Fill((*it_labr__).EnergyCal, (*it_hpge__).EnergyCal);
                            mGG_DC_long[cc_id]->Fill((*it_labr__).EnergyDC, (*it_hpge__).EnergyDC);
                            
                        }else{                       
                            mGG_time_diff[cc_id]->Fill((*it1__).domain,(*it1__).Time - (*it2__).Time -  GetCoincTimeCorrection((*it1__).domain,(*it2__).domain));
                            mGG[cc_id]->Fill((*it1__).EnergyCal, (*it2__).EnergyCal);
                            mGG_DC[cc_id]->Fill((*it1__).EnergyDC, (*it2__).EnergyDC);
                            mGG_long[cc_id]->Fill((*it1__).EnergyCal, (*it2__).EnergyCal);
                            mGG_DC_long[cc_id]->Fill((*it1__).EnergyDC, (*it2__).EnergyDC);
//                         if(cc_id == 13)std::cout<<"td hpge - labr "<<((*it1__).Time - (*it2__).Time -  GetCoincTimeCorrection((*it1__).domain,(*it2__).domain))<<"\n";
                        };
                    };
                  };
                 ggLaBr_HPGe_Qu.clear();
                 ggLaBr_HPGe_Qu.push_back(DelilaEvent);
              };
             return;
};


 int DelilaSelector::GetCoinc_det_def(int det_def1, int det_def2)
 {
     int id=0;
     if (det_def1<=det_def2) {id = det_def1*10+det_def2;}
         else {id = det_def2*10+det_def1; };
     return id;
 }


int DelilaSelector::GetCoincID(int dom1, int dom2){
    int id=0;
    if (dom1<=dom2) {id = dom1*100+dom2;}
        else {id = dom2*100+dom1; };
    return id;
}


void DelilaSelector::CheckPulserAllignement(int zero_dom)
{
   double time_diff_pulser;
   int cur_dom = DelilaEvent.domain;
   
   if (cur_dom != zero_dom){
       time_diff_pulser =  DelilaEvent.Time - PulserEvent.Time - GetCoincTimeCorrection(PulserEvent.domain, DelilaEvent.domain);;
       mPulser0TimeDiff->Fill(cur_dom, time_diff_pulser);
//        std::cout<<" CheckPulserAllignement "<<cur_dom <<"  DelilaEvent.domain "<<DelilaEvent.Time<<" PulserEvent.Time "<<PulserEvent.Time<<" dt ="<<time_diff_pulser<<" \n";
//        if (DelilaEvent.domain == 57) std::cout<<" time_diff_pulser #50 - #57 "<<time_diff_pulser<<" \n";
//        if (DelilaEvent.domain == 51) std::cout<<" time_diff_pulser #50 - #51 "<<time_diff_pulser<<" \n";
   }else PulserEvent = DelilaEvent;
   return;
}


