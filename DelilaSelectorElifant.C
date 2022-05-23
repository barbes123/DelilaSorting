#define DelilaSelectorElifant_cxx
// The class definition in data.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after hTimeDiffCoreCoreBegin(), when on PROOF called only on the
//                    slave servers.selected_trg_run_459_0.rselected_trg_run_459_0.roobunch_resettoot
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("DelilaSelectorEliade.C")
// root> T->Process("DelilaSelectorEliade.C","some options")
// root> T->Process("DelilaSelectorEliade.C+")
//


#include "DelilaSelectorElifant.h"
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TObjString.h>
#include <unordered_set>
#include <iomanip>      // sTreatDelilaEvent_td::setwsorted
using namespace std;


////////////////////////////////Please, modify if needed////////////////////////////////////////////
bool blGammaGamma       = true;
bool blCS               = false;
bool blOutTree          = true;
bool blFold             = false;
bool blTimeAlignement   = true;
////////////////////////////////Please, DO NOT modify ////////////////////////////////////////////
bool blIsTrigger            = false; //the trigger is open
bool blIsWindow             = false; //the trigger is open
bool blFirstTrigger         = false;
bool blAddTriggerToQueue    = false;
bool blFillAmaxEnergyDom    = false;

bool debug            = false;
bool blDebugElissa    = false;

ULong64_t trigger_cnt = 0;
ULong64_t trigger_events = 0;

// const int NumberOfClovers = 2;
const int max_domain = 500;
const int nbr_of_ch = 500;
// ULong64_t lastTimeStampTrigger = 0;

double beta = 0;

std::stringstream OutputFile;


void DelilaSelectorElifant::Read_ELIADE_LookUpTable() {
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
      int threshold = 1e6;
      std::istringstream is(oneline);
      if (debug) std::cout << is.str().c_str() << std::endl;
//       is >> curDet.ch >> curDet.dom >> curDet.theta >> curDet.phi >> curDet.TimeOffset >> curDet.threshold;
      is >> curDet.ch >> curDet.dom >> curDet.detType >> curDet.serial >> curDet.bgo_time_offset>> curDet.theta >> curDet.phi >> curDet.threshold >> curDet.cs_dom;
    //  std::cout<<" curDfalseet.ch  "<<curDet.ch <<" curDet.TimeOffset " <<curDet.TimeOffset<<std::endl;
      
      if (curDet.ch >= 0) {
          curDet.theta *= TMath::DegToRad();
          curDet.phi *= TMath::DegToRad();
	int pol_order = 0;
 	float offset_gate(0.),slope_gate(1.);
    is >> offset_gate;
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


void DelilaSelectorElifant::Read_TimeAlignment_LookUpTable() {
  std::cout << "I am Reading TimeAlignment LookUpTable ... ";
 
  char* pLUT_Path;
  pLUT_Path = getenv ("DELILA_LUT");
  if (pLUT_Path!=NULL)
    printf ("The LookUpTable path is: %s \n",pLUT_Path);


  std::stringstream LUTFile;
  LUTFile << pLUT_Path <<"/"<<"LUT_TA.dat";
//   const int nbr_of_ch = 200;
  std::ifstream lookuptable(LUTFile.str().c_str());

  if (!lookuptable.good()) {
    std::ostringstream os;
    os << "Could not open " << LUTFile.str().c_str()
    << " no time alignment will be done\n";
//        << " and I need it ;(\n";
//     Abort(os.str().c_str());
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

void DelilaSelectorElifant::Read_TimeAlignment_Trigger() {
  std::cout << "I am Reading Read_TimeAlignment_Trigger LookUpTable ... ";
 
  char* pLUT_Path;
  pLUT_Path = getenv ("DELILA_LUT");
  if (pLUT_Path!=NULL)
    printf ("The LookUpTable path is: %s \n",pLUT_Path);


  std::stringstream LUTFile;
  LUTFile << pLUT_Path <<"/"<<"LUT_TRIGGER.dat";
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

void DelilaSelectorElifant::Read_Confs() {
  std::cout << "I am Reading Your Configuration parameters ... ";
 
  char* pLUT_Path;
  pLUT_Path = getenv ("DELILA_LUT");
  if (pLUT_Path!=NULL)
    printf ("The LookUpTable path is: %s \n",pLUT_Path);


  std::stringstream LUTFile;
  LUTFile << pLUT_Path <<"/"<<"LUT_CONF.dat";
  std::ifstream lookuptable(LUTFile.str().c_str());
  
        
  event_length = 40000;
  pre_event_length = 0;
  bunch_reset = 0;
  det_def_trg = -1;
  channel_trg = -1;
  TriggerTimeFlag = 0;

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
      int coinc_id = 0; Float_t value = 0;
      is >> coinc_name>> coinc_id >> value;

      switch (coinc_id){
          case 1000: {
              event_length = value;
              std::cout<<"event_length "<<event_length<<" ps \n";
              break;
          }
          case 1001:
              {
              bunch_reset = value;
              std::cout<<"bunch_reset "<<bunch_reset<<" ps \n";
              break;
          }
          case 1002:
              {
              pre_event_length = value;
              std::cout<<"pre_event_length "<<pre_event_length<<" ps \n";
              break;
          }
          case 1111:{
               beta = value;
               std::cout<<"Beta is "<<beta<<" % \n";
              break;
          };
          case 9997:{
               std::cout<<" Enabled detectors \n";
               int det_list = value/1;
               while(det_list!= 0){
                   has_detector[detector_name[det_list%10]] = true;
                   std::cout<<" id " << det_list%10 <<" det " << detector_name[det_list%10] << " is "<<has_detector[detector_name[det_list%10]]<<"\n";
                   det_list = det_list/10;
            };
              break;
          }
          case 9998:{
               det_def_trg = value/1;
               std::cout<<"trg_det_type "<<det_def_trg<<" \n";
              break;
          }
          case 9999:{
               channel_trg = value/1;
               std::cout<<"channel_trg  "<<channel_trg<<" \n";
              break;
          }
          default:
          {
              //LUT_CONF[coinc_id] = value;
              coinc_gates[coinc_id] = value;
              //std::cout<<coinc_name<<" coin_id " << coinc_id <<" value "<<value <<" ps \n";
//               std::cout<<coinc_name<<" coin_id " << coinc_id <<"  coinc_gates[coinc_id] "<< coinc_gates[coinc_id] <<" ps \n";
              break;
          };
      };
  }
  lookuptable.close();
  }
  std::cout << " done" << std::endl;
}



void DelilaSelectorElifant::Print_ELIADE_LookUpTable()
{
    std::cout<<"Print_ELIADE_LookUpTable \n";		
    std::map<int, TDelilaDetector > ::iterator it__ = LUT_DELILA.begin();
    for (; it__ != LUT_DELILA.end(); ++it__) {
     // is >> curDet.ch >> curDet.dom >> theta >> phi >> curDet.TimeOffset >> curDet.threshold;
	std::cout<<" Ch "<<LUT_DELILA[it__->first].ch<<" Dom "<< LUT_DELILA[it__->first].dom<<" "<< LUT_DELILA[it__->first].theta<<" "<< LUT_DELILA[it__->first].phi <<" offset "<< LUT_DELILA[it__->first].TimeOffset<<" Thr "<< LUT_DELILA[it__->first].threshold<<" serial "<<LUT_DELILA[it__->first].serial<<" theta "<<LUT_DELILA[it__->first].theta<<" phi "<<LUT_DELILA[it__->first].phi <<" cs_dom: "<<LUT_DELILA[it__->first].cs_dom<<" pol_order: " <<LUT_DELILA[it__->first].pol_order <<std::endl;
    }
};

void DelilaSelectorElifant::Print_TimeAlignment_LookUpTable()
{
    std::cout<<"Print_TimeAlignment_LookUpTable \n";		
    std::map<int, int > ::iterator it__ = LUT_TA.begin();
    for (; it__ != LUT_TA.end(); ++it__) {
     // is >> curDet.ch >> curDet.dom >> theta >> phi >> curDet.TimeOffset >> curDet.threshold;
	std::cout<<" coinc_id "<<it__->first<<" time_corr "<< it__->second<<std::endl;
    }
};

void DelilaSelectorElifant::Print_TimeAlignment_Trigger_LookUpTable()
{
    std::cout<<"Print_TimeAlignment_trigger_LookUpTable \n";		
    std::map<int, double_t > ::iterator it__ = LUT_TA_TRG.begin();
    for (; it__ != LUT_TA_TRG.end(); ++it__) {
     // is >> curDet.ch >> curDet.dom >> theta >> phi >> curDet.TimeOffset >> curDet.threshold;
	std::cout<<" coinc_id "<<it__->first<<" time_corr "<< it__->second<<std::endl;
    }
};



float DelilaSelectorElifant::CalibDet(float val, int ch)
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

void DelilaSelectorElifant::Begin(TTree * tree)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
  
  detector_name[1]="HPGe";     has_detector[detector_name[1]] = false;
  detector_name[2]="SEG";      has_detector[detector_name[2]] = false;
  detector_name[3]="LaBr";     has_detector[detector_name[3]] = false;
  detector_name[4]="ACS";      has_detector[detector_name[4]] = false;
  detector_name[5]="BGOs";     has_detector[detector_name[5]] = false; //side
  detector_name[6]="BGOf";     has_detector[detector_name[6]] = false;//front
  detector_name[7]="Elissa";   has_detector[detector_name[7]] = false;
  detector_name[8]="neutron";  has_detector[detector_name[8]] = false;
  detector_name[9]="pulser";   has_detector[detector_name[9]] = false;  

  Read_Confs();

  std::map<UInt_t, Float_t>::iterator it_c_gates_ =  coinc_gates.begin();
  for(;it_c_gates_!=coinc_gates.end();++it_c_gates_){
       
      switch (it_c_gates_->first){
           case 11: {
               gg_coinc_id[it_c_gates_->first]="mgg_hpge_hpge";
               break;
         };
           case 12: {
               gg_coinc_id[it_c_gates_->first]="mgg_core_seg";
               break;
         };
          case 13: {
               gg_coinc_id[it_c_gates_->first]="mgg_labr_hpge";
               break;
         };
           case 33: {
               gg_coinc_id[it_c_gates_->first]="mgg_labr_labr";
               break;
         };
           case 37: {
               gg_coinc_id[it_c_gates_->first]="mgg_labr_elissa";
               break;
         };
         default: {
//               gg_coinc_id[it_c_gates_->first]=Form("mgg_not_defined_%i",it_c_gates_->first);
               break;
         };
       };
    };
    
   Read_ELIADE_LookUpTable();

}




void DelilaSelectorElifant::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
    
//    TString option = GetOption();
//    toks = option.Tokenize(",");
//    TString RunID = ((TObjString*) toks->At(0))->GetString();
//   double beta = ((TObjString*) toks->At(2))->GetString().Atof();
//    beta = ((TObjString*) toks->At(2))->GetString().Atof();
//    addBackMode = atoi(((TObjString*) toks->At(3))->GetString());
//    std::cout << "addBackMode  " << addBackMode <<std::endl;
   
   
   TString option = GetOption();
   toks = option.Tokenize(",");
   TString RunID = ((TObjString*) toks->At(0))->GetString();
   TString VolID = ((TObjString*) toks->At(1))->GetString();
   TString ServerID = ((TObjString*) toks->At(3))->GetString();
   
   

  
   nevents = 0;
//  nevents_reset=0;
  reset_counter = 0;

   hTimeSort = new TH1F("hTimeSort", "time_diff: current-last", 1e3, -1e5,1e5);
   hTimeSort->GetXaxis()->SetTitle("ps");
   fOutput->Add(hTimeSort);
   
   hTimeZero = new TH1F("hTimeZero", "Events with zero time", 500, -0.5, 499.5);
   hTimeZero->GetXaxis()->SetTitle("ch");
   hTimeZero->GetYaxis()->SetTitle("counts");   
   fOutput->Add(hTimeZero);
    
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
   
   
   mElissa = new TH2F("mElissa", "mElissa", max_domain, 0, max_domain, 4096, -0.5, 32770);
   mElissa->GetXaxis()->SetTitle("domain");
   mElissa->GetYaxis()->SetTitle("keV");
   fOutput->Add(mElissa);
   
   
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
   
   
   hdelilaQu_size = new TH1F("hdelilaQu_size", "hdelilaQu_size", 100,0,100);
   hdelilaQu_size->GetXaxis()->SetTitle("size");
   hdelilaQu_size->GetYaxis()->SetTitle("counts");
   fOutput->Add(hdelilaQu_size);
   
   
   
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
   
   
   
   
   if (has_detector["Elissa"]){
        std::map<int, TDelilaDetector > ::iterator it_lut_ = LUT_DELILA.begin();
        
        if (blFillAmaxEnergyDom){
            for (; it_lut_ != LUT_DELILA.end(); ++it_lut_) {
                if (LUT_DELILA[it_lut_->first].detType == 7){
                    int dom = LUT_DELILA[it_lut_->first].dom;
    //                 std::cout<<"dom "<<dom<<"\n";
                    mAmaxEnergyDom[dom] = new TH2F(Form("mAmaxEnergy_dom%i",dom), Form("mAmaxEnergy_dom%i",dom), 4096,0, 16384, 2e4,0,2e4);
                    mAmaxEnergyDom[dom] ->GetXaxis()->SetTitle("Energy, a.u.");
                    mAmaxEnergyDom[dom] ->GetYaxis()->SetTitle("rise time (Amax)");
                    fOutput->Add(mAmaxEnergyDom[dom]); 
                };
            };
        };
    
        hAmax = new TH1F("hAmax", "hAmax", 10e4,0,10e4);
        hAmax->GetXaxis()->SetTitle("rise time (Amax)");
        hAmax->GetYaxis()->SetTitle("counts");
        fOutput->Add(hAmax);
        
        mAmaxEnergy = new TH2F("mAmaxEnergy", "mAmaxEnergy", 2048,0, 16384, 10e2,0,10e4);
        mAmaxEnergy->GetXaxis()->SetTitle("Energy, a.u.");
        mAmaxEnergy->GetYaxis()->SetTitle("rise time (Amax)");
        fOutput->Add(mAmaxEnergy);
   };
//   std::map<UInt_t, Float_t>::iterator it_c_gates_ =  coinc_gates.begin();

//   for(;it_c_gates_!=coinc_gates.end();++it_c_gates_){
//       
//       switch (it_c_gates_->first){
//           case 11: {
//               gg_coinc_id[it_c_gates_->first]="mgg_hpge_hpge";
//               break;
//         };
//           case 12: {
//               gg_coinc_id[it_c_gates_->first]="mgg_core_seg";
//               break;
//         };
//          case 13: {
//               gg_coinc_id[it_c_gates_->first]="mgg_labr_hpge";
//               break;
//         };
//           case 33: {
//               gg_coinc_id[it_c_gates_->first]="mgg_labr_labr";
//               break;
//         };
//           case 37: {
//               gg_coinc_id[it_c_gates_->first]="mgg_labr_elissa";
//               break;
//         };
//         default: {
// //               gg_coinc_id[it_c_gates_->first]=Form("mgg_not_defined_%i",it_c_gates_->first);
//               break;
//         };
//       }
//   }
   

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
  
  hMult[itna->first] = new TH1F(Form("%s_hMult",itna->second.c_str()), Form("%s_hMult",itna->second.c_str()), 20,0,20);
  hMult[itna->first]->GetXaxis()->SetTitle("Multiplicity");
  hMult[itna->first]->GetYaxis()->SetTitle("Counts");
  fOutput->Add(hMult[itna->first]);

//   mGG_theta[itna->first] = new TH2F(Form("%s_mGG_theta",itna->second.c_str()), Form("%s",itna->second.c_str()),  360, -180, 180, 4096, -0.5, 16383.5);
  
//   if (itna->first == 33){
//   cGGTheta = new TH3F("cGGTheta", "cGGTheta",  6400, -3.2, 3.2, 4096, -0.5, 16383.5, 4096, -0.5, 16383.5);
//   cGGTheta->GetXaxis()->SetTitle("#delta#theta, rad");
//   cGGTheta->GetYaxis()->SetTitle("4 keV / bin");
//   cGGTheta->GetZaxis()->SetTitle("4 keV / bin");
//   fOutput->Add(cGGTheta);
// 
//   //cGGG = new TH3F("cGGG", "cGGG", 4096, -0.5, 16383.5, 4096, -0.5, 16383.5, 4096, -0.5, 16383.5);
//   cGGG = new TH3F("cGGG", "cGGG", 2048, -0.5, 16383.5, 2048, -0.5, 16383.5, 2048, -0.5, 16383.5);
//   cGGG->GetXaxis()->SetTitle("4 keV / bin");
//   cGGG->GetYaxis()->SetTitle("4 keV / bin");
//   cGGG->GetZaxis()->SetTitle("4 keV / bin");
//   fOutput->Add(cGGG);  
//   }
  
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
   
   if (itna->first == 37){
       mGG[itna->first]->GetXaxis()->SetTitle("LaBr, keV"); mGG[itna->first]->GetYaxis()->SetTitle("Solar, keV");
   };
   
   if (itna->first == 11){//for core-core
//     mGG_time_diff[itna->first] = new TH2F(Form("%s_time_diff",itna->second.c_str()), Form("%s_time_diff",itna->second.c_str()), 100, 0, 100, 10e3, -2e6, 2e6);//was tuned like that
       mGG_time_diff[itna->first] = new TH2F(Form("%s_time_diff",itna->second.c_str()), Form("%s_time_diff",itna->second.c_str()), max_domain, 0, max_domain, 4e2, -2e6, 2e6);

   }else{
     mGG_time_diff[itna->first] = new TH2F(Form("%s_time_diff",itna->second.c_str()), Form("%s_time_diff",itna->second.c_str()), max_domain, 0, max_domain, 4e4, -2e6, 2e6);
   };
     mGG_time_diff[itna->first]->GetXaxis()->SetTitle("domain"); mGG_DC[itna->first]->GetYaxis()->SetTitle("ps");
     fOutput->Add(mGG_time_diff[itna->first]);
//    };
//   std::cout<<itna->first<<" "<< Form("%s",itna->second.c_str())<<" Initialized \n" ;
};
   
  std::map<UInt_t,std::string>::iterator itna1 =  detector_name.begin();

  for(;itna1!=detector_name.end();++itna1){
   hDelila0[itna1->first] = new TH1F(Form("%s",itna1->second.c_str()), Form("%s before EB",itna1->second.c_str()), 4096, -0.5, 16383.5);
   hDelila0[itna1->first]->GetYaxis()->SetTitle("counts");
   hDelila0[itna1->first]->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelila0[itna1->first]);
   
   hDelila_single[itna1->first] = new TH1F(Form("%s_single",itna1->second.c_str()), Form("%s_single",itna1->second.c_str()), 4096, -0.5, 16383.5);
   hDelila_single[itna1->first]->GetYaxis()->SetTitle("counts");
   hDelila_single[itna1->first]->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelila_single[itna1->first]);
   
   
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
   
   mDomainTimeDiff_trigger = new TH2F("mDomainTimeDiff_trigger", "mDomainTimeDiff_trigger", max_domain, 0, max_domain, 10e3, -1e6, 9e6);
//    mDomainTimeDiff_trigger = new TH2F("mDomainTimeDiff_trigger", "mDomainTimeDiff_trigger", max_domain, 0, max_domain, 3e2, 0, 3e5);
   mDomainTimeDiff_trigger->GetXaxis()->SetTitle("domain");
   mDomainTimeDiff_trigger->GetYaxis()->SetTitle("Time diff, ps");
   fOutput->Add(mDomainTimeDiff_trigger);
   
   mDomainTimeDiff_bunch = new TH2F("mDomainTimeDiff_bunch", "mDomainTimeDiff_bunch", max_domain, 0, max_domain, 10e3, -1e6, 9e6);
//    mDomainTimeDiff_bunch = new TH2F("mDomainTimeDiff_bunch", "mDomainTimeDiff_bunch", max_domain, 0, max_domain, 3e2, 0, 3e5);
   mDomainTimeDiff_bunch->GetXaxis()->SetTitle("domain");
   mDomainTimeDiff_bunch->GetYaxis()->SetTitle("Time diff, ps");
   fOutput->Add(mDomainTimeDiff_bunch);
   
 //  hTriggerTrigger = new TH1F("hTriggerTrigger", "hTriggerTrigger", 100e4, 0, 100e6);
     hTriggerTrigger = new TH1F("hTriggerTrigger", "hTriggerTrigger", 1e4, 0, 1e10);
   hTriggerTrigger->GetYaxis()->SetTitle("Counts");
   hTriggerTrigger->GetXaxis()->SetTitle("Time diff, ps");
   hTriggerTrigger->SetTitle("Time between two trigger signals");
   fOutput->Add(hTriggerTrigger);
   
   hBunchBunch = new TH1F("hBunchBunch", "hBunchBunch", 10e4, -1e6, 9e6);
   hBunchBunch->GetYaxis()->SetTitle("Counts");
   hBunchBunch->GetXaxis()->SetTitle("Time diff, ps");
   hBunchBunch->SetTitle("Time between two Bunch signals");
   fOutput->Add(hBunchBunch);
   
   hCoincID= new TH1F("hCoincID", "hCoincID", 100,0,100);
   hCoincID->GetYaxis()->SetTitle("Counts");
   hCoincID->GetXaxis()->SetTitle("CoincID");
   hCoincID->SetTitle("CoincID");
   fOutput->Add(hCoincID);
   
   hTriggerDomain= new TH1F("hTriggerDomain", "hTriggerDomain", max_domain,-0.5,max_domain-0.5);
   hTriggerDomain->GetYaxis()->SetTitle("Counts");
   hTriggerDomain->GetXaxis()->SetTitle("Trigger Domain");
   hTriggerDomain->SetTitle("Trigger Domain");
   fOutput->Add(hTriggerDomain);
   
   hNBunch = new TH1F("hNBunch", "hNBunch", 100, 0, 100);
   hNBunch->GetYaxis()->SetTitle("Counts");
   hNBunch->GetXaxis()->SetTitle("N bunches after Trigger");
   hNBunch->SetTitle("N bunches");
   fOutput->Add(hNBunch);
   
   hTimeInBunch = new TH1F("hTimeInBunch", "hTimeInBunch", 10e3, -1e5, 9e5);
   hTimeInBunch->GetYaxis()->SetTitle("Counts");
   hTimeInBunch->GetXaxis()->SetTitle("hTimeInBunch");
   hTimeInBunch->SetTitle("time, ps");
   fOutput->Add(hTimeInBunch);
   
   hEventsPerTrigger = new TH1F("hEventsPerTrigger", "hEventsPerTrigger", 100, 0, 100);
   hEventsPerTrigger->GetYaxis()->SetTitle("Counts");
   hEventsPerTrigger->GetXaxis()->SetTitle("hEventsPerTrigger");
   hEventsPerTrigger->SetTitle("trigger");
   fOutput->Add(hEventsPerTrigger);
   
   hBunchFold = new TH1F("hBunchFold", "hBunchFold", 10, 0, 10);
   hBunchFold->GetYaxis()->SetTitle("Counts");
   hBunchFold->GetXaxis()->SetTitle("Det fold in bunch");
   hBunchFold->SetTitle("fold");
   fOutput->Add(hBunchFold);
   
   mFoldEnergy = new TH2F("mFoldEnergy", "mFoldEnergy", 10, 0, 10, 4096, -0.5, 65535.5);
   mFoldEnergy->GetXaxis()->SetTitle("fold");
   mFoldEnergy->GetYaxis()->SetTitle("Energy, 16 keV/pin");
   mFoldEnergy->SetTitle("Fold vs Energy");
   fOutput->Add(mFoldEnergy);
   
//    mTimeCalib = new TH2F("mTimeCalib", "mTimeCalib", 10000, 0, 10000, 2e3, -1e6, 1e6);
   mTimeCalib = new TH2F("mTimeCalib", "mTimeCalib", 10000, 0, 10000, 4e3, -2e5, 2e5);
   mTimeCalib->GetXaxis()->SetTitle("coinc ID");
   mTimeCalib->GetYaxis()->SetTitle("100 ps / bin");
   mTimeCalib->SetTitle("Sci time diff");
   fOutput->Add(mTimeCalib);
   
  mTimeCalibBGO = new TH2F("mTimeCalibBGO", "mTimeCalibBGO", max_domain, 0, max_domain, 2e3, -1e6, 1e6);
//    mTimeCalibBGO = new TH2F("mTimeCalibBGO", "mTimeCalibBGO", 10000, 0, 10000, 4e3, -2e5, 2e5);
   mTimeCalibBGO->GetXaxis()->SetTitle("coinc ID");
   mTimeCalibBGO->GetYaxis()->SetTitle("100 ps / bin");
   mTimeCalibBGO->SetTitle("BGO - HPGe/LaBr time diff");
   fOutput->Add(mTimeCalibBGO);
  /* 
   mTimeCalibBGO_cs_dom = new TH2F("mTimeCalibBGO_cs_dom", "mTimeCalibBGO_cs_dom", max_domain, 0, max_domain, 4e3, -2e5, 2e5);
   mTimeCalibBGO_cs_dom->GetXaxis()->SetTitle("coinc ID");
   mTimeCalibBGO_cs_dom->GetYaxis()->SetTitle("100 ps / bin");
   mTimeCalibBGO_cs_dom->SetTitle("BGO - HPGe/LaBr time diff");
   fOutput->Add(mTimeCalibBGO_cs_dom);
   */
   mTimeCalibDomain0 = new TH2F("mTimeCalibDomain0", "mTimeCalibDomain0", max_domain, -0.5, max_domain-0.5,  1e5, -1e6, 9e6);
   mTimeCalibDomain0->GetXaxis()->SetTitle("coinc ID");
   mTimeCalibDomain0->GetYaxis()->SetTitle("ps");
   mTimeCalibDomain0->SetTitle(Form("TimeDiff domain%i vs domain",channel_trg));
   fOutput->Add(mTimeCalibDomain0);
   
   
   mTimeDiffCS = new TH2F("mTimeDiffCS", "mTimeDiffCS", 100, 0, 100, 4e3, -2e6, 2e6);
   mTimeDiffCS->GetXaxis()->SetTitle("domain");
   mTimeDiffCS->GetYaxis()->SetTitle("ps");
   mTimeDiffCS->SetTitle("LaBr_BGO time diff");
   fOutput->Add(mTimeDiffCS);
   
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
   
   mNeutron = new TH2F("mNeutron", "mNeutron", 4096, 0, 16394, 4096, 0,1);
   mNeutron->GetXaxis()->SetTitle("Charge, a.u.");
   mNeutron->GetYaxis()->SetTitle("Short / Long");
   mNeutron->SetTitle("Neutron/#gamma separation");
   fOutput->Add(mNeutron);
   
   mShortLong = new TH2F("mShortLong", "mShortLong", 4096, 0, 16394, 4096, 0,16394);
   mShortLong->GetXaxis()->SetTitle("LongCharge, a.u.");
   mShortLong->GetYaxis()->SetTitle("ShortCharge, a.u.");
   mShortLong->SetTitle("Long vs short");
   fOutput->Add(mShortLong);
   
   
       
   std::map<int, TDelilaDetector > ::iterator it_lut_ = LUT_DELILA.begin();
    for (; it_lut_ != LUT_DELILA.end(); ++it_lut_) {
        if ((LUT_DELILA[it_lut_->first].detType == 1)||(LUT_DELILA[it_lut_->first].detType == 3)){
            int dom = LUT_DELILA[it_lut_->first].dom;
            mEnergy_time_diff[dom] = new TH2F(Form("mEnergy_time_diff_dom%i",dom), Form("mEnergy_time_diff_dom%i",dom), 4096, -0.5, 65535.5, 10e3, -1e5, 9e5);
            mEnergy_time_diff[dom] ->GetXaxis()->SetTitle("16 keV / bin");
            mEnergy_time_diff[dom] ->GetYaxis()->SetTitle("100 ps / bin");
            mEnergy_time_diff[dom]->SetTitle(" Time Diff Event.Time - Bunch.Time ");
            fOutput->Add(mEnergy_time_diff[dom]); 
        };
    };

//   TString option = GetOption();
//   toks = option.Tokenize(",");
//   TString RunID = ((TObjString*) toks->At(0))->GetString();
//   TString VolID = ((TObjString*) toks->At(1))->GetString();
//   TString ServerID = ((TObjString*) toks->At(3))->GetString();
    
  //std::stringstream OutputFile;
  OutputFile.str(std::string());
//   OutputFile << "selected_trg_run" << "_" << RunID <<"_"<<VolID;
    OutputFile << "selected_run" << "_" << RunID <<"_"<<VolID;
  if (atoi(ServerID) != 0) {OutputFile<<"_eliadeS"<<ServerID;};
  OutputFile << ".root";
//   std::cout <<"ServerID "<<ServerID<<" "<< OutputFile.str().c_str() <<std::endl;
   lastTime = 0;
   
  foutFile = new TFile (OutputFile.str().c_str(),"recreate"); 
   
   std::cout <<"ServerID "<<ServerID<<" "<< OutputFile.str().c_str() <<std::endl;
    start = std::clock();
    
    lastDelilaEvent.Time = 0;
    LastTriggerEvent.Time = 0;
    PulserEvent.Time = 0;
    DomainZeroEvent.Time = 0;
    LastBunchEvent.Time = 0;
//     nn_current_bunch = 0;
    delilaQu.clear();
    
    

 
    
    if (blTimeAlignement && blAddTriggerToQueue) {
        blAddTriggerToQueue = false;
        std::cout<<"blAddTriggerToQueue is set false because incompatibe with  blTimeAlignement true \n";
    }
    
    std::cout << "SlaveBegin() is done \n";

//    std::cout<<"event_length "<<event_length<<" ps \n";
//    std::cout<<"bunch_reset "<<bunch_reset <<" ps \n";
//    std::cout<<"Beta is "<<beta<<" % \n"; 

    
//     outputQu.clear();
}

Bool_t DelilaSelectorElifant::Process(Long64_t entry)
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
    nb_entries = GetEntries();

    int daq_ch = (fMod)*100+fChannel;
    DelilaEvent_.det_def = LUT_DELILA[daq_ch].detType;
    DelilaEvent_.channel = daq_ch;
 	hChannelHit->Fill(daq_ch);

    //Check that daq_ch is defined in LUT
      bool check_daq_ch = false;
      std::map<int, TDelilaDetector >::iterator it_daq_ch_ = LUT_DELILA.begin();
      for (; it_daq_ch_!= LUT_DELILA.end();++it_daq_ch_){
             if (LUT_DELILA[it_daq_ch_->first].ch == daq_ch){
               check_daq_ch = true;
//                std::cout<<"i am here daq_ch  "<< daq_ch<<" it_daq_ch_->first "<<LUT_DELILA[it_daq_ch_->first].ch<<"\n";
               continue;
             };
     };
     if (!check_daq_ch) return kTRUE;

     if (debug){std::cout<<"I am doing new entry, ch:"<< daq_ch << "\n";}

    if (LUT_DELILA.empty()){std::cout<<"LUT is empty \n"; return kTRUE;};//did not work well
    if (LUT_DELILA.find(daq_ch) == LUT_DELILA.end()){return kTRUE;};//did not work well
    
    if (debug){std::cout<<"I am doing new entry l.1043, ch:"<< daq_ch << "\n";}

    
    DelilaEvent_.domain = LUT_DELILA[daq_ch].dom;   
    int domain = DelilaEvent_.domain;

//     if (debug){std::cout<<"I am doing entry here, ch:"<< daq_ch << "\n";}
    
    DelilaEvent_.fEnergy = fEnergyLong;
    
    if ((DelilaEvent_.Energy_kev < LUT_DELILA[daq_ch].threshold)&&(DelilaEvent_.det_def < 9)) return kTRUE;
    
    DelilaEvent_.cs_domain = LUT_DELILA[daq_ch].cs_dom;
    DelilaEvent_.theta= LUT_DELILA[daq_ch].theta;
    DelilaEvent_.phi= LUT_DELILA[daq_ch].phi;
    
    
    hDomainHit->Fill(domain);
    hDetTypeHit->Fill(DelilaEvent_.det_def);   
    mDelila_raw->Fill(domain,DelilaEvent_.fEnergy);
    
    if (debug){std::cout<<"I am doing new entry l.1072, ch:"<< daq_ch << "\n";}

     //Check if the tree is time sorted
      DelilaEvent_.Time = fTimeStampFS;
//      DelilaEvent_.Time=fTimeStamp;
     double time_diff_last = DelilaEvent_.Time - lastDelilaTime;
     
     if (debug){std::cout<<"I am doing new entry l.1079, ch:"<< daq_ch << "\n";}
     
     //Check that the Tree is sorted in Time
     if (time_diff_last<0){std::cout<<"Warning time_diff_last: .Time  TTree may be not sorted by time"<< time_diff_last<<" \n";};
     if (DelilaEvent_.Time == 0) {hTimeZero->Fill(daq_ch);};
     hTimeSort->Fill(time_diff_last);
     
     if (debug){std::cout<<"I am doing new entry l.1084, ch:"<< daq_ch << "\n";}
     
     lastDelilaTime = DelilaEvent_.Time;     
     
     DelilaEvent_.Time-=LUT_TA[domain];
    
     if (DelilaEvent_.det_def == 9){//pulser
        CheckPulserAllignement(90);
        return kTRUE;
     }else if (DelilaEvent_.det_def == 98){
        return kTRUE;
     }else if (DelilaEvent_.det_def == 8){
          TreatNeutronSingle();
     }else if ((DelilaEvent_.det_def == 7) && has_detector["Elissa"]){
         TreatElissaSingle();
     }else if (((DelilaEvent_.det_def==4)||(DelilaEvent_.det_def==5)||(DelilaEvent_.det_def==6))&&has_detector["BGO"]){
         TreatBGOSingle();
     }else if ((DelilaEvent_.det_def == 3) && has_detector["LaBr"]) {
        TreatLaBrSingle();
     }else if ((DelilaEvent_.det_def == 1) && has_detector["HPGe"] ){
        TreatHpGeSingle();
     };
    
  if (debug){std::cout<<"I did TreatDelilaEvent_() \n";}
  
  EventBuilderPreTrigger();
  

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
    };
//    };
   
   nevents++;
//    nevents_reset++;
   return kTRUE;
}

void DelilaSelectorElifant::FillOutputTree(){
    
  std::deque<DelilaEvent>::iterator it_delila_ = delilaQu.begin();
  
  if (has_detector["HPGe"])   HPGeEvent   ->clear();
  if (has_detector["LaBr"])   LabrEvent   ->clear();
  if (has_detector["Elissa"]) ElissaEvent ->clear();;
  
  for (; it_delila_!= delilaQu.end();++it_delila_){
      switch (it_delila_->det_def){
          case 1:{
              if (has_detector["HPGe"]){
                  hpge_tree_event.SetEnergy(it_delila_->Energy_kev);
                  hpge_tree_event.SetEnergyDC(it_delila_->EnergyDC);
                  hpge_tree_event.SetChannel(it_delila_->channel);
                  hpge_tree_event.SetDomain(it_delila_->domain);
                  HPGeEvent->push_back(hpge_tree_event);
              };
              break;
          }
          case 3:{
              if (has_detector["LaBr"]){
                  labr_tree_event.SetEnergy(it_delila_->Energy_kev);
                  labr_tree_event.SetEnergyDC(it_delila_->EnergyDC);
                  labr_tree_event.SetChannel(it_delila_->channel);
                  labr_tree_event.SetDomain(it_delila_->domain);
                  LabrEvent->push_back(labr_tree_event);
              };
              break;
          }
          case 7:{//ELISSA
               if (has_detector["Elissa"]){
                  elissa_tree_event.SetEnergy(it_delila_->Energy_kev);
                  elissa_tree_event.SetAmax(it_delila_->Amax);
                  elissa_tree_event.SetChannel(it_delila_->channel);
                  elissa_tree_event.SetDomain(it_delila_->domain);
                  ElissaEvent->push_back(elissa_tree_event);
               };
              break;
          }
          default : {
           
              break;   
          }
      };
   };
   
   
   if((has_detector["HPGe"]   && HPGeEvent->size()>0)    ||
      (has_detector["LaBr"]   && LabrEvent->size()>0)    ||
      (has_detector["Elissa"] && ElissaEvent->size()>0))  {
        outputTree->Fill();
    };
   
   
   
  outputTree->Fill();
  return;
}


void DelilaSelectorElifant::TreatFold(int det)
{
   if (delilaQu.empty())return;
//      Int_t fold_size = foldQu.size();
     int nfold = 0;
     std::deque<DelilaEvent>::iterator it_delila_ = delilaQu.begin();

     for (; it_delila_!= delilaQu.end();++it_delila_){
         if (it_delila_->det_def == det) nfold++;
        };
         
    hBunchFold->Fill(nfold);
    it_delila_ = delilaQu.begin();
     
     for (; it_delila_!= delilaQu.end();++it_delila_){
         if (it_delila_->det_def == det) {
             it_delila_->fold = nfold;
//              mFoldEnergy->Fill(nfold,it_delila_->Energy_kev);
//              if (blOutTree) {
//                 DelilaEvent_Treated = *it_delila_;                
//                 outputTree->Fill();
//              };
            };
        };
    nfold = 0;
}

// bool check_gamma_det(std::deque<DelilaEvent>::iterator it_check_) {return ((it_check_->det_def != 1 )&&(it_check_->det_def != 2 )&&(it_check_->det_def != 3 )&&(it_check_->det_def != 7 ));};

void DelilaSelectorElifant::TreatGammaGammaCoinc()
{
   if (delilaQu.empty())return;
   std::map<int, int> nmult;
   std::map<UInt_t, std::string>  ::iterator it_mult_ =  gg_coinc_id.begin();

   for(;it_mult_!=gg_coinc_id.end();++it_mult_) {
        if ((it_mult_->first != 11 )&&(it_mult_->first != 13 )&&(it_mult_->first != 33 )) continue;
        nmult[it_mult_->first] = 0;   
   };
   
   std::deque<DelilaEvent>::iterator it1_= delilaQu.begin();
   std::deque<DelilaEvent>::iterator it2_= delilaQu.begin();
//    std::deque<DelilaEvent>::iterator it3_= delilaQu.begin();
   
   std::deque<DelilaEvent>::iterator it_tmp_= delilaQu.begin();
     for (; it1_!= delilaQu.end();++it1_){
          if ((it1_->det_def != 1 )&&(it1_->det_def != 3 )) continue;
           it2_= delilaQu.begin();
         for (; it2_!= delilaQu.end();++it2_){
             if (it1_ == it2_) continue;              
             if ((it2_->det_def != 1 )&&(it2_->det_def != 3 )) continue;
             int coinc_id = GetCoinc_det_def(it1_->det_def, it2_->det_def);

             //Check that daq_ch is defined in LUT
            std::map<UInt_t, Float_t> ::iterator it_c_gates_ = coinc_gates.find(coinc_id);
            if(it_c_gates_ == coinc_gates.end()) continue;

            if ((it1_->det_def == 1)&&(it2_->det_def == 3)){//ge-labr the time of labr is first
                it_tmp_ = it1_;
                it1_ = it2_;
                it2_ = it_tmp_;
            };
            
            if ((it1_->det_def == 7)&&(it2_->det_def == 3)){//ge-labr the time of labr is first
                it_tmp_ = it1_;
                it1_ = it2_;
                it2_ = it_tmp_;
            };
            double_t time_diff_gg = it2_->Time - it1_->Time;

            mGG_time_diff[coinc_id]->Fill(it2_->domain,time_diff_gg);
            hCoincID->Fill(coinc_id);
               double_t delta_theta = it1_->theta - it2_->theta;
            if (abs(time_diff_gg) < coinc_gates[coinc_id]){
                  mGG[coinc_id]->Fill(it1_->Energy_kev, it2_->Energy_kev);
                  nmult[coinc_id]++;
              };
        };
//                     std::cout<<"line 1268 \n";
     };
        
   it_mult_ =  gg_coinc_id.begin();
   for(;it_mult_!=gg_coinc_id.end();++it_mult_) hMult[it_mult_->first]->Fill(nmult[it_mult_->first]);   
//    std::cout<<"line 1275 \n";     
};


void DelilaSelectorElifant::TreatDelilaEvent()
{
    
    UShort_t daq_ch = DelilaEvent_.channel;
    UShort_t domain = DelilaEvent_.domain;
    
//     DelilaEvent_.Energy_kev = CalibDet(DelilaEvent_.fEnergy, daq_ch);
    double costheta = TMath::Cos(LUT_DELILA[daq_ch].theta);
    DelilaEvent_.EnergyDC = DelilaEvent_.Energy_kev*(1./sqrt(1 - beta*beta) * (1 - beta*costheta));
    
//     hDelila[DelilaEvent_.det_def]->Fill(DelilaEvent_.Energy_kev);
    mDelila->Fill(domain,DelilaEvent_.Energy_kev);
//     mDelila_long->Fill(domain,DelilaEvent_.Energy_kev);
    mDelilaDC->Fill(domain,DelilaEvent_.EnergyDC);
//     mDelilaDC_long->Fill(domain,DelilaEvent_.EnergyDC);
    
//     mThetaPhi->Fill(DelilaEvent_.theta, DelilaEvent_.phi);
    
    

    return;
}


void DelilaSelectorElifant::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void DelilaSelectorElifant::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the clieint coreID, nt, it can be used to present
   // the results graphically or save the results to file.
    
  std::cout<<"I  will terminate soon... "<<std::endl;  

  TIter iter(fOutput);
      

  std::cout << std::endl << "Finished after " << round(duration / 60) << ":"
	    << std::setw(2) << round(((int) duration) % 60) << std::setw(8)
	    << ". Write spectra to file \n" 
        << " A box of valpolicella can be sent to \n"
        << " Strada Reactorului 30, Magurele 077125 \n"
        << " office 420 and/or 418 \n" << std::endl;


      TObject *obj; 
      
      foutFile->mkdir("GammaGamma","GammaGamma");
      foutFile->mkdir("long","long");
      foutFile->mkdir("Energy_time_diff","Energy_time_diff");
      
      if (has_detector["Elissa"]) foutFile->mkdir("AmaxEnergy","AmaxEnergy");
      
      outputTree->Write();
      
       while ((obj = iter())) {
           TString name = obj->GetName();
//             std::cout<<"name "<<name<<std::endl;

           if(name.Contains("mDelila_long")){
           foutFile->cd(Form("%s:/", OutputFile.str().c_str()));
           }else if(name.Contains("mEnergy_time_diff")){
               foutFile->cd(Form("%s:/Energy_time_diff", OutputFile.str().c_str()));      
           }else if((name.Contains("mAmaxEnergy_dom")) && (has_detector["Elissa"]) ){
               foutFile->cd(Form("%s:/AmaxEnergy", OutputFile.str().c_str()));      
           } else if (name.Contains("mgg_")){
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
            }else if (obj->IsA()->InheritsFrom(TH3F::Class())){
            TH3 *h3 = (TH3*)obj;
                if (h3->GetEntries()>0) obj->Write();
            }
        };
        
//        outputTree->Write();
       foutFile->Close();
       delilaQu.clear();
       
//        std::cout<<"trigger_cnt "<<trigger_cnt<<std::endl;
   

}


int DelilaSelectorElifant::GetCoincTimeCorrection(int dom1, int dom2)
{
 int coin_id = GetCoincID(dom1, dom2);
//  if (dom1 >= dom2 ){coin_id = dom1*100+dom2;}else {coin_id = dom2*100+dom1;};
    
//  int coin_id = dom1*100+dom2;
 int time_corr = 0;
 std::map<int, int >::iterator it = LUT_TA.find(coin_id);
 if(it != LUT_TA.end()){time_corr = LUT_TA[coin_id];};
 return time_corr;
}

void DelilaSelectorElifant::cs()
{
    if (delilaQu.empty()) return;
    
    //BGO-GAMMA forward
    std::map<int, double> last_bgo_time;
    
    std::map<int, TDelilaDetector >::iterator it_lut_ = LUT_DELILA.begin();
    for (; it_lut_ != LUT_DELILA.end(); ++it_lut_){
        last_bgo_time[LUT_DELILA[it_lut_->first].cs_dom] = -1;
    }
    
    
    double time_diff_bgo = -1;
    std::deque<DelilaEvent>  ::iterator it_ev__ = delilaQu.begin();
    
    for (; it_ev__ != delilaQu.end(); ++it_ev__){
        int cs_dom = (*it_ev__).cs_domain;
        int det = (*it_ev__).det_def;

        switch (det)
        {
            case 5: {
                last_bgo_time[cs_dom] = (*it_ev__).Time; 
                break;
            };
            case 1: case 3: {
               if (last_bgo_time[cs_dom] != -1)
               {
                time_diff_bgo =  (*it_ev__).Time - last_bgo_time[cs_dom];// GetCoincTimeCorrection(it_ev__->domain,
//                  time_diff_bgo =  (*it_ev__).Time - last_bgo_time[cs_dom] - LUT_DELILA[it_ev__->channel].bgo_time_offset ;// GetCoincTimeCorrection(it_ev__->domain,
                 mTimeDiffCS ->Fill(cs_dom, time_diff_bgo);
//                  std::cout<<" coinc_gates "<<det<<" "<< coinc_gates[det*10+5] <<" "<<time_diff_bgo << " \n";
                 if (abs(time_diff_bgo) < coinc_gates[det*10+5]) //10000)
                 {
                    (*it_ev__).CS = 1;  std::cout<<" here cs \n";
                 }else{
                    last_bgo_time[cs_dom] = -1;
                 };                 
               };
            };
            default: break;
        };
        
    };

     //GAMMA-BGO backward
    time_diff_bgo = -1;
    it_lut_ = LUT_DELILA.begin();
    for (; it_lut_ != LUT_DELILA.end(); ++it_lut_){
        last_bgo_time[LUT_DELILA[it_lut_->first].cs_dom] = -1;
    }
    time_diff_bgo = 0;
    it_ev__ = delilaQu.end();
    
    for (; it_ev__ != delilaQu.begin(); --it_ev__){
        int cs_dom = (*it_ev__).cs_domain;
        int det = (*it_ev__).det_def;
        
        switch (det)
        {
            case 5: {
                last_bgo_time[cs_dom] = (*it_ev__).Time; 
                break;
            };
            case 1: case 3: {
               if ((last_bgo_time[cs_dom] == -1)||((*it_ev__).CS == 1))
               {
                 break;  
               }else{
                 time_diff_bgo =  (*it_ev__).Time - last_bgo_time[cs_dom];
//                  time_diff_bgo =  (*it_ev__).Time - last_bgo_time[cs_dom]  - LUT_DELILA[it_ev__->channel].bgo_time_offset;
                 mTimeDiffCS ->Fill(cs_dom, time_diff_bgo); 
//                  std::cout<<"time_diff_bgo back "<<time_diff_bgo<<" id "<<det*10+5 <<" "<<coinc_gates[det*10+5]<<"\n";
                 if (abs(time_diff_bgo) < coinc_gates[det*10+5]) //10000)
                 {
                      (*it_ev__).CS = 1;   std::cout<<" here cs \n";
                 }else{
                     last_bgo_time[cs_dom] = -1;
                 };
                 
               };
            };
            default: break;
        };
        
    };
    
    it_ev__ = delilaQu.begin();
    
    for (; it_ev__ != delilaQu.end(); ++it_ev__){
         if (((*it_ev__).CS == 0)&&((*it_ev__).det_def == 3)) hDelilaCS[(*it_ev__).det_def]->Fill((*it_ev__).Energy_kev);
         if (((*it_ev__).CS == 0)&&((*it_ev__).det_def == 1)) hDelilaCS[(*it_ev__).det_def]->Fill((*it_ev__).Energy_kev);
    };
    
};


 int DelilaSelectorElifant::GetCoinc_det_def(int det_def1, int det_def2)
 {
     int id=0;
     if (det_def1<=det_def2) {id = det_def1*10+det_def2;}
         else {id = det_def2*10+det_def1; };
     return id;
 }


int DelilaSelectorElifant::GetCoincID(int dom1, int dom2){
    int id=0;
    if (dom1<=dom2) {id = dom1*100+dom2;}
        else {id = dom2*100+dom1; };
    return id;
}


void DelilaSelectorElifant::CheckPulserAllignement(int zero_dom)
{
   double time_diff_pulser;
   int cur_dom = DelilaEvent_.domain;
   
   if (cur_dom != zero_dom){
       time_diff_pulser =  DelilaEvent_.Time - PulserEvent.Time - GetCoincTimeCorrection(PulserEvent.domain, DelilaEvent_.domain);;
       mPulser0TimeDiff->Fill(cur_dom, time_diff_pulser);
//        std::cout<<" CheckPulserAllignement "<<cur_dom <<"  DelilaEvent_.domain "<<DelilaEvent_.Time<<" PulserEvent.Time "<<PulserEvent.Time<<" dt ="<<time_diff_pulser<<" \n";
//        if (DelilaEvent_.domain == 57) std::cout<<" time_diff_pulser #50 - #57 "<<time_diff_pulser<<" \n";
//        if (DelilaEvent_.domain == 51) std::cout<<" time_diff_pulser #50 - #51 "<<time_diff_pulser<<" \n";
   }else PulserEvent = DelilaEvent_;
   return;
}

void DelilaSelectorElifant::PrintDelilaEvent(DelilaEvent &ev_)
{
    std::cout<<" === DelilaEvent_ === \n"<<
    "  ev_.domain "<<ev_.domain<<" \n"<<
    "  ev_.Energy_kev "<<ev_.Energy_kev<<" \n"<<
    "  ev_.fold "<<ev_.fold<<" \n";
//     "  ============= \n";
}

bool DelilaSelectorElifant::TriggerDecision()
{
   if (det_def_trg == -1) return false;
   if (channel_trg == -1) return false; 

   if (det_def_trg > 0) return (DelilaEvent_.det_def == det_def_trg/1);
   
   return (DelilaEvent_.domain == channel_trg/1);
};

void DelilaSelectorElifant::TimeAlignementTrigger()
{
     std::deque<DelilaEvent>::iterator it_= delilaQu.begin();
     double time_diff_temp = 0;

     for (; it_!= delilaQu.end();++it_){
           time_diff_temp = delilaQu.front().Time - TriggerTimeFlag;
           mTimeCalibDomain0->Fill((*it_).domain, time_diff_temp);
     };
         
    return;
};
void DelilaSelectorElifant::TreatLaBrSingle()
{
    UShort_t daq_ch = DelilaEvent_.channel;
    UShort_t domain = DelilaEvent_.domain;
    
    DelilaEvent_.Energy_kev = CalibDet(DelilaEvent_.fEnergy, daq_ch);

    
    double costheta = TMath::Cos(LUT_DELILA[daq_ch].theta);
    if (beta >0) DelilaEvent_.EnergyDC = DelilaEvent_.Energy_kev*(1./sqrt(1 - beta*beta) * (1 - beta*costheta));
    
    hDelila0[DelilaEvent_.det_def]->Fill(DelilaEvent_.Energy_kev); 
    mDelila->Fill(domain,DelilaEvent_.Energy_kev);
    mDelilaDC->Fill(domain,DelilaEvent_.EnergyDC);

    return;
}

void DelilaSelectorElifant::TreatHpGeSingle()
{
    UShort_t daq_ch = DelilaEvent_.channel;
    UShort_t domain = DelilaEvent_.domain;
    
    DelilaEvent_.Energy_kev = CalibDet(DelilaEvent_.fEnergy, daq_ch);

    double costheta = TMath::Cos(LUT_DELILA[daq_ch].theta);
    if (beta >0) DelilaEvent_.EnergyDC = DelilaEvent_.Energy_kev*(1./sqrt(1 - beta*beta) * (1 - beta*costheta));
    
    mDelila->Fill(domain,DelilaEvent_.Energy_kev);
    mDelilaDC->Fill(domain,DelilaEvent_.EnergyDC);
    hDelila0[DelilaEvent_.det_def]->Fill(DelilaEvent_.Energy_kev); 
     
}

void DelilaSelectorElifant::TreatNeutronSingle()
{
    DelilaEvent_.fEnergyShort = fEnergyShort;
    mShortLong->Fill(DelilaEvent_.fEnergy, DelilaEvent_.fEnergyShort);
//     std::cout<<DelilaEvent_.det_def<<" "<<" Long "<<DelilaEvent_.fEnergy<<" short "<< DelilaEvent_.fEnergyShort <<" n/g "<<(DelilaEvent_.fEnergy - DelilaEvent_.fEnergyShort)*1.0/DelilaEvent_.fEnergy<<"\n";
    mNeutron->Fill(DelilaEvent_.fEnergy*1.0, ((DelilaEvent_.fEnergy*1.0 - DelilaEvent_.fEnergyShort*1.0)/DelilaEvent_.fEnergy*1.0));
    hDelila0[DelilaEvent_.det_def]->Fill(DelilaEvent_.Energy_kev); 
}

void DelilaSelectorElifant::TreatElissaSingle()
{
    vector<float> data_fil=DelilaSelectorElifant::trapezoidal(fSignal, 1000, 20, 0);//to be checked if signals do not present
    float trap_max=*max_element(data_fil.begin(),data_fil.end());
    float trap_min=*min_element(data_fil.begin(),data_fil.end());
    DelilaEvent_.Amax=trap_max-trap_min;
    
    DelilaEvent_.Energy_kev = CalibDet(DelilaEvent_.fEnergy, daq_ch);

    
    hAmax->Fill(DelilaEvent_.Amax);
//     hAmax->Fill(DelilaEvent_.Amax/DelilaEvent_.Energy_kev);
    mAmaxEnergy->Fill(DelilaEvent_.Energy_kev,DelilaEvent_.Amax);
    hDelila0[DelilaEvent_.det_def]->Fill(DelilaEvent_.Energy_kev); 
    mElissa->Fill(DelilaEvent_.domain, DelilaEvent_.Energy_kev);

    
    if (blFillAmaxEnergyDom) mAmaxEnergyDom[DelilaEvent_.domain]->Fill(DelilaEvent_.Energy_kev,DelilaEvent_.Amax);
    
    if (blDebugElissa) cout << trap_max << " " << trap_min <<" Amax "<<DelilaEvent_.Amax<< endl;   
    
}

void DelilaSelectorElifant::TreatSolarLaBrCoinc()
{
    
}

void DelilaSelectorElifant::CheckPreQu()
{
    if (delilaPreQu.empty()) return;
    
    std::deque<DelilaEvent>::iterator it_= delilaPreQu.begin();
    double time_diff_temp;
    
//     for (; it_!= delilaPreQu.end();++it_){
    while ((it_!= delilaPreQu.end())&&(!delilaPreQu.empty())){
        time_diff_temp = it_->Time - DelilaEvent_.Time;  
       if (abs(time_diff_temp) > pre_event_length) delilaPreQu.erase(it_);
       ++it_;
//        if (delilaPreQu.empty()) break;
    };
}

void DelilaSelectorElifant::MovePreQu2Qu()
{
    if (delilaPreQu.empty()) return;
    
    if (!delilaQu.empty()) {
        std::cout<<"CheckPreQu delilaQu is not empty \n";
        delilaQu.clear();
    };
    
    std::deque<DelilaEvent>::iterator it_= delilaPreQu.begin();
    double time_diff_temp;
    
    for (; it_!= delilaPreQu.end();++it_){
           time_diff_temp = it_->Time - DelilaEvent_.Time;                      
           if (abs(time_diff_temp) < pre_event_length) delilaQu.push_back(*it_);
//            std::cout<<"time_diff_temp  "<< time_diff_temp << " \n";
    };
//     std::cout<<"PreQueSize is "<< delilaPreQu.size() << " \n";
//     std::cout<<"delilaQueSize is "<< delilaQu.size() << " \n";
    
    delilaPreQu.clear();
}

void DelilaSelectorElifant::EventBuilderSimple()
{
    
    if (!blIsTrigger) {
    if (TriggerDecision()) SetUpNewTrigger();
//        return kTRUE;
   }else {//if (blIsTrigger){
       
        double time_diff_trigger = DelilaEvent_.Time - LastTriggerEvent.Time;
        
        if (time_diff_trigger > event_length){//close event
            if (blTimeAlignement) TimeAlignementTrigger();
            if (blCS) cs();
            if (blGammaGamma) TreatGammaGammaCoinc();
            if (blFold) TreatFold(3);
            if (blOutTree) FillOutputTree();
            
            hdelilaQu_size->Fill(delilaQu.size());
//             outputTree->Fill();
            delilaQu.clear();
//             ElifantEvent->clear();
            blIsTrigger = false;
            
            if (TriggerDecision()) SetUpNewTrigger();          
            
        }else{
          
            delilaQu.push_back(DelilaEvent_);
            DelilaEvent_.trg = trigger_cnt;
            trigger_cnt++;
        };
        
    };
    
}

void DelilaSelectorElifant::EventBuilderPreTrigger()
{
//     std::cout<<" line 1685 \n";
    if (blIsWindow){//open
       double time_diff_trigger = DelilaEvent_.Time - TriggerTimeFlag;
       if (abs(time_diff_trigger) > event_length){//close event
//            std::cout<<" gates 33 "<<coinc_gates[33]<<"\n";
           if (blTimeAlignement)    TimeAlignementTrigger();
           if (blCS)                cs();
           if (blGammaGamma)        TreatGammaGammaCoinc();
           if (blFold)              TreatFold(3);
           if (blOutTree)           FillOutputTree();
           
           hdelilaQu_size->Fill(delilaQu.size());          
           
           delilaQu.clear();
           delilaPreQu.clear();
           blIsTrigger = false;
           blIsWindow = false;
           
           if (TriggerDecision()) SetUpNewTrigger();          
           
       }else{
//           std::cout<<" line 1705 \n";
           hDelila_single[DelilaEvent_.det_def]->Fill(DelilaEvent_.Energy_kev);

           delilaQu.push_back(DelilaEvent_);
           DelilaEvent_.trg = trigger_cnt;
       
           trigger_cnt++;
       };
       
   }else{//closed
       if (TriggerDecision()) {MovePreQu2Qu(); SetUpNewTrigger();}
       else{ 
//            CheckPreQu();
           DelilaEvent_.trg = trigger_cnt;
           delilaPreQu.push_back(DelilaEvent_);
    };
   };
}


void DelilaSelectorElifant::SetUpNewTrigger(){

//     MovePreQu2Qu();
    
    hTriggerTrigger->Fill(DelilaEvent_.Time - LastTriggerEvent.Time);
    hTriggerDomain->Fill(DelilaEvent_.domain);
    LastTriggerEvent = DelilaEvent_;
    
    if (blAddTriggerToQueue) delilaQu.push_back(DelilaEvent_);
    
    TriggerTimeFlag = DelilaEvent_.Time - pre_event_length;
    
//     std::cout<<" TriggerTimeFlag is "<<TriggerTimeFlag<<" DelilaEvent_.Time "<<DelilaEvent_.Time<<"  \n";  
    
    
    if (TriggerTimeFlag < 0){
      std::cout<<"SetUpNewTrigger() TriggerTimeFlag is < 0 \n";  
      TriggerTimeFlag = 0;
    };
    
    blIsTrigger = true;
    blIsWindow = true;

    trigger_cnt++;
    return;
}

void DelilaSelectorElifant::TreatACS()
{
    int gate = 20000;
    double_t time_diff_cs = 0;
    int cs_dom = DelilaEvent_.cs_domain;
    
    if (DelilaEvent_.det_def == 1) gate = 100000;

    if ((DelilaEvent_.det_def == 3)||(DelilaEvent_.det_def == 1)) {
        waitingQu_gamma[cs_dom].push_back(DelilaEvent_); 
        hDelila_single[DelilaEvent_.det_def]->Fill(DelilaEvent_.Energy_kev);
//         hDelilaDC[DelilaEvent_.det_def]->Fill(DelilaEvent_.EnergyDC);
//         hDelila_long[DelilaEvent_.det_def]->Fill(DelilaEvent_.EnergyCal);
//         hDelilaDC_long[DelilaEvent_.det_def]->Fill(DelilaEvent_.EnergyDC);
    }
    if (DelilaEvent_.det_def == 5) waitingQu_bgo[cs_dom].push_back(DelilaEvent_);
        
        if ((!waitingQu_gamma[cs_dom].empty())&&(!waitingQu_bgo.empty()))
        {
            std::deque<DelilaEvent>  ::iterator it_g__ = waitingQu_gamma[cs_dom].begin();
            std::deque<DelilaEvent>  ::iterator it_bgo__ = waitingQu_bgo[cs_dom].begin();
            for (; it_g__ != waitingQu_gamma[cs_dom].end();++it_g__){
                for (; it_bgo__ != waitingQu_bgo[cs_dom].end();++it_bgo__){   
                time_diff_cs = it_g__->Time - it_bgo__->Time ;//-  GetCoincTimeCorrection(it_g__->domain,it_bgo__->domain);
                mTimeCalibBGO->Fill(it_g__->domain, time_diff_cs);
                if (abs(time_diff_cs)<gate){
                    if (it_g__->CS == 1) continue;
                    it_g__->CS = 1;   
//                     it_g__->bgo_time_diff = time_diff_cs;
//                      if (it_g__->det_def == 1) {std::cout<<time_diff_cs<<" time_diff_cs \n";};
                    //mLaBr_BGO_time_diff->Fill(it_g__->domain, time_diff_cs);
                    }
                }
            }
        }
        
         if (!waitingQu_bgo[cs_dom].empty())
         {
             std::deque<DelilaEvent>  ::iterator it1_ = waitingQu_bgo[cs_dom].begin();
             for (; it1_ != waitingQu_bgo[cs_dom].end();)
             {
              if (abs(DelilaEvent_.Time - it1_->Time)>gate) {it1_=waitingQu_bgo[cs_dom].erase(it1_);}
                  else ++it1_;
             }
         };
         
         if (!waitingQu_gamma[cs_dom].empty())
         {
             std::deque<DelilaEvent>  ::iterator it2_ = waitingQu_gamma[cs_dom].begin();
             for (; it2_ != waitingQu_gamma[cs_dom].end();)
             {
              if (abs(DelilaEvent_.Time - it2_->Time)>gate) 
              {
                 // output_pQu.push(*it2_);                  
                  if (it2_->CS ==0){
                      hDelilaCS[it2_->det_def]->Fill(it2_->Energy_kev);
//                       hDelilaCS_long[it2_->det_def]->Fill(it2_->EnergyCal);
//                       mDelilaCS->Fill(it2_->cs_domain, it2_->EnergyCal);
//                       gamma_gamma_cs(*it2_);
//                       hDelilaCS_DC[it2_->det_def]->Fill(it2_->EnergyDC);
//                       hDelilaCS_DC_long[it2_->det_def]->Fill(it2_->EnergyDC);
//                       mDelilaCS_DC->Fill(it2_->cs_domain, it2_->EnergyDC);
                };
                  
                 if (blOutTree) {/*DelilaEventCS = *it2_; outputTree->Fill();*/};
                  
                  it2_=waitingQu_gamma[cs_dom].erase(it2_);
              }
              else ++it2_;
             };
         };

};


// add by saka
std::vector<float> DelilaSelectorElifant::trapezoidal(short wave[],int length, int L, int G)//L = 20; G = 0
{
  vector<float> ret;
  for(int j=0,n=length-(2*L+G); j<n; j++){
    float tmp=0;
    for(int k=0;k<L;k++){
      tmp+=(float)wave[j+k]*(-1/(float)L);
    }
    for(int k=0;k<G;k++){
      tmp+=wave[j+k+G]*0;
    }
    for(int k=0;k<L;k++){
      tmp+=(float)wave[j+k+G+L]*(1/(float)L);
    }
    ret.push_back(tmp);
  }
  return ret;

}
// add by saka

void DelilaSelectorElifant::TreatBGOSingle()
{
    hDelila0[DelilaEvent_.det_def]->Fill(DelilaEvent_.Energy_kev);
}










