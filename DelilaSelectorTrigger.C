#define DelilaSelectorTrigger_cxx
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
// root> T->Process("DelilaSelectorTrigger.C")
// root> T->Process("DelilaSelectorTrigger.C","some options")
// root> T->Process("DelilaSelectorTrigger.C+")
//


#include "DelilaSelectorTrigger.h"
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TObjString.h>
#include <unordered_set>
#include <iomanip>      // std::setwsorted
using namespace std;


////////////////////////////////Please, modify if needed////////////////////////////////////////////
bool blGammaGamma = true;
bool blGammaGammaCS;// = true;
bool blCS = true;
bool blIsTrigger = true; //trigger signal is present
// bool blTriggerMode = false;//to reset the queues after each trigger - to be implemented
bool blOutTree = true;
// bool blTimeEnergy = false;


////////////////////////////////Please, DO NOT modify ////////////////////////////////////////////
int addBackMode = 0; //0 - no addback; 1- addback;//not in use for ELIFANT
// bool blTrigger = false;
// double_t trigger_shift = 112969; //ps
// double_t bunch_distance = 400000;//ps
// bool blBunchIsValidated = true;
// bool blBunch = true;
bool blFirstTrigger = false;
// bool blPulserTimeAllignement = false; // for normal use should be false

bool debug = false;

ULong64_t trigger_cnt = 0;
ULong64_t trigger_events = 0;
// Int_t nfold = 0;
Int_t nn_current_bunch = 0;

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


void DelilaSelectorTrigger::Read_ELIADE_LookUpTable() {
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
      int threshold = 1;
      std::istringstream is(oneline);
      if (debug) std::cout << is.str().c_str() << std::endl;
//       is >> curDet.ch >> curDet.dom >> curDet.theta >> curDet.phi >> curDet.TimeOffset >> curDet.upperThreshold;
      is >> curDet.ch >> curDet.dom >> curDet.detType >> curDet.serial >> curDet.TimeOffset >> curDet.theta >> curDet.phi >> curDet.threshold >> curDet.cs_dom;
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


void DelilaSelectorTrigger::Read_TimeAlignment_LookUpTable() {
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

void DelilaSelectorTrigger::Read_TimeAlignment_Trigger() {
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

void DelilaSelectorTrigger::Read_Confs() {
  std::cout << "I am Reading Your Configuration parameters ... ";
 
  char* pLUT_Path;
  pLUT_Path = getenv ("DELILA_LUT");
  if (pLUT_Path!=NULL)
    printf ("The LookUpTable path is: %s \n",pLUT_Path);


  std::stringstream LUTFile;
  LUTFile << pLUT_Path <<"/"<<"LUT_CONF.dat";
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
      int coinc_id = 0; Float_t value = 0;
      is >> coinc_name>> coinc_id >> value;
      
      
      switch (coinc_id){
          case 1000: {
              bunch_length = value;
//               std::cout<<"bunch_length "<<bunch_length<<" ps \n";
              break;
          }
          case 1001:
              {
              bunch_reset = value;
//               std::cout<<"bunch_reset "<<bunch_reset<<" ps \n";
              break;
          }
          case 1111:{
              beta = value;
//               std::cout<<"Beta is "<<beta<<" % \n";
              break;
          }
          default:
          {
              //LUT_CONF[coinc_id] = value;
              coinc_gates[coinc_id] = value;
              std::cout<<coinc_name<<" coin_id " << coinc_id <<" value "<<value <<" ps \n";
              break;
          };
      };
      
      
      
  }
  lookuptable.close();
  }
  std::cout << " done" << std::endl;
}



void DelilaSelectorTrigger::Print_ELIADE_LookUpTable()
{
    std::cout<<"Print_ELIADE_LookUpTable \n";		
    std::map<unsigned int, TDelilaDetector > ::iterator it__ = LUT_DELILA.begin();
    for (; it__ != LUT_DELILA.end(); ++it__) {
     // is >> curDet.ch >> curDet.dom >> theta >> phi >> curDet.TimeOffset >> curDet.upperThreshold;
// 	std::cout<<" Ch "<<LUT_DELILA[it__->first].ch<<" Dom "<< LUT_DELILA[it__->first].dom<<" "<< LUT_DELILA[it__->first].theta<<" "<< LUT_DELILA[it__->first].phi <<" offset "<< LUT_DELILA[it__->first].TimeOffset<<" Thr "<< LUT_DELILA[it__->first].threshold<<" serial "<<LUT_DELILA[it__->first].serial<<" theta" <<LUT_DELILA[it__->first].theta<<" phi "<<LUT_DELILA[it__->first].phi <<" cs_dom: "<<LUT_DELILA[it__->first].cs_dom<<" pol_order: " <<LUT_DELILA[it__->first].pol_order <<std::endl;
        std::cout<<" Ch: "<<LUT_DELILA[it__->first].ch<<" Dom: "<< LUT_DELILA[it__->first].dom<<" type: "<< LUT_DELILA[it__->first].detType<<" serial: "<< LUT_DELILA[it__->first].serial <<" offset: "<< LUT_DELILA[it__->first].TimeOffset<<" Theta: "<< LUT_DELILA[it__->first].theta<<" phi: "<<LUT_DELILA[it__->first].phi<<" thr:" <<LUT_DELILA[it__->first].threshold<<" cs_dom "<<LUT_DELILA[it__->first].cs_dom<<" pol_order: "<<LUT_DELILA[it__->first].pol_order<<std::endl;  
    }
};

void DelilaSelectorTrigger::Print_TimeAlignment_LookUpTable()
{
    std::cout<<"Print_TimeAlignment_LookUpTable \n";		
    std::map<int, int > ::iterator it__ = LUT_TA.begin();
    for (; it__ != LUT_TA.end(); ++it__) {
     // is >> curDet.ch >> curDet.dom >> theta >> phi >> curDet.TimeOffset >> curDet.upperThreshold;
	std::cout<<" coinc_id "<<it__->first<<" time_corr "<< it__->second<<std::endl;
    }
};



float DelilaSelectorTrigger::CalibDet(float val, int ch)
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

void DelilaSelectorTrigger::Begin(TTree * tree)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

  
  TString option = GetOption();
  toks = option.Tokenize(",");
  TString RunID = ((TObjString*) toks->At(0))->GetString();
//   double beta = ((TObjString*) toks->At(2))->GetString().Atof();
  beta = ((TObjString*) toks->At(2))->GetString().Atof();
  addBackMode = atoi(((TObjString*) toks->At(3))->GetString());
  std::cout << "addBackMode  " << addBackMode <<std::endl;
  
//   std::cout<<"Beta is "<<beta<<" % \n";

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
//   fChain->SetBranchAddress("ChargeShort", 	&fEnergyShort, 	               &b_energ_short);
//   fChain->SetBranchAddress("ChargeShort", 	&DelilaEvent.fEnergy, 	       &b_energ_short);
  fChain->SetBranchAddress("FineTS", 	    &fTimeStampFS, 	               &b_tstmp_fine);
  fChain->SetBranchAddress("TimeStamp",     &fTimeStamp,                   &b_tstmp); 
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

  bunch_length = 40000;
  bunch_reset = 38000;
  
  Read_ELIADE_LookUpTable();
  Read_TimeAlignment_LookUpTable();
  Read_TimeAlignment_Trigger();
  Read_Confs();
//   Print_TimeAlignment_LookUpTable();

  
   
}




void DelilaSelectorTrigger::SlaveBegin(TTree * /*tree*/)
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
   
   
   hdelilaQu_size = new TH1F("hdelilaQu_size", "hdelilaQu_size", 20,0,20);
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
   
   coinc_gates[11]=600000;//hpge-hpge 40 ns
   coinc_gates[15]=100000;//hpge-bgo 40 ns
   
   coinc_gates[33]=10000;//labr-labr 40 ns
   coinc_gates[35]=10000;//labr-bgo 40 ns
   
   coinc_gates[13]=1000000;//hpge-bgo 40 ns
   
   gg_coinc_id[11]="mgg_hpge_hpge";
   gg_coinc_id[13]="mgg_labr_hpge";
   gg_coinc_id[33]="mgg_labr_labr";
   
   detector_name[1]="HPGe";
   detector_name[3]="LabBr";
   detector_name[5]="BGO";
   detector_name[9]="pulser";

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
   
//    mDomainTimeDiff_trigger = new TH2F("mDomainTimeDiff_trigger", "mDomainTimeDiff_trigger", max_domain, 0, max_domain, 100e3, -1e6, 99e6);
   mDomainTimeDiff_trigger = new TH2F("mDomainTimeDiff_trigger", "mDomainTimeDiff_trigger", max_domain, 0, max_domain, 100e3, 0, 1e9);
//    mDomainTimeDiff_trigger = new TH2F("mDomainTimeDiff_trigger", "mDomainTimeDiff_trigger", max_domain, 0, max_domain, 3e2, 0, 3e5);
   mDomainTimeDiff_trigger->GetXaxis()->SetTitle("domain");
   mDomainTimeDiff_trigger->GetYaxis()->SetTitle("Time diff, ps");
   fOutput->Add(mDomainTimeDiff_trigger);
   
   //     mDomainTimeDiff_trigger_TA = new TH2F("mDomainTimeDiff_trigger_TA", "mDomainTimeDiff_trigger_TA", max_domain, 0, max_domain, 10e3, -1e6, 9e6);
   mDomainTimeDiff_trigger_TA = new TH2F("mDomainTimeDiff_trigger_TA", "mDomainTimeDiff_trigger_TA", max_domain, 0, max_domain, 1e4, -1e5, 9e5);
   mDomainTimeDiff_trigger_TA->GetXaxis()->SetTitle("domain");
   mDomainTimeDiff_trigger_TA->GetYaxis()->SetTitle("Time diff, ps");
   fOutput->Add(mDomainTimeDiff_trigger_TA);
   
   //mDomainTimeDiff_bunch = new TH2F("mDomainTimeDiff_bunch", "mDomainTimeDiff_bunch", max_domain, 0, max_domain, 10e3, -1e6, 9e6);
   mDomainTimeDiff_bunch = new TH2F("mDomainTimeDiff_bunch", "mDomainTimeDiff_bunch", max_domain, 0, max_domain, 1e4, -1e5, 9e5);
//    mDomainTimeDiff_bunch = new TH2F("mDomainTimeDiff_bunch", "mDomainTimeDiff_bunch", max_domain, 0, max_domain, 3e2, 0, 3e5);
   mDomainTimeDiff_bunch->GetXaxis()->SetTitle("domain");
   mDomainTimeDiff_bunch->GetYaxis()->SetTitle("Time diff, ps");
   fOutput->Add(mDomainTimeDiff_bunch);
   
   hTriggerTrigger = new TH1F("hTriggerTrigger", "hTriggerTrigger", 100e4, 0, 100e6);
   hTriggerTrigger->GetYaxis()->SetTitle("Counts");
   hTriggerTrigger->GetXaxis()->SetTitle("Time diff, ps");
   hTriggerTrigger->SetTitle("Time between two trigger signals");
   fOutput->Add(hTriggerTrigger);
   
   hBunchBunch = new TH1F("hBunchBunch", "hBunchBunch", 10e4, -1e6, 9e6);
   hBunchBunch->GetYaxis()->SetTitle("Counts");
   hBunchBunch->GetXaxis()->SetTitle("Time diff, ps");
   hBunchBunch->SetTitle("Time between two Bunch signals");
   fOutput->Add(hBunchBunch);
   
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
   
   mTimeDiffCS = new TH2F("mTimeDiffCS", "mTimeDiffCS", 100, 0, 100, 4e4, -2e6, 2e6);
   mTimeDiffCS->GetXaxis()->SetTitle("domin");
   mTimeDiffCS->GetYaxis()->SetTitle("ps");
   mTimeDiffCS->SetTitle("LaBr_BGO time diff");
   fOutput->Add(mTimeDiffCS);
   
   mPulser0TimeDiff = new TH2F("mPulser0TimeDiff", "mPulser0TimeDiff", 100, 0.5, 100.5, 6e3, -3e6, 3e6);
   mPulser0TimeDiff->GetXaxis()->SetTitle("domain");
   mPulser0TimeDiff->GetYaxis()->SetTitle("ps");
   mPulser0TimeDiff->SetTitle("PulsePulser time diff");
   fOutput->Add(mPulser0TimeDiff);//time_diff relevant to the 1st channel (101), i.e. ch 101 is a trigger
   
/*   std::map<unsigned int, TDelilaDetector > ::iterator it_lut_ = LUT_DELILA.begin();
    for (; it_lut_ != LUT_DELILA.end(); ++it_lut_) {
        if ((LUT_DELILA[it_lut_->first].detType == 1)||(LUT_DELILA[it_lut_->first].detType == 3)){
            int dom = LUT_DELILA[it_lut_->first].dom;
            mEnergy_time_diff[dom] = new TH2F(Form("mEnergy_time_diff_dom%i",dom), Form("mEnergy_time_diff_dom%i",dom), 4096, -0.5, 65535.5, 10e3, -1e5, 9e5);
            mEnergy_time_diff[dom] ->GetXaxis()->SetTitle("16 keV / bin");
            mEnergy_time_diff[dom] ->GetYaxis()->SetTitle("100 ps / bin");
            mEnergy_time_diff[dom]->SetTitle(" Time Diff Event.Time - Bunch.Time ");
            fOutput->Add(mEnergy_time_diff[dom]); 
        };
    }*/;

  TString option = GetOption();
  toks = option.Tokenize(",");
  TString RunID = ((TObjString*) toks->At(0))->GetString();
  TString VolID = ((TObjString*) toks->At(1))->GetString();
  TString ServerID = ((TObjString*) toks->At(3))->GetString();
    
  //std::stringstream OutputFile;
  OutputFile.str(std::string());
  OutputFile << "selected_trg_run" << "_" << RunID <<"_"<<VolID;
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
    nn_current_bunch = 0;
    delilaQu.clear();
    
   std::cout<<"bunch_length "<<bunch_length<<" ps \n";
   std::cout<<"bunch_reset "<<bunch_reset <<" ps \n";
   std::cout<<"Beta is "<<beta<<" % \n"; 
//    Print_ELIADE_LookUpTable();
}

Bool_t DelilaSelectorTrigger::Process(Long64_t entry)
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

    int mod = DelilaEvent.fMod;
	int ch = DelilaEvent.fChannel;
    int daq_ch = (mod)*100+ch;
    DelilaEvent.det_def = LUT_DELILA[daq_ch].detType;
    DelilaEvent.channel = daq_ch;
 	hChannelHit->Fill(daq_ch);
    
     //Waiting for the first trigger to come
    if ((DelilaEvent.det_def == 99)&&(!blFirstTrigger)) {
       blFirstTrigger = true;
       LastTriggerEvent = DelilaEvent;
       LastBunchEvent = DelilaEvent;// + n_bunches * 400000;
       trigger_cnt++;
       std::cout<<"I found the first trigger signal! Cool! \n";//<<  LastTriggerEvent.Time <<" "<<LastBunchEvent.Time<< "\n";
    };     
    if (!blFirstTrigger) return kTRUE;

    
    //Check that daq_ch is defined in LUT
    std::map<unsigned int, TDelilaDetector >::iterator it = LUT_DELILA.find(daq_ch);
    if(it == LUT_DELILA.end()){return kTRUE;};
    
    
    DelilaEvent.domain = LUT_DELILA[daq_ch].dom;   
    int domain = DelilaEvent.domain;
    
    mDelila_raw->Fill(domain,DelilaEvent.fEnergy);
    hDomainHit->Fill(domain);
    hDetTypeHit->Fill(DelilaEvent.det_def);
    
    DelilaEvent.EnergyCal = CalibDet(DelilaEvent.fEnergy, daq_ch);
   
    if (DelilaEvent.EnergyCal < LUT_DELILA[daq_ch].threshold) return kTRUE;
    
    DelilaEvent.cs_domain = LUT_DELILA[daq_ch].cs_dom;
    DelilaEvent.theta= LUT_DELILA[daq_ch].theta;
    DelilaEvent.phi= LUT_DELILA[daq_ch].phi;

     //Check if the tree is time sorted
     DelilaEvent.Time = fTimeStampFS;
     double time_diff_last = DelilaEvent.Time - lastDelilaEvent.Time;
     //Check that the Tree is sorted in Time
     if (time_diff_last<0){std::cout<<"Warning time_diff_last: .Time  TTree may be not sorted by time"<< time_diff_last<<" \n";};
     if (DelilaEvent.Time == 0) {hTimeZero->Fill(DelilaEvent.fChannel);};
     hTimeSort->Fill(time_diff_last);
     DelilaEvent.Time-=LUT_TA_TRG[DelilaEvent.domain];
     lastDelilaEvent = DelilaEvent;
          
    if (DelilaEvent.det_def == 9){//pulser
        CheckPulserAllignement(90);
          return kTRUE;
    };
    
   if (DelilaEvent.domain == 99){//new trigger
        hTriggerTrigger->Fill(DelilaEvent.Time - LastTriggerEvent.Time);
        LastTriggerEvent = DelilaEvent;
        LastBunchEvent = DelilaEvent;
        if (blCS) cs();
        if (blGammaGamma) gamma_gamma();
        TreatFold();
        
        hdelilaQu_size->Fill(delilaQu.size());
        delilaQu.clear();
        
        nn_current_bunch = 0;
        hEventsPerTrigger->Fill(trigger_events);   
        trigger_events = 0;        
        trigger_cnt++;
        return kTRUE;        
    }else {
        trigger_events++;
        //DelilaEvent.TimeBunch = DelilaEvent.Time - LastBunchEvent.Time;
        DelilaEvent.TimeTrg = DelilaEvent.Time - LastTriggerEvent.Time;
        int nn =  DelilaEvent.TimeTrg  / bunch_length;
        DelilaEvent.TimeBunch = DelilaEvent.TimeTrg  - nn * bunch_length;


        if (DelilaEvent.TimeTrg < 1e6) mDomainTimeDiff_trigger_TA->Fill(DelilaEvent.domain,DelilaEvent.TimeTrg);
        
        mDomainTimeDiff_trigger->Fill(DelilaEvent.domain,DelilaEvent.TimeTrg);
        mDomainTimeDiff_bunch->Fill(DelilaEvent.domain,DelilaEvent.TimeBunch); 
        
        hNBunch->Fill(nn);
        
        TreatDelilaEvent();
        if (DelilaEvent.TimeBunch >= bunch_length){
               if (blCS) cs();
               if (blGammaGamma) gamma_gamma();
               TreatFold();
               hdelilaQu_size->Fill(delilaQu.size());
               delilaQu.clear();
               
               if (DelilaEvent.TimeBunch < 0) std::cout<<"Warning  DelilaEvent.TimeBunch < 0 "<< DelilaEvent.TimeBunch <<" \n";
                nn_current_bunch = nn;
            };
//             if ((DelilaEvent.det_def == 3)||(DelilaEvent.det_def == 1)){
                   // mEnergy_time_diff[domain]->Fill(DelilaEvent.EnergyCal,DelilaEvent.TimeBunch); //uncomment if needed
//             }; 
            
        if (DelilaEvent.TimeBunch < bunch_reset){
                hTimeInBunch->Fill(DelilaEvent.TimeBunch);         
                DelilaEvent.bunch = nn_current_bunch;
                DelilaEvent.trg = trigger_cnt;
                delilaQu.push_back(DelilaEvent);
//                 if (DelilaEvent.det_def == 3) {
//                     nfold++;
//                     foldQu.push_back(DelilaEvent);
//                 };
              };
             
//              mDomainTimeDiff_trigger->Fill(DelilaEvent.domain,DelilaEvent.TimeTrg);
//              mDomainTimeDiff_bunch->Fill(DelilaEvent.domain,DelilaEvent.TimeBunch);
                 
    };
    


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
   
   nevents++;
   nevents_reset++;
   return kTRUE;
}

void DelilaSelectorTrigger::TreatFold()
{
   if (delilaQu.empty())return;
//      Int_t fold_size = foldQu.size();
     int nfold = 0;
     std::deque<TDelilaEvent>::iterator it_delila_ = delilaQu.begin();

     for (; it_delila_!= delilaQu.end();++it_delila_){
         if (it_delila_->det_def == 3) nfold++;
        };
         
    hBunchFold->Fill(nfold);
    it_delila_ = delilaQu.begin();
     
     for (; it_delila_!= delilaQu.end();++it_delila_){
         if (it_delila_->det_def == 3) {
             it_delila_->fold = nfold;
//              mFoldEnergy->Fill(nfold,it_delila_->EnergyCal);//uncomment if want to fill
             if (blOutTree) {
                DelilaEventCS = *it_delila_;                
                outputTree->Fill();
             };
            };
        };
    nfold = 0;
}


void DelilaSelectorTrigger::gamma_gamma()//one gamma start the seocnd stop
{
   if (delilaQu.empty())return;
   
   std::map<int, int> nmult;
   std::map<UInt_t, std::string>  ::iterator it_mult_ =  gg_coinc_id.begin();

   for(;it_mult_!=gg_coinc_id.end();++it_mult_) nmult[it_mult_->first] = 0;   
   
   std::deque<TDelilaEvent>::iterator it1_= delilaQu.begin();
   std::deque<TDelilaEvent>::iterator it2_= delilaQu.begin();

     for (; it1_!= delilaQu.end();++it1_){
         
         if ((it1_->det_def != 1 )&&(it1_->det_def != 3 )) continue;
//          if (it1_->TimeBunch < 0 ) continue;//do not need because we shift the TimeBunch below zero
         
         for (; it2_!= delilaQu.end();++it2_){
             
             if (it1_ == it2_) continue;
//              if (it2_->TimeBunch < 0 ) continue;//do not need because we shift the TimeBunch below zero
             if ((it2_->det_def != 1 )&&(it2_->det_def != 3 )) continue;
             
             int coinc_id = GetCoincID(it1_->det_def, it2_->det_def);
             
             mGG_time_diff[coinc_id]->Fill(it1_->domain,it1_->TimeBunch );
             mGG_time_diff[coinc_id]->Fill(it2_->domain,it2_->TimeBunch );
            
             if (abs(it1_->TimeBunch) > coinc_gates[coinc_id]) continue;
             if (abs(it2_->TimeBunch) > coinc_gates[coinc_id]) continue;
             
//              double_t time_diff_coinc = it1_->TimeBunch - it2_->TimeBunch;
//              std::cout<<coinc_id<<"\n";
             nmult[coinc_id]++;
                 
                 if (coinc_id == 13) 
                 {
                      mGG_long[coinc_id]->Fill(it2_->EnergyCal, it1_->EnergyCal);
                      mGG_DC_long[coinc_id]->Fill(it2_->EnergyDC, it1_->EnergyDC);
                      
//                       mTimeDiff[it2_->domain]->Fill(it2_->domain,it2_->TimeBunch );
//                       mTimeDiff[coinc_id]->Fill(it1_->domain,it1_->TimeBunch );
                      
                      if ((it1_->CS == 0)&&(it2_->CS == 0))
                      {
                          mGG_CS_long[coinc_id]->Fill(it2_->EnergyCal, it1_->EnergyCal);
                          mGG_CS_DC_long[coinc_id]->Fill(it2_->EnergyDC, it1_->EnergyDC);
                        }
                      
                 } else
                 {
                      mGG_long[coinc_id]->Fill(it1_->EnergyCal, it2_->EnergyCal);
                      mGG_DC_long[coinc_id]->Fill(it1_->EnergyDC, it2_->EnergyDC);
                      
//                       mTimeDiff[coinc_id]->Fill(it2_->domain,it2_->TimeBunch );
//                       mTimeDiff[it1_->domain]->Fill(it1_->domain,it1_->TimeBunch );
                      
                      if ((it1_->CS == 0)&&(it2_->CS == 0))
                      {
                          mGG_CS_long[coinc_id]->Fill(it1_->EnergyCal, it2_->EnergyCal);
                          mGG_CS_DC_long[coinc_id]->Fill(it1_->EnergyDC, it2_->EnergyDC);
                    };
                 };
             };
        };
//         
   it_mult_ =  gg_coinc_id.begin();
   for(;it_mult_!=gg_coinc_id.end();++it_mult_) hMult[it_mult_->first]->Fill(nmult[it_mult_->first]);   
        
};

/*
void DelilaSelectorEliade::gamma_gamma()//one gamma start the seocnd stop
{
   if (delilaQu.empty())return;
   
   std::map<int, int> nmult;
   std::map<UInt_t, std::string>  ::iterator it_mult_ =  gg_coinc_id.begin();

   for(;it_mult_!=gg_coinc_id.end();++it_mult_) nmult[it_mult_->first] = 0;   
   
   std::deque<TDelilaEvent>::iterator it1_= delilaQu.begin();
   std::deque<TDelilaEvent>::iterator it2_= delilaQu.begin();

     for (; it1_!= delilaQu.end();++it1_){
         if ((it1_->det_def != 1 )||(it1_->det_def != 3 )) continue;
         for (; it2_!= delilaQu.end();++it2_){
             if (it1_ == it2_) continue;
             int coinc_id = GetCoincID(it1_->det_def, it2_->det_def);
             double_t time_diff_coinc = it1_->TimeBunch - it2_->TimeBunch;
             if (abs(time_diff_coinc) < coinc_gates[coinc_id])
             {
                  nmult[coinc_id]++;
                 
                 if (coinc_id == 13) 
                 {
                      mGG[coinc_id]->Fill(it2_->EnergyCal, it1_->EnergyCal);
                      mTimeDiff[coinc_id]->Fill(it2_->domain,time_diff_coinc);
                 } else
                 {
                      mGG[coinc_id]->Fill(it1_->EnergyCal, it2_->EnergyCal);
                      mTimeDiff[coinc_id]->Fill(it1_->domain,time_diff_coinc);
                 };
             };
        };
     };
        
   it_mult_ =  gg_coinc_id.begin();
   for(;it_mult_!=gg_coinc_id.end();++it_mult_) hMult[it_mult_->first]->Fill(nmult[it_mult_->first]);   
        
};*/

// void DelilaSelectorTrigger::TreatFold()//older version
// {
//    if (!foldQu.empty()){
//      Int_t fold_size = foldQu.size();
//      std::deque<TDelilaEvent>::iterator it_fold_ = foldQu.begin();
// 
//      for (; it_fold_ != foldQu.end();++it_fold_){
//           mFoldEnergy->Fill(fold_size,it_fold_->EnergyCal);
//            it_fold_->fold = fold_size;
//            if (blOutTree) {
//                 DelilaEventCS = *it_fold_;                
//                 outputTree->Fill();
//              };
//          };
//         
//      };
//     
//     hBunchFold->Fill(nfold);
//     nfold = 0;
//     foldQu.clear();
// }


void DelilaSelectorTrigger::TreatDelilaEvent()
{
    
    UShort_t daq_ch = DelilaEvent.channel;
    UShort_t domain = DelilaEvent.domain;
    
//     DelilaEvent.EnergyCal = CalibDet(DelilaEvent.fEnergy, daq_ch);
    double costheta = TMath::Cos(LUT_DELILA[daq_ch].theta);
    DelilaEvent.EnergyDC = DelilaEvent.EnergyCal*(1./sqrt(1 - beta*beta) * (1 - beta*costheta));
    
    hDelila[DelilaEvent.det_def]->Fill(DelilaEvent.EnergyCal);
    mDelila->Fill(domain,DelilaEvent.EnergyCal);
    mDelila_long->Fill(domain,DelilaEvent.EnergyCal);
    mDelilaDC->Fill(domain,DelilaEvent.EnergyDC);
    mDelilaDC_long->Fill(domain,DelilaEvent.EnergyDC);
    
//     mThetaPhi->Fill(DelilaEvent.theta, DelilaEvent.phi);

    return;
}


void DelilaSelectorTrigger::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void DelilaSelectorTrigger::Terminate()
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

      outputTree->Write();
      
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
        
//        outputTree->Write();
       foutFile->Close();
       delilaQu.clear();
       
//        std::cout<<"trigger_cnt "<<trigger_cnt<<std::endl;
   

}


int DelilaSelectorTrigger::GetCoincTimeCorrection(int dom1, int dom2)
{
 int coin_id = GetCoincID(dom1, dom2);
//  if (dom1 >= dom2 ){coin_id = dom1*100+dom2;}else {coin_id = dom2*100+dom1;};
    
//  int coin_id = dom1*100+dom2;
 int time_corr = 0;
 std::map<int, int >::iterator it = LUT_TA.find(coin_id);
 if(it != LUT_TA.end()){time_corr = LUT_TA[coin_id];};
 return time_corr;
}

void DelilaSelectorTrigger::cs()
{
    if (delilaQu.empty()) return;
    
    //BGO-GAMMA forward
    std::map<int, double> last_bgo_time;
    
    std::map<unsigned int, TDelilaDetector >::iterator it_lut_ = LUT_DELILA.begin();
    for (; it_lut_ != LUT_DELILA.end(); ++it_lut_){
        last_bgo_time[LUT_DELILA[it_lut_->first].cs_dom] = -1;
    }
    

//     double last_bgo_time = 0;
    double time_diff_bgo = -1;
    std::deque<TDelilaEvent>  ::iterator it_ev__ = delilaQu.begin();
    
    for (; it_ev__ != delilaQu.end(); ++it_ev__){
        int cs_dom = (*it_ev__).cs_domain;
        int det = (*it_ev__).det_def;
        
        switch (det)
        {
            case 5: {
                last_bgo_time[cs_dom] = (*it_ev__).TimeBunch; 
                break;
            };
            case 1: case 3: {
               if (last_bgo_time[cs_dom] == -1)
               {
                 break;  
               }else{
                 time_diff_bgo =  (*it_ev__).TimeBunch - last_bgo_time[cs_dom];
                 mTimeDiffCS ->Fill(cs_dom, time_diff_bgo);
//                  std::cout<<"time_diff_bgo forw "<<time_diff_bgo<<" bunch ID "<<(*it_ev__).bunch <<"\n";
                 if (abs(time_diff_bgo) < coinc_gates[det*10+5]) //10000)
                 {
                      (*it_ev__).CS = 1;  
                 }else{
                     last_bgo_time[cs_dom] = -1;
                 };
                 
               };
            };
            default: break;
        };
        
    };

     //GAMMA-BGO backward
//     last_bgo_time = 0;
    time_diff_bgo = -1;
    it_lut_ = LUT_DELILA.begin();
    for (; it_lut_ != LUT_DELILA.end(); ++it_lut_){
        last_bgo_time[LUT_DELILA[it_lut_->first].cs_dom] = -1;
    }
    time_diff_bgo = -1;
    it_ev__ = delilaQu.end();
    
    for (; it_ev__ != delilaQu.begin(); --it_ev__){
        int cs_dom = (*it_ev__).cs_domain;
        int det = (*it_ev__).det_def;
        
        switch (det)
        {
            case 5: {
                last_bgo_time[cs_dom] = (*it_ev__).TimeBunch; 
                break;
            };
            case 1: case 3: {
               if ((last_bgo_time[cs_dom] == -1)||((*it_ev__).CS == 1))
               {
                 break;  
               }else{
                 time_diff_bgo =  (*it_ev__).TimeBunch - last_bgo_time[cs_dom];
                 mTimeDiffCS ->Fill(cs_dom, time_diff_bgo);
//                  std::cout<<"time_diff_bgo back "<<time_diff_bgo<<" bunch ID "<<(*it_ev__).bunch <<"\n";
                if (abs(time_diff_bgo) < coinc_gates[det*10+5]) //10000)
                 {
                      (*it_ev__).CS = 1;  
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
         if (((*it_ev__).CS == 0)&&((*it_ev__).det_def == 3)) hDelilaCS[(*it_ev__).det_def]->Fill((*it_ev__).EnergyCal);
         if (((*it_ev__).CS == 0)&&((*it_ev__).det_def == 1)) hDelilaCS[(*it_ev__).det_def]->Fill((*it_ev__).EnergyCal);
    };
    
};


 int DelilaSelectorTrigger::GetCoinc_det_def(int det_def1, int det_def2)
 {
     int id=0;
     if (det_def1<=det_def2) {id = det_def1*10+det_def2;}
         else {id = det_def2*10+det_def1; };
     return id;
 }


int DelilaSelectorTrigger::GetCoincID(int dom1, int dom2){
    int id=0;
    if (dom1<=dom2) {id = dom1*10+dom2;}//was dom*100
        else {id = dom2*10+dom1; };
    return id;
}


void DelilaSelectorTrigger::CheckPulserAllignement(int zero_dom)
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

void DelilaSelectorTrigger::PrintDelilaEvent(TDelilaEvent &ev_)
{
    std::cout<<" === DelilaEvent === \n"<<
    "  ev_.domain "<<ev_.domain<<" \n"<<
    "  ev_.EnergyCal "<<ev_.EnergyCal<<" \n"<<
    "  ev_.fold "<<ev_.fold<<" \n";
//     "  ============= \n";
}



