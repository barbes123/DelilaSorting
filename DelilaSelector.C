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


int addBackMode = 0; //0 - no addback; 1- addback;//not in use for ELIFANT
bool blGammaGamma = false;
bool blGammaGammaCS = false;
bool blCS = true;


bool debug = false;
// bool doCS = false;

const int NumberOfClovers = 2;
const int max_domain = 400;
const int max_mod = 7;


// const int current_clover = 1;
const int nbr_of_boards = 8;
const int nbr_of_ch = 200;
// const int zero_channel = 101; //for time allignement
ULong64_t lastTime_pulser = 0;
ULong64_t lastTime_dom0 = 0;
ULong64_t lastTimeStamp = 0;

// std::unordered_set<int> cores = { 109,119,129,139,171,172};
std::unordered_set<int> pulsers = { 150,151,152,153};
std::unordered_set<int> crystal2mask = {1,3,5};
//std::unordered_set<int> frontsegments = {};
//std::unordered_set<int> cores = { 101,111,121};

// bool & DelilaSelector::operator<(const TDelilaEvent& ev1, const TDelilaEvent& ev2){return ev1.TimeStamp < ev2.TimeStamp;}


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
      TDelilaDetector curDet;
      Float_t theta(-1.), phi(-1.);
      int upperThreshold = 1e6;
      std::istringstream is(oneline);
      if (debug) std::cout << is.str().c_str() << std::endl;
//       is >> curDet.ch >> curDet.dom >> curDet.theta >> curDet.phi >> curDet.TimeOffset >> curDet.upperThreshold;
      is >> curDet.ch >> curDet.dom >> curDet.detType >> curDet.serial >> curDet.TimeOffset >> curDet.theta >> curDet.phi >> curDet.upperThreshold >> curDet.cs_dom;
    //  std::cout<<" curDfalseet.ch  "<<curDet.ch <<" curDet.TimeOffset " <<curDet.TimeOffset<<std::endl;
      
      if (curDet.ch >= 0) {
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
  addBackMode = atoi(((TObjString*) toks->At(2))->GetString());
  std::cout << "addBackMode  " << addBackMode <<std::endl;

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
   
        
   hDelila = new TH1F("hDelila", "hDelila", 4096, -0.5, 16383.5);
   hDelila->GetYaxis()->SetTitle("counts");
   hDelila->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelila);
   
   hDelilaCS = new TH1F("hDelilaCS", "hDelilaCS", 4096, -0.5, 16383.5);
   hDelilaCS->GetYaxis()->SetTitle("counts");
   hDelilaCS->GetXaxis()->SetTitle("keV");
   fOutput->Add(hDelilaCS);
   
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
   
   mThetaPhi = new TH2F("mThetaPhi", "mThetaPhi", 90,60,150,360,0,360);
   mThetaPhi->GetXaxis()->SetTitle("theta, degrees");
   mThetaPhi->GetYaxis()->SetTitle("phi, degrees");
   fOutput->Add(mThetaPhi);
   
   
   mDelilaCS = new TH2F("mDelilaCS", "mDelilaCS", max_domain, 0, max_domain, 16384, -0.5, 16383.5);
   mDelilaCS->GetXaxis()->SetTitle("domain");
   mDelilaCS->GetYaxis()->SetTitle("keV");
   fOutput->Add(mDelilaCS);
      
   mGammaGamma = new TH2F("mGammaGamma", "mGammaGamma", 4096, -0.5, 16383.5, 16384, -0.5, 16383.5);
   mGammaGamma->GetXaxis()->SetTitle("domain");
   mGammaGamma->GetYaxis()->SetTitle("keV");
   fOutput->Add(mGammaGamma);
   
   mGammaGammaCS = new TH2F("mGammaGammaCS", "mGammaGammaCS", 16384, -0.5, 16383.5, 16384, -0.5, 16383.5);
   mGammaGammaCS->GetXaxis()->SetTitle("domain");
   mGammaGammaCS->GetYaxis()->SetTitle("keV");
   fOutput->Add(mGammaGammaCS);
   
   mTimeDiffCS = new TH2F("mTimeDiffCS", "mTimeDiffCS", max_domain, 0, max_domain, 4e3, -2e5, 2e5);
   mTimeDiffCS->GetXaxis()->SetTitle("domain");
   mTimeDiffCS->GetYaxis()->SetTitle("ps/bin");
   fOutput->Add(mTimeDiffCS);
      
   mSegments = new TH2F("mSegments", "mSegments", max_domain, 0, max_domain, 16384, -0.5, 16383.5);
   mSegments->GetXaxis()->SetTitle("domain");
   mSegments->GetYaxis()->SetTitle("keV");
   fOutput->Add(mSegments);
   
   mCores = new TH2F("mCores", "mCores", max_domain, 0, max_domain, 16384, -0.5, 16383.5);
   mCores->GetXaxis()->SetTitle("domain");
   mCores->GetYaxis()->SetTitle("keV");
   fOutput->Add(mCores);
   
   mTimeCalib = new TH2F("mTimeCalib", "mTimeCalib", 10000, 0, 10000, 4e3, -2e5, 2e5);
   mTimeCalib->GetXaxis()->SetTitle("coinc ID");
   mTimeCalib->GetYaxis()->SetTitle("ps");
   mTimeCalib->SetTitle("Sci time diff");
   fOutput->Add(mTimeCalib);

   mPulser0TimeDiff = new TH2F("mPulser0TimeDiff", "mPulser0TimeDiff", 200, 0.5, 200.5, 500, -99.5, 899.5);
   fOutput->Add(mPulser0TimeDiff);//time_diff relevant to the 1st channel (101), i.e. ch 101 is a trigger
   
   mTimeDiff_gg = new TH2F("mTimeDiff_gg", "mTimeDiff_gg", 200, 0.5, 200.5, 500, -99.5, 399.5);
   fOutput->Add(mTimeDiff_gg);
   
   mTimeDiff_gg_CS = new TH2F("mTimeDiff_gg_CS", "mTimeDiff_gg_CS", 200, 0.5, 200.5, 500, -99.5, 399.5);
   fOutput->Add(mTimeDiff_gg_CS);
       
   mPulserPulser = new TH2F("mPulserPulser", "mPulserPulser",4096, -0.5, 8191.5, 4096, -0.5, 8195.5);
   fOutput->Add(mPulserPulser);
    
  TString option = GetOption();
  toks = option.Tokenize(",");
  TString RunID = ((TObjString*) toks->At(0))->GetString();
  TString VolID = ((TObjString*) toks->At(1))->GetString();
  TString ServerID = ((TObjString*) toks->At(3))->GetString();
  
  std::stringstream OutputFile;
  OutputFile << "selected_run" << "_" << RunID <<"_"<<VolID;
  if (atoi(ServerID) != 0) {OutputFile<<"_eliadeS"<<ServerID;};
  OutputFile << ".root";
//   std::cout <<"ServerID "<<ServerID<<" "<< OutputFile.str().c_str() <<std::endl;
   lastTime = 0;
   
  outputFile = new TFile (OutputFile.str().c_str(),"recreate"); 
   
   std::cout <<"ServerID "<<ServerID<<" "<< OutputFile.str().c_str() <<std::endl;
    start = std::clock();
    
    lastDelilaEvent.fTimeStamp = 0;
    lastEliadeZeroEvent.fTimeStamp = 0;
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
	DelilaEvent.EnergyCal = CalibDet(DelilaEvent.fEnergy, daq_ch);
	DelilaEvent.domain = LUT_DELILA[daq_ch].dom;
    DelilaEvent.det_def = LUT_DELILA[daq_ch].detType;
    DelilaEvent.cs_domain = LUT_DELILA[daq_ch].cs_dom;
    DelilaEvent.theta= LUT_DELILA[daq_ch].theta;
    DelilaEvent.phi= LUT_DELILA[daq_ch].phi;
        
    hDetTypeHit->Fill(DelilaEvent.det_def);
    mThetaPhi->Fill(DelilaEvent.theta, DelilaEvent.phi);
	
    int domain = DelilaEvent.domain;
    if ((DelilaEvent.fEnergy < LUT_DELILA[daq_ch].upperThreshold))  {/*std::cout<<DelilaEvent.fEnergy<< " "<< LUT_ELIADE[daq_ch].upperThreshold<<std::endl; */return kTRUE;}
            
  	hDomainHit->Fill(domain);
	mDelila_raw->Fill(domain,DelilaEvent.fEnergy);
	mDelila->Fill(domain,DelilaEvent.EnergyCal);
    //Check if the tree is time sorted
    if (lastDelilaEvent.fTimeStamp > DelilaEvent.fTimeStamp){std::cout<<"Warning: .fTimeStamp TTree may be not sorted by time \n";};
    
     DelilaEvent.fTimeStamp = DelilaEvent.fTimeStamp + LUT_DELILA[daq_ch].TimeOffset;
     
//      DelilaEvent.Time = 1000*DelilaEvent.fTimeStamp + DelilaEvent.fTimeStampFS;
//         DelilaEvent.Time = DelilaEvent.fTimeStamp*1.0;
         DelilaEvent.Time = DelilaEvent.fTimeStampFS;

     if (lastDelilaEvent.Time > DelilaEvent.Time){std::cout<<"Warning: .Time  TTree may be not sorted by time \n";};
        
        
     hTimeSort->Fill(DelilaEvent.Time - lastDelilaEvent.Time);
     
        
        ///To produce matrix for to Check Time Alignment
        if (coincQu_TA.empty()){coincQu_TA.push_back(DelilaEvent);}
         else
         {
             double time_diff = DelilaEvent.Time - coincQu_TA.front().Time - GetCoincTimeCorrection(DelilaEvent.domain,coincQu_TA.front().domain);

             if (std::abs(time_diff) < 10000000) //ps
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
         
         
         
     //Doing Gamma-Gamma matrix     
      if (blGammaGamma)
      {
            //coming soon
            
      };
         
          

     ///Doing CS rejection

     if (blCS){

         //coming soon
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
    
  //finish filling the Tree if the queue was not empty   
    while ((!output_pQu.empty())) {
//      cout << output_pQu.top().fTimeStamp<< " " << output_pQu.top().fTimeStamp<< " 111 \n";
     DelilaEventCS =  output_pQu.top();
     outputTree->Fill();
     output_pQu.pop();
  }; 
    

  std::cout<<"I  will terminate soon... "<<std::endl;  

  TIter iter(fOutput);
      

  std::cout << std::endl << "Finished after " << round(duration / 60) << ":"
	    << std::setw(2) << round(((int) duration) % 60) << std::setw(8)
	    << ". Write spectra to file" << std::endl;


      TObject *obj; 
       while ((obj = iter())) {
           if ( obj->IsA()->InheritsFrom(TH1F::Class())){
            TH1 *h1 = (TH1*)obj;
                if (h1->GetEntries()>0) obj->Write();
            }else if (obj->IsA()->InheritsFrom(TH2F::Class())){
            TH2 *h2 = (TH2*)obj;
                if (h2->GetEntries()>0) obj->Write();
            }
        };
        
       outputTree->Write();
   

}

int DelilaSelector::GetCoincTimeCorrection(int dom1, int dom2)
{
 int coin_id = dom1*100+dom2;
 int time_corr = 0;
 std::map<int, int >::iterator it = LUT_TA.find(coin_id);
 if(it != LUT_TA.end()){time_corr = LUT_TA[coin_id];};
 coin_id = dom2*100+dom1;
 if(it != LUT_TA.end()){time_corr = LUT_TA[coin_id];};
 return time_corr;
}
