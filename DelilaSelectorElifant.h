//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 23 13:02:03 2021 by ROOT version 6.22/06
// from TTree DelilaSelectorElifant/Energy Station
// found on file: run1005_03.root
//////////////////////////////////////////////////////////

#pragma once

#ifndef DelilaSelectorElifant_h
#define DelilaSelectorElifant_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH2.h>
#include <TString.h>
#include <TObjString.h>

#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>
#include <fstream>
#include <limits>
#include <csignal>
#include <ctime>
#include <TEntryList.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TVector3.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TCutG.h>
#include <TList.h>
#include <TKey.h>
#include <TSysEvtHandler.h>
#include <TSystem.h>
#include <TApplication.h>
#include <deque>
#include "TTreeIndex.h"
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* getenv */
#include <map>
#include <functional>
#include <queue>
#include <vector>

#include "DelilaEvent.h"
#include "HPGeTreeEvent.h"
#include "LaBrTreeEvent.h"
#include "ElissaTreeEvent.h"

//#include "TObjString.h."
// Headers needed by this particular selector


class DelilaSelectorElifant : public TSelector {
public :
    TTreeReader     fReader;  //!the tree reader
    TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
    TFile          *foutFile;
    TTree          *outputTree;
    ULong64_t	    fTimeStamp;
    double_t	    fTimeStampFS;//FineTS
    UShort_t	    fEnergyShort;
    UShort_t	    fEnergyLong;

  class TDelilaDetector { 
  public:
    Int_t	 dom;
    int	 ch;//ch daq
//     Int_t	 serial;
    TString	 serial;
    Float_t  theta;
    Float_t	 phi;  
    UShort_t detType;//0 - nothing; 1 - core; 2 - segment; 3 - CeBr; 4 - CsI; 5 - BGO1; 6 - BGO2; 9 - pulser
    Int_t	 TimeOffset;
    double_t bgo_time_offset;
    Int_t 	 threshold; 
    Int_t	 pol_order;
    Int_t    cs_dom;
    std::vector<float> calibE;
    TDelilaDetector(): dom(-1),phi(-1),theta(-1),TimeOffset(0),calibE(0),threshold(-1),ch(-1),pol_order(-1){};
 };
 
//   DelilaEvent     delila_tree_event;
  HPGeTreeEvent   hpge_tree_event;
  LaBrTreeEvent   labr_tree_event;
  ElissaTreeEvent elissa_tree_event;
 
  std::vector<HPGeTreeEvent>   *HPGeEvent;
  std::vector<LaBrTreeEvent>   *LabrEvent;
  std::vector<ElissaTreeEvent> *ElissaEvent;
  
  std::deque<DelilaEvent>      delilaQu;
  std::deque<DelilaEvent>      delilaPreQu;
  
  std::map<int,std::deque<DelilaEvent>> waitingQu_gamma; //for CeBr
  std::map<int,std::deque<DelilaEvent>> waitingQu_bgo; //for CS
  
  //Part for PHA
  short          fSignal[2000]; // add by saka //RecordLength branch is not read we assume rl = 2000
  
//   std::map<unsigned int, TDelilaDetector > LUT_DELILA;
  std::map<int, TDelilaDetector >       LUT_DELILA;    
  std::map<int, int >                   LUT_TA;
  std::map<int, double_t >              LUT_TA_TRG;

  DelilaEvent DelilaEvent_;  
  DelilaEvent DelilaEventTreated ;
  DelilaEvent lastDelilaEvent;  
  DelilaEvent lastEliadeZeroEvent;
  DelilaEvent LastTriggerEvent;
  DelilaEvent LastBunchEvent;    
  
  DelilaEvent PulserEvent;
  DelilaEvent DomainZeroEvent;    
  

  TBranch *b_channel;
  TBranch *b_tstmp;
  TBranch *b_tstmp_fine;
  TBranch *b_energ;  //ChargeLong
  TBranch *b_energ_short;  //ChargeShot
  TBranch *b_mod;
  TBranch *b_signal;
  
  Long64_t nb_entries;
  
  Long64_t event_length;
  Long64_t pre_event_length;
  Long64_t bunch_reset;
  
  UChar_t	        fMod; 
  UChar_t	        fChannel; 
 
  TH1F* hChannelHit;
  TH1F* hDomainHit;
  TH1F* hSegmentHit;
  TH1F* hDetTypeHit;

  std::map<int, TH1F*> hDelila0; //before event_builder
  std::map<int, TH1F*> hDelilaCS;
  std::map<int, TH1F*> hDelilaDC;
  std::map<int, TH1F*> hDelilaCS_DC;
  
  std::map<int, TH1F*> hDelila_single; //after event builder
  
  std::map<int, TH1F*> hDelila_long;
  std::map<int, TH1F*> hDelilaCS_long;
  std::map<int, TH1F*> hDelilaDC_long;
  std::map<int, TH1F*> hDelilaCS_DC_long;

  std::map<UInt_t, std::string> detector_name;

  TH1F* hTriggerTrigger;
  TH1F* hBunchBunch;
  TH1F* hNBunch;
  TH1F* hBunchFold;
  TH1F* hTimeInBunch;
  TH1F* hEventsPerTrigger;
  TH2F* mFoldEnergy;
  
  TH1F* hCoincID;
  TH1F* hTriggerDomain;
  
  TH1F* hdelilaQu_size;
  
  TH2F* mDelila;
  TH2F* mDelila_raw;
  TH2F* mDelilaCS;
  TH2F* mDelilaDC;//keV
  TH2F* mDelilaCS_DC;//keV
  
  TH2F* mDelila_long;
  TH2F* mDelila_raw_long;
  TH2F* mDelilaCS_long;
  TH2F* mDelilaDC_long;//keV
  TH2F* mDelilaCS_DC_long;//keV

  TH2F* mGammaGammaDC;
  TH2F* mGammaGammaCS_DC;
  
  TH2F* mEnergyTimeDiff_trigger;
  TH2F* mDomainTimeDiff_trigger;
  TH2F* mDomainTimeDiff_bunch;
  
  std::map<int, TH2F*> mGG;
  std::map<int, TH2F*> mGG_theta;
  std::map<int, TH2F*> mGG_CS;
  std::map<int, TH2F*> mGG_DC;
  std::map<int, TH2F*> mGG_CS_DC;
  std::map<int, TH2F*> mGG_time_diff;
  
  std::map<int, TH2F*> mGG_long;
  std::map<int, TH2F*> mGG_CS_long;
  std::map<int, TH2F*> mGG_DC_long;
  std::map<int, TH2F*> mGG_CS_DC_long;
  
//   std::map<int, TH2F*> mTimeDiffselected_run_2022_0.root;
  std::map<int, TH1F*> hMult;

  std::map<UInt_t, std::string> gg_coinc_id;
  std::map<UInt_t, Float_t> coinc_gates;//in ps

  
  std::map<int, std::string> domain_list;
  std::map<int, TH2F*> mEnergy_time_diff;
  
  TH2F* mTimeDiffCS;
  
  TH2F* mThetaPhi; 
  TH2F* mGammaGamma;
  TH2F* mTimeDiff_gg;
  TH1F* hMult_gg;
  TH2F* mGammaGammaCS;
  TH2F* mTimeDiff_gg_CS;
  TH1F* hMult_gg_CS;
  
  TH2F* mLaBr_LabBr_time_diff;

  
  TH2F *mPulser0TimeDiff;
  
  TH1F* hTimeDiffPulser;
  TH2F* mPulserPulser;
  
  TH1F* hTimeZero;
  TH1F* hTimeSort;
    
  TH2F* mTimeCalibDomain0;
  TH2F* mTimeCalib;
  TH2F* mTimeCalibBGO;
  TH2F* mTimeCalibBGO_cs_dom;
  
  TH2F* mNeutron;
  TH2F* mShortLong;
//   TH2F* mTimeCalibDomain0;
  
  TH1F *hAmax;
  TH2F *mAmaxEnergy;

  int det_def_trg;
  int channel_trg;
  double TriggerTimeFlag;
  double lastDelilaTime;

    
  std::clock_t start;
  double duration;
 
  ULong64_t nevents;
//   ULong64_t nevents_reset;
  int reset_counter;
 
  ULong64_t lastTime;
  
  TObjArray *toks;

   DelilaSelectorElifant(TTree * /*tree*/ =0) { }
   virtual ~DelilaSelectorElifant() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   virtual Long64_t GetEntries() { return fChain ? fChain->GetEntries() : 0;}

   virtual void  Read_ELIADE_LookUpTable();
   virtual void  Read_TimeAlignment_LookUpTable();
   virtual void  Read_TimeAlignment_Trigger();
   virtual void  Read_Confs();
   virtual void  Print_ELIADE_LookUpTable();
   virtual void  Print_TimeAlignment_LookUpTable();
   virtual void  Print_TimeAlignment_Trigger_LookUpTable();

   virtual float CalibDet(float,int);
   virtual int GetCoincTimeCorrection(int dom1, int dom2);
   virtual void cs();
   
   virtual void TreatDelilaEvent();
   virtual void TreatFold(int det); 
   virtual int  GetCoincID(int dom1, int dom2);
   virtual int  GetCoinc_det_def(int det_def1, int det_def2);
   virtual void CheckPulserAllignement(int zero_dom);
   virtual void PrintDelilaEvent(DelilaEvent &ev_);
   
   virtual void FillOutputTree();
   virtual void TimeAlignementTrigger();
   
   
   virtual void EventBuilderSimple();
   virtual void EventBuilderPreTrigger();
   virtual bool TriggerDecision();
   virtual void SetUpNewTrigger();
   virtual void CheckPreQu();
   virtual void MovePreQu2Qu();

   
   virtual void TreatLaBrSingle();
   virtual void TreatHpGeSingle();
   virtual void TreatNeutronSingle();
   virtual void TreatElissaSingle();
   virtual void TreatACS();
   
   virtual void TreatGammaGammaCoinc();
   virtual void TreatSolarLaBrCoinc();
   
   virtual std::vector<float> trapezoidal(short wave[],int length, int L, int G);//L = 20; G = 0
   
   std::map<std::string, bool> has_detector;


   ClassDef(DelilaSelectorElifant,0);

  
};

#endif

#ifdef DelilaSelectorElifant_cxx
void DelilaSelectorElifant::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
  
  if(!tree) {std::cout<<" TTree NOT found "<<std::endl; return;};
  
//   std::cout<<" Delila Sorting "<<std::endl;
  std::cout<<" TTree found "<<std::endl;
  fChain = tree;
  fChain->SetMakeClass(1);  
  fChain->SetBranchAddress("ChargeLong", 	&fEnergyLong,            	   &b_energ);
  fChain->SetBranchAddress("ChargeShort", 	&fEnergyShort, 	               &b_energ_short);
  fChain->SetBranchAddress("FineTS", 	    &fTimeStampFS, 	               &b_tstmp_fine);
  fChain->SetBranchAddress("TimeStamp",     &fTimeStamp,                   &b_tstmp); 
  fChain->SetBranchAddress("Ch", 	        &fChannel,                     &b_channel);
  fChain->SetBranchAddress("Mod", 	        &fMod,                         &b_mod);
  fChain->SetBranchAddress("Signal",        &fSignal,                      &b_signal);// add by saka 
  
  fReader.SetTree(tree);
       
   
  foutFile->cd();
  outputTree = new TTree("SelectedDelila","SelectedDelila");
  

  lastDelilaTime = 0;
  
  Read_ELIADE_LookUpTable();
  Read_TimeAlignment_LookUpTable();
  
//   Print_ELIADE_LookUpTable();
  
  if (has_detector["HPGe"]){
      HPGeEvent= new std::vector<HPGeTreeEvent>;
      outputTree->Branch("HPGeEvents",&HPGeEvent);
    };
  
  if (has_detector["LaBr"]){
      LabrEvent= new std::vector<LaBrTreeEvent>;
      outputTree->Branch("LaBrEvents",&LabrEvent);
    };
  
  if (has_detector["Elissa"]){
      ElissaEvent= new std::vector<ElissaTreeEvent>;
      outputTree->Branch("ElissaEvents",&ElissaEvent);
  };
  
  
  std::cout<<" === Present detectors === \n";
  std::cout<<" HPGe   " << has_detector["HPGe"] <<"  \n";
  std::cout<<" LaBr   " << has_detector["LaBr"] <<"  \n";
  std::cout<<" ACS    " << has_detector["ACS"] <<"  \n";
  std::cout<<" Elissa " << has_detector["Elissa "] <<" \n";
  std::cout<<" === Time settings ps === \n";
  std::map<UInt_t, Float_t> ::iterator itcc_ = coinc_gates.begin();
  for (; itcc_ != coinc_gates.end(); ++itcc_) {
     std::cout<<" coin_id " << itcc_->first <<"  coinc_gate "<< itcc_->second <<" ps \n";
  }
  std::cout<<" event_length     " << event_length      <<" ps \n";
  std::cout<<" pre_event_length " << pre_event_length  <<" ps \n";
  std::cout<<" ===                            === \n";
  
  std::cout << "Init() is done \n";

  
//   outputTree->Branch("fTEventTS",&DelilaEventTreated .fTimeStamp,"TimeStamp/l");
//   outputTree->Branch("fTEventFS",&DelilaEventTreated .fTimeStampFS,"TimeStamp/D");
//   outputTree->Branch("fTimeBunch",&DelilaEventTreated .TimeBunch,"TimeBunch/D");
//   outputTree->Branch("fTime",&DelilaEventTreated .Time,"Time/D");
//   outputTree->Branch("fEnergy",&DelilaEventTreated .fEnergy,"Energy/F");
//   outputTree->Branch("fEnergy_kev",&DelilaEventTreated .EnergyCal,"Energy_kev/F");
//   outputTree->Branch("fEnergyDC",&DelilaEventTreated .EnergyDC,"EnergyDC/F");
//   outputTree->Branch("fDomain",&DelilaEventTreated .domain,"Domain/b");
//   outputTree->Branch("fDetType",&DelilaEventTreated .det_def,"def/b");
//   outputTree->Branch("fCS",&DelilaEventTreated .CS,"CS/b");
//   outputTree->Branch("fTRG",&DelilaEventTreated .trg,"Trigger/b");
//   outputTree->Branch("fFold",&DelilaEventTreated .fold,"Fold/b");
  
//   LaBrEvent= new std::vector<LaBrEvent>;//dt_mon
//   outputTree->Branch("LaBrEvents",&LaBrEvent);//dt_mon
  
//   ElissaEvent= new std::vector<ElissaEvent>;//dt_mon
//   outputTree->Branch("ElissaEvents",&ElissaEvent);//dt_mon
}

Bool_t DelilaSelectorElifant::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef DelilaSelectorElifant_cxx
