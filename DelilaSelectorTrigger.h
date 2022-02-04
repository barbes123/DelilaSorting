//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 23 13:02:03 2021 by ROOT version 6.22/06
// from TTree DelilaSelectorTrigger/Energy Station
// found on file: run1005_03.root
//////////////////////////////////////////////////////////

#pragma once

#ifndef DelilaSelectorTrigger_h
#define DelilaSelectorTrigger_h

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

//#include "TObjString.h."
// Headers needed by this particular selector


class DelilaSelectorTrigger : public TSelector {
public :
    TTreeReader     fReader;  //!the tree reader
    TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
    TFile          *foutFile;
    TTree          *outputTree;


//   UChar_t	 uMod; 
//   UChar_t	 uChannel; 

 class TDelilaEvent {
 public:
    UChar_t	        fMod; 
    UChar_t	        fChannel;    
    ULong64_t	    fTimeStamp;
    double_t	    fTimeStampFS;//FineTS
    UShort_t	    fEnergy;//ChargeLong
    UShort_t	    fEnergyShort;//ChargeShot  
    UShort_t        det_def;//0 - nothing; 1 - core/single HPge; 2 - segment; 3 - CeBr; 4 - CsI; 5 - BGO1; 6 - BGO2; 7 -BGO - 3; 8 - solar cell; 9 - pulser
    float	        EnergyCal;
    float           EnergyDC;
    UShort_t        domain;
    UShort_t        cs_domain;
    UShort_t        channel;//ch daq
    UShort_t        CS;//0 - no; 1 - yes
    double_t        Time;
    double_t        TimeTrg;
    double_t        TimeBunch;
    Float_t         theta;
    Float_t	        phi;
    double_t        bgo_time_diff;
    ULong64_t       trg;
    UShort_t        bunch;
    UShort_t        fold;
   //int make_board_ID(){return fMod*100}
    TDelilaEvent(): domain(-1),channel(-1),fTimeStamp(0),fEnergy(-1),CS(0),cs_domain(0),Time(0),bgo_time_diff(-1),trg(0),bunch(0),fold(0){};
 };

  class TDelilaDetector { 
  public:
    Int_t	 dom;
    Int_t	 ch;//ch daq
    Int_t	 serial;
    Float_t  theta;
    Float_t	 phi;  
    UShort_t detType;//0 - nothing; 1 - core; 2 - segment; 3 - CeBr; 4 - CsI; 5 - BGO1; 6 - BGO2; 9 - pulser
    Int_t	 TimeOffset; 
    Int_t 	 upperThreshold; 
    Int_t	 pol_order;
    Int_t    cs_dom;
    std::vector<float> calibE;
    TDelilaDetector(): dom(-1),phi(-1),theta(-1),TimeOffset(0),calibE(0),upperThreshold(-1),ch(-1),pol_order(-1){};
 };

//   struct TDelilaEventCoinc {
//       TDelilaEvent LastEvent;
//       bool coinc;
//       int mult;
//  };
 
//   std::deque<TDelilaEvent> eliadeQu;
//   std::deque<TDelilaEvent> coincQu_pulser;
  std::deque<TDelilaEvent> coincQu_TA;//for time alignment
  std::deque<TDelilaEvent> outputQu;
//   std::deque<TDelilaEvent> gammagammaQu;
  std::deque<TDelilaEvent> gammagammaQu_CS;
  std::deque<TDelilaEvent> ggLaBr_HPGe_Qu;
  std::deque<TDelilaEvent> ggBGO_Qu;
  std::deque<TDelilaEvent> foldQu;

  
//   std::deque<TDelilaEvent> bgo_Qu;
//   std::deque<float> enrergyQu;
  
  std::map<int,std::deque<TDelilaEvent>> waitingQu_gamma; //for CeBr
  std::map<int,std::deque<TDelilaEvent>> waitingQu_bgo; //for CS
  
//   std::map<int,std::deque<TDelilaEvent>> waitingQu_gg_; //for LaBr-HPGe 1 - HPGe 3 - LaBr



//   std::deque<TDelilaEvent> eliadeQu_sorted;
  std::map<unsigned int, TDelilaDetector > LUT_DELILA;  
  std::map<int, int > LUT_TA;
  std::map<int, double_t > LUT_TA_TRG;
  std::map<int, int > LUT_CONF;

  TDelilaEvent DelilaEvent;  
  TDelilaEvent DelilaEventCS;
  TDelilaEvent lastDelilaEvent;  
  TDelilaEvent lastEliadeZeroEvent;
  TDelilaEvent LastTriggerEvent;
  TDelilaEvent LastBunchEvent;    
  
  TDelilaEvent PulserEvent;
  TDelilaEvent DomainZeroEvent;    

  
//   TDelilaEventCoinc EliadeCoincEvent[4];
  
  TDelilaEvent startEventCore;  
  TDelilaEvent startEventSegment;  

  TBranch *b_channel;
  TBranch *b_tstmp;
  TBranch *b_tstmp_fine;
  TBranch *b_energ;  //ChargeLong
  TBranch *b_energ_short;  //ChargeShot
  TBranch *b_mod;  
  
  Long64_t nb_entries;
  
  Long64_t bunch_length;
  Long64_t bunch_reset;
  double   beta;
  
 // std::map<UInt_t,TH1F*> hEnergy_raw;
  //std::map<UInt_t,TH1F*> hEnergy_cal;
 
//   std::map<UInt_t, TDelilaEvent> last_board_event;
 
  TH1F* hChannelHit;
  TH1F* hDomainHit;
  TH1F* hSegmentHit;
  TH1F* hDetTypeHit;
//   TH2F* mDelila;//keV
//   TH2F* mDelilaCS;//keV
  //TH1F* hDelila;//keV
  //TH1F* hDelilaCS;//keV
  std::map<int, TH1F*> hDelila;
  std::map<int, TH1F*> hDelilaCS;
  std::map<int, TH1F*> hDelilaDC;
  std::map<int, TH1F*> hDelilaCS_DC;
  
  std::map<int, TH1F*> hDelila_long;
  std::map<int, TH1F*> hDelilaCS_long;
  std::map<int, TH1F*> hDelilaDC_long;
  std::map<int, TH1F*> hDelilaCS_DC_long;

  std::map<UInt_t, std::string> detector_name;

  TH1F* hTriggerTrigger;
  TH1F* hBunchBunch;
  TH1F* hNBunch;
  TH1F* hBunchFold;
  TH2F* mFoldEnergy;
  
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

  
  
//   TH1F* hDelilaDC;//keV
//   TH1F* hDelilaCS_DC;//keV
  TH2F* mGammaGammaDC;
  TH2F* mGammaGammaCS_DC;
  
  TH2F* mEnergyTimeDiff_trigger;
  TH2F* mDomainTimeDiff_trigger;
  TH2F* mDomainTimeDiff_bunch;
  
//   std::map<int, TH2F*> mapTimeEnergy;

//    std::map<UInt_t,std::map<UInt_t,TH1F*> > hIntegral_spectra_us;//4 beta, 5 gamma, 6 nim, 7 num

  
//   TH2F* mgg_labr_hpge;
//   TH2F* mgg_labr_hpgeCS;
//   TH2F* mgg_labr_hpgeDC;
//   TH2F* mgg_labr_hpgeCS_DC;  
//   TH2F* mgg_labr_time_diff;
  
//   TH2F* mgg_hpge_hpge;
//   TH2F* mgg_hpge_hpgeCS;
//   TH2F* mgg_hpge_hpgeDC;
//   TH2F* mgg_hpge_hpgeCS_DC;  
//   TH2F* mgg_hpge_time_diff;

  
  std::map<int, TH2F*> mGG;
  std::map<int, TH2F*> mGG_CS;
  std::map<int, TH2F*> mGG_DC;
  std::map<int, TH2F*> mGG_CS_DC;
  std::map<int, TH2F*> mGG_time_diff;
  
  std::map<int, TH2F*> mGG_long;
  std::map<int, TH2F*> mGG_CS_long;
  std::map<int, TH2F*> mGG_DC_long;
  std::map<int, TH2F*> mGG_CS_DC_long;
  
  std::map<int, TH2F*> mTimeDiff;

//   TH2F* mTimeDiff_trg_laBr;
//   TH2F* mTimeDiff_trg_bgo;
//   TH2F* mTimeDiff_trg_hpge;

  std::map<UInt_t, std::string> gg_coinc_id;
  std::map<UInt_t, UInt_t> coinc_gates;//in ps

  
  std::map<int, std::string> domain_list;
  std::map<int, TH2F*> mEnergy_time_diff;


  TH2F* mDelilaTD;
  
  
  TH2F* mThetaPhi; 
  TH2F* mGammaGamma;
  TH2F* mTimeDiff_gg;
  TH1F* hMult_gg;
  TH2F* mGammaGammaCS;
  TH2F* mTimeDiff_gg_CS;
  TH1F* hMult_gg_CS;
  
  TH2F* mLaBr_BGO_time_diff;
  TH2F* mLaBr_LabBr_time_diff;

  
//  TH2F* mDelilaSegEnergy;
//   TH2F *mDomTimeDiff;
  TH2F *mPulser0TimeDiff;
//   TH2F *mDom0TimeDiff;
//   TH2F *mDom0TimeDiffEnergy;
  
  TH1F* hTimeDiffPulser;
  TH2F* mPulserPulser;
  
  TH2F* mTimeDiffCS;
  
  TH1F* hTimeZero;
  TH1F* hTimeSort;
    
  TH2F* mTimeCalib;
  TH2F* mTimeCalibDomain0;


    
  std::clock_t start;
  double duration;
 
  ULong64_t nevents;
  ULong64_t nevents_reset;
  int reset_counter;
 
  ULong64_t lastTime;
  int n_bunches;
  
  TObjArray *toks;

   DelilaSelectorTrigger(TTree * /*tree*/ =0) { }
   virtual ~DelilaSelectorTrigger() { }
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
   virtual float CalibDet(float,int);
//    virtual int CheckTimeAlignment(int to_domain);
   virtual int GetCoincTimeCorrection(int dom1, int dom2);
   virtual void cs();
   virtual void gamma_gamma();
   virtual void TreatDelilaEvent();
   virtual void TreatFold(); 
   virtual int GetCoincID(int dom1, int dom2);
   virtual int GetCoinc_det_def(int det_def1, int det_def2);
   virtual void CheckPulserAllignement(int zero_dom);
   virtual void PrintDelilaEvent(TDelilaEvent &ev_);

   ClassDef(DelilaSelectorTrigger,0);
   
   
   struct CompareTimeStamp{
   bool operator()(TDelilaEvent const& ev1, TDelilaEvent const& ev2)
    {
        // return "true" if "p1" is ordered
        // before "p2", for example:
        return ev1.fTimeStamp> ev2.fTimeStamp;
        }
    };
    
    priority_queue<TDelilaEvent, vector<TDelilaEvent>, CompareTimeStamp> output_pQu;
 
// bool operator<(const TDelilaEvent& ev1, const TDelilaEvent& ev2)
// {
//  
//     // this will return true when second person
//     // has greater height. Suppose we have p1.height=5
//     // and p2.height=5.5 then the object which
//     // have max height will be at the top(or
//     // max priority)
//     return ev1.fTimeStamp < ev2.fTimeStamp;
// }   
   
   
//    DelilaSelectorTrigger& operator<(const TDelilaEvent& ev1, const TDelilaEvent& ev2);
//     bool friend operator<(const TDelilaEvent& ev1, const TDelilaEvent& ev2);//{return ev1.TimeStamp < ev2.TimeStamp;}
     // this will return true when second person
    // has greater height. Suppose we have p1.height=5
    // and p2.height=5.5 then the object which
    // have max height will be at the top(or
    // max priority)
 

};

#endif

#ifdef DelilaSelectorTrigger_cxx
void DelilaSelectorTrigger::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
   
   
  foutFile->cd();
  outputTree = new TTree("SelectedDelila","SelectedDelila");
//   outputTree->Branch("fTEventTS",&DelilaEventCS.fTimeStamp,"TimeStamp/l");
  outputTree->Branch("fTEventFS",&DelilaEventCS.fTimeStampFS,"TimeStamp/D");
  outputTree->Branch("fTimeBunch",&DelilaEventCS.TimeBunch,"fTimeBunch/D");
  outputTree->Branch("fEnergy",&DelilaEventCS.fEnergy,"Energy/F");
  outputTree->Branch("fEnergy_kev",&DelilaEventCS.EnergyCal,"Energy_kev/F");
//   outputTree->Branch("fEnergyDC_kev",&DelilaEventCS.fEnergyDC_kev,"Energy_kev/F");
  outputTree->Branch("fDomain",&DelilaEventCS.domain,"Domain/b");
  outputTree->Branch("fDetType",&DelilaEventCS.det_def,"def/b");
  outputTree->Branch("fCS",&DelilaEventCS.CS,"CS/b");
  outputTree->Branch("fTRG",&DelilaEventCS.trg,"Trigger/b");
  outputTree->Branch("fFold",&DelilaEventCS.fold,"Fold/b");
}

Bool_t DelilaSelectorTrigger::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef DelilaSelectorTrigger_cxx
