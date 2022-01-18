//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 23 13:02:03 2021 by ROOT version 6.22/06
// from TTree DelilaSelector/Energy Station
// found on file: run1005_03.root
//////////////////////////////////////////////////////////

#pragma once

#ifndef DelilaSelector_h
#define DelilaSelector_h

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


class DelilaSelector : public TSelector {
public :
    TTreeReader     fReader;  //!the tree reader
    TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
    TFile          *outputFile;
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
    UShort_t        det_def;//0 - nothing; 1 - core; 2 - segment; 3 - CeBr; 4 - CsI; 5 - BGO1; 6 - BGO2; 9 - pulser
    float	        EnergyCal;
    float           EnergyDC;
    UShort_t        domain;
    int             cs_domain;
    UShort_t        channel;//ch daq
    UShort_t        core;
    UShort_t        segment;
    UShort_t        CS;//0 - no; 1 - yes
    double_t        Time;
    Float_t         theta;
    Float_t	        phi;
    double_t        bgo_time_diff;
   //int make_board_ID(){return fMod*100}
    TDelilaEvent(): domain(-1),channel(-1),fTimeStamp(0),fEnergy(-1),CS(0),cs_domain(0),Time(0),bgo_time_diff(-1){};
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
  std::deque<TDelilaEvent> gammagammaQu;
  std::deque<TDelilaEvent> gammagammaQu_CS;
  
//   std::deque<TDelilaEvent> bgo_Qu;
//   std::deque<float> enrergyQu;
  
  std::map<int,std::deque<TDelilaEvent>> waitingQu_gamma; //for CeBr
  std::map<int,std::deque<TDelilaEvent>> waitingQu_bgo; //for CS

//   std::deque<TDelilaEvent> eliadeQu_sorted;
  std::map<unsigned int, TDelilaDetector > LUT_DELILA;
  std::map<int, int > LUT_TA;

 
 
  TDelilaEvent DelilaEvent;  
  TDelilaEvent DelilaEventCS;
  TDelilaEvent lastDelilaEvent;  
  TDelilaEvent lastEliadeZeroEvent;  
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
  
 // std::map<UInt_t,TH1F*> hEnergy_raw;
  //std::map<UInt_t,TH1F*> hEnergy_cal;
 
//   std::map<UInt_t, TDelilaEvent> last_board_event;
 
  TH1F* hChannelHit;
  TH1F* hDomainHit;
  TH1F* hSegmentHit;
  TH1F* hDetTypeHit;
  TH2F* mDelila;//keV
  TH2F* mDelilaCS;//keV
  TH1F* hDelila;//keV
  TH1F* hDelilaCS;//keV
  TH2F* mDelila_raw;
  
  TH2F* mDelilaDC;//keV
  TH2F* mDelilaCS_DC;//keV
  TH1F* hDelilaDC;//keV
  TH1F* hDelilaCS_DC;//keV
  TH2F* mGammaGammaDC;
  TH2F* mGammaGammaCS_DC;


  
  TH2F* mDelilaTD;
  
  
  TH2F* mThetaPhi; 
  TH2F* mGammaGamma;
  TH2F* mTimeDiff_gg;
  TH1F* hMult_gg;
  TH2F* mGammaGammaCS;
  TH2F* mTimeDiff_gg_CS;
  TH1F* hMult_gg_CS;

  
  //Part For RoSPHERE
//   TH2F* mLaBr_raw;  
//   TH2F* mLaBr_kev;
//   TH1F* hLaBr_kev;
//   TH1F* hLaBrCS_kev;
  
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

    
  std::clock_t start;
  double duration;
 
  ULong64_t nevents;
  ULong64_t nevents_reset;
  int reset_counter;
 
  ULong64_t lastTime;
  
  TObjArray *toks;

   DelilaSelector(TTree * /*tree*/ =0) { }
   virtual ~DelilaSelector() { }
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
   virtual void  Print_ELIADE_LookUpTable();
   virtual void  Print_TimeAlignment_LookUpTable();
   virtual float CalibDet(float,int);
//    virtual int CheckTimeAlignment(int to_domain);
   virtual int GetCoincTimeCorrection(int dom1, int dom2);
   virtual void cs();
   virtual void gamma_gamma();
   virtual void gamma_gamma_cs(TDelilaEvent &ev_);
   virtual void time_alignment();
   
   

   ClassDef(DelilaSelector,0);
   
   
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
   
   
//    DelilaSelector& operator<(const TDelilaEvent& ev1, const TDelilaEvent& ev2);
//     bool friend operator<(const TDelilaEvent& ev1, const TDelilaEvent& ev2);//{return ev1.TimeStamp < ev2.TimeStamp;}
     // this will return true when second person
    // has greater height. Suppose we have p1.height=5
    // and p2.height=5.5 then the object which
    // have max height will be at the top(or
    // max priority)
 

};

#endif

#ifdef DelilaSelector_cxx
void DelilaSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
   
   
  outputFile->cd();
  outputTree = new TTree("SelectedDelila","SelectedDelila");
  outputTree->Branch("fTEventTS",&DelilaEventCS.fTimeStamp,"TimeStamp/l");
  outputTree->Branch("fTEventFS",&DelilaEventCS.fTimeStampFS,"TimeStamp/D");
  outputTree->Branch("fEnergy",&DelilaEventCS.EnergyCal,"Energy/F");
  outputTree->Branch("fDomain",&DelilaEventCS.domain,"Domain/b");
  outputTree->Branch("fDetType",&DelilaEventCS.det_def,"def/b");
  outputTree->Branch("fCS",&DelilaEventCS.CS,"CS/b");
}

Bool_t DelilaSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef DelilaSelector_cxx
// R G // 2 1
// W B // 3 0
// B0 G1 R2 W3

