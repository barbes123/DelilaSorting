#ifndef LABR_TREE_EVENT_H_INCLUDED
#define LABR_TREE_EVENT_H_INCLUDED

#include <map>
#include <vector>
#include <climits>
#include <TObject.h>


class LaBrTreeEvent : public TObject{
 public:
  LaBrTreeEvent(){;}
  ~LaBrTreeEvent(){;}


  virtual void SetEnergy       (float v) {fEnergy       =v;}
  virtual void SetEnergyDC     (float v) {fEnergyDC     =v;}
  virtual void SetTheta        (float v) {fTheta        =v;}
  virtual void SetTimeDiff     (float v) {fTimeDiff     =v;}
  virtual void SetDomain       (int   v) {fDomain       =v;}
  virtual void SetChannel      (int   v) {fChannel      =v;}
  virtual void SetTS           (UInt_t v){fTimeStamp    =v;}
//   virtual void SetRing         (UInt_t v){fRing         =v;}
//   virtual void SetPartCoinc    (Int_t v) {fPartCoinc    =v;}
//   virtual void SetGammaCoinc   (Int_t v) {fGammaCoinc   =v;}
  
  virtual Float_t GetEnergy       () {return fEnergy       ;}
  virtual Float_t GetEnergyDC     () {return fEnergyDC     ;}
  virtual Float_t GetTheta        () {return fTheta        ;}
  virtual Int_t   GetDomain       () {return fDomain       ;}
  virtual Int_t   GetChannel      () {return fChannel       ;}
  virtual UInt_t  GetTS           () {return fTimeStamp    ;}
//   virtual Int_t   GetRing         () {return fRing         ;}
//   virtual Int_t   GetPartCoinc    () {return fPartCoinc    ;}
//   virtual Int_t   GetGammaCoinc   () {return fGammaCoinc   ;}

	
 private:
  Float_t fEnergy       ;
  Float_t fEnergyDC     ;
  Float_t fTheta        ;
  Float_t fTimeDiff     ;
  Int_t   fDomain       ;
  Int_t   fChannel      ;
  UInt_t  fTimeStamp    ;
//   Int_t   fRing         ;
//   Int_t   fPartCoinc    ;
//   Int_t   fGammaCoinc   ;
  ClassDef(LaBrTreeEvent,1)
};

#endif
