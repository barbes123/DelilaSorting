#ifndef ELISSA_TREE_EVENT_H_INCLUDED
#define ELISSA_TREE_EVENT_H_INCLUDED

#include <map>
#include <vector>
#include <climits>
#include <TObject.h>


class ElissaTreeEvent : public TObject{
 public:
  ElissaTreeEvent(){;}
  ~ElissaTreeEvent(){;}

  virtual void SetEnergy       (float v) {fEnergy       =v;}
  virtual void SetTheta        (float v) {fTheta        =v;}
  virtual void SetTimeDiff     (float v) {fTimeDiff     =v;}
  virtual void SetDomain       (int   v) {fDomain       =v;}
  virtual void SetChannel      (int   v) {fChannel      =v;}
  virtual void SetTS           (UInt_t v){fTimeStamp    =v;}
  virtual void SetAmax         (float  v){fAmax         =v;}
  
  virtual Float_t GetEnergy       () {return fEnergy       ;}
  virtual Float_t GetTheta        () {return fTheta        ;}
  virtual Int_t   GetDomain       () {return fDomain       ;}
  virtual Int_t   GetChannel      () {return fChannel       ;}
  virtual UInt_t  GetTS           () {return fTimeStamp    ;}
	
 private:
  Float_t fEnergy       ;
  Float_t fTheta        ;
  Float_t fTimeDiff     ;
  Int_t   fDomain       ;
  Int_t   fChannel      ;
  UInt_t  fTimeStamp    ;
  float   fAmax; //add by saka

  ClassDef(ElissaTreeEvent,1)
};

#endif
