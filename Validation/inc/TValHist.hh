#ifndef Tools_TValHist_hh
#define Tools_TValHist_hh

#include "Validation/inc/TValPar.hh"
#include "TH1.h"
#include "TString.h"


class TValHist: public TObject {

public:

  TValHist() {
    ClearB();
  }

  virtual ~TValHist() {}

  Double_t   GetKsProb() { return fKsProb; }
  Double_t   GetFrProb() { return fFrProb; }
  Bool_t     GetDiff()   { return fDiff; }
  Int_t      GetStatus() { return fStatus; }
  TValPar&   GetPar()    { return fPar;}
  TString&   GetTag()    { return fTag;}
  Float_t    GetFontScale() { return fFontScale; }

  virtual const char* GetName() const=0;
  virtual const char* GetTitle() const=0;

  void SetPar(TValPar& x){ fPar = x; }
  void SetTag(TString& x){ fTag = x; }
  void SetFontScale(Float_t x){ fFontScale = x; }

  virtual Int_t Analyze(Option_t* Opt="")=0;
  virtual void  Summary(Option_t* Opt="")=0;
  virtual void  Draw(Option_t* Opt="")=0;
  virtual void  Dump() const=0;

  virtual void  ClearB(Option_t* Opt="");
  virtual void  ClearProb() { fKsProb=0.0; fFrProb=0.0; fDiff=true; }
  virtual void  Clear(Option_t* Opt="")=0;

protected:
  Double_t fKsProb;  // KS test, ~1 is good
  Double_t fFrProb;  // fraction difference test, ~1 is good
  Bool_t   fDiff;    // true if any difference at all

  Int_t fStatus;

  TValPar fPar;
  TString fTag;

  Float_t fFontScale;

  ClassDef(TValHist,1)

};

#endif  /* Tools_TValHist_hh*/
