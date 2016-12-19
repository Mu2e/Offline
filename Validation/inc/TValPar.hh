#ifndef Tools_TValPar_hh
#define Tools_TValPar_hh

#include "TObject.h"

class TValPar: public TObject {

public:

  TValPar() {
    Clear();
  }

  ~TValPar() {}

  Int_t      GetMode()   const { return fMode; }
  Int_t      GetIndependent()   const { return fIndep; }
  Double_t   GetScale1() const { return fScale1; }
  Double_t   GetScale2() const { return fScale2; }
  Int_t      GetUnder()  const { return fUnder; }
  Int_t      GetOver()   const { return fOver; }
  Double_t   GetLoose()  const { return fLoose; }
  Double_t   GetTight()  const { return fTight; }


    // 0 = no scaling
    // 1 = scale hist2 to integral of hist1
    // 2 = scale hist1 to scale1 and hist2 to scale2
  void SetMode(Int_t x)      { fMode = x; }
    // how to interpret the KS stat
    // 0 = samples should be identical
    // 1 = samples are from same parent distribution, but statistically indep
  void SetIndependent(Int_t x);
  void SetScale1(Double_t x) { fScale1 = x; }
  void SetScale2(Double_t x) { fScale2 = x; }
  void SetUnder(Int_t x)     { fUnder = x; }
  void SetOver(Int_t x)      { fOver = x; }
  void SetLoose(Double_t x)  { fLoose = x; }
  void SetTight(Double_t x)  { fTight = x; }

  virtual void  Clear(Option_t* Opt="");

protected:
  Int_t    fMode;
  Int_t    fIndep;
  Double_t fScale1;
  Double_t fScale2;
  Int_t    fUnder; //0=include underflow, 1=ignore it
  Int_t    fOver;  //0=include  overflow, 1=ignore it
  Float_t  fLoose; // loose threshold (0.99 or 0.001)
  Float_t  fTight; // tight threshold (0.999 or 0.01)
  Bool_t   fChanged; // set if loose and tight were set by user

  ClassDef(TValPar,1)

};

#endif  /* Tools_TValPar_hh*/
