#ifndef Validation_TValHistE_hh
#define Validation_TValHistE_hh

#include "Validation/inc/TValHist.hh"
#include "TEfficiency.h"
#include "TString.h"


class TValHistE: public TValHist {

public:

  TValHistE(TEfficiency* Hist1=0, TEfficiency* Hist2=0) {
    Clear();
    fEff1  = Hist1;
    fEff2  = Hist2;
  }

  virtual ~TValHistE() {}

  TEfficiency*       GetHist1 () { return fEff1;  }
  TEfficiency*       GetHist2 () { return fEff2;  }

  virtual const char* GetName() const {
    static const char* name = "TValHistE name";
    if (fEff1) { return fEff1->GetName();
    } else { return name; }
  }
  virtual const char* GetTitle() const {
    static const char* name = "TValHistE title";
    if (fEff1) { return fEff1->GetTitle();
    } else { return name; }
  }

  void SetHist1(TEfficiency* x)  { fEff1 = x; ClearProb(); }
  void SetHist2(TEfficiency* x)  { fEff2 = x; ClearProb(); }

  virtual Int_t Analyze(Option_t* Opt="");
  virtual void  Summary(Option_t* Opt="");
  virtual void  Draw(Option_t* Opt="");
  virtual void  Dump() const;

  virtual void  Clear(Option_t* Opt="");

protected:
  TEfficiency*     fEff1;
  TEfficiency*     fEff2;

  Double_t fSum1;
  Double_t fSum2;

  ClassDef(TValHistE,1)

};

#endif  /* Validation_TValHistE_hh*/
