#ifndef Validation_TValHistP_hh
#define Validation_TValHistP_hh

#include "Validation/inc/TValHist.hh"
#include "TProfile.h"
#include "TString.h"


class TValHistP: public TValHist {

public:

  TValHistP(TProfile* Prof1=0, TProfile* Prof2=0) {
    Clear();
    fProf1  = Prof1;
    fProf2  = Prof2;
  }

  virtual ~TValHistP() {}

  TProfile*       GetProf1 () { return fProf1;  }
  TProfile*       GetProf2 () { return fProf2;  }

  virtual const char* GetName() const {
    static const char* name = "TValHistP name";
    if (fProf1) { return fProf1->GetName();
    } else { return name; }
  }
  virtual const char* GetTitle() const {
    static const char* name = "TValHistP title";
    if (fProf1) { return fProf1->GetTitle();
    } else { return name; }
  }

  void SetProf1(TProfile* x)  { fProf1 = x; ClearProb(); }
  void SetProf2(TProfile* x)  { fProf2 = x; ClearProb(); }

  virtual Int_t Analyze(Option_t* Opt="");
  virtual void  Summary(Option_t* Opt="");
  virtual void  Draw(Option_t* Opt="");
  virtual void  Dump() const;

  virtual void  Clear(Option_t* Opt="");

protected:
  TProfile*     fProf1;
  TProfile*     fProf2;

  Double_t fSum1;
  Double_t fSum2;

  ClassDef(TValHistP,1)

};

#endif  /* Validation_TValHistP_hh*/
