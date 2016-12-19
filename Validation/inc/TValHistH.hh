#ifndef Validation_TValHistH_hh
#define Validation_TValHistH_hh

#include "Validation/inc/TValHist.hh"
#include "TH1.h"
#include "TString.h"


class TValHistH: public TValHist {

public:

  TValHistH(TH1* Hist1=0, TH1* Hist2=0) {
    Clear();
    fHist1  = Hist1;
    fHist2  = Hist2;
  }

  virtual ~TValHistH() {}

  TH1*       GetHist1 () { return fHist1;  }
  TH1*       GetHist2 () { return fHist2;  }

  virtual const char* GetName() const {
    static const char* name = "TValHistH name";
    if (fHist1) { return fHist1->GetName();
    } else { return name; }
  }
  virtual const char* GetTitle() const {
    static const char* name = "TValHistH title";
    if (fHist1) { return fHist1->GetTitle();
    } else { return name; }
  }

  void SetHist1(TH1* x)  { fHist1 = x; ClearProb(); }
  void SetHist2(TH1* x)  { fHist2 = x; ClearProb(); }

  virtual Int_t Analyze(Option_t* Opt="");
  virtual void  Summary(Option_t* Opt="");
  virtual void  Draw(Option_t* Opt="");
  virtual void  Dump() const;

  virtual void  Clear(Option_t* Opt="");

protected:
  TH1*     fHist1;
  TH1*     fHist2;

  Double_t fSum1;
  Double_t fSum2;
  Double_t fNorm2;

  ClassDef(TValHistH,1)

};

#endif  /* Validation_TValHistH_hh*/
