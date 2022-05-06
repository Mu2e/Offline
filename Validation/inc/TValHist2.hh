#ifndef Validation_TValHist2_hh
#define Validation_TValHist2_hh

#include "Offline/Validation/inc/TValHist.hh"
#include "TH2.h"
#include "TString.h"

class TValHist2 : public TValHist {
 public:
  TValHist2(TH2* Hist1 = 0, TH2* Hist2 = 0) {
    Clear();
    fHist1 = Hist1;
    fHist2 = Hist2;
  }

  virtual ~TValHist2() {}

  TH2* GetHist1() { return fHist1; }
  TH2* GetHist2() { return fHist2; }

  virtual const char* GetName() const {
    static const char* name = "TValHist2 name";
    if (fHist1) {
      return fHist1->GetName();
    } else {
      return name;
    }
  }
  virtual const char* GetTitle() const {
    static const char* name = "TValHist2 title";
    if (fHist1) {
      return fHist1->GetTitle();
    } else {
      return name;
    }
  }

  void SetHist1(TH2* x) {
    fHist1 = x;
    ClearProb();
  }
  void SetHist2(TH2* x) {
    fHist2 = x;
    ClearProb();
  }

  virtual Int_t Analyze(Option_t* Opt = "");
  virtual void Summary(Option_t* Opt = "");
  virtual void Draw(Option_t* Opt = "");
  virtual void Dump() const;

  virtual void Clear(Option_t* Opt = "");

 protected:
  TH2* fHist1;
  TH2* fHist2;

  Double_t fSum1;
  Double_t fSum2;
  Double_t fNorm2;

  ClassDef(TValHist2, 1)
};

#endif /* Validation_TValHist2_hh*/
