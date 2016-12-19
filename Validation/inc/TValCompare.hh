#ifndef Tools_TValCompare_hh
#define Tools_TValCompare_hh

#include "TObject.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TObjArray.h"
#include "Validation/inc/TValPar.hh"
#include "Validation/inc/TValHistH.hh"
#include "Validation/inc/TValHistP.hh"
#include "Validation/inc/TValHistE.hh"


class TValCompare: public TObject {

public:

  TValCompare(TString File1="", TString File2="") {
    fFileN1 = File1;
    fFileN2 = File2;
    fFile1 = nullptr;
    fFile2 = nullptr;
    fVerbose = 1;
    fMinStat = -1;
    fMaxStat = 999;
  }

  ~TValCompare() {}

  TValPar&   GetPar()   { return fPar; }
  Int_t      GetVerbose() { return fVerbose; }
  Int_t      GetMinStat() { return fMinStat; }
  Int_t      GetMaxStat() { return fMaxStat; }

  void SetFile1(TString x) { fFileN1 = x; }
  void SetFile2(TString x) { fFileN2 = x; }

  virtual Int_t Analyze(Option_t* Opt="");
  // get a histogram based on searching name, path and title
  TValHist* GetHist(TString str);
  virtual void Report(Option_t* Opt="");
  virtual void Summary(Option_t* Opt="");
  virtual void Display(Option_t* Opt="");
  virtual void SaveAs(const char *filename="",Option_t *option="") const;
  //Int_t Write(const char *name=0, Int_t option=0, Int_t bufsize=0);
  //Int_t Write(const char *name=0, Int_t option=0, Int_t bufsize=0) const;

  virtual void  Delete(Option_t* Opt="");
  void SetVerbose(Int_t x) { fVerbose = x; }
  void SetMinStat(Int_t x) { fMinStat = x; }
  void SetMaxStat(Int_t x) { fMaxStat = x; }

protected:
  TString  fFileN1;
  TString  fFileN2;
  TFile*  fFile1;
  TFile*  fFile2;
  TValPar fPar;

  TObjArray fList; // list of histogram objects being compared

  Int_t fVerbose;
  Int_t fMinStat;
  Int_t fMaxStat;

  ClassDef(TValCompare,1)

};

#endif  /* Tools_TValCompare_hh*/
