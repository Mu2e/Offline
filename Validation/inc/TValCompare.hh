#ifndef Tools_TValCompare_hh
#define Tools_TValCompare_hh

#include "Offline/Validation/inc/TValHist2.hh"
#include "Offline/Validation/inc/TValHistE.hh"
#include "Offline/Validation/inc/TValHistH.hh"
#include "Offline/Validation/inc/TValHistP.hh"
#include "Offline/Validation/inc/TValPar.hh"
#include "TFile.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TString.h"

class TValCompare : public TObject {
 public:
  TValCompare(TString File1 = "", TString File2 = "") {
    fFileN1 = File1;
    fFileN2 = File2;
    fFile1 = nullptr;
    fFile2 = nullptr;
    fVerbose = 1;
    fMinStat = -1;
    fMaxStat = 999;
  }

  ~TValCompare() {}

  TValPar& GetPar() { return fPar; }
  Int_t GetVerbose() { return fVerbose; }
  Int_t GetMinStat() { return fMinStat; }
  Int_t GetMaxStat() { return fMaxStat; }

  void SetFile1(TString x) { fFileN1 = x; }
  void SetFile2(TString x) { fFileN2 = x; }

  // open and traverse the first file
  virtual Int_t GetDirs();
  virtual Int_t OneFile(Option_t* Opt = "");
  virtual Int_t Analyze(Option_t* Opt = "");
  // get a histogram based on searching name, path and title
  TValHist* GetHist(TString str);
  virtual void Report(Option_t* Opt = "");
  virtual void Summary(Option_t* Opt = "");
  virtual void Display(Option_t* Opt = "");
  virtual void SaveAs(const char* filename = "", Option_t* option = "") const;
  // save for one file option
  virtual void SaveAs1(const char* filename = "", Option_t* option = "") const;

  virtual void Delete(Option_t* Opt = "");
  void SetVerbose(Int_t x) { fVerbose = x; }
  void SetMinStat(Int_t x) { fMinStat = x; }
  void SetMaxStat(Int_t x) { fMaxStat = x; }

 protected:
  TString fFileN1;
  TString fFileN2;
  TFile* fFile1;
  TFile* fFile2;
  TValPar fPar;

  TObjArray fList;  // list of histogram objects being compared
  TObjArray fDirs;  // list of directories to analyze

  Int_t fVerbose;
  Int_t fMinStat;
  Int_t fMaxStat;

  ClassDef(TValCompare, 1)
};

#endif /* Tools_TValCompare_hh*/
