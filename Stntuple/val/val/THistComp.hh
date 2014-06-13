#ifndef THistComp_hh
#define THistComp_hh

#include "TH1.h"

class THistComp: public TObject {
protected:
  TH1*     fHist1;
  TH1*     fHist2;
  Double_t fKsProb;
  Double_t fNorm;
  TString  fHistory;
public:

  THistComp(TH1* Hist1=0, TH1* Hist2=0, Double_t KsProb=0) {
    fHist1  = Hist1;
    fHist2  = Hist2;
    fKsProb = KsProb;
    fNorm = 0.0;
  }

  ~THistComp() {}
//-----------------------------------------------------------------------------
//  accessors
//-----------------------------------------------------------------------------
  // 0 = no norm, <0=by area, >0 expected hist1 area/hist2 area
  void SetNorm(Double_t n = 0.0) { fNorm = n; }

  TH1*       GetHist1 () { return fHist1;  }
  TH1*       GetHist2 () { return fHist2;  }
  Double_t   GetKsProb() { return fKsProb; }
  Double_t   GetNorm()   { return fNorm; }
//-----------------------------------------------------------------------------
//  overloaded functions of TObject
//-----------------------------------------------------------------------------
  const char* GetName() const {
    static const char* name = "no_name";
    if (fHist1) {
      return fHist1->GetName();
    }
    else {
      return name;
    }
  }

  void SetHistory(TString h) { fHistory=h;}
  TString& GetHistory() {return fHistory;}

  virtual void Browse(TBrowser* b) { Draw("ep"); }

  virtual void Draw(Option_t* Opt="");  // *MENU*;
  void DrawEP() { Draw("ep"); }         // *MENU*;
  virtual void        Dump() const;    // *MENU*

  ClassDef(THistComp,2)
};



class TGoodHistComp: public THistComp {
public:
  TGoodHistComp(TH1* Hist1=0, TH1* Hist2=0, Double_t KsProb=0): 
    THistComp(Hist1,Hist2,KsProb) {
  }

  virtual ~TGoodHistComp(){}


  ClassDef(TGoodHistComp,1)
};


class TBadHistComp: public THistComp {
public:
  TBadHistComp(TH1* Hist1=0, TH1* Hist2=0, Double_t KsProb=0): 
    THistComp(Hist1,Hist2,KsProb) {
  }

  virtual ~TBadHistComp(){}


  ClassDef(TBadHistComp,1)
};

#endif
