///////////////////////////////////////////////////////////////////////////////
// 2013-05-31: P.Murat this is a 2D calorimeter-based likelihood which 
//             distinguishes electrons and muons
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_inc_TEmuLogLH_hh
#define murat_inc_TEmuLogLH_hh

// class smooth_new;

#include "TObject.h"
#include "TH1.h"

class TEmuLogLH: public TObject {
public:

  struct CalData_t {
    double fDt;				// delta(T) = T(cluster)-T(track)
    double fEp;				// E(cluster)/P(track)
  };
					// electron histograms
  struct TrkData_t {
    double fXs;				// slope/sig(slope)
    double fDe;				// DeDx probability
  };
					// electron histograms
  TH1*          fEleDtHist;
  TH1*          fEleEpHist;
  TH1*          fEleXsHist;
  TH1*          fEleDeHist;		// DeDx probability
					// muon histograms
  TH1*          fMuoDtHist;
  TH1*          fMuoEpHist;
  TH1*          fMuoXsHist;
  TH1*          fMuoDeHist;		// DeDx probability
					// this part may not be needed
  // smooth_new*   fEleDtFunc;
  // smooth_new*   fEleEpFunc;

  // smooth_new*   fMuoDtFunc;
  // smooth_new*   fMuoEpFunc;
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  TEmuLogLH();
  TEmuLogLH(const char* FnEle, const char* FnMuo);

  ~TEmuLogLH();
//-----------------------------------------------------------------------------
// initialization
//-----------------------------------------------------------------------------
  void InitEleDtHist(TH1*& Hist);
  void InitMuoDtHist(TH1*& Hist);
  void InitEleEpHist(TH1*& Hist);
  void InitMuoEpHist(TH1*& Hist);
  int  Init();
//-----------------------------------------------------------------------------
// log(LH) of a given hypothesis is normally negative. 
// If the calculated likelihood is zero, the returned value of Log(LH) 
// is set to 999.
//-----------------------------------------------------------------------------
					// assume that there is no under/overflows
  void   SetEleDtHist(TH1* Hist);
  void   SetEleEpHist(TH1* Hist);
  void   SetEleXsHist(TH1* Hist);
  //  void   SetEleDeHist(TH1* Hist);

  void   SetMuoDtHist(TH1* Hist);
  void   SetMuoEpHist(TH1* Hist);
  void   SetMuoXsHist(TH1* Hist);
  //  void   SetMuoDeHist(TH1* Hist);

  double LogLHCal(CalData_t* Data, int PdgCode);
  //  double LogLHTrk(TrkData_t* Data, int PdgCode);

  double LogLHDt  (double Ep, int PdgCode);
  double LogLHEp  (double Ep, int PdgCode);
  double LogLHXs  (double Ep, int PdgCode);
  //  double LogLHDe  (double Ep, int PdgCode);

  double LogLHREp (double Ep);
  double LogLHRDt (double Dt);
  double LogLHRXs (double Xs);
  //  double LogLHRDe (double De);
					// log_lhr = log_lh(ele)-log_lh(muo)
  double LogLHRCal(CalData_t* Data);
  //  double LogLHRTrk(TrkData_t* Data);

  ClassDef (TEmuLogLH,0)
};

#endif
