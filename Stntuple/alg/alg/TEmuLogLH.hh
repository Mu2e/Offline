///////////////////////////////////////////////////////////////////////////////
// 2013-05-31: P.Murat this is a 2D calorimeter-based likelihood which 
//             distinguishes electrons and muons
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_inc_TEmuLogLH_hh
#define murat_inc_TEmuLogLH_hh

// class smooth_new;

#include "TObject.h"
#include "TH1.h"
#include "TH2.h"

class TEmuLogLH: public TObject {
public:

  struct PidData_t {
    double fPath;			// trajectory length inside the disk
    double fDt;				// delta(T) = T(cluster)-T(track)
    double fEp;				// E(cluster)/P(track)
    double fFrE2;			// cluster (e1+e2)/etot
    double fXs;				// slope/sig(slope)
    double fDe;				// DeDx probability
  };

  int            fNEpSlices;		// assume < 10, currently - 7
  float          fPath[10];
					// electron histograms
  TH2F*          fEleEpVsPath;
  TH1F*          fEleEpHist[10];
  TH1F*          fEleDtHist;
  TH1F*          fEleXsHist;
  TH1F*          fEleDeHist;		// DeDx probability
					// muon histograms
  TH2F*          fMuoEpVsPath;
  TH1F*          fMuoEpHist[10];
  TH1F*          fMuoDtHist;
  TH1F*          fMuoXsHist;
  TH1F*          fMuoDeHist;		// DeDx probability
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
  void InitEleDtHist(const char* Fn);
  void InitMuoDtHist(const char* Fn);

  void InitEleEpHist(const char* Fn);
  void InitMuoEpHist(const char* Fn);

  void InitEleDtHist(const TH1F* Hist);
  void InitMuoDtHist(const TH1F* Hist);

  void InitEleEpHist(const TH2F* Hist);
  void InitMuoEpHist(const TH2F* Hist);

  void InitEleXsHist(const char* Fn);
  void InitMuoXsHist(const char* Fn);
					// versions: "v4_2_4"
  int  Init(const char* Version);

  int  Init_v4_2_4();
//-----------------------------------------------------------------------------
// log(LH) of a given hypothesis is normally negative. 
// If the calculated likelihood is zero, the returned value of Log(LH) 
// is set to 999.
//-----------------------------------------------------------------------------
					// assume that there is no under/overflows
  void   SetEleXsHist(TH1* Hist);

  void   SetMuoXsHist(TH1* Hist);

  double LogLHCal(PidData_t* Data, int PdgCode);

  //  double LogLHTrk(TrkData_t* Data, int PdgCode);

  double LogLHDt  (PidData_t* Data, int PdgCode);
  double LogLHEp  (PidData_t* Data, int PdgCode);

  double LogLHXs  (double Xs, int PdgCode);

  double LogLHREp (PidData_t* Data);
  double LogLHRDt (PidData_t* Data);
  double LogLHRXs (double Xs);
					// log_lhr = log_lh(ele)-log_lh(muo)
  double LogLHRCal(PidData_t* Data);
  //  double LogLHRTrk(TrkData_t* Data);

  int ReadHistogram1D(const char* Fn, TH1F** Hist);
  int ReadHistogram2D(const char* Fn, TH2F** Hist);

  ClassDef (TEmuLogLH,0)
};

#endif
