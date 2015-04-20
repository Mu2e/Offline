///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_TValidationModule_hh
#define Stntuple_ana_TValidationModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/obj/TDiskCalorimeter.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

class TValidationModule: public TStnModule {
public:

  struct TrackPar_t {
    int     fNHPl;
    int     fNEPl;
    int     fNDPl;
    float   fDpF ;    // tracker-only resolution
    float   fDp0 ;
    float   fDp2 ;
    float   fDpFSt;
    double  fDioWt;

    double  fEcl;
    double  fEp;
    double  fDx;
    double  fDy;
    double  fDz;
    double  fDt;
    double  fDu;			// rotated residuals
    double  fDv;
    double  fChi2Match;
    double  fPath;
  };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct CaloHist_t {
    TH1F*    fVaneID;		       // per crystal hit
    TH1F*    fEnergy  [4];
    TH1F*    fTime    [4];
    TH1F*    fNHits   [4];
    TH1F*    fRadius  [4];
    TH1F*    fRadiusWE[4];
    TH1F*    fE700    [4];
    TH1F*    fT700    [4];
    TH1F*    fN700    [4];
    TH1F*    fR700    [4];
    TH1F*    fRWE700  [4];
  };

  struct ClusterHist_t {
    TH1F*    fVaneID;
    TH1F*    fEnergy;
    TH1F*    fT0;
    TH1F*    fRow;
    TH1F*    fCol;
    TH1F*    fX;
    TH1F*    fY;
    TH1F*    fZ;
    TH1F*    fR;
    TH1F*    fNCr0;			// all clustered
    TH1F*    fNCr1;			// above 1MeV
    TH1F*    fYMean;
    TH1F*    fZMean;
    TH1F*    fSigY;
    TH1F*    fSigZ;
    TH1F*    fSigR;
    TH1F*    fFrE1;
    TH1F*    fFrE2;
    TH1F*    fSigE1;
    TH1F*    fSigE2;
  };

  struct EventHist_t {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleMom;
    TH1D*    fDioMom;
    TH1F*    fEleCosTh;
    TH1F*    fNClusters;
    TH1F*    fNTracks;
    TH1F*    fNStrawHits[2];
    TH1F*    fNGoodSH;
    TH1F*    fDtClT;
    TH1F*    fEMax;			// energy of the first reco cluster
    TH1F*    fDtClS;
    TH1F*    fSHTime;
    TH1F*    fNHyp;
    TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency
    TH1F*    fNGenp;                    // N(particles in GENP block)

    TH1F*    fNCaloCrystalHits[2];
    TH2F*    fNCaloHitsVsVane[2];
    TH2F*    fNCaloHitsVsRow[2];
    TH2F*    fNCaloHitsVsCol[2];
    // calorimeter hit histograms

    TH1F*    fETot        [4];            // total energy/event 
    TH2F*    fECrVsR      [4];            // total energy_per_crystal/event vs radius
    TH2F*    fNCrVsR      [4];            // total energy_per_crystal/event vs radius

    TH2F*    fNCrystalHitsVsR[4];            //
    TH2F*    fNHitCrystalsVsR[4];            //

    TH1F*    fNHitCrystalsTot;
    TH1F*    fECal;
    TH1F*    fECalOverEKin;
      
  };

  struct TrackHist_t {
    TH1F*    fP[3];			// total momentum, 3 hists with different binning
    TH1F*    fP0;
    TH1F*    fP2;
    TH1F*    fPt;
    TH1D*    fPDio;                     // momentum dist weighted with the DIO weight
    TH1F*    fFitMomErr;
    TH1F*    fPFront;
    TH1F*    fDpFront;
    TH1F*    fDpFront0;
    TH1F*    fDpFront2;
    TH2F*    fDpFVsZ1;
    TH1F*    fPStOut;
    TH1F*    fDpFSt;			// P(TT_Hollow) - P(ST_Out)
    TH1F*    fCosTh;
    TH1F*    fChi2;
    TH1F*    fNDof;
    TH1F*    fChi2Dof;
    TH1F*    fChi2DofC;
    TH1F*    fNActive;
    TH1F*    fT0;
    TH1F*    fT0Err;
    TH1F*    fQ;
    TH1F*    fFitCons[2];		// fit consistency (0 to 1)
    TH1F*    fD0;
    TH1F*    fZ0;
    TH1F*    fTanDip;
    TH1F*    fResid;
    TH1F*    fAlgMask;
					// matching histograms
    TH1F*    fNClusters;
    TH1F*    fVaneID;
    TH1F*    fXCal;
    TH1F*    fYCal;
    TH1F*    fZCal;
    TH1F*    fXTrk;
    TH1F*    fYTrk;
    TH1F*    fZTrk;
    TH1F*    fRTrk;
    TH1F*    fDt;			// track-cluster residuals
    TH1F*    fChi2Match;
    TH1F*    fDt_eMinus;
    TH1F*    fDt_ePlus;
    TH1F*    fDt_muMinus;
    TH1F*    fDt_muPlus;
    TH1F*    fDx;
    TH1F*    fDy;
    TH1F*    fDz;
    TH1F*    fDu;
    TH1F*    fDv;
    TH2F*    fDvVsDu;
    TH1F*    fPath;
    TH2F*    fDuVsPath;
    TH2F*    fDucVsPath;
    TH2F*    fDvVsPath;
    TH2F*    fDvcVsPath;
    TH2F*    fDtVsPath;
    TH2F*    fDuVsTDip;
    TH2F*    fDvVsTDip;
    TH1F*    fZ1;
    TH1F*    fECl;
    TH1F*    fEClEKin;
    TH1F*    fEp;
    TH2F*    fEpVsPath;
    TH1F*    fEp_eMinus;
    TH1F*    fEp_ePlus;
    TH1F*    fEp_muMinus;
    TH1F*    fEp_muPlus;
    TH2F*    fNHVsStation;
    TH2F*    fNHVsNSt;

    TH1F*    fRSlope;
    TH1F*    fXSlope;
					// likelihoods
    TH2F*    fEpVsDt;
    TH1F*    fEleLogLHCal;
    TH1F*    fMuoLogLHCal;
    TH1F*    fLogLHRCal;
    TH1F*    fLogLHRDeDx;
    TH1F*    fLogLHRXs;
    TH1F*    fLogLHRTrk;
    TH1F*    fLogLHR;
					// MC truth
    TH1F*    fPdgCode;	                // PDG code of the particle produced most hits
    TH1F*    fFrGH;			// fraction of hits produced by the particle

    TH2F*    fNEPlVsNHPl;
    TH2F*    fNDPlVsNHPl;
    TH2F*    fChi2dVsNDPl;
    TH2F*    fDpFVsNDPl;

    TH1F*    fFrE1;
    TH1F*    fFrE2;
  };

  struct GenpHist_t {
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fGenID;			// 
    TH1F*    fZ0;			// 
    TH1F*    fT0;			// 
    TH1F*    fR0;			// 
    TH1F*    fP;			// 
    TH1F*    fCosTh;			// 
  };
					// histograms for the simulated CE
  struct SimpHist_t {
    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
    TH1F*    fNStrawHits;
  };

  struct TrackEffHist_t {
    TH1F*    fPtMc;			// denominator
    TH1F*    fPtReco;			// numerator
  };

//-----------------------------------------------------------------------------
//  fTrackHist[ 0]: all the tracks
//  fTrackHist[ 1]: all the tracks Pt > PtMin and within the fiducial
//  fTrackHist[ 2]: [1]+ matched to MC
//  fTrackHist[ 3]: [1]+ not matched to MC
//  fTrackHist[ 4]: [2]+ inside  the jet
//  fTrackHist[ 5]: [3]+ inside  the jet
//  fTrackHist[ 6]: [2]+ outside the jet
//  fTrackHist[ 7]: [3]+ outside the jet
//  fTrackHist[ 8]:
//  fTrackHist[ 9]:
//  fTrackHist[10]:
//  fTrackHist[11]: tracks with pt.10 inside the COT
//
//  fTrackEffHist[0]
//  fTrackEffHist[1]
//  fTrackEffHist[2]
//  fTrackEffHist[3]
//-----------------------------------------------------------------------------
  enum { kNEventHistSets   = 100 };
  enum { kNTrackHistSets   = 400 };
  enum { kNClusterHistSets = 100 };
  enum { kNCaloHistSets    = 100 };
  enum { kNGenpHistSets    = 100 };
  enum { kNSimpHistSets    = 100 };

  struct Hist_t {
    TH1F*          fCrystalR[2];	          // crystal radius
    EventHist_t*   fEvent  [kNEventHistSets];
    TrackHist_t*   fTrack  [kNTrackHistSets];
    ClusterHist_t* fCluster[kNClusterHistSets];
    CaloHist_t*    fCalo   [kNCaloHistSets];
    GenpHist_t*    fGenp   [kNGenpHistSets];
    SimpHist_t*    fSimp   [kNSimpHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*   fTrackBlock;
  TStnClusterBlock* fClusterBlock;
  TCalDataBlock*    fCalDataBlock;
  TStrawDataBlock*  fStrawDataBlock;
  TGenpBlock*       fGenpBlock;
  TSimpBlock*       fSimpBlock;
					// additional track parameters (assume ntracks < 20)
  TrackPar_t        fTrackPar[20];
					// histograms filled
  Hist_t            fHist;
					// cut values
  double            fPtMin;

  TGenParticle*     fParticle;		// electron or muon
  int               fPdgCode;		// determines which one
  int               fGeneratorCode;      

  TSimParticle*     fSimp;
  double            fEleE;		// electron energy

  int               fCalorimeterType;

  int               fNClusters;
  int               fNTracks[10];
  int               fNGoodTracks;
  int               fNMatchedTracks;
  int               fNStrawHits;
  int               fNCalHits;
  int               fNGenp;		// N(generated particles)

  int               fNHyp;
  int               fBestHyp[10];
  int               fFillDioHist;
					// fTrackNumber[i]: track number, 
					// corresponding to OBSP particle #i
					// or -1
  TStnArrayI        fTrackNumber;

  TStnTrack*        fTrack;
  TStnCluster*      fCluster;

  TDiskCalorimeter* fDiskCalorimeter;

  TStnTrackID*      fTrackID;
  TEmuLogLH*        fLogLH;

  double            fMinT0;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TValidationModule(const char* name="Validation", const char* title="Validation");
  ~TValidationModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
  TStnTrackBlock*    GetTrackBlock  () { return fTrackBlock;   }
  TStnClusterBlock*  GetClusterBlock() { return fClusterBlock; }

  TStnTrackID*       GetTrackID     () { return fTrackID; }
  TEmuLogLH*         GetLogLH       () { return fLogLH; }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void               SetFillDioHist  (int YesNo) { fFillDioHist   = YesNo; }
  void               SetPdgCode      (int Code ) { fPdgCode       = Code ; }
  void               SetGeneratorCode(int Code ) { fGeneratorCode = Code ; }
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int     BeginJob();
  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookCaloHistograms    (CaloHist_t*    Hist, const char* Folder);
  void    BookClusterHistograms (ClusterHist_t* Hist, const char* Folder);
  void    BookGenpHistograms    (GenpHist_t*    Hist, const char* Folder);
  void    BookEventHistograms   (EventHist_t*   Hist, const char* Folder);
  void    BookSimpHistograms    (SimpHist_t*    Hist, const char* Folder);
  void    BookTrackHistograms   (TrackHist_t*   Hist, const char* Folder);

  void    FillEventHistograms    (EventHist_t* Hist);
  void    FillCaloHistograms     (CaloHist_t*    Hist, TStnCrystal*  Crystal);
  void    FillClusterHistograms  (ClusterHist_t* Hist, TStnCluster*  Cluster);
  void    FillGenpHistograms     (GenpHist_t*    Hist, TGenParticle* Genp   );
  void    FillSimpHistograms     (SimpHist_t*    Hist, TSimParticle* Simp   );
  void    FillTrackHistograms    (TrackHist_t*   Hist, TStnTrack*    Trk    );

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TValidationModule,0)
};

#endif
