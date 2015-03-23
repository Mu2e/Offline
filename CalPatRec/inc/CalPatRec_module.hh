///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_CalPatRec_module
#define CalPatRec_CalPatRec_module

#ifdef __GCCXML__A
namespace art {
  class EDProducer;
  class Run;
  class Event;
};
#else
#  include "art/Framework/Core/EDProducer.h"
#  include "art/Framework/Principal/Event.h"
#endif

// data
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// BaBar
#include "BaBar/BaBar.hh"
#include "BaBar/BbrStringUtils.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "TrkBase/HelixParams.hh"
#include "TrkBase/TrkPoca.hh"
// #include "KalmanTests/inc/KalFitMC.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"
#include "CalPatRec/inc/CalTimePeak.hh"
#include "CalPatRec/inc/Doublet.hh"

#include "TROOT.h"
#include "TFolder.h"
#include "CalPatRec/inc/KalFitHack.hh"
#include "CalPatRec/inc/HelixFitHack.hh"
#include "CalPatRec/inc/THackData.hh"

// Mu2e
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
#include "TrkPatRec/inc/PayloadSaver.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"
#include "TFolder.h"

//#include "TStopwatch.h"
// #include "TSpectrum.h"
// #include "TSpectrum2.h"
// #include "TSpectrum3.h"
// #include "TMVA/Reader.h"
// boost
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>

class Ref;
class THackData;

namespace fhicl {
  class ParameterSet;
}

namespace mu2e {

  class CalPatRec : public art::EDProducer {
  public:
    struct Hist_t {
      TH1F* _cutflow;	     // diagnostic flow
      TH1F* _hTpeaks;	     // number of time-peaks found per event
      TH1F* _hNfitIter;	     // number of call to kalman addHits
			     //     TH1F* _hTfit[2];     //time spent on kalman filter per event
			     //     TH1F* _hTtot;     //total time spent per event
      TH1F* _hdfdzmode;	     // distribution of the loop index where the findTrack search converged
      TH1F* _hradius;        // radius of the helix used by findTrack
      TH1F* _hphi0;          // phi0 of the helix used by findTrack
      TH1F* _htanlambda;     // tanLambda (tangent of the pitch angle) of the helix
      TH1F* _hdfdz;          // dfdz of the theretical helix. dfdz = tanLambda/radius
      TH1F* _hdist;
      TH1F* _hdz;
      TH1F* _hNpoints;
      TH1F* _hchi2;
      TH2F* _hdistvsdz;
      
      TH1F* _hkdfdzmode;     // distribution of the loop index where the findTrack search converged
      TH1F* _hkradius[2];    // radius of the helix used by findTrack
      TH1F* _hkphi0;         // phi0 of the helix used by findTrack
      TH1F* _hktanlambda;    // tanLambda (tangent of the pitch angle) of the theretical helix
      TH1F* _hkdfdz[2];      // dfdz of the theretical helix. dfdz = tanLambda/radius
      TH1F* _h0mode;
      TH1F* _hk0mode;
      TH1F* _hkdist;
      TH1F* _hkdz;
      TH1F* _hkNpoints;
      TH1F* _hkchi2;
      TH2F* _hkdistvsdz[2];
      
      TH1F* _hdrw  [2];
      TH1F* _hchi2w[2];
      
      TH1F* _hchi2zphi[2];
      
      TH1F* _hseeddoca[3];
      TH1F* _hseeddr  [2];
      TH1F* _hseeddfdz[2];
      TH1F* _hNpointsSeed[2]; //
      TH1F* _hkaldoca[2];
//-----------------------------------------------------------------------------
// histograms for doublets 0:all, 1:OS, 2:SS
//-----------------------------------------------------------------------------
      TH1F* _ndoublets[3];      // N(doublets)
      TH1F* _ndstraws [3];      // N(straws) per "doublet"
      TH1F* _tslope   [3];      // track slope dx/dz wrt the plane
      TH1F* _dsbest   [3];      // delta (slope) best
      TH1F* _dsnext   [3];      // delta (slope) next to the best

      TH1F* _ntracks;
//-----------------------------------------------------------------------------
// time diag histograms (inherited)
//-----------------------------------------------------------------------------
      TH1F *ctsp, *rtsp, *ttsp, *ltsp, *tdtsp;
    };

    Ref*    _ref;

  protected:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    //    TStopwatch*   fStopwatch;

    unsigned     _iev;
					// configuration parameters
    int          _diag; 
    int          _debug;
    int          _printfreq;
    bool         _addhits; 
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string  _shLabel;
    std::string  _shpLabel;
    std::string  _shfLabel;
    std::string  _ccmLabel; // caloClusterModuleLabel

    std::string  _dtspecpar;

    StrawHitFlag     _tsel, _hsel, _addsel, _ksel;
    StrawHitFlag     _bkgsel, _addbkg;
    double           _maxedep;
    int              fUseDoublets;
    double           _mindt;
    double           _maxdt;
    double           _maxdtmiss;
    double           _fbf;
					// time spectrum parameters
    //    bool             _findtpeak;
    unsigned         _maxnpeak;
    int              _minnhits;
    double           _tmin;
    double           _tmax;
    double           _tbin;
    unsigned         _nbins;
    double           _minClusterEnergy;	// min seed energy
    int              _minClusterSize;   // min size of the seeding cluster
    double           _ymin;
    double           _1dthresh;
    double           _pitchAngle;
					// outlier cuts
    double           _maxseeddoca;
    double           _maxhelixdoca;
    double           _maxadddoca;
    double           _maxaddchi;
    TrkParticle      _tpart;	        // particle type being searched for
    TrkFitDirection  _fdir;		// fit direction in search
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;
    const CaloClusterCollection*          _ccCollection;
    const StepPointMCCollection*          _stepcol;
    const PtrStepPointMCVectorCollection* _listOfMCStrawHits;

    StrawHitFlagCollection*               _flags;    // flags are not const - for a reason ?
					
    HelixFitHack             _hfit;	// robust helix fitter
    KalFitHack               _seedfit;  // Kalman filter config for the Seed fit ( fit using hit wires)
    KalFitHack               _kfit;     // full-blown Kalman filter
    CalTimePeakCollection*   _tpeaks;   // cache of time peaks
    std::string              _iname;	// data instance name

    PayloadSaver             _payloadSaver;

    XYZPHackVector           _index;
    int                      _nindex;

    const TTracker*          _tracker;  // straw tracker geometry
    TFolder*                 _folder;
					// MC tools
    //    KalFitMC _kfitmc;
//-----------------------------------------------------------------------------
// strawhit tuple variables
//-----------------------------------------------------------------------------
    TTree*   _shdiag;
    Int_t    _eventid;
    //    threevec _shp;
    Float_t  _edep;
    Float_t  _time, _rho;
    Int_t    _nmcsteps;
    Int_t    _mcnunique,_mcnmax;
    Int_t    _mcpdg,_mcgen,_mcproc;
    //    threevec _mcshp, _mcop;
    Float_t  _mcshlen;
    Float_t  _mcedep,_mcemax;
    Float_t  _pdist,_pperp,_pmom;
    Float_t  _mctime;
    Int_t    _esel,_rsel, _timesel,  _delta, _stereo, _isolated;
    Int_t    _device, _sector, _layer, _straw;
    Int_t    _ishpeak, _ntpeak, _nshtpeak;
    //    Float_t  _shtpeak;
    Float_t  _shpres, _shrres, _shchisq, _shdt, _shdist;
    Float_t  _shmct0, _shmcmom, _shmctd;
					// fit tuple variables
    Int_t    _nadd,_ipeak;
    //    Float_t  _hcx, _hcy, _hr, _hdfdz, _hfz0;
    Float_t  _mccx, _mccy, _mcr, _mcdfdz, _mcfz0;
    Int_t    _helixfail,_seedfail,_kalfail;
    //    helixpar _hpar,_spar;
    //    helixpar _hparerr,_sparerr;
    Int_t    _snhits, _snactive, _sniter, _sndof, _snweediter;
    Float_t  _schisq, _st0;
    Int_t    _nchit;
    Int_t    _npeak, _nmc;
    Float_t  _tpeak;
    int      _ntracks;
//-----------------------------------------------------------------------------
// hit filtering tuple variables
//-----------------------------------------------------------------------------
    std::vector<TrkHitFilter> _sfilt, _hfilt;
//-----------------------------------------------------------------------------
// diagnostics histograms
//-----------------------------------------------------------------------------
    Hist_t     _hist;

    THackData* fHackData;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    enum fitType {helixFit=0,seedFit,kalFit};
    explicit CalPatRec(const fhicl::ParameterSet& PSet);
    virtual ~CalPatRec();
    
    virtual void beginJob();
    virtual void beginRun(art::Run&);
    virtual void produce (art::Event& event ); 
    virtual void endJob();
//-----------------------------------------------------------------------------
// helper functions
//-----------------------------------------------------------------------------
    bool findData         (const art::Event& e);
    void findTimePeaks    (CalTimePeakCollection* TimePeakColl);
    void createTimePeak   (CalTimePeakCollection* TimePeakColl);
    void filterOutliers   (TrkDef& mytrk,Trajectory const& traj,double maxdoca,std::vector<TrkHitFilter>& thfvec);

//----------------------------------------------------------------------
// 2015 - 02 - 16 Gianipez added the two following functions
//----------------------------------------------------------------------
    void findDoublets     (KalRep* krep, DoubletCollection *dcol);//search doublets in a giventimepeak
    void findLoopApex     (){}//search the straw hits closer to the apexes of the helix loops

    void findMissingHits  (KalFitResult& kalfit, std::vector<hitIndex>& indices);
    void createDiagnostics();
    void fillStrawDiag    ();
    void fillTimeDiag     ();
    void fillFitDiag      (art::Event&               Event   ,
			   int                       ipeak   , 
			   HelixFitHackResult const& helixfit,
			   KalFitResult       const& seedfit ,
			   KalFitResult       const& kalfit  );
  };
}
#endif
