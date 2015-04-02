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

  class Calorimeter;
  class TTracker;

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
      TH1F* _hNpointsSeed   [2]; //
      TH1F* _hdoca          [4];
      TH1F* _hkaldoca       [2];
      TH1F* _hNpointsRescued[2];
      TH1F* _hPhiResid      [2];
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

    std::string             _iname_seed;// data instance name for the output 
					 // of seed fit (used for diagnostics only)

    XYZPHackVector           _index;
    int                      _nindex;

    const TTracker*          _tracker;     // straw tracker geometry
    const Calorimeter*       _calorimeter; // cached pointer to the calorimeter geometry

    TFolder*                 _folder;
// //-----------------------------------------------------------------------------
// // strawhit tuple variables
// //-----------------------------------------------------------------------------
    int    _eventid;
    int    _ntracks;
    int    _ipeak;
    int    _helixfail,_seedfail,_kalfail;
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
