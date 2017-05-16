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
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"

// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/BbrStringUtils.hh"
#include "CalPatRec/inc/TrkDefHack.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "CalPatRec/inc/CalTimePeak.hh"
#include "RecoDataProducts/inc/Doublet.hh"

#include "TROOT.h"
#include "TFolder.h"
#include "CalPatRec/inc/KalFitHack.hh"
#include "CalPatRec/inc/HelixFitHack.hh"
#include "CalPatRec/inc/THackData.hh"

// Mu2e
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

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
    struct HelixFitHist_t {
      TH1F*  nhits;           // number of hits on a helix  
    };

    struct Hist_t {
      HelixFitHist_t  helixFit;  // helix fit histograms

      TH1F* _cutflow[2];       // diagnostic flow
      TH1F* _dt     [2];       // distribution in Dt = t_calo-cluster - t_straw-hit
      TH1F* _Tpeaks;	       // number of time-peaks found per event
      TH1F* _NfitIter;	       // number of call to kalman addHits
	     		       //     TH1F* _hTfit[2];     //time spent on kalman filter per event
	     		       //     TH1F* _hTtot;     //total time spent per event
      TH1F* _radius;           // radius of the helix used by findTrack
      TH1F* _phi0;             // phi0 of the helix used by findTrack
      TH1F* _tanlambda;        // tanLambda (tangent of the pitch angle) of the helix
      TH1F* _dfdz;             // dfdz of the theretical helix. dfdz = tanLambda/radius
      TH2F* _distvsdz;
      TH1F* _dphidz[4];
      TH1F* _dphi0 [4];
      TH1F* _dr[2];

      TH1F* _kradius[2];    // radius of the helix used by findTrack

      //      TH1F* _kdist;
      TH1F* _kdz;
      TH1F* _kNpoints;
      TH1F* _kchi2;
      TH2F* _kdistvsdz[2];
      TH1F* _kdphidz[4];
      TH1F* _kdphi0[4];
      TH1F* _kdr[2];

      TH1F* _drw  [2];
      TH1F* _chi2w[2];
      
      TH1F* _chi2zphi[2];
      
      TH1F* _seeddoca[3];
      TH1F* _seeddr  [2];
      TH1F* _seeddfdz[2];
      TH1F* _NpointsSeed   [2]; //
      TH1F* _doca          [4];
      TH1F* _kaldoca       [2];
      TH1F* _NpointsRescued[2];
      TH1F* _PhiResid      [3];
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
      TH1F*  _nhits[2];         // same distribution with different limits
      TH1F   *ctsp, *rtsp, *ttsp, *ltsp, *tdtsp;
    };

    Ref*    _ref;

  protected:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    //    TStopwatch*   fStopwatch;

    unsigned     _iev;
					// configuration parameters
    int          _diagLevel; 
    int          _debugLevel;
    int          _printfreq;
    bool         _addhits; 
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string      _shLabel ; // MakeStrawHit label (makeSH)
    std::string      _shDigiLabel ;
    std::string      _shpLabel;
    std::string      _shfLabel;
    std::string      _ccmLabel; // caloClusterModuleLabel
    std::string      _crmLabel;
    std::string      _chmccpLabel;

    //    std::string      _dtspecpar;

    StrawHitFlag     _tsel, _hsel, _addsel, _ksel;
    StrawHitFlag     _bkgsel, _addbkg;
    double           _maxedep;
    //    int              _useDoublets;
    double           _mindt;
    double           _maxdt;
    double           _maxdtmiss;
    int              _final;         // 1: make final ambig resolution decision
    double           _fbf;
					// time spectrum parameters
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

    int              _nhits_from_gen;
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;
    const CaloClusterCollection*          _ccCollection;
    const StepPointMCCollection*          _stepcol;
    const PtrStepPointMCVectorCollection* _listOfMCStrawHits;

    const CaloHitCollection*              _chcol;
    const CaloHitMCTruthCollection*       _chmccol;
    const PtrStepPointMCVectorCollection* _listOfMCCrystals;
    const CaloHitSimPartMCCollection*     _chsmccol;
    CaloHitMCNavigator*                   _caloHitNavigator;

    StrawHitFlagCollection*               _flags;    // flags are not const - for a reason ?
					
    HelixFitHack                          _hfit;	// robust helix fitter
    KalFitHack                            _seedfit;  // Kalman filter config for the Seed fit ( fit using hit wires)
    KalFitHack                            _kfit;     // full-blown src/Kalman filter

    KalFitResult*                         _sfresult; // seed fit result
    KalFitResult*                         _kfresult; // full fit result

    CalTimePeakCollection*                _tpeaks;   // cache of time peaks

    XYZPHackVector                        _index;
    int                                   _nindex;
    int                                   _nrescued;    // by the seed fit

    const TTracker*                       _tracker;     // straw tracker geometry
    const Calorimeter*                    _calorimeter; // cached pointer to the calorimeter geometry

    const TrackerCalibrations*            _trackerCalib;

    TFolder*                              _folder;
    int                                   _eventid;
    int                                   _ntracks;
//-----------------------------------------------------------------------------
// diagnostics histograms
//-----------------------------------------------------------------------------
    Hist_t                                _hist;

    THackData*                            fHackData;

    int                                   _minNMCHits;
    int                                   fQualityTrack;
    int                                   fCaloClusterFromCE;
    int                                   fNCaloEnergyCut, fNCaloSizeCut, fNHitsTimePeakCut, fNTimeWindow;
    int                                   fSHSel[3], fSHBkg[3];
    bool                                  _clCE;

    double                                _mbtime;               // period of 1 microbunch
    SimParticleTimeOffset*                fgTimeOffsets;

    HelixTraj*                            _helTraj;
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
    void filterOutliers   (TrkDefHack& mytrk,Trajectory const& traj,double maxdoca,std::vector<TrkHitFilter>& thfvec);
//----------------------------------------------------------------------
// 2015 - 02 - 16 Gianipez added the two following functions
//----------------------------------------------------------------------
    void findDoublets     (KalRep* krep, DoubletCollection *dcol);//search doublets in a giventimepeak
    void findLoopApex     (){}//search the straw hits src/closer to the apexes of the helix loops

    void findMissingHits  (KalFitResult& kalfit, std::vector<StrawHitIndex>& indices);
    void bookHistograms   ();
    void fillStrawDiag    ();
    void fillTimeDiag     ();
    void fillFitDiag      (art::Event&       Event   ,
			   int               ipeak   , 
			   HelixFitHackResult&  helixfit,
			   KalFitResult      & seedfit ,
			   KalFitResult      & kalfit  );

    void fillSeedFitHistograms(KalFitResult& SFResult);

    void init             (KalFitResult*&  KRes, TrkDefHack* TDef);

  };
}
#endif

