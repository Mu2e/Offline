//
// Object to perform helix fit to straw hits
//
// $Id: CalHelixFinderAlg.hh,v 1.8 2014/05/18 13:56:50 murat Exp $
// $Author: murat $
// $Date: 2014/05/18 13:56:50 $
//
#ifndef CalHelixFinderAlg_HH
#define CalHelixFinderAlg_HH

// data
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
// BaBar
#include "BTrk/TrkBase/TrkErrCode.hh"
//root
#include "TString.h"

#include "CalPatRec/inc/LsqSums4.hh"
#include "CalPatRec/inc/CalTimePeak.hh"
#include "CalPatRec/inc/CalHelixPoint.hh"
#include "CalPatRec/inc/CalHelixFinderData.hh"

class TH1F;

namespace fhicl {
  class ParameterSet;
}

namespace mu2e {
  class Calorimeter;
  class TTracker;
//-----------------------------------------------------------------------------
// output struct
//-----------------------------------------------------------------------------
  class CalHelixFinderAlg {
  public:
    enum { kMaxNHits = 10000 } ;
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    const TTracker*            _tracker;
    const Calorimeter*         _calorimeter;

    const CalTimePeak*         fTimePeak;
    const TimeCluster*         fTimeCluster; //needed for debugging
    
    double                     fCaloTime;
    double                     fCaloX;   
    double                     fCaloY;   
    double                     fCaloZ;   
    
    std::vector<CalHelixPoint> _xyzp;        // normally includes only hits from the time peak
//-----------------------------------------------------------------------------
// for diagnostics purposes save several states of _xyzp (only if _diag > 0)
//
// [0]: copy of initial _xyzp state
// [1]: after filterDist
// [2]: doPatternRecognition: after rescueHitsBeforeSeed
// [3]: doPatternRecognition: after rescueHits  (inside doPatternRecognition)
// [4]: doPatternRecognition: after final findZ (inside doPatternRecognition)
// [5]: doPatternRecognition: after doPatternRecognition
//
// each XYZPHack point has an 'outlier' flag telling if it belongs to the track
//-----------------------------------------------------------------------------
    struct SaveResults_t {
      std::vector<CalHelixPoint> _xyzp;    // points with used flags
      CalHelixFinderData         _helix;   // helix, found at this stage
    };

    SaveResults_t        _results[6];

    int                  fSeedIndex;
    int                  fCandidateIndex;
    int                  fLastIndex;
    int                  fUseDefaultDfDz;
    double               fHitChi2Max;

    int                  _diag;
    int                  _debug;
    int                  _debug2;
    StrawHitFlag         _hsel;         // good hit selection
    StrawHitFlag         _bkgsel;       // background hit selection
    int                  _minnhit;      // minimum # of hits to work with
    int                  _minnstereo;   // minimum # of stereo hits, I believe should  not be used
                                        // parameters for AGE center determination
    double               _minzsep;
    double               _maxzsep;      // Z separation of points for pitch estimate
    double               _maxdz;        // stereo selection parameters
    double               _maxdot;
                                        // 2014-03-10 Gianipez and P. Murat: limit
                                        // the dfdz value in the pattern-recognition stage
    double               _mpDfDz;
    double               _maxDfDz;
    double               _minDfDz;
    double               _sigmaPhi;     // hit phi resolution (wrt the trajectory axis, assume R=25-30cm)
    double               _weightXY;     //scale factor for makeing the xy-chi2 with a mean close to 1
    double               _weightZPhi;
    double               _maxXDPhi;     // max normalized hit residual in phi (findRZ)

    double               _hdfdz;        // estimated d(phi)/dz value
    double               _sdfdz;        // estimated d(phi)/dz error
    double               _hphi0;
					// 201-03-31 Gianipez added for changing the value of the
					// squared distance requed bewtween a straw hit and its predicted 
					// position used in the patter recognition procedure
    double               _distPatRec;

    double               _rbias;         // robust fit parameter bias
    double               _efac;          // error factor
    double               _rhomin;
    double               _rhomax;        // crude cuts on tranvservse radius for stereo hits
    double               _mindist;       // minimum distance between points used in circle initialization
    double               _maxdist;       // maximum distance in hits
    double               _pmin, _pmax;   // range of total momentum
    double               _tdmin, _tdmax; // range of abs(tan(dip)
    double               _rcmin,_rcmax;  // maximum transverse radius of circle
    double               _sfactor;       // stereo hit error factor
    bool                 _xyweights;
    bool                 _zweights;      // weight points by estimated errors
    bool                 _filter;        // filter hits
    bool                 _plotall;       // plot also failed fits
    bool                 _usetarget;     // constrain to target when initializing
    bool                 _allstereo;     // Use all stereo combinations (true) or only use each hit once
    mutable double       _bz;            // cached value of Field Z component at the tracker origin
//-----------------------------------------------------------------------------
// cached value of radius and pitch sign: these depend on the particle type
// and direction
//-----------------------------------------------------------------------------
    double    _rmin, _rmax, _smin, _smax, _dfdzsign;
//-----------------------------------------------------------------------------//
// store the paramters value of the most reliable track candidate
//-----------------------------------------------------------------------------//
    double    _x0, _y0, _phi0, _radius, _dfdz;
    int       _goodPointsTrkCandidate;
    int       _minPointsTrkCandidate;
    double    _chi2TrkCandidate;
    double    _maxChi2TrkCandidate;
    int       _markCandidateHits;
                                        // thresholds for XY and ZPhi chi2 fits
    double    _chi2xyMax;
    double    _chi2zphiMax;

    // indices, distance from prediction and distance along z axis from the seeding hit
    // of the hits found in the pattern recognition

    int       _indicesTrkCandidate[kMaxNHits];
    double    _distTrkCandidate   [kMaxNHits];
    double    _dzTrkCandidate     [kMaxNHits];

    double    _phiCorrected       [kMaxNHits];
    int       _phiCorrectedDefined;


    double    _dfdzErr;                 // error on dfdz by ::findDfDz

    TH1F*     _hDist;
    double    _chi2nFindZ;
    double    _eventToLook;
    TH1F*     _hDfDzRes;
    TH1F*     _hPhi0Res;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
                                        // parameter set should be passed in on construction

    explicit CalHelixFinderAlg(fhicl::ParameterSet const&);
    virtual ~CalHelixFinderAlg();
                                        // cached bfield accessor
    double bz() const;

    void   calculateDfDz(double &phi0, double &phi1, double &z0,  double &z1, double &dfdz);

    //projects the straw hit error along the radial direction of the circle-helix
    double calculateWeight     (CLHEP::Hep3Vector& HitPos, 
				CLHEP::Hep3Vector& StrawDir, 
				CLHEP::Hep3Vector& HelCenter, 
				double             Radius,
                                int                Print, 
				const char*        Banner);

    double calculatePhiWeight  (CLHEP::Hep3Vector HitPos, CLHEP::Hep3Vector StrawDir, CLHEP::Hep3Vector HelCenter, 
				double Radius, int Print, TString Banner);

    //calculates the residual along the radial direction of the helix-circle
    double calculateRadialDist (const CLHEP::Hep3Vector& HitPos, 
				const CLHEP::Hep3Vector& HelCenter, 
				double                   Radius);

    void   calculateTrackParameters(CLHEP::Hep3Vector& p0, double&radius,
                                    double& phi0, double& tanLambda,
                                    CLHEP::Hep3Vector p1, CLHEP::Hep3Vector p2,
                                    CLHEP::Hep3Vector p3, CalHelixFinderData& mytrk,
                                    bool cleanPattern=false);

    static double deltaPhi   (double phi1, double phi2);

   // returns the index of the hit which provides the highest contribute to the chi2
    void   doCleanUpWeightedCircleFit(::LsqSums4     &TrkSxy,
                                       int            SeedIndex,
                                       int            *IdVec,
                                       CLHEP::Hep3Vector     &HelCenter,
                                       double         &Radius,
                                       double         *Weights,
                                       int            &Iworst);

    bool   doLinearFitPhiZ     (CalHelixFinderData& Helix, int SeedIndex, int *indexVec,
				int UseInteligentWeight=0);

   //perfoms the weighted circle fit, update the helix parameters (HelicCenter, Radius) and
    // fills the vector Weights which holds the calculated weights of the hits
    void   doWeightedCircleFit (::LsqSums4 &TrkSxy, int SeedIndex,int *IdVec,
                                CLHEP::Hep3Vector &HelCenter, double &Radius, double *Weights,
                                int Print=0, TString Banner="");

    void doPatternRecognition(CalHelixFinderData& mytrk);
  
    void defineHelixParams(CalHelixFinderData& Helix) const;

    // TH1F* hDist() {return _hDist;}

    int   isHitUsed(int index);

    void fillXYZP                     (CalHelixFinderData& Helix);
    void filterDist                   ();
    void filterUsingPatternRecognition(CalHelixFinderData& Helix);
    bool findHelix                    (CalHelixFinderData& Helix);
    bool findHelix                    (CalHelixFinderData& Helix, const CalTimePeak* TimePeak);
    bool findHelix                    (CalHelixFinderData& Helix, const TimeCluster* TimePeak );
    int  findDfDz                     (CalHelixFinderData& Helix, int SeedIndex, int *indexVec);
    void findTrack                    (int                  seedIndex,
			               double&              chi2,
			               int&                 countGoodPoint,
			               CalHelixFinderData&  mytrk,
			               int&                 mode,
			               bool                 useDefaultDfDz=false,
			               int                  useMPVdfdz = 0);
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
    void  setTracker    (const TTracker*    Tracker) { _tracker     = Tracker; }
    void  setCalorimeter(const Calorimeter* Cal    ) { _calorimeter = Cal    ; }
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
    void   plotXY               (int ISet);
    void   plotZPhi             (int ISet);
    void   printInfo            (CalHelixFinderData& myhel);
    void   printXYZP            (const char* Title);

    int    refineHelixParameters(CalHelixFinderData& Trk,
				 int seedIndex, int *indexVec,
				 int Print=0, TString Banner="");

                                        // 12-10-2013 Gianipez: new pattern recognition functions
    void   rescueHitsBeforeSeed (CalHelixFinderData&  mytrk);

    void   rescueHits           (CalHelixFinderData&  mytrk, int seedIndex       ,
				 int *indexVec             , int UsePhiResiduals = 0);


    void   resetTrackParamters  ();
//-----------------------------------------------------------------------------
// save intermediate results in diagnostics mode
//-----------------------------------------------------------------------------
    void   saveResults                    (std::vector<CalHelixPoint>& Xyzp, 
					   CalHelixFinderData&         Helix, 
					   int                         Index);

    void   searchWorstHitWeightedCircleFit(int               SeedIndex,
                                           int               *IdVec,
                                           CLHEP::Hep3Vector &HelCenter,
                                           double            &Radius,
                                           double            *Weights,
                                           int               &Iworst ,
                                           double            &HitChi2Worst);

  };
}
#endif
