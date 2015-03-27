//
// Object to perform helix fit to straw hits
//
// $Id: HelixFitHack.hh,v 1.8 2014/05/18 13:56:50 murat Exp $
// $Author: murat $ 
// $Date: 2014/05/18 13:56:50 $
//
#ifndef HelixFitHack_HH
#define HelixFitHack_HH

// data
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"

// HelixFitHack objects
#include "KalmanTests/inc/TrkDef.hh"
// BaBar
#include "TrkBase/TrkErrCode.hh"
//CLHEP
//root
#include "TString.h"

#include "LsqSums4.hh"

#include "CalPatRec/inc/HelixDefHack.hh"
#include "CalPatRec/inc/CalTimePeak.hh"
#include "CalPatRec/inc/XYZPHack.hh"
#include "CalPatRec/inc/HelixFitHackResult.hh"

namespace fhicl {
  class ParameterSet;
}

class TH1F;
class TFile;

namespace mu2e {
//-----------------------------------------------------------------------------
// output struct
//-----------------------------------------------------------------------------
  class HelixFitHack {
  public:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    const CalTimePeak*   fTimePeak;
    XYZPHackVector       _xyzp;        // normally includes only hits from the time peak
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
      XYZPHackVector     _xyzp;    // points with used flags
      HelixFitHackResult _helix;   // helix, found at this stage
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
    double               _mindelta;     // minimum slope difference to use a triple in circle center initialization
    int                  _minnhit;      // minimum # of hits to work with
    int                  _minnstereo;   // minimum # of stereo hits, I believe should  not be used
					// parameters for AGE center determination
    double               _lambda0;
    double               _lstep;
    double               _minlambda; 
    int                  _maxniter;     // maxium # of iterations to global minimum
    double               _nsigma;       // # of sigma for filtering outlyers
    double               _nssigma;      // # of sigma for filtering stereo time division
    double               _minzsep;
    double               _maxzsep;      // Z separation of points for pitch estimate
    double               _maxdz;	// stereo selection parameters
    double               _maxdot; 
					// 2014-03-10 Gianipez and P. Murat: limit 
					// the dfdz value in the pattern-recognition stage
    double               _mpDfDz;
    double               _maxDfDz;
    double               _minDfDz;
    double               _sigmaPhi;     // hit phi resolution (wrt the trajectory axis, assume R=25-30cm)
    double               _maxXDPhi;     // max normalized hit residual in phi (findRZ)
    
    double               _hdfdz;	// estimated d(phi)/dz value
    double               _sdfdz;	// estimated d(phi)/dz error
    double               _hphi0;

    //201-03-31 Gianipez added th following parameter for changing the value of the 
    // squared distance requed bewtween the strawhits and the theretical position
    // used in the patter recognition procedure
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
    bool                 _forcep;    // force p/pt to be in range (true), or exclude fits outside that range (false)
    bool                 _xyweights;
    bool                 _zweights;  // weight points by estimated errors 
    bool                 _filter;    // filter hits
    bool                 _plotall;   // plot also failed fits
    bool                 _usetarget; // constrain to target when initializing
    bool                 _allstereo; // Use all stereo combinations (true) or only use each hit once
    mutable double       _bz;        // cached value of Field Z component at the tracker origin
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
    int       _indicesTrkCandidate[400];
    double    _distTrkCandidate   [400];
    double    _dzTrkCandidate     [400];

    
    double    _dfdzErr;                 // error on dfdz by ::findDfDz

    TH1F*     _hDist;
    double    _chi2nFindZ;
    double    _eventToLook;
    TH1F*     _hDfDzRes;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
					// parameter set should be passed in on construction

    explicit HelixFitHack(fhicl::ParameterSet const&);
    virtual ~HelixFitHack();
					// cached bfield accessor
    double bz() const;

    void   calculateDfDz(double &phi0, double &phi1, double &z0,  double &z1, double &dfdz);

    //projects the straw hit error along the radial direction of the circle-helix
    double calculateWeight     (Hep3Vector HitPos, Hep3Vector StrawDir, Hep3Vector HelCenter, double Radius,
				int Print, TString Banner);

    //calculates the residual along the radial direction of the helix-circle
    double calculateRadialDist (Hep3Vector HitPos, Hep3Vector HelCenter, double Radius);

    void   calculateTrackParameters(Hep3Vector& p0, double&radius,
				    double& phi0, double& tanLambda,
				    Hep3Vector p1, Hep3Vector p2,
				    Hep3Vector p3, HelixFitHackResult& mytrk,
				    bool cleanPattern=false);

    // returns the index of the hit which provides the highest contribute to the chi2 
    void    doCleanUpWeightedCircleFit(::LsqSums4     &TrkSxy, 
				       int            SeedIndex, 
				       int            *IdVec,
				       Hep3Vector     &HelCenter, 
				       double         &Radius, 
				       double         *Weights,
				       int            &Iworst);

    void   doPatternRecognition(XYZPHackVector& xyzp, HelixFitHackResult& mytrk);

    //perfoms the weighted circle fit, update the helix parameters (HelicCenter, Radius) and 
    // fills the vector Weights which holds the calculated weights of the hits
    void   doWeightedCircleFit (::LsqSums4 &TrkSxy, int SeedIndex,int *IdVec,
				Hep3Vector &HelCenter, double &Radius, double *Weights,
				int Print=0, TString Banner="");
        
					// convert to BaBar helix parameters.  
					// Also return an error estimate

    void  helixParams (HelixFitHackResult const& helix, 
		       CLHEP::HepVector&         pvec , 
		       CLHEP::HepVector&         perr ) const;
        
    TH1F* hDist() {return _hDist;}

    int   isHitUsed(int index);
    
    //    void printInfo(HelixFitHackResult& myhel);

    XYZPHackVector  fxyzp;
    double          fPhiCorrected[1000];
    
					// allow passing in the struct by hand

    bool findHelix      (HelixFitHackResult& Helix);
    bool findHelix      (HelixFitHackResult& Helix, const CalTimePeak* TimePeak );
    bool doLinearFitPhiZ(HelixFitHackResult& Helix, int SeedIndex, int *indexVec);
    int  findDfDz       (HelixFitHackResult& Helix, int SeedIndex, int *indexVec);
 
    void findTrack(int                  seedIndex,
		   double&              chi2,
		   int&                 countGoodPoint,
		   HelixFitHackResult&  mytrk, 
		   int&                 mode,
		   bool                 useDefaultDfDz=false,
		   int                  useMPVdfdz = 0);

    void fillXYZP(HelixDefHack const& mytrk);

    static double deltaPhi(double phi1, double phi2);

    void          doPatternRecognition(HelixFitHackResult& mytrk);

// find the Absolute Geometric Error.  Returns the median radius as well.

// utility function to resolve phi wrapping    

//    static double deltaPhi(double phi1, double phi2);

// find the Absolute Geometric Error.  Returns the median radius as well.
    void filterDist();
    void filterUsingPatternRecognition(HelixFitHackResult&  mytrk);
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
    void   plotXY  (int ISet);
    void   plotZPhi(int ISet);
    void   printInfo(HelixFitHackResult& myhel);

    int  refineHelixParameters(HelixFitHackResult& Trk,
			       int seedIndex, int *indexVec, 
			       int Print=0, TString Banner="");
    
					// 12-10-2013 Gianipez: new pattern recognition functions
    void rescueHitsBeforeSeed(HelixFitHackResult&  mytrk);

    void rescueHits(HelixFitHackResult&  mytrk, int seedIndex       , 
		    int *indexVec             , int UsePhiResiduals = 0);

    
    void resetTrackParamters();
//-----------------------------------------------------------------------------
// save intermediate results in diagnostics mode
//-----------------------------------------------------------------------------
    void   saveResults(XYZPHackVector& Xyzp, HelixFitHackResult&  Helix, int Index); 

    void   searchWorstHitWeightedCircleFit(int             SeedIndex,
					   int            *IdVec,
					   Hep3Vector     &HelCenter, 
					   double         &Radius, 
					   double         *Weights,
					   int            &Iworst , 
					   double         &HitChi2Worst);
    
  };
}
#endif
