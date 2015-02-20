//
// Object to perform helix fit to straw hits
//
// $Id: HelixFitHack.hh,v 1.8 2014/05/18 13:56:50 murat Exp $
// $Author: murat $ 
// $Date: 2014/05/18 13:56:50 $
//
#ifndef HelixFitHack_HH
#define HelixFitHack_HH

// framework
#include "fhiclcpp/ParameterSet.h"
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

#include "LsqSums4.hh"

#include "HelixDefHack.hh"
#include "CalPatRec/inc/CalTimePeak.hh"

class TH1F;
class TFile;

namespace mu2e {
//-----------------------------------------------------------------------------
// output struct
//-----------------------------------------------------------------------------
  struct HelixFitHackResult {
    HelixDefHack       _hdef;         // must copy by value as references can't be re-assigned
    TrkErrCode         _fit;	      // fit status code from last fit
//-----------------------------------------------------------------------------
// circle parameters; the z center is ignored.
//-----------------------------------------------------------------------------
    ::LsqSums4         _sxy;
    ::LsqSums4         _srphi;
    double             _chi2;

    CLHEP::Hep3Vector  _center;
    double             _radius;
//-----------------------------------------------------------------------------
// 2015-02-06 P.Murat: fit with non-equal weights - XY -only
//-----------------------------------------------------------------------------
    ::LsqSums4         _sxyw;
    CLHEP::Hep3Vector  _cw;
    double             _rw;
    double             _chi2w;
//-----------------------------------------------------------------------------
// Z parameters; dfdz is the slope of phi vs z (=-sign(1.0,qBzdir)/(R*tandip)), 
// fz0 is the phi value of the particle where it goes through z=0
// note that dfdz has a physical ambiguity in q*zdir.
//-----------------------------------------------------------------------------
    double             _dfdz;
    double             _fz0;

    HelixFitHackResult(TrkDef const& tdef) : 
      _hdef(tdef),  
      _fit(TrkErrCode::fail),
      _radius(-1.0),
      _dfdz(0.0),
      _fz0(0.0) {
    }

    HelixFitHackResult(HelixDefHack const& hdef) : 
      _hdef(hdef),  
      _fit(TrkErrCode::fail),
      _radius(-1.0),
      _dfdz(0.0),
      _fz0(0.0) {
    }

    HelixFitHackResult(HelixFitHackResult & hfitRes) :
      _hdef(hfitRes._hdef),
      _fit(hfitRes._fit),
      _sxy(hfitRes._sxy),
      _srphi(hfitRes._srphi),
      _chi2(hfitRes._chi2),
      _center(hfitRes._center),
      _radius(hfitRes._radius),
      _dfdz(hfitRes._dfdz),
      _fz0(hfitRes._fz0){
    
    }

    HelixFitHackResult& operator =(HelixFitHackResult const& other);
 };

//-----------------------------------------------------------------------------
// utility struct; value plus error
//-----------------------------------------------------------------------------
  struct VALERR {
    double _val;
    double _err;
  };

  struct FZ {
    VALERR _phi;
    double _z;
  };
//-----------------------------------------------------------------------------
// utility struct
//-----------------------------------------------------------------------------
  struct XYZPHack {
    size_t             _ind;		// straw hit index
    CLHEP::Hep3Vector  _pos;		// position
    double             _phi;	        // ambiguity-resolved phi angle
    StrawHitFlag       _flag;		// flag
    int                _used;           // = 1 if the strawhit is used by another track-candidate
    CLHEP::Hep3Vector  _wdir;		// wire direction
    CLHEP::Hep3Vector  _sdir;           // straw radial direction, perp to Z and wire direction
					// errors are asymmetric; along the wire is given by time division, 
					// perp to the wire by the straw size/sqrt(12)
    const mu2e::Straw* _straw;          // pointer to the straw
    const mu2e::StrawHit* _strawhit;          // pointer to the strawHit
    double             _perr;
    double             _rerr;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
					// initialize some variables on construction
    XYZPHack():_phi(0.0) {
    }

    XYZPHack(size_t ind, StrawHit const& sh, StrawHitPosition const& shp, Straw const& straw, StrawHitFlag const& flag);

    XYZPHack(size_t ind, CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& wdir, double werr, double serr);
 
// radial position information
    virtual void rinfo     (CLHEP::Hep3Vector const& center, VALERR& rad) const;
    virtual void finfo     (CLHEP::Hep3Vector const& center, VALERR& phi) const;
    bool         use       () const;
    int          isUsed    () { return _used;}
    bool         stereo    () const;
    bool         isOutlier () const;
    bool         isCalosel () const;
    void         setOutlier();
    void         setUse    (bool use);

    static double       _efac;
    static StrawHitFlag _useflag;	// flag bits to define use
  };

  typedef std::vector<XYZPHack> XYZPHackVector;

//-----------------------------------------------------------------------------
// struct to hold AGE sums, contains weighted (s)ums of (c)osine and (s)in for 
// points on (c)ircumference, (o)utside the median radius, or (i)nside the median radius
//-----------------------------------------------------------------------------
  struct SUMS {
    double   _scc, _ssc, _sco, _sso, _sci, _ssi;
    unsigned _nc, _no, _ni;
    SUMS() : _scc(0.0),_ssc(0.0),_sco(0.0),_sso(0.0),_sci(0.0),_ssi(0.0){}
    void clear() { _scc = _ssc = _sco = _sso = _sci = _ssi = 0.0;
      _nc = _no = _ni = 0; }
  };


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  class HelixFitHack {
//-----------------------------------------------------------------------------
// configuration parameters
//-----------------------------------------------------------------------------
    const CalTimePeak*   fTimePeak;
    
    int  fSeedIndex;   
    int  fCandidateIndex;
    int  fLastIndex;   
    int  fUseDefaultDfDz;

    int                 _diag;
    int                 _debug;
    double    _mindelta; // minimum slope difference to use a triple in circle center initialization
    unsigned  _minnhit; // minimum # of hits to work with
    unsigned  _minnstereo; // minimum # of stereo hits
    double    _lambda0,_lstep,_minlambda; // parameters for AGE center determination
    unsigned  _maxniter; // maxium # of iterations to global minimum
    double    _nsigma; // # of sigma for filtering outlyers
    double    _nssigma; // # of sigma for filtering stereo time division
    double    _minzsep, _maxzsep; // Z separation of points for pitch estimate
    double    _maxdz;			// stereo selection parameters
    double    _maxdot; 

					//2014-03-10 Gianipez and P. Murat introduced the following paramter to limit 
					// the dfdz value in the pattern-recognition stage
    double   _mpDfDz, _maxDfDz, _minDfDz;
    
    double   _hdfdz;		        // estimated d(phi)/dz value
    double   _sdfdz;		        // estimated d(phi)/dz error
    double   _hphi0;

    //201-03-31 Gianipez added th following parameter for changing the value of the 
    // squared distance requed bewtween the strawhits and the theretical position
    // used in the patter recognition procedure
    double _distPatRec;

    double _rbias;  // robust fit parameter bias
    double _efac; // error factor
    double _rhomin, _rhomax; // crude cuts on tranvservse radius for stereo hits
    double _mindist; // minimum distance between points used in circle initialization
    double _maxdist; // maximum distance in hits
    double _pmin, _pmax; // range of total momentum
    double _tdmin, _tdmax; // range of abs(tan(dip)
    double _rcmin,_rcmax; // maximum transverse radius of circle
    double _sfactor; // stereo hit error factor
    bool   _forcep; // force the p/pt to be in range (true), or exclude fits outside that range (false)
    bool   _xyweights,_zweights; // weight points by estimated errors 
    bool   _filter; // filter hits
    bool   _plotall; // plot also failed fits
    bool   _usetarget; // constrain to target when initializing
    bool   _allstereo; // Use all stereo combinations (true) or only use each hit once
    mutable double _bz; // cached value of Field Z component at the tracker origin
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

    //thresholds for the xy and zphi chi2 fit
    double    _chi2xyMax;
    double    _chi2zphiMax;

    //indices, distance from prediction and distance along z axis from the seeding hit
    // of the hits found in the pattern recognition
    int       _indicesTrkCandidate[400];
    double    _distTrkCandidate[400];
    double    _dzTrkCandidate[400];

    //error on dfdz resulting from the proceedure ::findDfDz
    double    _dfdzErr;

    TH1F*     _hDist;
    double    _chi2nFindZ;
    double    _eventToLook;
    TH1F*     _hDfDzRes;

  public:
					// parameter set should be passed in on construction

    explicit HelixFitHack(fhicl::ParameterSet const&);
    virtual ~HelixFitHack();
					// main function: given a track definition, 
					// find the helix parameters

    bool findHelix(HelixFitHackResult& myfit, const CalTimePeak* TimePeak);

					// convert to BaBar helix parameters.  
					// Also return an error estimate

    void helixParams (HelixFitHackResult const& helix, 
		      CLHEP::HepVector&         pvec , 
		      CLHEP::HepVector&         perr ) const;
    
    
    TH1F* hDist() {return _hDist;}

    int   isHitUsed(int index);
    
    void printInfo(HelixFitHackResult& myhel);

    XYZPHackVector  _xyzp;
    
  protected:
//-----------------------------------------------------------------------------
// utlity functions
//-----------------------------------------------------------------------------
					// allow passing in the struct by hand

    bool findHelix     (XYZPHackVector& xyzp, HelixFitHackResult& myfit);
    bool findXY        (XYZPHackVector& xyzp, HelixFitHackResult& myhel);
    bool findXY_new    (XYZPHackVector& xyzp, HelixFitHackResult& myhel,
			int seedIndex, int *indexVec);
    bool findZ         (XYZPHackVector& xyzp, HelixFitHackResult& myhel, 
			int seedIndex, int *indexVec);
    void findDfDz      (XYZPHackVector& xyzp, HelixFitHackResult& myhel, 
			int seedIndex, int *indexVec);
    

    bool initCircle    (XYZPHackVector const& xyzp,HelixFitHackResult& myhel);
    bool initCircle_new(XYZPHackVector const& xyzp,HelixFitHackResult& myhel);

  private:

    //    void fillXYZP(HelixDefHack const& mytrk, XYZPHackVector& xyzp);
    void fillXYZP(HelixDefHack const& mytrk);

// cached bfield accessor

    double bz() const;

// utility function to resolve phi wrapping    

    static double deltaPhi(double phi1, double phi2);

// diagnostics

    void plotXY(HelixDefHack const& mytrk, XYZPHackVector const& xyzp, HelixFitHackResult const& myhel) const;

    void plotZ (HelixDefHack const& mytrk, XYZPHackVector const& xyzp, HelixFitHackResult const& myhel) const;

// find the Absolute Geometric Error.  Returns the median radius as well.

    bool findCenterAGE(XYZPHackVector const& xyzp, Hep3Vector& center , double& rmed, double& age, bool useweights=false);

    void findAGE (XYZPHackVector const& xyzp, Hep3Vector const& center, double& rmed, double& age, bool useweights=false);

    void fillSums(XYZPHackVector const& xyzp, Hep3Vector const& center, double rmed, SUMS&   sums, bool useweights=false);

    void filterXY(XYZPHackVector& xyzp      , Hep3Vector const& center, double rmed, bool& changed);

    void filterDist(XYZPHackVector& xyzp);
					// 12-10-2013 Gianipez: new pattern recognition functions
    void rescueHits(XYZPHackVector& xyzp, HelixFitHackResult&  mytrk);

    void filterUsingPatternRecognition(XYZPHackVector& xyzp, HelixFitHackResult&  mytrk);
    
    void resetTrackParamters();

    void doPatternRecognition(XYZPHackVector& xyzp, HelixFitHackResult& mytrk);

    void findTrack(XYZPHackVector&      xyzp, 
		   int                  seedIndex,
		   double&              chi2,
		   int&                 countGoodPoint,
		   HelixFitHackResult&  mytrk, 
		   int&                 mode,
		   bool                 useDefaultDfDz=false,
		   int                  useMPVdfdz = 0);

    void calculateTrackParameters(Hep3Vector& p0, double&radius,
				  double& phi0, double& tanLambda,
				  Hep3Vector p1, Hep3Vector p2,
				  Hep3Vector p3, HelixFitHackResult& mytrk,
				  bool cleanPattern=false);

    void calculateDfDz(double &phi0, double &phi1, double &z0,  double &z1, double &dfdz);

    int  refineHelixParameters(XYZPHackVector& Xyzp, HelixFitHackResult& Trk,
			       int seedIndex, int *indexVec);
  };
}
#endif
