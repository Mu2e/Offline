//
// Object to perform helix fit to straw hits
//
// $Id: HelixFitHack.hh,v 1.6 2014/04/08 04:25:46 murat Exp $
// $Author: murat $ 
// $Date: 2014/04/08 04:25:46 $
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
				      // circle parameters; the z center is ignored.
    ::LsqSums4         _sxy;
    ::LsqSums4         _srphi;
    double             _chi2;

    CLHEP::Hep3Vector  _center;
    double             _radius;
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
  struct XYZP {
    size_t            _ind;		// straw hit index
    CLHEP::Hep3Vector _pos;		// position
    double            _phi;	        // ambiguity-resolved phi angle
    StrawHitFlag      _flag;		// flag
    CLHEP::Hep3Vector _wdir;		// wire direction
    CLHEP::Hep3Vector _sdir;            // straw radial direction, perp to Z and wire direction
					// errors are asymmetric; along the wire is given by time division, 
					// perp to the wire by the straw size/sqrt(12)
    double _perr,_rerr;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
					// initialize some variables on construction
    XYZP():_phi(0.0) {
    }

    XYZP(size_t index, StrawHit const& sh, StrawHitPosition const& shp, Straw const& straw);

    XYZP(size_t ind, CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& wdir, double werr, double serr) :
      _ind(ind),
      _pos(pos),
      _phi(_pos.phi()),
      _wdir(wdir),
      _sdir(wdir.y(),-wdir.x(),0.0),
      _perr(_efac*werr),
      _rerr(_efac*serr) {}
 
// radial position information
    virtual void rinfo     (CLHEP::Hep3Vector const& center, VALERR& rad) const;
    virtual void finfo     (CLHEP::Hep3Vector const& center, VALERR& phi) const;
    bool         use       () const;
    bool         stereo    () const;
    bool         isOutlier ();
    void         setOutlier();
    void         setUse    (bool use);

    static double       _efac;
    static StrawHitFlag _useflag;	// flag bits to define use
  };

  typedef std::vector<XYZP> XYZPVector;

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
    
    int  fSeedIndex;   // hack->fData[2]
    int  fLastIndex;   // hack->fData[3]

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
    double   _maxDfDz;

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
    double _rmin, _rmax, _smin, _smax, _dfdzsign;
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

    TH1F*     _hDist;
    double    _chi2nFindZ;
    double    _eventToLook;

  protected:
//-----------------------------------------------------------------------------
// utlity functions
//-----------------------------------------------------------------------------
					// allow passing in the struct by hand

    bool findHelix     (XYZPVector& xyzp, HelixFitHackResult& myfit);
    bool findXY        (XYZPVector& xyzp, HelixFitHackResult& myhel);
    bool findXY_new    (XYZPVector& xyzp, HelixFitHackResult& myhel);
    bool findZ         (XYZPVector& xyzp, HelixFitHackResult& myhel);

    bool initCircle    (XYZPVector const& xyzp,HelixFitHackResult& myhel);
    bool initCircle_new(XYZPVector const& xyzp,HelixFitHackResult& myhel);

  private:

    void fillXYZP(HelixDefHack const& mytrk, XYZPVector& xyzp);

// cached bfield accessor

    double bz() const;

// utility function to resolve phi wrapping    

    static double deltaPhi(double phi1, double phi2);

// diagnostics

    void plotXY(HelixDefHack const& mytrk, XYZPVector const& xyzp, HelixFitHackResult const& myhel) const;

    void plotZ (HelixDefHack const& mytrk, XYZPVector const& xyzp, HelixFitHackResult const& myhel) const;

// find the Absolute Geometric Error.  Returns the median radius as well.

    bool findCenterAGE(XYZPVector const& xyzp, Hep3Vector& center , double& rmed, double& age, bool useweights=false);

    void findAGE (XYZPVector const& xyzp, Hep3Vector const& center, double& rmed, double& age, bool useweights=false);

    void fillSums(XYZPVector const& xyzp, Hep3Vector const& center, double rmed, SUMS&   sums, bool useweights=false);

    void filterXY(XYZPVector& xyzp      , Hep3Vector const& center, double rmed, bool& changed);

    void filterDist(XYZPVector& xyzp);
					// 12-10-2013 Gianipez: new pattern recognition functions
    
    void doPatternRecognition(XYZPVector& xyzp, HelixFitHackResult& mytrk);

    void findTrack(XYZPVector&          xyzp, 
		   int                  seedIndex,
		   double&              chi2,
		   int&                 countGoodPoint,
		   HelixFitHackResult&  mytrk, 
		   bool                 cleanPatetrn=false);

    void calculateTrackParameters(Hep3Vector& p0, double&radius,
				  double& phi0, double& tanLambda,
				  Hep3Vector p1, Hep3Vector p2,
				  Hep3Vector p3, HelixFitHackResult& mytrk,
				  bool cleanPattern=false);

    void calculateDfDz(double &phi0, double &phi1, double &z0,  double &z1, double &dfdz);
  };
}
#endif
