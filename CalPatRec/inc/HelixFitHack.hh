//
// Object to perform helix fit to straw hits
//
// $Id: HelixFitHack.hh,v 1.2 2014/03/04 20:42:02 gianipez Exp $
// $Author: gianipez $ 
// $Date: 2014/03/04 20:42:02 $
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
//#include "CLHEP/Matrix/Vector.h"
//root

#include "LsqSums4.hh"
class TH1F;
class TFile;
// C+

namespace mu2e 
{


// add the StrawHitPosition collection to TrkDef
  class HelixDef : public TrkDef {
    public:
      HelixDef(TrkDef const& tdef) : TrkDef(tdef), _shpos(0) {}
      HelixDef(const StrawHitCollection* strawcollection,
	       const StrawHitPositionCollection* shposcollection, 
	       const std::vector<hitIndex>& strawhits,
	       TrkParticle const& tpart=_eminus, 
	       TrkFitDirection const& fdir=_downstream) : TrkDef(strawcollection,strawhits,tpart,fdir), _shpos(shposcollection) {}
      HelixDef() {}
      const StrawHitPositionCollection* strawHitPositionCollection() const { return _shpos; }
    private:
      const StrawHitPositionCollection* _shpos;
  };
// output struct
  struct HelixFitHackResult {
    HelixDef _hdef; // must copy by value as references can't be re-assigned
// fit status
    TrkErrCode _fit; // error code from last fit
// circle parameters; the z center is ignored.

    ::LsqSums4          _sxy;
    ::LsqSums4          _srphi;
    double            _chi2;

    CLHEP::Hep3Vector _center;
    double            _radius;
// Z parameters; dfdz is the slope of phi vs z (=-sign(1.0,qBzdir)/(R*tandip)), fz0 is the phi value of the particle where it goes through z=0
// note that dfdz has a physical ambiguity in q*zdir.
    double _dfdz, _fz0;

    HelixFitHackResult(TrkDef const& tdef) : _hdef(tdef),  _fit(TrkErrCode::fail),_radius(-1.0),_dfdz(0.0),_fz0(0.0) {}
    HelixFitHackResult(HelixDef const& hdef) : _hdef(hdef),  _fit(TrkErrCode::fail),_radius(-1.0),_dfdz(0.0),_fz0(0.0) {}
    HelixFitHackResult& operator =(HelixFitHackResult const& other);
 };

// utility struct; value plus error
  struct VALERR {
    double _val;
    double _err;
  };
  struct FZ {
    VALERR _phi;
    double _z;
  };

// utility struct
  struct XYZP {
// straw hit index
    size_t _ind;
// position
    CLHEP::Hep3Vector _pos;
// ambiguity-resolved phi angle
    double _phi;
// flag
    StrawHitFlag _flag;
// wire direction
    CLHEP::Hep3Vector _wdir;
// straw radial direction, perp to Z and wire direction
    CLHEP::Hep3Vector _sdir;
// errors are asymmetric; along the wire is given by time division, perp to the wire by the straw size/sqrt(12)
    double _perr,_rerr;
// initialize some variables on construction
    XYZP():_phi(0.0){}
    XYZP(size_t index,StrawHit const& sh, StrawHitPosition const& shp,Straw const& straw);
    XYZP(size_t ind,CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& wdir, double werr, double serr) :
      _ind(ind),_pos(pos),_phi(_pos.phi()),_wdir(wdir),_sdir(wdir.y(),-wdir.x(),0.0),_perr(_efac*werr),_rerr(_efac*serr) {}
 
// radial position information
    virtual void rinfo(CLHEP::Hep3Vector const& center, VALERR& rad) const;
    virtual void finfo(CLHEP::Hep3Vector const& center, VALERR& phi) const;
    bool use       () const;
    bool stereo    () const;
    bool isOutlier ();
    void setOutlier();
    void setUse(bool use);
    static double _efac;
// flag bits to define use
    static StrawHitFlag _useflag;
  };

  typedef std::vector<XYZP> XYZPVector;
// struct to hold AGE sums
  struct SUMS {
// weighted (s)ums of (c)osine and (s)in for points on (c)ircumference, (o)utside the median radius, or (i)nside the median radius
    double _scc, _ssc, _sco, _sso, _sci, _ssi;
    unsigned _nc, _no, _ni;
    SUMS() : _scc(0.0),_ssc(0.0),_sco(0.0),_sso(0.0),_sci(0.0),_ssi(0.0){}
    void clear() { _scc = _ssc = _sco = _sso = _sci = _ssi = 0.0;
      _nc = _no = _ni = 0; }
  };

  class HelixFitHack
  {
  public:
// parameter set should be passed in on construction
    explicit HelixFitHack(fhicl::ParameterSet const&);
    virtual ~HelixFitHack();
// main function: given a track definition, find the helix parameters
    bool findHelix(HelixFitHackResult& myfit);
// allow passing in the struct by hand
    bool findHelix(XYZPVector& xyzp, HelixFitHackResult& myfit);
// convert to BaBar helix parameters.  Also return an error estimate
    void helixParams (HelixFitHackResult const& helix,CLHEP::HepVector& pvec,CLHEP::HepVector& perr) const;
    TH1F *_hDist;
    double _chi2nFindZ, _eventToLook;
  protected:
// utlity functions
    bool findXY        (XYZPVector& xyzp,HelixFitHackResult& myhel);
    bool findXY_new    (XYZPVector& xyzp,HelixFitHackResult& myhel);
    bool findZ         (XYZPVector& xyzp,HelixFitHackResult& myhel);

    bool initCircle    (XYZPVector const& xyzp,HelixFitHackResult& myhel);
    bool initCircle_new(XYZPVector const& xyzp,HelixFitHackResult& myhel);
  private:
    void fillXYZP(HelixDef const& mytrk, XYZPVector& xyzp);
// cached bfield accessor
    double bz() const;
// utility function to resolve phi wrapping    
    static double deltaPhi(double phi1, double phi2);
// diagnostics
    void plotXY(HelixDef const& mytrk, XYZPVector const& xyzp, HelixFitHackResult const& myhel) const;
    void plotZ(HelixDef const& mytrk, XYZPVector const& xyzp, HelixFitHackResult const& myhel) const;
// find the Absolute Geometric Error.  Returns the median radius as well.
    bool findCenterAGE(XYZPVector const& xyzp,Hep3Vector& center, double& rmed, double& age,bool useweights=false);
    void findAGE(XYZPVector const& xyzp, Hep3Vector const& center,double& rmed, double& age,bool useweights=false);
    void fillSums(XYZPVector const& xyzp, Hep3Vector const& center,double rmed,SUMS& sums,bool useweights=false);
    void filterXY(XYZPVector& xyzp, Hep3Vector const& center,double rmed,bool& changed);
    void filterDist(XYZPVector& xyzp);
//12 - 10 - 2013 Gianipez added few functions for integrating a new
//pattern recognition
    
    void doPatternRecognition(XYZPVector& xyzp, HelixFitHackResult& mytrk);
    void findTrack(XYZPVector& xyzp, int seedIndex,
		   double&chi2,int &countGoodPoint,
		   HelixFitHackResult& mytrk, bool cleanPatetrn=false);
    void calculateTrackParameters(Hep3Vector& p0, double&radius,
				  double& phi0, double& tanLambda,
				  Hep3Vector p1, Hep3Vector p2,
				  Hep3Vector p3, HelixFitHackResult& mytrk,
				  bool cleanPattern=false);

    void calculateDfDz(double &phi0, double &phi1, 
		       double &z0,  double &z1,
		       double &dfdz);
//--------------------------------------------------------------------------------//



// configuration parameters
    int _diag,_debug;
    double _mindelta; // minimum slope difference to use a triple in circle center initialization
    unsigned _minnhit; // minimum # of hits to work with
    unsigned _minnstereo; // minimum # of stereo hits
    double _lambda0,_lstep,_minlambda; // parameters for AGE center determination
    unsigned _maxniter; // maxium # of iterations to global minimum
    double _nsigma; // # of sigma for filtering outlyers
    double _nssigma; // # of sigma for filtering stereo time division
    double _minzsep, _maxzsep; // Z separation of points for pitch estimate
    double _maxdz, _maxdot; // stereo selection parameters
    double _rbias;  // robust fit parameter bias
    double _efac; // error factor
    double _rhomin, _rhomax; // crude cuts on tranvservse radius for stereo hits
    double _mindist; // minimum distance between points used in circle initialization
    double _maxdist; // maximum distance in hits
    double _pmin, _pmax; // range of total momentum
    double _tdmin, _tdmax; // range of abs(tan(dip)
    double _rcmin,_rcmax; // maximum transverse radius of circle
    double _sfactor; // stereo hit error factor
    bool _forcep; // force the p/pt to be in range (true), or exclude fits outside that range (false)
    bool _xyweights,_zweights; // weight points by estimated errors 
    bool _filter; // filter hits
    bool _plotall; // plot also failed fits
    bool _usetarget; // constrain to target when initializing
    bool _allstereo; // Use all stereo combinations (true) or only use each hit once
    mutable double _bz; // cached value of Field Z component at the tracker origin
// cached value of radius and pitch sign: these depend on the particle type
// and direction
    double _rmin, _rmax, _smin, _smax, _dfdzsign;

 //ROOT file for collecting histograms//
    // TFile *_fOut;
//-    -------------------------------------------------//
  };
}
#endif
