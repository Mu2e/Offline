//
// Object to perform helix fit to straw hits
//
// $Id: HelixFit.hh,v 1.8 2014/07/10 14:47:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/07/10 14:47:26 $
//
#ifndef HelixFit_HH
#define HelixFit_HH

// framework
#include "fhiclcpp/ParameterSet.h"
// data
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// HelixFit objects
#include "TrkReco/inc/HelixDef.hh"
#include "TrkReco/inc/HelixFitResult.hh"
// BaBar
#include "BTrk/TrkBase/TrkErrCode.hh"
//root
class TH1F;
// C+

namespace mu2e 
{

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
    int _ind;
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
// if this hit is from a MC conversion
    bool _conversion;
// initialize some variables on construction
    XYZP():_ind(-1),_phi(0.0),_perr(1000.0),_rerr(1000.0),_conversion(false){}
    XYZP(size_t index,StrawHit const& sh, StrawHitPosition const& shp,Straw const& straw);
    XYZP(size_t ind,CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& wdir, double werr, double serr);
     XYZP(CLHEP::Hep3Vector const& pos, double size);
// radial position information
    virtual void rinfo(CLHEP::Hep3Vector const& center, VALERR& rad) const;
    virtual void finfo(CLHEP::Hep3Vector const& center, VALERR& phi) const;
    bool use() const;
    bool stereo() const;
    void setOutlier();
    void setUse(bool use);
    static double _efac;
// flag bits to define use
    static StrawHitFlag _useflag, _dontuseflag;
    static int _debug;
// for checking if it's a conversion hit
    bool conversion() const { return _conversion; }
    void setConversion(bool conv) { _conversion = conv; }
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
  
  class RobustHelixFit
  {
  public:
// parameter set should be passed in on construction
    explicit RobustHelixFit(fhicl::ParameterSet const&);
    virtual ~RobustHelixFit();
// main function: given a track definition, find the helix parameters
    bool findHelix(HelixFitResult& myfit,bool plothelix=false);
// allow passing in the struct by hand
    bool findHelix(XYZPVector& xyzp, HelixFitResult& myfit);
// convert to BaBar helix parameters.  Also return an error estimate
    void helixParams (HelixFitResult const& helix,CLHEP::HepVector& pvec,CLHEP::HepVector& perr) const;
  protected:
// utlity functions
    bool findXY(XYZPVector& xyzp,HelixFitResult& myhel);
    bool findZ(XYZPVector& xyzp,HelixFitResult& myhel);
    bool initCircle(XYZPVector const& xyzp,HelixFitResult& myhel);
  private:
    void fillXYZP(HelixDef const& mytrk, XYZPVector& xyzp);
// cached bfield accessor
    double bz() const;
// utility function to resolve phi wrapping    
    static double deltaPhi(double phi1, double phi2);
// diagnostics
    void plotXY(HelixDef const& mytrk, XYZPVector const& xyzp, HelixFitResult const& myhel) const;
    void plotZ(HelixDef const& mytrk, XYZPVector const& xyzp, HelixFitResult const& myhel) const;
// find the Absolute Geometric Error.  Returns the median radius as well.
    bool findCenterAGE(XYZPVector const& xyzp,CLHEP::Hep3Vector& center, double& rmed, double& age,bool useweights=false);
    void findAGE(XYZPVector const& xyzp, CLHEP::Hep3Vector const& center,double& rmed, double& age,bool useweights=false);
    void fillSums(XYZPVector const& xyzp, CLHEP::Hep3Vector const& center,double rmed,SUMS& sums,bool useweights=false);
    void filterXY(XYZPVector& xyzp, CLHEP::Hep3Vector const& center,double rmed,bool& changed);
    void filterDist(XYZPVector& xyzp);
// configuration parameters
    int _diag,_debug;
    double _mindelta; // minimum slope difference to use a triple in circle center initialization
    unsigned _minnhit; // minimum # of hits to work with
    unsigned _minnstereo; // minimum # of stereo hits
    double _lambda0,_lstep,_minlambda; // parameters for AGE center determination
    unsigned _maxniter; // maxium # of iterations to global minimum
    double _nsigma; // # of sigma for filtering outlyers
    double _minzsep, _maxzsep; // Z separation of points for pitch estimate
    double _rbias;  // robust fit parameter bias
    double _efac; // error factor
    double _mindist; // minimum distance between points used in circle initialization
    double _maxdist; // maximum distance in hits
    double _pmin, _pmax; // range of total momentum
    double _tdmin, _tdmax; // range of abs(tan(dip)
    double _rcmin,_rcmax; // maximum transverse radius of circle
    double _sfactor; // stereo hit error factor
    bool _forcep; // force the p/pt to be in range (true), or exclude fits outside that range (false)
    bool _xyweights,_zweights; // weight points by estimated errors 
    bool _filter; // filter hits
    bool _stereoinit; // require stereo hits to initialize
    bool _stereofit; // require stereo hits 
    bool _targetpoint; // use target as a point in the circle fit
    bool _targetinit; // require consistency with target when initializing circle
    bool _targetinter; // require fit to intersect the target
    double _targetradius; // target size to use in constraint or init
    double _trackerradius; // tracker radius to use in init
    mutable double _bz; // cached value of Field Z component at the tracker origin
// cached value of radius and pitch sign: these depend on the particle type
// and direction
    double _rmin, _rmax, _smin, _smax, _dfdzsign;
// diagnostic histograms
    TH1F *_rdiff, *_fdiff;
    TH1F *_rpull, *_fpull;
  };
}
#endif
