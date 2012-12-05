//
// Object to perform helix fit to straw hits
//
// $Id: HelixFit.hh,v 1.2 2012/12/05 18:48:21 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/12/05 18:48:21 $
//
#ifndef HelixFit_HH
#define HelixFit_HH

// framework
#include "fhiclcpp/ParameterSet.h"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// HelixFit objects
#include "KalmanTests/inc/TrkDef.hh"
// BaBar
#include "TrkBase/TrkErrCode.hh"
//CLHEP
//#include "CLHEP/Matrix/Vector.h"
//root
class TH1F;
// C+

namespace mu2e 
{
// output struct
  struct HelixFitResult {
    TrkDef _tdef; // must copy by value as references can't be re-assigned
// fit status
    TrkErrCode _fit; // error code from last fit
// circle parameters; the z center is ignored.
    CLHEP::Hep3Vector _center;
    double _radius;
// Z parameters; dfdz is the slope of phi vs z (=-sign(1.0,qBzdir)/(R*tandip)), fz0 is the phi value of the particle where it goes through z=0
// note that dfdz has a physical ambiguity in q*zdir.
    double _dfdz, _fz0;
    HelixFitResult(TrkDef const& tdef) : _tdef(tdef),  _fit(TrkErrCode::fail),_radius(-1.0),_dfdz(0.0),_fz0(0.0) {}
    HelixFitResult& operator =(HelixFitResult const& other);
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
    enum useBit{good=0,stereohit=1,stereopoint=2,outlier=3};
// position
    CLHEP::Hep3Vector _pos;
// ambiguity-resolved phi angle
    double _phi;
// use flag
    unsigned _use;
// subhits
    std::vector<XYZP*> _subhits;
// initialize some variables on construction
    XYZP():_phi(0.0),_use(1<<good){}
    XYZP(CLHEP::Hep3Vector const& pos) : _pos(pos),_phi(pos.phi()),_use(1<<good) {} 
// radial position information
    virtual void rinfo(CLHEP::Hep3Vector const& center, VALERR& rad) const = 0;
    virtual void finfo(CLHEP::Hep3Vector const& center, VALERR& phi) const = 0;
    virtual void setUse(bool use,useBit ibit=good) { if(use) _use |= (1<<ibit); else _use &= !(1<<ibit); }
    virtual bool use(useBit ibit=good) const {return (_use&(1<<ibit)) == (1<<ibit); }
  };
// struct for hits
  struct HitXYZP : public XYZP {
// straw hit index
    size_t _ind;
// wire direction
    CLHEP::Hep3Vector _wdir;
// straw radial direction, perp to Z and wire direction
    CLHEP::Hep3Vector _sdir;
// errors are asymmetric; along the wire is given by time division, perp to the wire by the straw size/sqrt(12)
    double _werr, _serr;
// constructor
    HitXYZP(size_t ind,CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& wdir, double werr, double serr) :
    XYZP(pos),_ind(ind),_wdir(wdir),_sdir(wdir.y(),-wdir.x(),0.0),_werr(werr),_serr(serr) {}
// compute radius and radial error, given the circle center
    virtual void rinfo(CLHEP::Hep3Vector const& center, VALERR& rad) const;
    virtual void finfo(CLHEP::Hep3Vector const& center, VALERR& phi) const;
    virtual bool use(useBit ibit=good) const;
  };

// struct for stereo pair  
  struct StereoXYZP : public XYZP {
    StereoXYZP(HitXYZP& h1,HitXYZP& h2);
    StereoXYZP(); // needed for std::vector
    HitXYZP* _h1;
    HitXYZP* _h2;
    double _d1, _d2;
    double _perr; // pixel size, assumed constant for now
    virtual void rinfo(CLHEP::Hep3Vector const& center, VALERR& rad) const;
    virtual void finfo(CLHEP::Hep3Vector const& center, VALERR& phi) const;
// override setting use flag to keep bookkeepping right.
    virtual void setUse(bool use,useBit ibit=good);
    virtual bool use(useBit ibit=good) const;
  };

// container struct to keep data coherent
  struct XYZPVector {
    std::vector<HitXYZP> _hxyzp;
    std::vector<StereoXYZP> _sxyzp;
    std::vector<XYZP*> _xyzp; // note: these objects are owned by the vectors above!!!
// insure consistency on construction
    XYZPVector(std::vector<HitXYZP> const& hits);
    XYZPVector(std::vector<HitXYZP> const& hits,std::vector<StereoXYZP> const& shits);
  };

// struct to hold AGE sums
  struct SUMS {
// weighted (s)ums of (c)osine and (s)in for points on (c)ircumference, (o)utside the median radius, or (i)nside the median radius
    double _scc, _ssc, _sco, _sso, _sci, _ssi;
    unsigned _nc, _no, _ni;
    SUMS() : _scc(0.0),_ssc(0.0),_sco(0.0),_sso(0.0),_sci(0.0),_ssi(0.0){}
    void clear() { _scc = _ssc = _sco = _sso = _sci = _ssi = 0.0;
      _nc = _no = _ni = 0; }
  };
  
  class HelixFit
  {
  public:
// parameter set should be passed in on construction
    explicit HelixFit(fhicl::ParameterSet const&);
    virtual ~HelixFit();
// main function: given a track definition, find the helix parameters
    bool findHelix(HelixFitResult& myfit);
// allow passing in the struct by hand
    bool findHelix(XYZPVector& xyzp, HelixFitResult& myfit);
// convert to BaBar helix parameters.  Also return an error estimate
    void helixParams (HelixFitResult const& helix,CLHEP::HepVector& pvec,CLHEP::HepVector& perr) const;
  protected:
// utlity functions
    bool findXY(std::vector<XYZP*>& xyzp,HelixFitResult& myhel);
    bool findZ(std::vector<XYZP*>& xyzp,HelixFitResult& myhel);
    bool initCircle(std::vector<XYZP*> const& xyzp,HelixFitResult& myhel);
  private:
    void fillHitXYZP(TrkDef const& mytrk, std::vector<HitXYZP>& xyzp);
    void findStereoPairs(std::vector<HitXYZP>& xyzp, std::vector<StereoXYZP>& spairs);
// cached bfield accessor
    double bz() const;
// utility function to resolve phi wrapping    
    static double deltaPhi(double phi1, double phi2);
// diagnostics
    void plotXY(TrkDef const& mytrk, std::vector<HitXYZP>const& xyzp, std::vector<StereoXYZP> const& pairs, HelixFitResult const& myhel) const;
    void plotZ(TrkDef const& mytrk, std::vector<HitXYZP> const& xyzp, std::vector<StereoXYZP> const& pairs, HelixFitResult const& myhel) const;
// find the Absolute Geometric Error.  Returns the median radius as well.
    bool findCenterAGE(std::vector<XYZP*> const& xyzp,Hep3Vector& center, double& rmed, double& age,bool useweights=false);
    void findAGE(std::vector<XYZP*> const& xyzp, Hep3Vector const& center,double& rmed, double& age,bool useweights=false);
    void fillSums(std::vector<XYZP*> const& xyzp, Hep3Vector const& center,double rmed,SUMS& sums,bool useweights=false);
    void filterXY(std::vector<XYZP*>& xyzp, Hep3Vector const& center,double rmed,bool& changed);
    void filterDist(std::vector<XYZP*>& xyzp);
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
    double _sfac; // error factor  straw position perp to wire direction
    double _rhomin, _rhomax; // crude cuts on tranvservse radius for stereo hits
    double _mindist; // minimum distance between points used in circle initialization
    double _maxdist; // maximum distance in hits
    double _pmin, _pmax; // range of total momentum
    double _tdmin, _tdmax; // range of abs(tan(dip)
    double _rcmin,_rcmax; // maximum transverse radius of circle
    bool _forcep; // force the p/pt to be in range (true), or exclude fits outside that range (false)
    bool _useweights; // weight points by estimated errors 
    bool _filter; // filter hits
    bool _plotall; // plot also failed fits
    bool _usetarget; // constrain to target when initializing
    bool _allstereo; // Use all stereo combinations (true) or only use each hit once
    mutable double _bz; // cached value of Field Z component at the tracker origin
// cached value of radius and pitch sign: these depend on the particle type
// and direction
    double _rmin, _rmax, _smin, _smax, _dfdzsign;
    mutable TH1F *_sdist, *_spull;
  };
}
#endif
