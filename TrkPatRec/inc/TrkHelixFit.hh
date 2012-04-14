//
// Object to perform helix fit to straw hits
//
// $Id: TrkHelixFit.hh,v 1.5 2012/04/14 17:54:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/04/14 17:54:02 $
//
#ifndef TrkHelixFit_HH
#define TrkHelixFit_HH

// framework
#include "fhiclcpp/ParameterSet.h"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// TrkHelixFit objects
#include "KalmanTests/inc/TrkDef.hh"
// BaBar
#include "TrkBase/TrkErrCode.hh"
//CLHEP
//#include "CLHEP/Matrix/Vector.h"
// C++
namespace mu2e 
{
// output struct
  struct TrkHelix {
// fit status
    TrkErrCode _fit; // error code from last fit
// circle parameters; the z center is ignored.
    CLHEP::Hep3Vector _center;
    double _radius;
// Z parameters; dfdz is the slope of phi vs z (=1/(R*tandip)), fz0 is the phi value of the particle where it goes through z=0
    double _dfdz, _fz0;
    TrkHelix() : _fit(TrkErrCode::fail),_radius(-1.0),_dfdz(0.0),_fz0(0.0) {}
 };
  
// utility struct
  struct RAD {
    double _radius;
    double _rerr;
  };
  
// utility struct
  struct XYZP {
    CLHEP::Hep3Vector _pos;
// wire direction
    CLHEP::Hep3Vector _wdir;
// straw radial direction, perp to Z and wire direction
    CLHEP::Hep3Vector _sdir;
// errors are asymmetric; along the wire is given by time division, perp to the wire by the straw size/sqrt(12)
    double _werr, _serr;
// use this hit or not
    bool _use;
// constructor
    XYZP(CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& wdir, double werr, double serr) :
    _pos(pos),_wdir(wdir),_sdir(wdir.y(),-wdir.x(),0.0),_werr(werr),_serr(serr), _use(true) {}
// compute radius and radial error, given the circle center
    void rinfo(CLHEP::Hep3Vector const& center, RAD& rad) const;
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
  
  class TrkHelixFit
  {
  public:
// parameter set should be passed in on construction
    explicit TrkHelixFit(fhicl::ParameterSet const&);
    virtual ~TrkHelixFit();
// main function: given a track definition, find the helix parameters
    bool findHelix(TrkDef const& mytrk,TrkHelix& myfit);
// convert to BaBar helix parameters.  Also return an error estimate
    void helixParams (TrkHelix const& helix,CLHEP::HepVector& pvec,CLHEP::HepVector& perr) const;
  private:
// utlity functions
    bool findXY(std::vector<XYZP>& xyzp,TrkHelix& myhel);
    bool findZ(std::vector<XYZP> const& xyzp,TrkHelix& myhel);
    bool initCircle(std::vector<XYZP> const& xyzp,TrkHelix& myhel);
    void fillXYZP(TrkDef const& mytrk, std::vector<XYZP>& xyzp);
    
// find the Absolute Geometric Error.  Returns the median radius as well.
    bool findCenterAGE(std::vector<XYZP> const& xyzp,Hep3Vector& center, double& rmed, double& age,bool useweights=false);
    void findAGE(std::vector<XYZP> const& xyzp, Hep3Vector const& center,double& rmed, double& age,bool useweights=false);
    void fillSums(std::vector<XYZP> const& xyzp, Hep3Vector const& center,double rmed,SUMS& sums,bool useweights=false);
    void filter(std::vector<XYZP>& xyzp, Hep3Vector const& center,double rmed,bool& changed);
// configuration parameters
    int _diag,_debug;
    double _mindelta; // minimum slope difference to use a triple in circle center initialization
    unsigned _minnhit; // minimum # of hits to work with
    unsigned _maxnhit; // maximum # of hits to use for initial center finding
    double _rfactor,_lambda0,_lstep,_minlambda; // parameters for AGE center determination
    unsigned _maxniter; // maxium # of iterations to global minimum
    double _nsigma; // # of sigma for filtering outlyers
    double _minzsep, _maxzsep; // Z separation of points for tandip estimate
    bool _target; // include target constraint
    double _tsig; // sigma of target constraint
    double _rbias;  // robust fit parameter bias
    double _sfac; // error factor on straw position perp to wire direction
    static double _targetz; // z position of target hit
// 
  };
}
#endif
