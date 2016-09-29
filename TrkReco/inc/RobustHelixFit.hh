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
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/HelixHit.hh"
// HelixFit objects
#include "RecoDataProducts/inc/HelixSeed.hh"
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

// struct to hold AGE sums
  struct AGESums {
// (s)ums of (c)osine and (s)in for points on (c)ircumference, (o)utside the median radius, or (i)nside the median radius
    double _scc, _ssc, _sco, _sso, _sci, _ssi;
    unsigned _nc, _no, _ni;
    AGESums() : _scc(0.0),_ssc(0.0),_sco(0.0),_sso(0.0),_sci(0.0),_ssi(0.0){}
    void clear() { _scc = _ssc = _sco = _sso = _sci = _ssi = 0.0;
      _nc = _no = _ni = 0; }
  };

  class RobustHelixFit
  {
  public:
// parameter set should be passed in on construction
    explicit RobustHelixFit(fhicl::ParameterSet const&);
    virtual ~RobustHelixFit();
// main function: given a helix seed and the hits, fill in the helix parameters.  Note
// the HelixSeed object is both input and output
    void findHelix(HelixSeed& myfit);
// Helicity of this fit
    Helicity const& helicity() const { return _helicity; }
  protected:
// utlity functions
    bool findXY(HelixHitCollection& hhits,RobustHelix& myhel);
    bool findZ(HelixHitCollection& hhits,RobustHelix& myhel);
    bool initCircle(HelixHitCollection const& hhits,RobustHelix& myhel);
  private:
// utility function to resolve phi wrapping    
    static double deltaPhi(double phi1, double phi2);
// find the Absolute Geometric Error.  Returns the median radius as well.
    bool findCenterAGE(HelixHitCollection const& hhits,CLHEP::Hep3Vector& center, double& rmed, double& age);
    void findAGE(HelixHitCollection const& hhits, CLHEP::Hep3Vector const& center,double& rmed, double& age);
    void fillSums(HelixHitCollection const& hhits, CLHEP::Hep3Vector const& center,double rmed,AGESums& sums);
    void filterXY(HelixHitCollection& hhits, CLHEP::Hep3Vector const& center,double rmed,bool& changed);
    void filterDist(HelixHitCollection& hhits);
    unsigned hitCount(HelixHitCollection const& hhits) const; // count good hits
// interact with HelixHits
    bool use(HelixHit const&) const;
    bool stereo(HelixHit const&) const;
    void setOutlier(HelixHit&) const;
    void setResolvedPhi(HelixHit&) const;
    void radInfo(CLHEP::Hep3Vector const& center, HelixHit const& hhit, VALERR& rad) const;
    void phiInfo(CLHEP::Hep3Vector const& center, HelixHit const& hhit, VALERR& phi) const;

// configuration parameters
    int _debug;
    StrawHitFlag _useflag, _dontuseflag;
    double _mindelta; // minimum slope difference to use a triple in circle center initialization
    unsigned _minnhit; // minimum # of hits to work with
    unsigned _minnstereo; // minimum # of stereo hits
    double _lambda0,_lstep,_minlambda; // parameters for AGE center determination
    unsigned _maxniter; // maxium # of iterations to global minimum
    double _minzsep, _maxzsep; // Z separation of points for pitch estimate
    double _mindphi, _maxdphi; // phi separation of points for pitch estimate
    double _rbias;  // robust fit parameter bias
    double _efac; // error factor
    double _mindist; // minimum distance between points used in circle initialization
    double _maxdist; // maximum distance in hits
    double _rmin,_rmax; // circle radius range
    double _tdmin, _tdmax; // range of abs(tan(dip)
    double _sfactor; // stereo hit error factor
    bool _force; // force the fit values to be in range
    bool _filterxy, _filterz; // filter hits
    bool _stereoinit; // require stereo hits to initialize
    bool _stereofit; // require stereo hits 
    bool _targetpoint; // use target as a point in the circle fit
    bool _targetinit; // require consistency with target when initializing circle
    bool _targetinter; // require fit to intersect the target
    double _targetradius; // target size to use in constraint or init
    double _trackerradius; // tracker radius to use in init
    double _rwind; // raidus window for defining points to be 'on' the helix
    double _rout; // radius difference for a hit to be an xy outlier
    double _pout; // phi difference for a hit to be a z outlier
    Helicity _helicity; // helicity value to look for.  This defines the sign of dphi/dz
    double _smin, _smax;
 };
}
#endif
