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
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// HelixFit objects
#include "TrkReco/inc/XYZP.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
// BaBar
#include "BTrk/TrkBase/TrkErrCode.hh"
//root
class TH1F;
// C+

namespace mu2e 
{

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
// main function: given a helix seed and the hits, fill in the helix parameters.  Note
// the HelixSeed object is both input and output
    void findHelix(StrawHitCollection const& shcol,
	StrawHitPositionCollection const& shpos, HelixSeed& myfit);
  protected:
// utlity functions
    bool findXY(XYZPVector& xyzp,RobustHelix& myhel);
    bool findZ(XYZPVector& xyzp,RobustHelix& myhel);
    bool initCircle(XYZPVector const& xyzp,RobustHelix& myhel);
  private:
// utility function to resolve phi wrapping    
    static double deltaPhi(double phi1, double phi2);
// find the Absolute Geometric Error.  Returns the median radius as well.
    bool findCenterAGE(XYZPVector const& xyzp,CLHEP::Hep3Vector& center, double& rmed, double& age,bool useweights=false);
    void findAGE(XYZPVector const& xyzp, CLHEP::Hep3Vector const& center,double& rmed, double& age,bool useweights=false);
    void fillSums(XYZPVector const& xyzp, CLHEP::Hep3Vector const& center,double rmed,SUMS& sums,bool useweights=false);
    void filterXY(XYZPVector& xyzp, CLHEP::Hep3Vector const& center,double rmed,bool& changed);
    void filterDist(XYZPVector& xyzp);
// configuration parameters
    int _debug;
    double _mindelta; // minimum slope difference to use a triple in circle center initialization
    unsigned _minnhit; // minimum # of hits to work with
    unsigned _minnstereo; // minimum # of stereo hits
    double _lambda0,_lstep,_minlambda; // parameters for AGE center determination
    unsigned _maxniter; // maxium # of iterations to global minimum
    double _nsigma; // # of sigma for filtering outlyers
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
    bool _xyweights,_zweights; // weight points by estimated errors 
    bool _filterxy, _filterz; // filter hits
    bool _stereoinit; // require stereo hits to initialize
    bool _stereofit; // require stereo hits 
    bool _targetpoint; // use target as a point in the circle fit
    bool _targetinit; // require consistency with target when initializing circle
    bool _targetinter; // require fit to intersect the target
    double _targetradius; // target size to use in constraint or init
    double _trackerradius; // tracker radius to use in init
    Helicity _helicity; // helicity value to look for.  This defines the sign of dphi/dz
// cached value of radius and pitch sign: these depend on the particle type
// and direction
    double _smin, _smax;
 };
}
#endif
