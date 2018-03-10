//
// Object to perform helix fit to straw hits
//
// $Id: HelixFit.hh,v 1.8 2014/07/10 14:47:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/07/10 14:47:26 $
//
#ifndef TrkReco_RobustHelixFit_HH
#define TrkReco_RobustHelixFit_HH

#include "fhiclcpp/ParameterSet.h"
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "TH1F.h"


namespace mu2e 
{
  // simple struct to keep track of azimuth/radius projection
  struct FZ 
  {
    float _phi;
    float _z;
    FZ(XYZVec const& hpos, XYZVec const& center);
  };

  // struct to hold AGE sums
  struct AGESums {
    // (s)ums of (c)osine and (s)in for points on (c)ircumference, (o)utside the median radius, or (i)nside the median radius
    float _scc, _ssc, _sco, _sso, _sci, _ssi;
    unsigned _nc, _no, _ni;
    AGESums() : _scc(0.0),_ssc(0.0),_sco(0.0),_sso(0.0),_sci(0.0),_ssi(0.0){}
    void clear() { _scc = _ssc = _sco = _sso = _sci = _ssi = 0.0;
      _nc = _no = _ni = 0; }
  };


  class RobustHelixFit
  {
    public:
      
      enum CircleFit {median=0, AGE=1 };

      explicit RobustHelixFit(fhicl::ParameterSet const&);
      virtual ~RobustHelixFit();

      void fitCircle(HelixSeed& hseed);
      void fitFZ(HelixSeed& hseed);
      bool goodHelix(RobustHelix const& rhel);
      Helicity const& helicity() const { return _helicity; }

    private:

      void fitHelix(HelixSeed& hseed);
      void fitCircleMedian(HelixSeed& hseed);
      void fitCircleAGE(HelixSeed& hseed);
      bool initFZ(ComboHitCollection& hhits,RobustHelix& myhel);
      void findAGE(HelixSeed const& hseed, XYZVec const& center,float& rmed, float& age);
      void fillSums(HelixSeed const& hseed, XYZVec const& center,float rmed,AGESums& sums);

      bool goodCircle(RobustHelix const& rhel);
      bool goodFZ(RobustHelix const& rhel);

      void forceTargetInter(XYZVec& center, float& radius);

      bool use(ComboHit const&) const;
      bool stereo(ComboHit const&) const;
      void setOutlier(ComboHit&) const;

      static float deltaPhi(float phi1, float phi2);
      void initPhi(ComboHit& hh, RobustHelix const& myhel) const;
      bool resolvePhi(ComboHit& hh, RobustHelix const& myhel) const;
      float hitWeight(ComboHit const& hhit) const;
      bool goodLambda(Helicity const& h, float lambda) const;

      int _debug;
      CircleFit _cfit; // type of circle fit
      StrawHitFlag _useflag, _dontuseflag;
      unsigned _minnhit; // minimum # of hits to work with
      float _lambda0,_lstep,_minlambda; // parameters for AGE center determination
      unsigned _nphibins; // # of bins in histogram for phi at z intercept
      float _phifactor; // range factr for phi z intercept histogram 
      unsigned _minnphi; // minimum # of entries in max bin of phi intercept histogram 
      unsigned _maxniter; // maxium # of iterations to global minimum
      float _minzsep, _maxzsep; // Z separation of points for pitch estimate
      float _mindphi, _maxdphi; // phi separation of points for pitch estimate
      float _mindist; // minimum distance between points used in circle initialization
      float _maxdist; // maximum distance in hits
      float _rmin,_rmax; // circle radius range
//      float _mindelta; // minimum slope difference to use a triple in circle center initialization
      float _minarea2; // minimum triangle area for triple (squared)
      float _lmin, _lmax; // range of lambda = dz/dphi
      bool _targetpoint; // use target as a point in the circle fit
      bool _targetcon; // require consistency with target
      bool _targetinter; // require fit to intersect the target
      bool _tripler; // use triples to compute r
      bool _errrwt; // use hit errors to weight radius calculation 
      bool _usecc; // use the calorimeter cluster in the fit (transverse only)
      float _ccwt; // weight of a calorimeter cluster in non-stereo hit units
      float _targetradius; // target size to use in constraint or init
      float _trackerradius; // tracker radius to use in init
      float _rwind; // raidus window for defining points to be 'on' the helix
      Helicity _helicity; // helicity value to look for.  This defines the sign of dphi/dz
      TH1F _hphi;
      unsigned _ntripleMax;
 };
}
#endif
