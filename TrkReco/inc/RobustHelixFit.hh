//
// Object to perform helix fit to straw hits
//
//
#ifndef TrkReco_RobustHelixFit_HH
#define TrkReco_RobustHelixFit_HH

#include "fhiclcpp/ParameterSet.h"
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"
//#include "Mu2eUtilities/inc/LsqSums4.hh"
#include "TrkReco/inc/RobustHelixFinderData.hh"

#include "Mu2eUtilities/inc/MedianCalculator.hh"

//using namespace ROOT::Math::VectorUtil;

namespace mu2e 
{

  class Calorimeter;
  class Tracker;

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
      
    enum CircleFit {median=0, mean, AGE };

    explicit RobustHelixFit(fhicl::ParameterSet const&);
    virtual ~RobustHelixFit();

    bool initCircle(RobustHelixFinderData& helixData, bool forceTargetCon, bool useTripleAreaWt=false);
    void fitCircle(RobustHelixFinderData& helixData, bool forceTargetCon, bool useTripleAreaWt=false);
    bool initFZ(RobustHelixFinderData& helixData, int initHitPhi=1);
    bool initFZ_2(RobustHelixFinderData& helixData);
    bool initFZ_from_dzFrequency(RobustHelixFinderData& helixData, int initHitPhi=1);
    bool fillArrayDz(RobustHelixFinderData& HelixData, std::vector<int> &v,float &bin_size, float& startDz);
    bool extractFZ0(RobustHelixFinderData& HelixData, float& fz0);
    bool extractLambdaFromDzHist(int *hist_sum, float& lambda);
    void findHistPeaks(std::vector<int> &input, int bin_size, 
		       float &start_dz,
		       std::vector<float> &xPeak, std::vector<float> &xSigma, std::vector<float>&swmax, std::vector<int> &iPeak, int &first_peak, int &peaks_found);
    void fitFZ(RobustHelixFinderData& helixData);
    void fitFZ_2(RobustHelixFinderData& helixData, int weightMode=1);
    bool goodHelix(RobustHelix const& rhel);
    Helicity const& helicity() const { return _helicity; }

    bool goodCircle(RobustHelix const& rhel);
    bool goodFZ(RobustHelix const& rhel);

    //function used to evaluate the hit weight used in the XY fit
    float evalWeightXY  (const ComboHit& Hit, XYVec& Center);
    float evalWeightZPhi(const ComboHit& Hit, XYVec& Center, float Radius);

    void  setTracker    (const Tracker*    Tracker) { _tracker     = Tracker; }
    void  setCalorimeter(const Calorimeter* Cal    ) { _calorimeter = Cal    ; }

    //    bool  targetcon()   {return _targetcon; }

    const Tracker*            _tracker;
    const Calorimeter*         _calorimeter;

    void fitCircleMedian(RobustHelixFinderData& helixData, bool forceTargetCon, bool useTripleAreaWt=false);
    
    float lambdaMin()  { return _lmin; }
    float lambdaMax()  { return _lmax; }

  private:

    void fitHelix(RobustHelixFinderData& helixData, bool forceTargetCon, bool useTripletAreaWt=false);
    void fitCircleAGE(RobustHelixFinderData& helixData);
    void fitCircleMean(RobustHelixFinderData& helixData);
    void findAGE(RobustHelixFinderData  const& helixData, XYZVec const& center,float& rmed, float& age);
    void fillSums(RobustHelixFinderData const& helixData, XYZVec const& center,float rmed,AGESums& sums);
    void forceTargetInter(XYZVec& center, float& radius);

    bool use(ComboHit const&) const;
    bool stereo(ComboHit const&) const;
    void setOutlier(ComboHit&) const;

    static float deltaPhi(float phi1, float phi2);
    void initPhi(ComboHit& hh, RobustHelix const& myhel) const;
    bool resolvePhi(ComboHit& hh, RobustHelix const& myhel) const;
    float hitWeight(ComboHit const& hhit) const;
    bool goodLambda(Helicity const& h, float lambda) const;

   

    int _diag;
    int _debug;
    CircleFit _cinit, _cfit; // type of circle fit
    StrawHitFlag _useflag, _dontuseflag;
    int      _minnsh;  // minimum # of StrawHits
    unsigned _minnhit; // minimum # of hits to work with
    float _minxyresid; // minimum distance used in the circle fit to be clusterized. units are mm
    float _lambda0,_lstep,_minlambda; // parameters for AGE center determination
    float _mindfdz, _maxdfdz;//parameters use for findDfDz function
    int   _nLoopsdfdz;//parameter for number of loops included in DfDz fit
    unsigned _nphibins; // # of bins in histogram for phi at z intercept
    float _phifactor; // range factr for phi z intercept histogram 
    unsigned _minnphi; // minimum # of entries in max bin of phi intercept histogram 
    unsigned _maxniter; // maxium # of iterations to global minimum
    float _minzsep, _maxzsep; // Z separation of points for pitch estimate
    float _mindphi, _maxdphi; // phi separation of points for pitch estimate
    float _sigmaPhi; //approximated uncertanty on the ComboHit helix-phi coordinate
    float _mindist; // minimum distance between points used in circle initialization
    float _maxdist; // maximum distance in hits
    float _maxdxy; // maximum distance in hits after the triplet loop in fitCiircleMedian
    float _maxXDPhi;//maximum normalized residual for a hit in the z-phi fit
    float _rmin,_rmax; // circle radius range
    float _rcmin,_rcmax; // circle centerradius range
    //      float _mindelta; // minimum slope difference to use a triple in circle center initialization
    float _minarea2; // minimum triangle area for triple (squared)
    float _lmin, _lmax; // range of lambda = dz/dphi
    bool _targetpoint; // use target as a point in the circle fit
    //    bool _targetcon; // require consistency with target
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
    unsigned _ntripleMin, _ntripleMax;
    bool     _use_initFZ_from_dzFrequency;
    float    _initFZFrequencyNSigma;
    int      _initFZFrequencyBinsToIntegrate;
    int      _initFZFrequencyArraySize;
    int      _initFZFrequencyNMaxPeaks;
    float    _initFZFrequencyTolerance;
    unsigned _initFZNBins;
    float    _initFZMinL, _initFZMaxL, _initFZStepL;
    unsigned _fitFZNBins;
    float    _fitFZMinL, _fitFZMaxL, _fitFZStepL;
    MedianCalculator  _medianCalculator;
  };
}
#endif
