//
// Object to perform helix fit to straw hits
//
//
#ifndef TrkReco_RobustHelixFit_HH
#define TrkReco_RobustHelixFit_HH

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"
//#include "Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/TrkReco/inc/RobustHelixFinderData.hh"

#include "Offline/Mu2eUtilities/inc/MedianCalculator.hh"

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
    FZ(XYZVectorF const& hpos, XYZVectorF const& center);
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

      struct Config {
        fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
        fhicl::Atom<int> debugLevel{fhicl::Name("debugLevel"), fhicl::Comment("debug level")};
        fhicl::Atom<int> CircleInitType{fhicl::Name("CircleInitType"), fhicl::Comment("initialize circle")};
        fhicl::Atom<int> CircleFitType{fhicl::Name("CircleFitType"), fhicl::Comment("fit circle")};
        fhicl::Sequence<std::string> DontUseFlag{fhicl::Name("DontUseFlag"), fhicl::Comment("don't use flag")};
        fhicl::Atom<unsigned> minNStrawHits{fhicl::Name("minNStrawHits"), fhicl::Comment("minimum number of straw hits")};
        fhicl::Atom<unsigned> minNHit{fhicl::Name("minNHit"), fhicl::Comment("minimum number of hits")};
        fhicl::Atom<float> minXYResid{fhicl::Name("minXYResid"), fhicl::Comment("minimum XY residual")};
        fhicl::Atom<float> lambda0{fhicl::Name("lambda0"), fhicl::Comment("initial lambda")};
        fhicl::Atom<float> lstep{fhicl::Name("lstep"), fhicl::Comment("lambda step")};
        fhicl::Atom<float> minlambda{fhicl::Name("minlambda"), fhicl::Comment("minimum lambda")};
        fhicl::Atom<float> minDfDz{fhicl::Name("minDfDz"), fhicl::Comment("minimum Df Dz")};
        fhicl::Atom<float> maxDfDz{fhicl::Name("maxDfDz"), fhicl::Comment("maximum Df Dz")};
        fhicl::Atom<int> nLoopsdfdz{fhicl::Name("nLoopsdfdz"), fhicl::Comment("number of loops df dz")};
        fhicl::Atom<unsigned> NPhiHistBins{fhicl::Name("NPhiHistBins"), fhicl::Comment("number of phi hist bins")};
        fhicl::Atom<float> PhiHistRangeFactor{fhicl::Name("PhiHistRangeFactor"), fhicl::Comment("phi hist range factor")};
        fhicl::Atom<unsigned> MinNPhi{fhicl::Name("MinNPhi"), fhicl::Comment("minimum number phi")};
        fhicl::Atom<unsigned> maxniter{fhicl::Name("maxniter"), fhicl::Comment("maximum number of iterations")};
        fhicl::Atom<float> minzsep{fhicl::Name("minzsep"), fhicl::Comment("minimum z separation")};
        fhicl::Atom<float> maxzsep{fhicl::Name("maxzsep"), fhicl::Comment("maximum z separation")};
        fhicl::Atom<float> mindphi{fhicl::Name("mindphi"), fhicl::Comment("minimum change in phi")};
        fhicl::Atom<float> sigmaPhi{fhicl::Name("sigmaPhi"), fhicl::Comment("sigma phi")};
        fhicl::Atom<float> mindist{fhicl::Name("mindist"), fhicl::Comment("minimum distance of combo hits in triplet in mm")};
        fhicl::Atom<float> maxdist{fhicl::Name("maxdist"), fhicl::Comment("maximum distance of combo hits in triplet in mm")};
        fhicl::Atom<float> maxdxy{fhicl::Name("maxdxy"), fhicl::Comment("maximum dxy")};
        fhicl::Atom<float> maxXDPhi{fhicl::Name("maxXDPhi"), fhicl::Comment("maximum XDPhi")};
        fhicl::Atom<float> minR{fhicl::Name("minR"), fhicl::Comment("minimum helix radius in mm")};
        fhicl::Atom<float> maxR{fhicl::Name("maxR"), fhicl::Comment("maximum helix radius in mm")};
        fhicl::Atom<float> minCenterR{fhicl::Name("minCenterR"), fhicl::Comment("minimum center radius in mm")};
        fhicl::Atom<float> maxCenterR{fhicl::Name("maxCenterR"), fhicl::Comment("maximum center radius in mm")};
        fhicl::Atom<float> minAbsLambda{fhicl::Name("minAbsLambda"), fhicl::Comment("minimum absolute lambda")};
        fhicl::Atom<float> maxAbsLambda{fhicl::Name("maxAbsLambda"), fhicl::Comment("maximum absolute lambda")};
        fhicl::Atom<bool> targetintersect{fhicl::Name("targetintersect"), fhicl::Comment("target intersect")};
        fhicl::Atom<bool> TripleRadius{fhicl::Name("TripleRadius"), fhicl::Comment("triple radius")};
        fhicl::Atom<bool> HitErrorWeight{fhicl::Name("HitErrorWeight"), fhicl::Comment("hit error weight")};
        fhicl::Atom<bool> UseCaloCluster{fhicl::Name("UseCaloCluster"), fhicl::Comment("use calo cluster bool")};
        fhicl::Atom<float> CaloClusterWeight{fhicl::Name("CaloClusterWeight"), fhicl::Comment("calorimeter cluster weight in units of non-stereo hits")};
        fhicl::Atom<float> targetradius{fhicl::Name("targetradius"), fhicl::Comment("radius of stopping target constraint, effective target radius in mm")};
        fhicl::Atom<float> trackerradius{fhicl::Name("trackerradius"), fhicl::Comment("tracker outer radius in mm")};
        fhicl::Atom<float> RadiusWindow{fhicl::Name("RadiusWindow"), fhicl::Comment("radius window for calling a point to be 'on' the helix in AGG fit in mm")};

        fhicl::Atom<unsigned> ntripleMin{fhicl::Name("ntripleMin"), fhicl::Comment("minimum number of triplets")};
        fhicl::Atom<unsigned> ntripleMax{fhicl::Name("ntripleMax"), fhicl::Comment("maximum number of triplets")};
        fhicl::Atom<bool> use_initFZ_from_dzFrequency{fhicl::Name("use_initFZ_from_dzFrequency"), fhicl::Comment("use initFZ from dz frequency")};
        fhicl::Atom<float> initFZMinLambda{fhicl::Name("initFZMinLambda"), fhicl::Comment("init FZ minimum lambda")};
        fhicl::Atom<float> initFZMaxLambda{fhicl::Name("initFZMaxLambda"), fhicl::Comment("init FZ maximum lambda")};
        fhicl::Atom<float> initFZStepLambda{fhicl::Name("initFZStepLambda"), fhicl::Comment("init FZ lambda step")};
        fhicl::Atom<float> fitFZMinLambda{fhicl::Name("fitFZMinLambda"), fhicl::Comment("fit FZ minimum lambda")};
        fhicl::Atom<float> fitFZMaxLambda{fhicl::Name("fitFZMaxLambda"), fhicl::Comment("fit FZ maximum lambda")};
        fhicl::Atom<float> fitFZStepLambda{fhicl::Name("fitFZStepLambda"), fhicl::Comment("fit FZ lambda step")};
        fhicl::Atom<float> minArea{fhicl::Name("minArea"), fhicl::Comment("minimum triplet area")};
        fhicl::Atom<float> initFZFrequencyNSigma{fhicl::Name("initFZFrequencyNSigma"), fhicl::Comment("init FZ N Sigma frequency")};
        fhicl::Atom<int> initFZFrequencyBinsToIntegrate{fhicl::Name("initFZFrequencyBinsToIntegrate"), fhicl::Comment("init FZ frequency number of bins to integrate")};
        fhicl::Atom<int> initFZFrequencyArraySize{fhicl::Name("initFZFrequencyArraySize"), fhicl::Comment("init FZ frequency array size")};
        fhicl::Atom<int> initFZFrequencyNMaxPeaks{fhicl::Name("initFZFrequencyNMaxPeaks"), fhicl::Comment("init FZ frequency number of max peaks")};
        fhicl::Atom<float> initFZFrequencyTolerance{fhicl::Name("initFZFrequencyTolerance"), fhicl::Comment("init FZ frequency tolerance")};
      };

      explicit RobustHelixFit(const Config& config);
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
      void findAGE(RobustHelixFinderData  const& helixData, XYZVectorF const& center,float& rmed, float& age);
      void fillSums(RobustHelixFinderData const& helixData, XYZVectorF const& center,float rmed,AGESums& sums);
      void forceTargetInter(XYZVectorF& center, float& radius);

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
