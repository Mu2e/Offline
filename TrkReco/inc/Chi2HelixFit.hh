//
// Object to perform helix fit to straw hits
//
//
#ifndef TrkReco_Chi2HelixFit_HH
#define TrkReco_Chi2HelixFit_HH

#include "fhiclcpp/ParameterSet.h"
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"
#include "Mu2eUtilities/inc/LsqSums4.hh"
#include "TrkReco/inc/RobustHelixFit.hh"
#include "TrkReco/inc/RobustHelixFinderData.hh"

using namespace ROOT::Math::VectorUtil;

namespace mu2e 
{

  class Calorimeter;
  class Tracker;

  class Chi2HelixFit
  {
  public:

    explicit Chi2HelixFit(fhicl::ParameterSet const&);
    virtual ~Chi2HelixFit();

    bool initChi2Circle(RobustHelixFinderData& helixData, bool TargetCon);
    void fitChi2Circle(RobustHelixFinderData& helixData, bool TargetCon);

    bool initFZ_2(RobustHelixFinderData& helixData);
    void initFitChi2FZ(RobustHelixFinderData& helixData);
    void fitChi2FZ(RobustHelixFinderData& helixData, int weightMode=1);

    bool goodHelixChi2(RobustHelixFinderData& helixData);
    void defineHelixParams(RobustHelixFinderData& helixData);

    bool goodZPhiFit(LsqSums4& lsqsum);

    //function used to evaluate the hit weight used in the XY fit
    float evalWeightXY  (const ComboHit& Hit, XYVec& Center);
    float evalWeightZPhi(const ComboHit& Hit, XYVec& Center, float Radius);

    //function to perfrom the XY and ZPhi fit using the Lsqsum4 class
    void  refineFitXY  (RobustHelixFinderData& helixData, bool TargetCon, int weightMode=1);
    void  refineFitZPhi(RobustHelixFinderData& helixData, int weightMode=1);

    void  setTracker    (const Tracker*    Tracker) { _tracker     = Tracker; }
    void  setCalorimeter(const Calorimeter* Cal    ) { _calorimeter = Cal    ; }
    
    void  setRobustHelixFitter(RobustHelixFit *RHFit) { _rhfit = RHFit; }

    const Tracker*            _tracker;
    const Calorimeter*         _calorimeter;


  private:

    void fitChi2CircleMedian(RobustHelixFinderData& helixData, bool TargetCon);

    bool use(ComboHit const&) const;
    void setOutlier(ComboHit&) const;

    static float deltaPhi(float phi1, float phi2);
    bool resolvePhi(ComboHit& hh, RobustHelix const& myhel) const;
    float hitWeight(ComboHit const& hhit) const;

    int _diag;
    int _debug;
    StrawHitFlag _dontuseflag;
    int      _minnsh;  // minimum # of StrawHits
    float _minxyresid; // minimum distance used in the circle fit to be clusterized. units are mm
    float _chi2xymax, _chi2zphimax; //maximum chi2 allowed for the XY and ZPhi fit respectively
    float _mindfdz, _maxdfdz;//paramters use for findDfDz function
    float _sigmaPhi; //approximated uncertanty on the ComboHit helix-phi coordinate
    float _maxdxy; // maximum distance in hits after the triplet loop in fitCiircleMedian
    float _maxXDPhi;//maximum normalized residual for a hit in the z-phi fit
    float _rwind; // raidus window for defining points to be 'on' the helix

    RobustHelixFit* _rhfit;

  };
}
#endif
