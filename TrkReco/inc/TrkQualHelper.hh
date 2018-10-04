//
// a class to calculate all the variables that go into TrkQual
// 
// FIXME! Need to find the correct package to put this in
// -- it is used by TrkDiag/src/KalDiag.cc and TrkPatRec/src/KalFinalFit_module.cc
#ifndef TrkQualHelper_HH
#define TrkQualHelper_HH
#include "RecoDataProducts/inc/TrkQual.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "Rtypes.h"
namespace mu2e
{
// general information about a track
  class TrkQualHelper {

  private:
    MVATools* _trkqualmva;

  public:
    TrkQualHelper(const fhicl::ParameterSet& pset)
    {
      _trkqualmva = new MVATools(pset.get<fhicl::ParameterSet>("TrkQualMVA", fhicl::ParameterSet()));
      _trkqualmva->initMVA();
    }
    void fillTrkQual(const KalRep* krep, TrkQual& trkqual);
    void fillTrkQual(const KalSeed& kseed, TrkQual& trkqual);

    const MVATools* MVATool() { return _trkqualmva; }

  private:
    inline bool checkStatus(const int& status);
    inline bool checkStatus(const TrkFitFlag& status);

    void findBestKSeg(const KalRep* krep, KalSegment& kseg);
    bool findBestKSeg(const KalSeed& kseed, KalSegment& kseg);

    void fillHitCountVariables(const std::vector<TrkStrawHitSeed>& hits, const std::vector<TrkStraw>& straws, TrkQual& trkqual);
    void fillTrkVariables(const double& fit_con, const TrkT0& trkt0, const KalSegment& bestkseg, TrkQual& trkqual);
  };
}
#endif

