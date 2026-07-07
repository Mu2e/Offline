//
// Functor to calculate a track t0 from a TimeCluster.  This takes into account the specific offsets of each subsystem
//

#ifndef TrkReco_TrkTimeCalculator_HH
#define TrkReco_TrkTimeCalculator_HH

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"


namespace mu2e
{
  class TrkTimeCalculator {

    public:

      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<int>     debug{          Name("debugLevel"),          Comment("Debug flag"),0 };
        fhicl::Atom<double>  avgDriftTime{   Name("AverageDriftTime"),    Comment("Average time offset for straw hits"),22.5 };
        fhicl::Atom<bool>    useTOTdrift{    Name("UseTOTDrift"),         Comment("Switch to use the time drift"),true };
        fhicl::Atom<double>  beta{           Name("StrawHitBeta"),        Comment("Beta of the particle-hypothesis used"),1.0 };
        fhicl::Atom<double>  shErr{          Name("StrawHitTimeErr"),     Comment("Effective hit time res. without TOT"),9.7 };
        fhicl::Atom<double>  caloTimeOffset{ Name("TrkToCaloTimeOffset"), Comment("Time offsets for downstream particles in the calorimeter"),-0.4 };
        fhicl::Atom<double>  caloZOffset{    Name("CaloShowermaxZ"),      Comment("Location of shower max WRT COG Z position"),-120.0 };
        fhicl::Atom<double>  caloTimeErr{    Name("CaloTimeErr"),         Comment("Time resolution for calorimeter"),0.5 };
      };

      // parameter set must be passed in on construction
      explicit TrkTimeCalculator(fhicl::ParameterSet const&);
      explicit TrkTimeCalculator(const Config& conf);

      virtual ~TrkTimeCalculator() {};

      // update the t0 value inside a TimeCluster
      void updateT0(TimeCluster& tc, StrawHitCollection const& shcol);
      // same, taking a HelixSeed.  This uses the pitch to make a more sophisticated z correction
      void updateT0(HelixSeed& hs, StrawHitCollection const& shcol);

      double timeOfFlightTimeOffset(double hitz,double pitch) const; // z position in the tracker frame!!!
      double strawHitTimeErr() const { return _shErr; }
      double trkToCaloTimeOffset() const { return _caloTimeOffset; }
      double caloClusterTimeErr() const { return _caloTimeErr; }
      // same for a ComboHit
      double comboHitTime(ComboHit const& ch,double pitch);
      // calculate the t0 for a calo cluster.
      double caloClusterTime(CaloCluster const& cc,double pitch) const;


    private:

      int _debug;
      //StrawHitFlag _useflag, _dontuseflag;
      double _avgDriftTime;
      bool   _useTOTdrift;
      double _beta;
      double _shErr;
      double _caloZOffset;
      double _caloTimeOffset;
      double _caloTimeErr;
  };
}
#endif
