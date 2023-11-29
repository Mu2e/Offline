//
// Compton electron finder algorithm
//
//
#ifndef CalPatRec_PhiZSeedFinderAlg_hh
#define CalPatRec_PhiZSeedFinderAlg_hh

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"

#include "Offline/CalPatRec/inc/PhiZSeedFinder_types.hh"

namespace fhicl {
  class ParameterSet;
}

namespace mu2e {
  class Calorimeter;
  class Tracker;

  using PhiZSeedFinderTypes::Data_t;

  class PhiZSeedFinderAlg {
  public:
//-----------------------------------------------------------------------------
// algorithm talk-to parameters
//-----------------------------------------------------------------------------
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>             debugLevel        {Name("debugLevel"        ), Comment("debug level"                 ) };
      fhicl::Atom<int>             diagLevel         {Name("diagLevel"         ), Comment("diag level"                  ) };
      fhicl::Atom<int>             printErrors       {Name("printErrors"       ), Comment("print errors"                ) };
      fhicl::Atom<int>             testOrder         {Name("testOrder"         ), Comment("if 1, test order"            ) };
    };
  public:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    Data_t*         _data;

    int             _debugLevel;
    int             _diagLevel;
    int             _printErrors;
    int             _testOrder;
    float           _timeBin;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    PhiZSeedFinderAlg() {}

    explicit     PhiZSeedFinderAlg(const fhicl::Table<PhiZSeedFinderAlg::Config>& config, Data_t* _data);

    int          orderHits           (const TimeCluster* Tc);
    void         run                 (const TimeCluster* Tc);
  };
}
#endif
