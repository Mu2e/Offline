// This module "mixes" requested data products from a secondary input
// file into the current event.
//
// Andrei Gaponenko, 2018

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art_root_io/RootIOPolicy.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TupleAs.h"
#include "canvas/Utilities/InputTag.h"

#include "EventMixing/inc/Mu2eProductMixer.hh"

//================================================================
namespace mu2e {

  namespace {
    using namespace fhicl;
    struct MyTopConfig {
      Table<Mu2eProductMixer::Config> products { Name("products") };

      Atom<unsigned> nSecondaries { Name("nSecondaries"),
          Comment("Number of secondary events per single primary, a positive integer."),
          1u };
    };

    // The following hack will hopefully go away after
    // https://cdcvs.fnal.gov/redmine/issues/19970
    // is resolved.
    MyTopConfig
    retrieveConfiguration(const std::string& subTableName, const fhicl::ParameterSet& pset)
    {
      std::set<std::string> ignorable_keys;

      // Ignore everything but the subtable
      const auto& allnames = pset.get_names();
      for(const auto& i: allnames) {
        if(i != subTableName) {
          ignorable_keys.insert(i);
        }
      }

      return fhicl::Table<MyTopConfig>(pset.get<fhicl::ParameterSet>(subTableName),
                                       ignorable_keys )();
    }
  }

  //----------------------------------------------------------------
  // Our "detail" class for art/Framework/Modules/MixFilter.h
  class ResamplingMixerDetail {
    Mu2eProductMixer spm_;
    const unsigned nSecondaries_;
  public:
    ResamplingMixerDetail(const fhicl::ParameterSet& pset, art::MixHelper &helper);
    size_t nSecondaries() const { return nSecondaries_; }
  };

  ResamplingMixerDetail::ResamplingMixerDetail(const fhicl::ParameterSet& pset, art::MixHelper& helper)
    : spm_{ retrieveConfiguration("mu2e", pset).products(), helper }
    , nSecondaries_{ retrieveConfiguration("mu2e", pset).nSecondaries() }
  {}

}

//================================================================
namespace mu2e {
  // This is the module class.
  typedef art::MixFilter<ResamplingMixerDetail,art::RootIOPolicy> ResamplingMixer;

}

DEFINE_ART_MODULE(mu2e::ResamplingMixer);
