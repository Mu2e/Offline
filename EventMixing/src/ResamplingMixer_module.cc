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

  //----------------------------------------------------------------
  // Our "detail" class for art/Framework/Modules/MixFilter.h
  class ResamplingMixerDetail {
    Mu2eProductMixer spm_;
    const unsigned nSecondaries_;
  public:

    struct Mu2eConfig {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<Mu2eProductMixer::Config> products { Name("products") };

      fhicl::Atom<unsigned> nSecondaries { Name("nSecondaries"),
          Comment("Number of secondary events per single primary, a positive integer."),
          1u };
    };

    struct Config {
      fhicl::Table<Mu2eConfig> mu2e { fhicl::Name("mu2e") };
    };

    using Parameters = art::MixFilterTable<Config>;
    ResamplingMixerDetail(const Parameters& pset, art::MixHelper &helper);
    size_t nSecondaries() const { return nSecondaries_; }
  };

  ResamplingMixerDetail::ResamplingMixerDetail(const Parameters& pars, art::MixHelper& helper)
    : spm_{ pars().mu2e().products(), helper }
    , nSecondaries_{ pars().mu2e().nSecondaries() }
  {}

}

//================================================================
namespace mu2e {
  // This is the module class.
  typedef art::MixFilter<ResamplingMixerDetail,art::RootIOPolicy> ResamplingMixer;

}

DEFINE_ART_MODULE(mu2e::ResamplingMixer);
