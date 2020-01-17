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

    bool writeEventIDs_;
    art::EventIDSequence idseq_;

  public:

    struct Mu2eConfig {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<Mu2eProductMixer::Config> products { Name("products"),
          Comment("A table specifying products to be mixed.  For each supported data type\n"
                  "there is a mixingMap sequence that defines mapping of inputs to outputs.\n"
                  "Each entry in the top-level mixingMap sequence is a sequence of two strings:\n"
                  "    [ \"InputTag\", \"outputInstanceName\" ]\n"
                  "The output instance name colon \":\" is special: it means take instance name from the input tag.\n"
                  "For example, with this config:\n"
                  "   mixingMap: [ [ \"detectorFilter:tracker\", \"tracker\" ], [ \"detectorFilter:virtualdetector\", \":\" ] ]\n"
                  "the outputs will be named \"tracker\" and \"virtualdetector\"\n"
                  )
          };

      fhicl::Atom<unsigned> nSecondaries { Name("nSecondaries"),
          Comment("Number of secondary events per single primary, a positive integer."),
          1u };

      fhicl::Atom<bool> writeEventIDs { Name("writeEventIDs"),
          Comment("Write out IDs of events on the secondary input stream."),
          true
          };
    };

    struct Config {
      fhicl::Table<Mu2eConfig> mu2e { fhicl::Name("mu2e") };
    };

    using Parameters = art::MixFilterTable<Config>;
    ResamplingMixerDetail(const Parameters& pset, art::MixHelper &helper);

    size_t nSecondaries() const { return nSecondaries_; }

    void processEventIDs(const art::EventIDSequence& seq);

    void finalizeEvent(art::Event& e);
  };

  ResamplingMixerDetail::ResamplingMixerDetail(const Parameters& pars, art::MixHelper& helper)
    : spm_{ pars().mu2e().products(), helper }
    , nSecondaries_{ pars().mu2e().nSecondaries() }
    , writeEventIDs_{ pars().mu2e().writeEventIDs() }
  {
    if(writeEventIDs_) {
      helper.produces<art::EventIDSequence>();
    }
  }

  void ResamplingMixerDetail::processEventIDs(const art::EventIDSequence& seq) {
    if(writeEventIDs_) {
      idseq_ = seq;
    }
  }

  void ResamplingMixerDetail::finalizeEvent(art::Event& e) {
    if(writeEventIDs_) {
      auto o = std::make_unique<art::EventIDSequence>();
      o->swap(idseq_);
      e.put(std::move(o));
    }
  }

}

//================================================================
namespace mu2e {
  // This is the module class.
  typedef art::MixFilter<ResamplingMixerDetail,art::RootIOPolicy> ResamplingMixer;

}

DEFINE_ART_MODULE(mu2e::ResamplingMixer);
