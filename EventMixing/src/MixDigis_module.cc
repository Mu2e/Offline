// Mix digi-level data products, as opposed to simulation-level data products,
// into art events
//
// Ed Callaghan, 2023

// stl

// art
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art_root_io/RootIOPolicy.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TupleAs.h"

// cetlib_except
#include "cetlib_except/exception.h"

// mu2e
#include "Offline/EventMixing/inc/Mu2eProductMixer.hh"

namespace mu2e{
  // detail class for art::MixFilter<>
  class MixDigisDetail{
    public:
      struct Mu2eConfig{
        // configuration is simple --- which data products to mix in
        fhicl::Table<Mu2eProductMixer::Config> products{
          fhicl::Name("products"),
          fhicl::Comment("Products to be mixed, as form of a mixingMap, i.e. tuples of InputTags to output instance names.")

        };
      };

      struct Config{
        fhicl::Table<Mu2eConfig> mu2e{ fhicl::Name("mu2e") };
      };

      using Parameters = art::MixFilterTable<Config>;
      explicit MixDigisDetail(const Parameters& parameters,
                                  art::MixHelper& helper);
      size_t nSecondaries();

      void processEventIDs(const art::EventIDSequence& seq);
      void beginSubRun(const art::SubRun& subrun);
      void endSubRun(art::SubRun& subrun);
      void startEvent(const art::Event& event);
      void finalizeEvent(art::Event& event);

    protected:
      /**/

    private:
      Mu2eProductMixer _mixer;
      art::EventIDSequence evids;
  };

  // implementation
  MixDigisDetail::MixDigisDetail(const Parameters& parameters,
                                         art::MixHelper& helper)
      : _mixer{parameters().mu2e().products(), helper}{
    helper.produces<art::EventIDSequence>();
  }

  // always overlay simulated event onto one event window
  size_t MixDigisDetail::nSecondaries(){
    size_t rv = 1;
    return rv;
  }

  // always keep track of secondary EventIDs
  void MixDigisDetail::processEventIDs(const art::EventIDSequence& seq){
    // intent is to overlay signals onto a nominally-complete event,
    // so we should only mix in one event
    if (seq.size() != 1){
      throw cet::exception("MIX") << "mu2e::DixDigisDetail: mixing more than one digi frame" << std::endl;
    }
    _mixer.processEventIDs(seq);
    evids = seq;
  }

  // forward beginSubRun
  void MixDigisDetail::beginSubRun(const art::SubRun& subrun){
    _mixer.beginSubRun(subrun);
  }

  // forward endSubRun
  void MixDigisDetail::endSubRun(art::SubRun& subrun){
    _mixer.endSubRun(subrun);
  }

  // forward startEvent
  void MixDigisDetail::startEvent(const art::Event& event){
    _mixer.startEvent(event);
  }

  void MixDigisDetail::finalizeEvent(art::Event& event){
    // "manually" mix in EventIDs, as these do not preexist as data products
    auto seq = std::make_unique<art::EventIDSequence>();
    seq->swap(evids);
    event.put(std::move(seq));
  }

  // define the module
  typedef art::MixFilter<MixDigisDetail,art::RootIOPolicy> MixDigis;
}

DEFINE_ART_MODULE(mu2e::MixDigis);
