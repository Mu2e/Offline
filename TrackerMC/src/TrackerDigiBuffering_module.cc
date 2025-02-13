// Ed Callaghan
// Implement per-panel event-level fifo
// February 2025

// stl
#include <map>
#include <string>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

// cetlib_except
#include "cetlib_except/exception.h"

// mu2e
#include "Offline/Mu2eUtilities/inc/StrawDigiBundle.hh"
#include "Offline/Mu2eUtilities/inc/StrawDigiBundleCollection.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
//#include "Offline/TrackerConditions/inc/PanelElectronics.hh"

namespace mu2e{
  class TrackerDigiBuffering: public art::EDProducer{
    public:
      struct Config{
        fhicl::Atom<art::InputTag> digi_tag{
          fhicl::Name("StrawDigiCollection"),
          fhicl::Comment("StrawDigiCollection to buffer")
        };
        fhicl::Atom<art::InputTag> adcs_tag{
          fhicl::Name("StrawDigiADCWaveformCollection"),
          fhicl::Comment("StrawDigiADCWaveformCollection to buffer")
        };
        fhicl::OptionalAtom<art::InputTag> dgmc_tag{
          fhicl::Name("StrawDigiMCCollection"),
          fhicl::Comment("StrawDigiMCCollection to buffer")
        };
        // TODO rm in favor of PanelElectronics in proditions
        fhicl::Atom<unsigned int> fifo_size{
          fhicl::Name("fifo_size"),
          fhicl::Comment("Length of per-panel fifo")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      TrackerDigiBuffering(const Parameters&);
     ~TrackerDigiBuffering() = default;

      void produce(art::Event&) override;

    protected:
      art::InputTag _digi_tag;
      art::InputTag _adcs_tag;
      std::optional<art::InputTag> _dgmc_tag;
      // TODO rm in favor of PanelElectronics in proditions
      unsigned int _fifo_size;

      /*
      void sort_into(const StrawDigiBundleCollection&,
                           StrawDigiBundleCollection&);
      */

    private:
      /**/
  };

  TrackerDigiBuffering::TrackerDigiBuffering(const Parameters& config):
      art::EDProducer(config),
      _digi_tag(config().digi_tag()),
      _adcs_tag(config().adcs_tag()),
      _dgmc_tag(config().dgmc_tag()),
      // TODO rm in favor of PanelElectronics in proditions
      _fifo_size(config().fifo_size()){
    // framework hooks
    this->consumes<StrawDigiCollection>(_digi_tag);
    this->produces<StrawDigiCollection>();
    this->consumes<StrawDigiADCWaveformCollection>(_adcs_tag);
    this->produces<StrawDigiADCWaveformCollection>();

    if (_dgmc_tag){
      this->consumes<StrawDigiMCCollection>(_dgmc_tag.value());
      this->produces<StrawDigiMCCollection>();
    }
  }

  // auxiliary comparison to facilitate time-sorting
  bool compare_tdcs(const StrawDigiBundle& lhs, const StrawDigiBundle& rhs){
    bool rv = (lhs.GetStrawDigi().TDC() < rhs.GetStrawDigi().TDC());
    return rv;
  }

  void TrackerDigiBuffering::produce(art::Event& event){
    // TODO rm in favor of PanelElectronics in proditions
    unsigned int fifo_size = _fifo_size;

    // fetch digis and waveforms, and apply sanity check
    auto digi_h = event.getValidHandle<StrawDigiCollection>(_digi_tag);
    auto adcs_h = event.getValidHandle<StrawDigiADCWaveformCollection>(_digi_tag);
    if (digi_h->size() != adcs_h->size()){
      std::string msg = "mismatched digi and waveform collections, lengths";
      msg += " " + std::to_string(digi_h->size());
      msg += " !=";
      msg += " " + std::to_string(adcs_h->size());
      throw cet::exception("TrackerDigiBuffering") << msg << std::endl;
    }

    // optionally fetch mc truth, and sanity check
    art::Handle<StrawDigiMCCollection> dgmc_h;
    if (_dgmc_tag){
      dgmc_h = event.getHandle<StrawDigiMCCollection>(_dgmc_tag.value());
      if (digi_h->size() != dgmc_h->size()){
        std::string msg = "mismatched digi and waveform collections, lengths";
        msg += " " + std::to_string(digi_h->size());
        msg += " !=";
        msg += " " + std::to_string(adcs_h->size());
        throw cet::exception("TrackerDigiBuffering") << msg << std::endl;
      }
    }

    // construct flat list of triplets, optionally including mc truth
    StrawDigiBundleCollection unbuffered;
    if (dgmc_h.isValid()){
      unbuffered.Append(*digi_h, *adcs_h, *dgmc_h);
    }
    else{
      unbuffered.Append(*digi_h, *adcs_h);
    }

    // construct aggregates to be truncated, in two stages:
    std::map<StrawId,StrawDigiBundleCollection> aggregates;
    // first, partition based on panel
    for (const auto& bundle: unbuffered){
      const auto& digi = bundle.GetStrawDigi();
      StrawId panel = digi.strawId().getPanelId();
      if (aggregates.count(panel) < 1){
        aggregates[panel] = StrawDigiBundleCollection();
      }
      aggregates[panel].Append(bundle);
    }
    // then, order each partition based on tdc
    for (auto& pair: aggregates){
      //const auto& key = pair.first;
      auto& unsorted = pair.second;
      //StrawDigiBundleCollection sorted;
      //this->sort_into(unsorted, sorted);
      //aggregates[key] = sorted;
      std::sort(unsorted.begin(), unsorted.end(), compare_tdcs);
    }

    // now, truncate each partition
    for (const auto& pair: aggregates){
      const auto& key = pair.first;
      const auto& overfilled = pair.second;
      if (fifo_size < overfilled.size()){
        StrawDigiBundleCollection truncated;
        for (size_t i = 0 ; i < fifo_size ; i++){
          truncated.Append(overfilled[i]);
        }
        aggregates[key] = truncated;
      }
    }

    // next, flatten the partitions into a single collection
    StrawDigiBundleCollection buffered;
    for (const auto& pair: aggregates){
      buffered += pair.second;
    }

    // finally, sort the flattened collection by tdc
    //StrawDigiBundleCollection buffered;
    //this->sort_into(flattened, buffered);
    std::sort(buffered.begin(), buffered.end(), compare_tdcs);

    // store in event
    auto digi_buffered = buffered.GetStrawDigiPtrs();
    event.put(std::move(digi_buffered));
    auto adcs_buffered = buffered.GetStrawDigiADCWaveformPtrs();
    event.put(std::move(adcs_buffered));
    if (_dgmc_tag){
      auto dgmc_buffered = buffered.GetStrawDigiMCPtrs();
      event.put(std::move(dgmc_buffered));
    }
  }

  /*
  // indirect sort, shamefully using bare pointers
  void TrackerDigiBuffering::sort_into(const StrawDigiBundleCollection& in,
                                             StrawDigiBundleCollection& out){
    // recast as ptrs
    std::vector<StrawDigiBundle*> sortable(in.size());
    for (size_t i = 0 ; i < in.size() ; i++){
      sortable[i] = const_cast<StrawDigiBundle*>(&in[i]);
    }

    // in-place sort
    //std::sort(sortable.begin(), sortable.end(), compare_tdcs);

    // recast as objects
    for (const auto& ptr: sortable){
      out.Append(*ptr);
    }
  }
  */
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TrackerDigiBuffering)
