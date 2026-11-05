//
// Filter out events with high occupancy that could disrupt DAQ processing
// Original author: Michael MacKenzie, 2026

// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "cetlib_except/exception.h"

// artdaq
#include <artdaq-core/Data/Fragment.hh>
#include <artdaq-core/Data/ContainerFragment.hh>
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"

// Offline
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"

// C++
#include <ostream>
#include <string>
#include <typeinfo>

// TRACE
#include "TRACE/tracemf.h"
#define TRACE_NAME "OccupancyFilter"

namespace mu2e {
  class OccupancyFilter : public art::EDFilter {
  public:

    // Table for a specific collection selection
    struct CutConfig {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<std::string>  tag    {Name("tag")    , Comment("Collection tag")};
      fhicl::Atom<int>          minSize{Name("minSize"), Comment("Minimum collection size (<= 0 = no minimum)"), -1};
      fhicl::Atom<int>          maxSize{Name("maxSize"), Comment("Maximum collection size (<  0 = no maximum)"), -1};
    };

    struct FragmentCutConfig { // no tag needed here
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>          minSize{Name("minSize"), Comment("Minimum number of fragments (<= 0 = no minimum)"), -1};
      fhicl::Atom<int>          maxSize{Name("maxSize"), Comment("Maximum number of fragments (<  0 = no maximum)"), -1};
    };

    // Main configuration
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::OptionalTable<CutConfig>          comboHits {Name("comboHits") , Comment("Tracker combo hits selection")};
      fhicl::OptionalTable<CutConfig>          strawHits {Name("strawHits") , Comment("Tracker straw hits selection")};
      fhicl::OptionalTable<CutConfig>          strawDigis{Name("strawDigis"), Comment("Tracker straw digis selection")};
      fhicl::OptionalTable<CutConfig>          caloHits  {Name("caloHits")  , Comment("Calorimeter hits selection")};
      fhicl::OptionalTable<CutConfig>          caloDigis {Name("caloDigis") , Comment("Calorimeter digis selection")};
      fhicl::OptionalTable<FragmentCutConfig>  fragments {Name("fragments") , Comment("Fragments selection")};
    };

    // Selection on a collection size: a negative bound is not applied
    struct Cuts {
      std::string tag = ""; // not used for fragments
      int min = -1;
      int max = -1;
      Cuts() {}
      Cuts(const CutConfig& cuts) : tag(cuts.tag()), min(cuts.minSize()), max(cuts.maxSize()) { validate(tag); }
      Cuts(const FragmentCutConfig& cuts) : min(cuts.minSize()), max(cuts.maxSize()) { validate("fragments"); }
      // A maximum, if set, must not be below the minimum
      void validate(const std::string& name) const {
        if(max >= 0 && max < min) {
          throw cet::exception("CONFIG") << "OccupancyFilter: invalid cut for " << name
                                         << ": maxSize (" << max << ") < minSize (" << min << ")\n";
        }
      }
      bool passes(const size_t n) const {
        return (min < 0 || n >= size_t(min)) && (max < 0 || n <= size_t(max));
      }
      friend std::ostream& operator<<(std::ostream& os, const Cuts& c) {
        return os << " (tag = '" << c.tag << "', min = " << c.min << ", max = " << c.max << ")";
      }
    };

    using Parameters = art::EDFilter::Table<Config>;

    explicit OccupancyFilter(const Parameters& config);


  private:
    bool filter  (art::Event& event) override;
    bool endRun  (art::Run&   run  ) override;

    // Check if a collection passes the selection
    template <typename T> bool check_collection(const art::Event& event, const Cuts& cuts);
    bool check_fragments(const art::Event& event, const Cuts& cuts);

    // Inputs
    bool                     _filterComboHits;
    bool                     _filterStrawHits;
    bool                     _filterStrawDigis;
    bool                     _filterCaloHits;
    bool                     _filterCaloDigis;
    bool                     _filterFragments;

    // Selections
    Cuts                     _comboHitCuts;
    Cuts                     _strawHitCuts;
    Cuts                     _strawDigiCuts;
    Cuts                     _caloHitCuts;
    Cuts                     _caloDigiCuts;
    Cuts                     _fragmentCuts;

    // Data
    int                      _nevt, _npass;

  };

  //-----------------------------------------------------------------------------
  OccupancyFilter::OccupancyFilter(const Parameters& conf)
    : art::EDFilter{conf}
    , _filterComboHits (conf().comboHits ())
    , _filterStrawHits (conf().strawHits ())
    , _filterStrawDigis(conf().strawDigis())
    , _filterCaloHits  (conf().caloHits  ())
    , _filterCaloDigis (conf().caloDigis ())
    , _filterFragments (conf().fragments ())
    , _nevt(0)
    , _npass(0)
  {
    // Check each collection option if it's set, getting its info if available
    if(_filterComboHits ) _comboHitCuts  = Cuts(*(conf().comboHits ()));
    if(_filterStrawHits ) _strawHitCuts  = Cuts(*(conf().strawHits ()));
    if(_filterStrawDigis) _strawDigiCuts = Cuts(*(conf().strawDigis()));
    if(_filterCaloHits  ) _caloHitCuts   = Cuts(*(conf().caloHits  ()));
    if(_filterCaloDigis ) _caloDigiCuts  = Cuts(*(conf().caloDigis ()));
    if(_filterFragments ) _fragmentCuts  = Cuts(*(conf().fragments ()));

    TLOG(TLVL_DEBUG + 1) << ":"
                         << " filter combo hits = "   << _filterComboHits  <<  _comboHitCuts
                         << " filter straw hits = "   << _filterStrawHits  <<  _strawHitCuts
                         << " filter straw digis = "  << _filterStrawDigis <<  _strawDigiCuts
                         << " filter calo hits =   "  << _filterCaloHits   <<  _caloHitCuts
                         << " filter calo digis = "   << _filterCaloDigis  <<  _caloDigiCuts
                         << " filter fragments = "    << _filterFragments  <<  _fragmentCuts ;
  }

  //-----------------------------------------------------------------------------
  // Check if a collection passes the selection
  template <typename T>
  bool OccupancyFilter::check_collection(const art::Event& event, const Cuts& cuts) {
    const std::string& tag = cuts.tag;
    // Attempt to get the collection
    art::Handle<T> handle;
    if(!event.getByLabel(tag, handle)) {
      TLOG(TLVL_WARNING) << ": Unable to find handle for " << typeid(T).name() << " with tag " << tag;
      return cuts.min <= 0; // a missing collection can't satisfy a minimum, but passes a max-only cut
    }

    // Check the collection
    const T* collection = handle.product();
    const size_t nobj = collection->size();
    const bool passed = cuts.passes(nobj);
    TLOG(TLVL_DEBUG + 3) << ":  handle " << typeid(T).name() << " with tag " << tag
                         << " has size " << nobj << " and result " << passed;
    if(!passed) TLOG(TLVL_DEBUG + 4) << ":  handle " << typeid(T).name() << " with tag " << tag
                                     << " has size " << nobj << " and fails the check (min = " << cuts.min << ", max = " << cuts.max << ")";
    return passed;
  }

  //-----------------------------------------------------------------------------
  // Check if the fragment collection passes the selection
  bool OccupancyFilter::check_fragments(const art::Event& event, const Cuts& cuts) {
    const int max_size = cuts.max;

    // Get all fragment handles
    std::vector<art::Handle<artdaq::Fragments>> fragmentHandles = event.getMany<std::vector<artdaq::Fragment>>();

    // Count all DTCEVT fragments
    size_t nfragments(0);
    for (const auto& handle : fragmentHandles) {
      if (!handle.isValid() || handle->empty()) {
        continue;
      }

      // Container type --> count fragments in each block
      if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
        for (const auto& cont : *handle) {
          artdaq::ContainerFragment contf(cont);
          if (contf.fragment_type() != mu2e::FragmentType::DTCEVT) {
            continue;
          }

          nfragments += contf.block_count();
          if(max_size >= 0 && nfragments > size_t(max_size)) {
            TLOG(TLVL_DEBUG + 10) << ": Reached maximum fragments at " << nfragments;
            return false; // stop processing if we ever fail
          }
        }
      } else { // directly a list of fragments
        if (handle->front().type() == mu2e::FragmentType::DTCEVT) {
          nfragments += handle->size();
          if(max_size >= 0 && nfragments > size_t(max_size)) {
            TLOG(TLVL_DEBUG + 10) << ": Reached maximum fragments at " << nfragments;
            return false; // stop processing if we ever fail
          }
        }
      }
    }
    const bool passed = cuts.passes(nfragments);
    TLOG(TLVL_DEBUG + 11) << ": N(fragments) = " << nfragments << " result = " << passed;

    return passed;
  }

  //-----------------------------------------------------------------------------
  bool OccupancyFilter::filter(art::Event& event) {

    // Count total events seen
    ++_nevt;

    // Filter flag
    bool passed = true;

    // Filter on each requested collection
    if(_filterComboHits ) passed &= check_collection<mu2e::ComboHitCollection >(event, _comboHitCuts);
    if(_filterStrawHits ) passed &= check_collection<mu2e::StrawHitCollection >(event, _strawHitCuts);
    if(_filterStrawDigis) passed &= check_collection<mu2e::StrawDigiCollection>(event, _strawDigiCuts);
    if(_filterCaloHits  ) passed &= check_collection<mu2e::CaloHitCollection  >(event, _caloHitCuts);
    if(_filterCaloDigis ) passed &= check_collection<mu2e::CaloDigiCollection >(event, _caloDigiCuts);
    if(_filterFragments ) passed &= check_fragments(event, _fragmentCuts);
    TLOG(TLVL_DEBUG + 5) << ": Event has status " << passed;

    // Count accepted events
    if(passed) ++_npass;
    else TLOG(TLVL_DEBUG + 6) << ": Event failed occupancy selection";

    // Return the result
    return passed;
  }

  //-----------------------------------------------------------------------------
  bool OccupancyFilter::endRun(art::Run& run) {
    // Print a summary of the filter results
    const float rate = (_nevt > 0) ? float(_npass)/float(_nevt) : 0.f;
    TLOG(TLVL_DEBUG + 2) << "passed " << _npass << " events out of " << _nevt << " for a ratio of " << rate;
    _nevt = 0; _npass = 0; // reset for the next run
    return true;
  }
}

DEFINE_ART_MODULE(mu2e::OccupancyFilter)
