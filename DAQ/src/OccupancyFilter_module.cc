//
// Filter out events with high occupancy that could disrupt DAQ processing
// Original author: Michael MacKenzie, 2026

// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalTable.h"

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
#include <iostream>
#include <memory>
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
      fhicl::Atom<int>          maxSize{Name("maxSize"), Comment("Maximum collection size")};
    };

    // Main configuration
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::OptionalTable<CutConfig>  comboHits {Name("comboHits") , Comment("Tracker combo hits selection")};
      fhicl::OptionalTable<CutConfig>  strawHits {Name("strawHits") , Comment("Tracker straw hits selection")};
      fhicl::OptionalTable<CutConfig>  strawDigis{Name("strawDigis"), Comment("Tracker straw digis selection")};
      fhicl::OptionalTable<CutConfig>  caloHits  {Name("caloHits")  , Comment("Calorimeter hits selection")};
      fhicl::OptionalTable<CutConfig>  caloDigis {Name("caloDigis") , Comment("Calorimeter digis selection")};
      fhicl::OptionalTable<CutConfig>  fragments {Name("fragments") , Comment("Fragments selection")};
      fhicl::Atom<int>                 diagLevel {Name("diagLevel") , Comment("Diagnostic printout level"), 0};
    };

    using Parameters = art::EDFilter::Table<Config>;

    explicit OccupancyFilter(const Parameters& config);


  private:
    bool filter  (art::Event& event) override;
    bool endRun  (art::Run&   run  ) override;

    // Check if a collection passes the selection
    template <typename T> bool check_collection(const art::Event& event, const std::string& tag, const size_t max_size);
    bool check_fragments(const art::Event& event, const size_t max_size);

    // Inputs
    bool                     _filterComboHits;
    bool                     _filterStrawHits;
    bool                     _filterStrawDigis;
    bool                     _filterCaloHits;
    bool                     _filterCaloDigis;
    bool                     _filterFragments;
    int                      _diagLevel;

    int                      _maxComboHits;
    std::string              _tagComboHits;
    int                      _maxStrawHits;
    std::string              _tagStrawHits;
    int                      _maxStrawDigis;
    std::string              _tagStrawDigis;
    int                      _maxCaloHits;
    std::string              _tagCaloHits;
    int                      _maxCaloDigis;
    std::string              _tagCaloDigis;
    int                      _maxFragments;

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
    , _diagLevel(conf().diagLevel())
    , _nevt(0)
    , _npass(0)
  {
    // Check each collection option if it's set, getting its info if available
    if(_filterComboHits) {
      _maxComboHits = conf().comboHits()->maxSize();
      _tagComboHits = conf().comboHits()->tag    ();
    }
    if(_filterStrawHits) {
      _maxStrawHits = conf().strawHits()->maxSize();
      _tagStrawHits = conf().strawHits()->tag    ();
    }
    if(_filterStrawDigis) {
      _maxStrawDigis = conf().strawDigis()->maxSize();
      _tagStrawDigis = conf().strawDigis()->tag    ();
    }
    if(_filterCaloHits) {
      _maxCaloHits = conf().caloHits()->maxSize();
      _tagCaloHits = conf().caloHits()->tag    ();
    }
    if(_filterCaloDigis) {
      _maxCaloDigis = conf().caloDigis()->maxSize();
      _tagCaloDigis = conf().caloDigis()->tag    ();
    }
    if(_filterFragments) {
      _maxFragments = conf().fragments()->maxSize();
    }

    TLOG(TLVL_DEBUG + 1) << ":"
                         << " filter combo hits = "   << _filterComboHits << " max = " << _maxComboHits
                         << " filter straw hits = "   << _filterStrawHits  << " max = " << _maxStrawHits
                         << " filter straw digis = "  << _filterStrawDigis << " max = " << _maxStrawDigis
                         << " filter calo hits =   "  << _filterCaloHits   << " max = " << _maxCaloHits
                         << " filter calo digis = "   << _filterCaloDigis  << " max = " << _maxCaloDigis
                         << " filter fragments = "    << _filterFragments  << " max = " << _maxFragments;
  }

  //-----------------------------------------------------------------------------
  // Check if a collection passes the selection
  template <typename T>
  bool OccupancyFilter::check_collection(const art::Event& event, const std::string& tag, const size_t max_size) {
    // Attempt to get the collection
    art::Handle<T> handle;
    if(!event.getByLabel(tag, handle)) {
      TLOG(TLVL_WARNING) << ": Unable to find handle for " << typeid(T).name() << " with tag " << tag;
      return true; // don't filter if it's missing, as it's not high occupancy
    }

    // Check the collection
    const T* collection = handle.product();
    const size_t nobj = collection->size();
    const bool passed = nobj < max_size;
    TLOG(TLVL_DEBUG + 3) << ":  handle " << typeid(T).name() << " with tag " << tag
                         << " has size " << nobj << " and result " << passed;
    if(!passed) TLOG(TLVL_DEBUG + 4) << ":  handle " << typeid(T).name() << " with tag " << tag
                                     << " has size " << nobj << " and fails the check (max = " << max_size << ")";
    return passed;
  }

  //-----------------------------------------------------------------------------
  // Check if the fragment collection passes the selection
  bool OccupancyFilter::check_fragments(const art::Event& event, const size_t max_size) {
    if(max_size == 0) return false; // impossible to pass

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
            break;
          }

          for (size_t ii = 0; ii < contf.block_count(); ++ii) {
            nfragments += contf[ii]->size();
            if(nfragments >= max_size) {
              TLOG(TLVL_DEBUG + 10) << ": Reached maximum fragments at " << nfragments;
              return false; // stop processing if we ever fail
            }
          }
        }
      } else { // directly a list of fragments
        if (handle->front().type() == mu2e::FragmentType::DTCEVT) {
          nfragments += handle->size();
          if(nfragments >= max_size) {
            TLOG(TLVL_DEBUG + 10) << ": Reached maximum fragments at " << nfragments;
            return false; // stop processing if we ever fail
          }
        }
      }
    }
    TLOG(TLVL_DEBUG + 11) << ": Passed event, N(fragments) = " << nfragments;

    return true; // can only reach here if we don't fail the checks
  }

  //-----------------------------------------------------------------------------
  bool OccupancyFilter::filter(art::Event& event) {

    // Count total events seen
    ++_nevt;

    // Filter flag
    bool passed = true;

    // Filter on each requested collection
    if(_filterComboHits ) passed &= check_collection<mu2e::ComboHitCollection >(event, _tagComboHits , _maxComboHits );
    if(_filterStrawHits ) passed &= check_collection<mu2e::StrawHitCollection >(event, _tagStrawHits , _maxStrawHits );
    if(_filterStrawDigis) passed &= check_collection<mu2e::StrawDigiCollection>(event, _tagStrawDigis, _maxStrawDigis);
    if(_filterCaloHits  ) passed &= check_collection<mu2e::CaloHitCollection  >(event, _tagCaloHits  , _maxCaloHits  );
    if(_filterCaloDigis ) passed &= check_collection<mu2e::CaloDigiCollection >(event, _tagCaloDigis , _maxCaloDigis );
    if(_filterFragments ) passed &= check_fragments(event, _maxFragments);
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
    return true;
  }
}

DEFINE_ART_MODULE(mu2e::OccupancyFilter)
