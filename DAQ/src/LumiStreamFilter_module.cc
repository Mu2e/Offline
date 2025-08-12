// Collect lumi stream information and write it out at lower frequency in subruns (and potentially
// events)

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

#include <artdaq-core-mu2e/Data/EventHeader.hh>
#ifdef ARTDAQ_DAQDATA_GLOBALS_HH
#include "artdaq/DAQdata/Globals.hh"
#endif

#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"

#include <algorithm>
#include <iostream>
#include <string>

#include "trace.h"
#define TRACE_NAME "LumiSteamFilter"

namespace mu2e {

  class LumiStreamFilter : public art::EDFilter {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int>                   diagLevel     {Name("diagLevel")     , Comment("diagnostic Level"), 0};
      fhicl::OptionalAtom<art::InputTag> caloTag       {Name("caloTag")       , Comment("Calo intensity info Tag")};
      fhicl::OptionalAtom<art::InputTag> timeClusterTag{Name("timeClusterTag"), Comment("Time cluster intensity info Tag")};
      fhicl::OptionalAtom<art::InputTag> trackerTag    {Name("trackerTag")    , Comment("Tracker intensity info Tag")};
      fhicl::OptionalAtom<art::InputTag> headerTag     {Name("headerTag")     , Comment("Event header Tag")};
      fhicl::OptionalAtom<int>           eventFreq     {Name("eventFreq")     , Comment("Frequency of event data product writing, defaulting to sub-runs if not set")};
      fhicl::Atom<bool>                  passFirst     {Name("passFirst")     , Comment("Pass first event per sub-run to ensure the output file is prepared"), false};
      fhicl::Atom<int>                   duplications  {Name("duplications")  , Comment("Duplicate input data N times"), 0};
    };

    explicit LumiStreamFilter(const art::EDFilter::Table<Config>& config);
    virtual bool filter(art::Event& event) override;
    virtual bool endSubRun(art::SubRun& sr) override;

  private:
    int           _diagLevel;
    art::InputTag _caloTag;
    bool          _useCalo;
    art::InputTag _timeClusterTag;
    bool          _useTimeCluster;
    art::InputTag _trackerTag;
    bool          _useTracker;
    art::InputTag _headerTag;
    bool          _useHeader;
    int           _eventFreq;
    bool          _useSubruns;
    bool          _passFirst;
    int           _duplications;
    std::string   _outName;
    long          _eventCount;
    bool          _firstInSubrun;

    std::unique_ptr<mu2e::IntensityInfosCalo>        _caloInfos;
    std::unique_ptr<mu2e::IntensityInfosTimeCluster> _timeClusterInfos;
    std::unique_ptr<mu2e::IntensityInfosTrackerHits> _trackerInfos;
    std::unique_ptr<mu2e::EventHeaders>              _headers;
  };

  LumiStreamFilter::LumiStreamFilter(const art::EDFilter::Table<Config>& config)
    : art::EDFilter{config}
    , _diagLevel(config().diagLevel())
    , _useCalo(config().caloTag(_caloTag))
    , _useTimeCluster(config().timeClusterTag(_timeClusterTag))
    , _useTracker(config().trackerTag(_trackerTag))
    , _useHeader(config().headerTag(_headerTag))
    , _useSubruns(!config().eventFreq(_eventFreq))
    , _passFirst(config().passFirst())
    , _duplications(config().duplications())
#ifdef ARTDAQ_DAQDATA_GLOBALS_HH
    , _outName(app_name)
#else
    , _outName("")
#endif
    , _eventCount(0)
    , _firstInSubrun(true)
    , _caloInfos(nullptr)
    , _timeClusterInfos(nullptr)
    , _trackerInfos(nullptr)
    , _headers(nullptr)
  {
    std::replace(_outName.begin(), _outName.end(), '_', 'U');  // protect against forbidden branch name characters
    std::replace(_outName.begin(), _outName.end(), '-', 'D');
    if(_diagLevel > 0) {
      std::cout << "LumiStreamFilter::" << __func__
                << ": Configured with:\n eventFreq = " << _eventFreq << std::endl;
      std::cout << " useCalo = " << _useCalo;
      if(_useCalo) std::cout << " (" << _caloTag.encode().c_str() << ")";
      std::cout << std::endl;
      std::cout << " useTimeCluster = " << _useTimeCluster;
      if(_useTimeCluster) std::cout << " (" << _timeClusterTag.encode().c_str() << ")";
      std::cout << std::endl;
      std::cout << " useTracker = " << _useTracker;
      if(_useTracker) std::cout << " (" << _trackerTag.encode().c_str() << ")";
      std::cout << " useHeader = " << _useHeader;
      if(_useHeader) std::cout << " (" << _headerTag.encode().c_str() << ")";
      std::cout << std::endl;
      std::cout << " passFirst = " << _passFirst << std::endl;
      std::cout << " duplications = " << _duplications << std::endl;
      std::cout << " outName = " << _outName.c_str() << std::endl;
    }

    if(_useCalo) {
      _caloInfos = std::unique_ptr<mu2e::IntensityInfosCalo>(new mu2e::IntensityInfosCalo);
      if(!_useSubruns) produces<mu2e::IntensityInfosCalo>();
      produces<mu2e::IntensityInfosCalo, art::InSubRun>(_outName);
    }
    if(_useTimeCluster) {
      _timeClusterInfos = std::unique_ptr<mu2e::IntensityInfosTimeCluster>(new mu2e::IntensityInfosTimeCluster);
      if(!_useSubruns) produces<mu2e::IntensityInfosTimeCluster>();
      produces<mu2e::IntensityInfosTimeCluster, art::InSubRun>(_outName);
    }
    if(_useTracker) {
      _trackerInfos = std::unique_ptr<mu2e::IntensityInfosTrackerHits>(new mu2e::IntensityInfosTrackerHits);
      if(!_useSubruns) produces<mu2e::IntensityInfosTrackerHits>();
      produces<mu2e::IntensityInfosTrackerHits, art::InSubRun>(_outName);
    }
    if(_useHeader) {
      _headers = std::unique_ptr<mu2e::EventHeaders>(new mu2e::EventHeaders);
      if(!_useSubruns) produces<mu2e::EventHeaders>();
      produces<mu2e::EventHeaders, art::InSubRun>(_outName);
    }
    if(_useSubruns) _eventFreq = 0;
    else {
      if(_caloInfos)        _caloInfos       ->reserve(_eventFreq);
      if(_timeClusterInfos) _timeClusterInfos->reserve(_eventFreq);
      if(_trackerInfos)     _trackerInfos    ->reserve(_eventFreq);
    }
  }

  bool LumiStreamFilter::endSubRun(art::SubRun& sr) {
    art::SubRunNumber_t subrunNumber = sr.subRun();
    art::RunNumber_t    runNumber    = sr.run();
    if(_diagLevel > 0)
      std::cout << "LumiStreamFilter::" << __func__ << ": Run:SubRun " << runNumber
                << ":" << subrunNumber << " out_name = " << _outName
                << ": Writing out the accumulated intensity info collections from "
                << _eventCount << " events\n";
    TLOG(TLVL_DEBUG + 5)
      << "LumiStreamFilter::" << __func__ << ": Run:SubRun " << runNumber << ":"
      << subrunNumber << " out_name = " << _outName
      << ": Writing out the accumulated intensity info collections from " << _eventCount
      << " events\n";

    // Add the data to the subrun
    if(_useCalo)        sr.put(std::move(_caloInfos)       , _outName, art::fullSubRun());
    if(_useTimeCluster) sr.put(std::move(_timeClusterInfos), _outName, art::fullSubRun());
    if(_useTracker)     sr.put(std::move(_trackerInfos)    , _outName, art::fullSubRun());
    if(_useHeader)      sr.put(std::move(_headers)         , _outName, art::fullSubRun());

    // Initialize new collections
    if(_useCalo)         _caloInfos = std::unique_ptr<mu2e::IntensityInfosCalo>(new mu2e::IntensityInfosCalo);
    if(_useTimeCluster) _timeClusterInfos = std::unique_ptr<mu2e::IntensityInfosTimeCluster>(new mu2e::IntensityInfosTimeCluster);
    if(_useTracker)     _trackerInfos = std::unique_ptr<mu2e::IntensityInfosTrackerHits>(new mu2e::IntensityInfosTrackerHits);
    if(_useHeader)      _headers = std::unique_ptr<mu2e::EventHeaders>(new mu2e::EventHeaders);

    _firstInSubrun = true;
    _eventCount    = 0;
    return true;
  }

  bool LumiStreamFilter::filter(art::Event& event) {
    ++_eventCount;  // increment the event counter

    // for printout use
    const art::EventNumber_t  eventNumber  = event.event();
    const art::SubRunNumber_t subrunNumber = event.subRun();
    const art::RunNumber_t    runNumber    = event.run();

    if(_diagLevel > 2)
      std::cout << "LumiStreamFilter::" << __func__ << ": Begin processing Event "
                << runNumber << ":" << subrunNumber << ":" << eventNumber << std::endl;
    //---------------------------------------
    // Retrieve the data

    if(_useCalo) {
      art::Handle<mu2e::IntensityInfoCalo> caloH;
      if(!event.getByLabel(_caloTag, caloH) || !caloH.product()) {
        std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber
                  << ":" << subrunNumber << ":" << eventNumber
                  << ": Calo intensity object not found!\n";
        _caloInfos->push_back(mu2e::IntensityInfoCalo());  // add an empty one to maintain the list alignment
      } else {
        const auto caloInfo = caloH.product();
        _caloInfos->push_back(mu2e::IntensityInfoCalo(*caloInfo));
        if(_diagLevel > 2)
          std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber
                    << ":" << subrunNumber << ":" << eventNumber
                    << ": Retrieved calo information\n";
      }
    }

    if(_useTimeCluster) {
      art::Handle<mu2e::IntensityInfoTimeCluster> timeClusterH;
      if(!event.getByLabel(_timeClusterTag, timeClusterH) || !timeClusterH.product()) {
        std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber
                  << ":" << subrunNumber << ":" << eventNumber
                  << ": Time cluster intensity object not found!\n";
        _timeClusterInfos->push_back(mu2e::IntensityInfoTimeCluster());  // add an empty one to maintain the list alignment
      } else {
        const auto timeClusterInfo = timeClusterH.product();
        _timeClusterInfos->push_back(mu2e::IntensityInfoTimeCluster(*timeClusterInfo));
        if(_diagLevel > 2)
          std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber
                    << ":" << subrunNumber << ":" << eventNumber
                    << ": Retrieved time cluster information\n";
      }
    }

    if(_useTracker) {
      art::Handle<mu2e::IntensityInfoTrackerHits> trackerH;
      if(!event.getByLabel(_trackerTag, trackerH) || !trackerH.product()) {
        std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber
                  << ":" << subrunNumber << ":" << eventNumber
                  << ": Tracker intensity object not found!\n";
        _trackerInfos->push_back(mu2e::IntensityInfoTrackerHits());  // add an empty one to maintain the list alignment
      } else {
        const auto trackerInfo = trackerH.product();
        _trackerInfos->push_back(mu2e::IntensityInfoTrackerHits(*trackerInfo));
        if(_diagLevel > 2)
          std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber
                    << ":" << subrunNumber << ":" << eventNumber
                    << ": Retrieved tracker information\n";
      }
    }

    if(_useHeader) {
      art::Handle<mu2e::EventHeader> headerH;
      if(!event.getByLabel(_headerTag, headerH) || !headerH.product()) {
        std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber
                  << ":" << subrunNumber << ":" << eventNumber
                  << ": Event Header not found!\n";
        _headers->push_back(mu2e::EventHeader());  // add an empty one to maintain the list alignment
      } else {
        const auto header = headerH.product();
        _headers->push_back(mu2e::EventHeader(*header));
        if(_diagLevel > 2)
          std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber
                    << ":" << subrunNumber << ":" << eventNumber
                    << ": Retrieved Event Header\n";
      }
      TLOG(TLVL_DEBUG + 10) << "Event " << runNumber << ":" << subrunNumber << ":"
                            << eventNumber << "Event header information: "
                            << "EWT = " << _headers->back().ewt
                            << " mode = " << _headers->back().mode;
    }

    // Check if data duplication is requested
    if(_duplications > 0) {
      for(int icopy = 0; icopy < _duplications; ++icopy) {
        if(_useCalo) _caloInfos->push_back(_caloInfos->back());
        if(_useTimeCluster) _timeClusterInfos->push_back(_timeClusterInfos->back());
        if(_useTracker) _trackerInfos->push_back(_trackerInfos->back());
        if(_useHeader) _headers->push_back(_headers->back());
      }
    }

    //---------------------------------------
    // Write out the data if requested

    bool retval = _eventFreq > 0 && (_eventCount % _eventFreq == 0);
    if(retval) {
      if(_diagLevel > 0)
        std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber
                  << ":" << subrunNumber << ":" << eventNumber
                  << ": Writing out the accumulated intensity info collections from "
                  << _eventCount << " events\n";

      // Add the data to the event
      if(_useCalo)        event.put(std::move(_caloInfos));
      if(_useTimeCluster) event.put(std::move(_timeClusterInfos));
      if(_useTracker)     event.put(std::move(_trackerInfos));
      if(_useHeader)      event.put(std::move(_headers));

      // Initialize new collections
      if(_useCalo)        _caloInfos        = std::unique_ptr<mu2e::IntensityInfosCalo>(new mu2e::IntensityInfosCalo);
      if(_useTimeCluster) _timeClusterInfos = std::unique_ptr<mu2e::IntensityInfosTimeCluster>(new mu2e::IntensityInfosTimeCluster);
      if(_useTracker)     _trackerInfos     = std::unique_ptr<mu2e::IntensityInfosTrackerHits>(new mu2e::IntensityInfosTrackerHits);
      if(_useHeader)      _headers          = std::unique_ptr<mu2e::EventHeaders>(new mu2e::EventHeaders);
      if(_useCalo)        _caloInfos       ->reserve(_eventFreq);
      if(_useTimeCluster) _timeClusterInfos->reserve(_eventFreq);
      if(_useTracker)     _trackerInfos    ->reserve(_eventFreq);
      if(_useHeader)      _headers         ->reserve(_eventFreq);

      // Reset the counter
      _eventCount = 0;
    } else if(!_useSubruns) {  // add empty information in failed events if writing into events
      if(_diagLevel > 2)
        std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber
                  << ":" << subrunNumber << ":" << eventNumber
                  << ": Adding empty intensity information to the event\n";
      if(_useCalo) {
        auto tmp_info = std::unique_ptr<mu2e::IntensityInfosCalo>(new mu2e::IntensityInfosCalo);
        event.put(std::move(tmp_info));
      }
      if(_useTimeCluster) {
        auto tmp_info = std::unique_ptr<mu2e::IntensityInfosTimeCluster>(new mu2e::IntensityInfosTimeCluster);
        event.put(std::move(tmp_info));
      }
      if(_useTracker) {
        auto tmp_info = std::unique_ptr<mu2e::IntensityInfosTrackerHits>(new mu2e::IntensityInfosTrackerHits);
        event.put(std::move(tmp_info));
      }
      if(_useHeader) {
        auto tmp_info = std::unique_ptr<mu2e::EventHeaders>(new mu2e::EventHeaders);
        event.put(std::move(tmp_info));
      }
    }

    // Force the output file to be created by passing the first event seen in a subrun, still passing
    // empty containers if it's not a selected event
    retval |= _passFirst && _firstInSubrun;
    _firstInSubrun = false;

    if(_diagLevel > 1)
      std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber << ":"
                << subrunNumber << ":" << eventNumber << ": Return value = " << retval
                << " for event count = " << _eventCount << " events and event freq "
                << _eventFreq << std::endl;
    return retval;
  }  // filter

}  // namespace mu2e

using mu2e::LumiStreamFilter;
DEFINE_ART_MODULE(LumiStreamFilter)
