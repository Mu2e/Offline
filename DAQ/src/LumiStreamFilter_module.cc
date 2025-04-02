// Collect lumi stream information and write it out at lower frequency in subruns (and potentially events)

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"

#include <iostream>

namespace mu2e
{

  class LumiStreamFilter : public art::EDFilter
  {
  public:
    struct Config
    {
      fhicl::Atom<int>                   diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level"), 0};
      fhicl::OptionalAtom<art::InputTag> caloTag{fhicl::Name("caloTag"), fhicl::Comment("Calo intensity info Tag")};
      fhicl::OptionalAtom<art::InputTag> timeClusterTag{fhicl::Name("timeClusterTag"), fhicl::Comment("Time cluster intensity info Tag")};
      fhicl::OptionalAtom<art::InputTag> trackerTag{fhicl::Name("trackerTag"), fhicl::Comment("Tracker intensity info Tag")};
      fhicl::OptionalAtom<int>           eventFreq{fhicl::Name("eventFreq"), fhicl::Comment("Frequency of event data product writing, defaulting to sub-runs if not set")};
    };

    explicit LumiStreamFilter(const art::EDFilter::Table<Config>& config);
    virtual bool filter(art::Event& event) override;
    virtual bool endSubRun(art::SubRun& sr ) override;

  private:
    int                                      _diagLevel;
    art::InputTag                            _caloTag;
    bool                                     _useCalo;
    art::InputTag                            _timeClusterTag;
    bool                                     _useTimeCluster;
    art::InputTag                            _trackerTag;
    bool                                     _useTracker;
    int                                      _eventFreq;
    bool                                     _useSubruns;
    long                                     _eventCount;

    std::unique_ptr<mu2e::IntensityInfosCalo>         _caloInfos;
    std::unique_ptr<mu2e::IntensityInfosTimeCluster>  _timeClusterInfos;
    std::unique_ptr<mu2e::IntensityInfosTrackerHits>  _trackerInfos;

  };

  LumiStreamFilter::LumiStreamFilter(const art::EDFilter::Table<Config>& config) :
    art::EDFilter{config}
    , _diagLevel(config().diagLevel())
    , _useCalo(config().caloTag(_caloTag))
    , _useTimeCluster(config().timeClusterTag(_timeClusterTag))
    , _useTracker(config().trackerTag(_trackerTag))
    , _useSubruns(config().eventFreq(_eventFreq))
    , _eventCount(0)
    , _caloInfos(nullptr)
    , _timeClusterInfos(nullptr)
    , _trackerInfos(nullptr)
  {
    if(_diagLevel > 0) {
      std::cout << "LumiStreamFilter::" << __func__ << ": Configured with:\n eventFreq = " << _eventFreq << std::endl;
      std::cout << " useCalo = " << _useCalo;
      if(_useCalo) std::cout << " (" << _caloTag.encode().c_str() << ")";
      std::cout << std::endl;
      std::cout << " useTimeCluster = " << _useTimeCluster;
      if(_useTimeCluster) std::cout << " (" << _timeClusterTag.encode().c_str() << ")";
      std::cout << std::endl;
      std::cout << " useTracker = " << _useTracker;
      if(_useTracker) std::cout << " (" << _trackerTag.encode().c_str() << ")";
      std::cout << std::endl;
    }			 

    if(_useCalo) {
      _caloInfos = std::unique_ptr<mu2e::IntensityInfosCalo>(new mu2e::IntensityInfosCalo);
      produces<mu2e::IntensityInfosCalo>();
      produces<mu2e::IntensityInfosCalo, art::InSubRun>();
    }
    if(_useTimeCluster) {
      _timeClusterInfos = std::unique_ptr<mu2e::IntensityInfosTimeCluster>(new mu2e::IntensityInfosTimeCluster);
      produces<mu2e::IntensityInfosTimeCluster>();
      produces<mu2e::IntensityInfosTimeCluster, art::InSubRun>();
    }
    if(_useTracker) {
      _trackerInfos = std::unique_ptr<mu2e::IntensityInfosTrackerHits>(new mu2e::IntensityInfosTrackerHits);
      produces<mu2e::IntensityInfosTrackerHits>();
      produces<mu2e::IntensityInfosTrackerHits, art::InSubRun>();
    }
    if(_useSubruns) _eventFreq = 0;
    else {
      if(_caloInfos       ) _caloInfos       ->reserve(_eventFreq);
      if(_timeClusterInfos) _timeClusterInfos->reserve(_eventFreq);
      if(_trackerInfos    ) _trackerInfos    ->reserve(_eventFreq);
    }
  }

  bool LumiStreamFilter::endSubRun(art::SubRun& sr ) {
    art::SubRunNumber_t subrunNumber = sr.subRun();
    art::RunNumber_t    runNumber    = sr.run   ();
    if(_diagLevel > 0) std::cout << "LumiStreamFilter::" << __func__ << ": Subrun " << runNumber << ":" << subrunNumber
				  << ": Writing out the accumulated intensity info collections from " << _eventCount << " events\n";

    // Add the data to the subrun
    if(_caloInfos       ) sr.put(std::move(_caloInfos       ), art::fullSubRun());
    if(_timeClusterInfos) sr.put(std::move(_timeClusterInfos), art::fullSubRun());
    if(_trackerInfos    ) sr.put(std::move(_trackerInfos    ), art::fullSubRun());

    // Initialize new collections
    if(_caloInfos       ) _caloInfos        = std::unique_ptr<mu2e::IntensityInfosCalo       >(new mu2e::IntensityInfosCalo       );
    if(_timeClusterInfos) _timeClusterInfos = std::unique_ptr<mu2e::IntensityInfosTimeCluster>(new mu2e::IntensityInfosTimeCluster);
    if(_trackerInfos    ) _trackerInfos     = std::unique_ptr<mu2e::IntensityInfosTrackerHits>(new mu2e::IntensityInfosTrackerHits);

    return true;
  }

  bool LumiStreamFilter::filter(art::Event& event)
  {
    ++_eventCount; // increment the event counter

    // for printout use
    const art::EventNumber_t  eventNumber  = event.event ();
    const art::SubRunNumber_t subrunNumber = event.subRun();
    const art::RunNumber_t    runNumber    = event.run   ();

    if(_diagLevel > 2) std::cout << "LumiStreamFilter::" << __func__ << ": Begin processing Event " << runNumber << ":" << subrunNumber << ":" << eventNumber
				 << std::endl;
    //---------------------------------------
    // Retrieve the data

    if(_useCalo) {
      art::Handle<mu2e::IntensityInfoCalo> caloH;
      if(!event.getByLabel(_caloTag, caloH) || !caloH.product()) {
	std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber << ":" << subrunNumber << ":" << eventNumber
		  << ": Calo intensity object not found!\n";
	_caloInfos->push_back(mu2e::IntensityInfoCalo()); // add an empty one to maintain the list alignment
      } else {
	const auto caloInfo = caloH.product();
	_caloInfos->push_back(mu2e::IntensityInfoCalo(*caloInfo));
	if(_diagLevel > 2) std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber << ":" << subrunNumber << ":" << eventNumber
				     << ": Retrieved calo information\n";
      }
    }

    if(_useTimeCluster) {
      art::Handle<mu2e::IntensityInfoTimeCluster> timeClusterH;
      if(!event.getByLabel(_timeClusterTag, timeClusterH) || !timeClusterH.product()) {
	std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber << ":" << subrunNumber << ":" << eventNumber
		  << ": Time cluster intensity object not found!\n";
	_timeClusterInfos->push_back(mu2e::IntensityInfoTimeCluster()); // add an empty one to maintain the list alignment
      } else {
	const auto timeClusterInfo = timeClusterH.product();
	_timeClusterInfos->push_back(mu2e::IntensityInfoTimeCluster(*timeClusterInfo));
	if(_diagLevel > 2) std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber << ":" << subrunNumber << ":" << eventNumber
				     << ": Retrieved time cluster information\n";
      }
    }

    if(_useTracker) {
      art::Handle<mu2e::IntensityInfoTrackerHits> trackerH;
      if(!event.getByLabel(_trackerTag, trackerH) || !trackerH.product()) {
	std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber << ":" << subrunNumber << ":" << eventNumber
		  << ": Tracker intensity object not found!\n";
	_trackerInfos->push_back(mu2e::IntensityInfoTrackerHits()); // add an empty one to maintain the list alignment
      } else {
	const auto trackerInfo = trackerH.product();
	_trackerInfos->push_back(mu2e::IntensityInfoTrackerHits(*trackerInfo));
	if(_diagLevel > 2) std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber << ":" << subrunNumber << ":" << eventNumber
				     << ": Retrieved tracker information\n";
      }
    }

    //---------------------------------------
    // Write out the data if requested

    const bool retval = _eventFreq > 0 && (_eventCount % _eventFreq == 0);
    if(retval) {
      if(_diagLevel > 0) std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber << ":" << subrunNumber << ":" << eventNumber
				   << ": Writing out the accumulated intensity info collections from " << _eventCount << " events\n";

      // Add the data to the event
      if(_useCalo       ) event.put(std::move(_caloInfos       ));
      if(_useTimeCluster) event.put(std::move(_timeClusterInfos));
      if(_useTracker    ) event.put(std::move(_trackerInfos    ));

      // Initialize new collections
      if(_useCalo       ) _caloInfos        = std::unique_ptr<mu2e::IntensityInfosCalo       >(new mu2e::IntensityInfosCalo       );
      if(_useTimeCluster) _timeClusterInfos = std::unique_ptr<mu2e::IntensityInfosTimeCluster>(new mu2e::IntensityInfosTimeCluster);
      if(_useTracker    ) _trackerInfos     = std::unique_ptr<mu2e::IntensityInfosTrackerHits>(new mu2e::IntensityInfosTrackerHits);
      if(!_useSubruns && _useCalo       ) _caloInfos       ->reserve(_eventFreq);
      if(!_useSubruns && _useTimeCluster) _timeClusterInfos->reserve(_eventFreq);
      if(!_useSubruns && _useTracker    ) _trackerInfos    ->reserve(_eventFreq);

      // Reset the counter
      _eventCount = 0;
    } else { // add empty information in failed events
      if(_diagLevel > 2) std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber << ":" << subrunNumber << ":" << eventNumber
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
    }

    if(_diagLevel > 1) std::cout << "LumiStreamFilter::" << __func__ << ": Event " << runNumber << ":" << subrunNumber << ":" << eventNumber
				 << ": Return value = " << retval << " for event count = " << _eventCount << " events\n";
    return retval;
  } //filter

} //namespace mu2e

using mu2e::LumiStreamFilter;
DEFINE_ART_MODULE(LumiStreamFilter)
