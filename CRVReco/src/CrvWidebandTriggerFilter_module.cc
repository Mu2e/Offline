#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"

#include <string>

using namespace std;

namespace mu2e
{

  class CrvWidebandTriggerFilter : public art::EDFilter
  {
    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{Name("crvRecoPulsesModuleLabel"), Comment("module label of the input CrvRecoPulses")};
      fhicl::Sequence<int> triggerSectorTypes{Name("triggerSectorTypes"), Comment("CRV sector types used for triggers")};
      fhicl::Atom<int> PEthreshold{Name("PEthreshold"), Comment("PE threshold required at triggers")};
      fhicl::Atom<double> maxTimeDifference{Name("maxTimeDifference"), Comment("maximum time difference between pulses at triggers")};
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit CrvWidebandTriggerFilter(const Parameters& conf);
    bool filter(art::Event& e);

    private:
    std::string       _crvRecoPulsesModuleLabel;
    std::vector<int>  _triggerSectorTypes;
    int               _PEthreshold;
    double            _maxTimeDifference;

  };

  CrvWidebandTriggerFilter::CrvWidebandTriggerFilter(const Parameters& conf) :
    art::EDFilter(conf),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _triggerSectorTypes(conf().triggerSectorTypes()),
    _PEthreshold(conf().PEthreshold()),
    _maxTimeDifference(conf().maxTimeDifference())
  {
  }

  bool CrvWidebandTriggerFilter::filter(art::Event& event)
  {
    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    if(!event.getByLabel(_crvRecoPulsesModuleLabel,crvRecoPulseCollection)) return(false);

    size_t nSectorTypes=_triggerSectorTypes.size();
    std::vector<std::vector<double> > triggerPulseTimes(nSectorTypes);  //organized by sector types
    for(size_t recoPulseIndex=0; recoPulseIndex<crvRecoPulseCollection->size(); ++recoPulseIndex)
    {
      const CrvRecoPulse &crvRecoPulse = crvRecoPulseCollection->at(recoPulseIndex);
      const CRSScintillatorBarIndex &barIndex = crvRecoPulse.GetScintillatorBarIndex();
      const CRSScintillatorBarId &crvCounterId = CRS->getBar(barIndex).id();
      int sectorNumber = crvCounterId.getShieldNumber();
      int sectorType = CRS->getCRSScintillatorShields().at(sectorNumber).getSectorType();
      if(std::find(_triggerSectorTypes.begin(),_triggerSectorTypes.end(), sectorType) == _triggerSectorTypes.end()) continue;

      if(crvRecoPulse.GetPEs()<_PEthreshold) continue;

      size_t vectorIndex=std::distance(_triggerSectorTypes.begin(),std::find(_triggerSectorTypes.begin(),_triggerSectorTypes.end(),sectorType));
      triggerPulseTimes[vectorIndex].push_back(crvRecoPulse.GetPulseTime());
    }

    for(size_t i=0; i<triggerPulseTimes.size(); ++i)
    {
      if(triggerPulseTimes[i].empty()) return(false);
      std::sort(triggerPulseTimes[i].begin(),triggerPulseTimes[i].end());
    }

    std::vector<size_t> currentIndices(triggerPulseTimes.size());  //starts at 0th time of all sector types
    std::vector<double> currentElements(triggerPulseTimes.size());
    for(size_t i=0; i<triggerPulseTimes.size(); ++i) currentElements[i]=triggerPulseTimes[i].at(0);  //0th time entry of all sector types
    while(true)
    {
      auto minTime=std::min_element(currentElements.begin(),currentElements.end());
      auto maxTime=std::max_element(currentElements.begin(),currentElements.end());
      double diff=*maxTime-*minTime;
      if(diff<=_maxTimeDifference) return(true);

      size_t elementIndexOfMinTime=std::distance(currentElements.begin(),minTime);  //which sector type had the earliest time for the current combination of time entries
      ++currentIndices[elementIndexOfMinTime];  //go to the next time entry of the sector type that had the earliest time
      if(currentIndices[elementIndexOfMinTime]>=triggerPulseTimes.at(elementIndexOfMinTime).size()) return(false);  //no more times left for this sector type
                                                                                                                    // --> won't find smaller time difference
      currentElements[elementIndexOfMinTime]=triggerPulseTimes[elementIndexOfMinTime].at(currentIndices[elementIndexOfMinTime]); //update the current time entry for this sector type
    }

    return(false);
  }
}

DEFINE_ART_MODULE(mu2e::CrvWidebandTriggerFilter)
