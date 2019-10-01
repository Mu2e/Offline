#include "CLHEP/Units/SystemOfUnits.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

using namespace CLHEP;

#include <string>
#include <vector>

using namespace std;

namespace mu2e 
{
  class TrackerStepPointFilter : public art::EDFilter 
  {
    public:
    explicit TrackerStepPointFilter(fhicl::ParameterSet const& pset);
    virtual ~TrackerStepPointFilter() { }
    virtual bool filter(art::Event& e);

    private:
    std::string _g4ModuleLabel;
    std::string _stepPointInstance;
    size_t _minStepPoints;
  };

  TrackerStepPointFilter::TrackerStepPointFilter(fhicl::ParameterSet const& pset) :
    art::EDFilter{pset},
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel","g4run")),
    _stepPointInstance(pset.get<string>("stepPointInstance","tracker")),
    _minStepPoints(pset.get<size_t>("minStepPoints",15))
  {
  }

  bool TrackerStepPointFilter::filter(art::Event& event) 
  {
    art::Handle<StepPointMCCollection> stepPoints;
    if(event.getByLabel(_g4ModuleLabel, _stepPointInstance, stepPoints))
    {
      if(stepPoints->size()>=_minStepPoints) return(true);
    }
    return(false);
  }
}

DEFINE_ART_MODULE(mu2e::TrackerStepPointFilter);
