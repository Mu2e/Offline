#include "CLHEP/Units/SystemOfUnits.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "fhiclcpp/ParameterSet.h"

using namespace CLHEP;
#include "BTrk/KalmanTrack/KalRep.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include <string>
#include <vector>

using namespace std;

namespace mu2e 
{
  class TrkPatRecFilter : public art::EDFilter 
  {
    public:
    explicit TrkPatRecFilter(fhicl::ParameterSet const& pset);
    virtual ~TrkPatRecFilter() { }
    virtual bool filter(art::Event& e);

    private:
    std::vector<std::string> _fitterModuleLabels;
    std::vector<std::string> _trkPatRecInstances;
    double _minMomentum;
    double _maxMomentum;
    double _minCETime; // minimum time of CE
    bool _reverse; // reverse the logic of the final decision
    std::string _mcdigislabel;
    int _minCEHits; // minimum number of true CE hits to require
  };

  TrkPatRecFilter::TrkPatRecFilter(fhicl::ParameterSet const& pset) :
    art::EDFilter{pset},
    _fitterModuleLabels(pset.get<std::vector<std::string> >("fitterModuleLabels")),
    _trkPatRecInstances(pset.get<std::vector<std::string> >("trkPatRecInstances")),
    _minMomentum(pset.get<double>("minMomentum", -1.0)),
    _maxMomentum(pset.get<double>("maxMomentum", 1.0e6)),
    _minCETime(pset.get<double>("minCETime", 550.0)),
    _reverse(pset.get<bool>("ReverseTrackSelection",false)),
    _mcdigislabel(pset.get<string>("StrawHitMCLabel","makeSH")),
    _minCEHits(pset.get<int>("MinCEHits",-1))
  {
    if(_fitterModuleLabels.size() != _trkPatRecInstances.size()) 
    {
        throw cet::exception("G4CONTROL")
          << "Sizes of fitterModuleLabels and trkPatRecInstances do not match.\n";
    }
  }

  bool TrkPatRecFilter::filter(art::Event& event) 
  {
    bool retval(false);
    for(unsigned int i=0; i<_fitterModuleLabels.size(); i++)
    {
      art::Handle<KalRepPtrCollection> trksHandle;
      if(event.getByLabel(_fitterModuleLabels[i],_trkPatRecInstances[i],trksHandle))
      {
	KalRepPtrCollection const& kalReps = *trksHandle;
	for(unsigned int j=0; j<kalReps.size(); j++)
        {
          KalRep const* kalrep = kalReps.at(j).get();
          if(kalrep!=NULL)
          {
            double fltLMin=kalrep->startFoundRange();
            double fltLMax=kalrep->endFoundRange();
            double p1=kalrep->momentum(fltLMin).mag();
            double p2=kalrep->momentum(fltLMax).mag();
            if(p1/CLHEP::MeV<_minMomentum ) continue;
            if(p2/CLHEP::MeV>_maxMomentum ) continue;
            retval = true;
	    break;
          }
        }
      }
    }
    if(_reverse) retval = !retval;
    // check for MC truth
    if(retval && _minCEHits > 0){
      int ncehits(0);
      art::Handle<StrawDigiMCCollection> mcdigisHandle;
      if(event.getByLabel(_mcdigislabel,mcdigisHandle)){
	const StrawDigiMCCollection* mcdigis = mcdigisHandle.product();
	for(auto imcdigi = mcdigis->begin(); imcdigi != mcdigis->end(); ++imcdigi){
	  if( imcdigi->wireEndTime(StrawEnd::cal) > _minCETime ) {
	    if(imcdigi->strawGasStep(StrawEnd::cal)->simParticle()->genParticle().isNonnull() &&
		imcdigi->strawGasStep(StrawEnd::cal)->simParticle()->genParticle()->generatorId().isConversion()) ++ncehits;
	  }
	}
      }
      retval &= ncehits >= _minCEHits;
    }
    return(retval);
  }
}

DEFINE_ART_MODULE(mu2e::TrkPatRecFilter);
