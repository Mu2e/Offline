#include "CLHEP/Units/SystemOfUnits.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

using namespace CLHEP;
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTrack/KalRep.hh"

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
  };

  TrkPatRecFilter::TrkPatRecFilter(fhicl::ParameterSet const& pset) :
    _fitterModuleLabels(pset.get<std::vector<std::string> >("fitterModuleLabels")),
    _trkPatRecInstances(pset.get<std::vector<std::string> >("trkPatRecInstances")),
    _minMomentum(pset.get<double>("minMomentum", -1.0)),
    _maxMomentum(pset.get<double>("maxMomentum", -1.0))
  {
    if(_fitterModuleLabels.size() != _trkPatRecInstances.size()) 
    {
        throw cet::exception("G4CONTROL")
          << "Sizes of fitterModuleLabels and trkPatRecInstances do not match.\n";
    }
  }

  bool TrkPatRecFilter::filter(art::Event& event) 
  {
    art::Handle<KalRepCollection> kalReps;
    for(unsigned int i=0; i<_fitterModuleLabels.size(); i++)
    {
      if(event.getByLabel(_fitterModuleLabels[i],_trkPatRecInstances[i],kalReps))
      {
        for(unsigned int j=0; j<kalReps->size(); j++)
        {
          KalRep const* kalrep = kalReps->at(j);
          if(kalrep!=NULL)
          {
            double fltLMin=kalrep->startValidRange();
            double fltLMax=kalrep->endValidRange();
            double p1=kalrep->momentum(fltLMin).mag();
            double p2=kalrep->momentum(fltLMax).mag();
            if(p1/CLHEP::MeV<_minMomentum && _minMomentum!=-1.0) continue;
            if(p2/CLHEP::MeV>_maxMomentum && _maxMomentum!=-1.0) continue;
            return(true);
          }
        }
      }
    }
    return(false);
  }
}

DEFINE_ART_MODULE(mu2e::TrkPatRecFilter);
