//
//  Filter for selecting good helices (pat. rec. output): this is part of the track trigger
//  Original author: Dave Brown (LBNL) 3/1/2017
//
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"
#include "fhiclcpp/ParameterSet.h"
//#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/GeomHandle.hh"
//#include "GeometryService/inc/DetectorSystem.hh"
// data
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

using namespace CLHEP;
// c++
#include <string>
#include <vector>
#include <iostream>
#include <memory>

using namespace std;

namespace mu2e
{
  class CaloCosmicCalib : public art::EDFilter
  {
  public:
    explicit CaloCosmicCalib(fhicl::ParameterSet const& pset);
    virtual bool filter(art::Event& event) override;
    virtual bool beginRun(art::Run&   run   );
    virtual bool endRun( art::Run& run ) override;

  private:
    art::InputTag _clTag;
    int           _minncrystalhits;
    double        _minenergy, _maxenergy;
    std::string   _trigPath;
    int           _debug;
    // counters
    unsigned _nevt, _npass;
  };

  CaloCosmicCalib::CaloCosmicCalib(fhicl::ParameterSet const& pset) :
    art::EDFilter{pset},
    _clTag          (pset.get<art::InputTag>("CaloClusterCollection")),
    _minncrystalhits(pset.get<int>          ("MinNCrystalHits")),
    _minenergy      (pset.get<double>       ("MinEnergy")), //MeV
    _maxenergy      (pset.get<double>       ("MaxEnergy")), //MeV
    _trigPath       (pset.get<std::string>  ("triggerPath")),
    _debug          (pset.get<int>          ("debugLevel")),
    _nevt(0), _npass(0)
  {
    produces<TriggerInfo>();
  }

  bool CaloCosmicCalib::beginRun(art::Run & run){
    // get bfield
    // GeomHandle<BFieldManager> bfmgr;
    // GeomHandle<DetectorSystem> det;
    // Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    // _bz0 = bfmgr->getBField(vpoint_mu2e).z();
    return true;
  }

  bool CaloCosmicCalib::filter(art::Event& evt){
    // create output
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false); // preset to fail
    // find the collection
    auto clH = evt.getValidHandle<CaloClusterCollection>(_clTag);
    const CaloClusterCollection* clcol = clH.product();
    size_t trig_ind(0);
    // loop over the collection: if any pass the selection, pass this event
    for(auto icl = clcol->begin();icl != clcol->end(); ++icl) {
      auto const& cl = *icl;
      
      // get energy and cluster size
      float energy     = cl.energyDep();
      int   clsize     = cl.size();
      
      if(_debug > 2){
        std::cout << moduleDescription().moduleLabel() << " nhits = " << cl.size() << " energy = " << energy << std::endl;
      }
      if( (energy >= _minenergy) && 
	  (energy <= _maxenergy) && 
	  (clsize >= _minncrystalhits) ) {
        retval = true;
        ++_npass;
        // Fill the trigger info object
        if (trig_ind == 0){
	  triginfo->_triggerBits.merge(TriggerFlag::caloCalib);
	  triginfo->_triggerPath = _trigPath;
	}
        // associate to the caloCluster which triggers.  Note there may be other caloClusters which also pass the filter
        // but filtering is by event!
        size_t index = std::distance(clcol->begin(),icl);
	triginfo->_caloClusters.push_back(art::Ptr<CaloCluster>(clH,index));
	++trig_ind;
        if(_debug > 1){
          std::cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << std::endl;
        }
      }
    }
    evt.put(std::move(triginfo));
    return retval;
  }

  bool CaloCosmicCalib::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      std::cout << moduleDescription().moduleLabel() << " passed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << std::endl;
    }
    return true;
  }
}
using mu2e::CaloCosmicCalib;
DEFINE_ART_MODULE(CaloCosmicCalib);
