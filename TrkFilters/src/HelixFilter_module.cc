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
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
// data
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "DataProducts/inc/Helicity.hh"
// mu2e
#include "Mu2eUtilities/inc/HelixTool.hh"
using namespace CLHEP;
// c++
#include <string>
#include <vector>
#include <iostream>
#include <memory>

using namespace std;

namespace mu2e
{
  class HelixFilter : public art::EDFilter
  {
  public:
    explicit HelixFilter(fhicl::ParameterSet const& pset);
    virtual bool filter(art::Event& event) override;
    virtual bool beginRun(art::Run&   run   );
    virtual bool endRun( art::Run& run ) override;

  private:
    art::InputTag _hsTag;
    bool          _hascc; // Calo Cluster
    int           _hel;
    int           _minnstrawhits;
    double        _minmom, _maxmom;
    double        _maxpT;
    double        _minpT;
    double        _maxchi2XY;
    double        _maxchi2PhiZ;
    double        _maxd0;
    double        _mind0;
    double        _maxlambda;
    double        _minlambda;
    int           _maxnloops;
    double        _minnloops;
    double        _bz0;
    TrkFitFlag    _goodh; // helix fit flag
    std::string   _trigPath;
    int           _debug;
    // counters
    unsigned      _nevt, _npass;
  };

  HelixFilter::HelixFilter(fhicl::ParameterSet const& pset) :
    _hsTag        (pset.get<art::InputTag>("helixSeedCollection","PosHelixFinder")),
    _hascc        (pset.get<bool>  ("requireCaloCluster",false)),
    _hel          (pset.get<int>   ("helicity")),
    _minnstrawhits(pset.get<int>   ("minNStrawHits",15)),
    _minmom       (pset.get<double>("minMomentum",70.0)),
    _maxmom       (pset.get<double>("maxMomentum",120.0)),
    _minpT        (pset.get<double>("minPt", 0.)),
    _maxchi2XY    (pset.get<double>("maxChi2XY", 8.)),
    _maxchi2PhiZ  (pset.get<double>("maxChi2PhiZ", 8.)),
    _maxd0        (pset.get<double>("maxD0", 200.)),
    _mind0        (pset.get<double>("minD0", -200.)),
    _maxlambda    (pset.get<double>("maxAbsLambda",350.)),
    _minlambda    (pset.get<double>("minAbsLambda",150.)),
    _maxnloops    (pset.get<double>("maxNLoops",30.)),
    _minnloops    (pset.get<double>("minNLoops",0.)),
    _goodh        (pset.get<vector<string> >("helixFitFlag",vector<string>{"HelixOK"})),
    _trigPath     (pset.get<std::string>("triggerPath")),
    _debug        (pset.get<int>   ("debugLevel",0)),
    _nevt(0), _npass(0)
  {
    produces<TriggerInfo>();
  }

  bool HelixFilter::beginRun(art::Run & run){
    // get bfield
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();
    return true;
  }

  bool HelixFilter::filter(art::Event& evt){
    // create output
    unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false); // preset to fail
    // find the collection
    auto hsH = evt.getValidHandle<HelixSeedCollection>(_hsTag);
    const HelixSeedCollection* hscol = hsH.product();
    float mm2MeV = 3./10.*_bz0;
    // loop over the collection: if any pass the selection, pass this event
    for(auto ihs = hscol->begin();ihs != hscol->end(); ++ihs) {
      auto const& hs = *ihs;
      
      //check the helicity
      if (!(hs.helix().helicity() == Helicity(_hel)))        continue;

      HelixTool helTool(&hs, 3);
      // compute the helix momentum.  Note this is in units of mm!!!
      float hmom       = hs.helix().momentum()*mm2MeV;
      int   nstrawhits = helTool.nstrawhits();
      float hpT        = hs.helix().radius()*mm2MeV;
      float chi2XY     = hs.helix().chi2dXY();
      float chi2PhiZ   = hs.helix().chi2dZPhi();
      float d0         = hs.helix().rcent() - hs.helix().radius();
      float lambda     = std::fabs(hs.helix().lambda());
      int   nLoops     = helTool.nLoops();
      
      if(_debug > 2){
        cout << *currentContext()->moduleLabel() << "status = " << hs.status() << " nhits = " << hs.hits().size() << " mom = " << hmom << endl;
      }
      if( hs.status().hasAllProperties(_goodh) &&
          (!_hascc || hs.caloCluster().isNonnull()) &&
	  nstrawhits >= _minnstrawhits &&
          hpT        >= _minpT &&
	  chi2XY     <= _maxchi2XY &&
	  chi2PhiZ   <= _maxchi2PhiZ &&
	  d0         <= _maxd0 &&
	  d0         >= _mind0 &&
	  lambda     <= _maxlambda &&
	  lambda     >= _minlambda &&
	  nLoops     <= _maxnloops &&
	  nLoops     >= _minnloops &&
          hmom       >= _minmom    && 
	  hmom       <= _maxmom) {
        retval = true;
        ++_npass;
        // Fill the trigger info object
        triginfo->_triggerBits.merge(TriggerFlag::helix);
	triginfo->_triggerPath = _trigPath;
        // associate to the helix which triggers.  Note there may be other helices which also pass the filter
        // but filtering is by event!
        size_t index = std::distance(hscol->begin(),ihs);
        triginfo->_helix = art::Ptr<HelixSeed>(hsH,index);
        if(_debug > 1){
          cout << *currentContext()->moduleLabel() << " passed event " << evt.id() << endl;
        }
        break;
      }
    }
    evt.put(std::move(triginfo));
    return retval;
  }

  bool HelixFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      cout << *currentContext()->moduleLabel() << " paassed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << endl;
    }
    return true;
  }
}
using mu2e::HelixFilter;
DEFINE_ART_MODULE(HelixFilter);
