/*Author: S Middleton
//Date: Dec 2019
//Purpose:  Preselector for Straight Tracks
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>

using namespace std;

namespace mu2e {

  class CosmicTriggerPreselect : public art::EDFilter {
  public:
    explicit CosmicTriggerPreselect(fhicl::ParameterSet const& pset);
    virtual ~CosmicTriggerPreselect() { }

    bool filter( art::Event& event);
    virtual void beginJob() override;
    virtual void endJob() override;

  private:

    int _diagLevel;
    std::string _g4ModuleLabel;
    std::string _virtualDetectorLabel;
    art::InputTag _shTag;   
    art::InputTag _shpTag;


    std::string _trkPatRecModuleLabel;
    art::InputTag _tqTag;
    std::string _instanceName;
    TrkParticle _tpart;
    TrkFitDirection _fdir;
    SimParticleTimeOffset _toff;  // time offset smearing
    int nch; // number of combo hits
    const ComboHit*     comboHit;
    const ComboHitCollection* _chcol;
 
    int _virtualhit_good;
    int _trk_good;
   
    float _VDPID;
    float _PMIN;
    float _T0MIN;
    float _T0MAX;
    
    int    _nProcessed;
   
  };

  CosmicTriggerPreselect::CosmicTriggerPreselect(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _diagLevel(pset.get<int>("diagLevel",0)),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel","g4run")),
    _virtualDetectorLabel(pset.get<std::string>("virtualDetectorName","virtualdetector")),
    _shTag(pset.get<art::InputTag>("StrawHitCollectionTag","makePH")),
    _trkPatRecModuleLabel(pset.get<std::string>("trkPatRecModuleLabel")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
    _virtualhit_good(pset.get<int>("virtualhit",1)),
    _trk_good(pset.get<int>("trkgood",1)),
    _VDPID                    (pset.get<int>("VDPID",13)),
    _PMIN                     (pset.get<float>("PMIN",150.)),//1GeV maybe? i.e. high mom selection
    _T0MIN                    (pset.get<float>("T0MIN",700.)),
    _T0MAX                    (pset.get<float>("T0MAX",1720.)),
    _nProcessed(0){
  }

  void CosmicTriggerPreselect::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("diag");
  }

  bool CosmicTriggerPreselect::filter(art::Event& event) {
    ++_nProcessed;
    if (_nProcessed%100==0 && _diagLevel > 0) std::cout<<"Processing event from CosmicTriggerPreselect =  "<<_nProcessed <<std::endl;

    //Handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Tracker>() ) return false;

    //load the timeoffset
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    double _mbtime = accPar->deBuncherPeriod;
    _toff.updateMap(event);

    //Fetch ComboHits:
    auto chcolH = event.getValidHandle<mu2e::ComboHitCollection>(_shTag);

    if (chcolH.product() != 0){
      _chcol= chcolH.product();
      nch=_chcol->size();
      for(int istr=0; istr<nch;++istr) {
	comboHit=&_chcol->at(istr);
      }
    }
    else{
      std::cout << "CosmicTriggerPreselect: ComboHitCollection Handle NOT VALID" << std::endl;
      return false;
    }
    if (_diagLevel > 1) std::cout << "Combo Hits Number: "  <<  nch <<endl;
    // Handle ECAL virtual hits
    art::Handle<StepPointMCCollection> vdhits;
    event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);
        
    // Handle tracks
    art::Handle<CosmicTrackSeedCollection> trksHandle;
    event.getByLabel(_trkPatRecModuleLabel, trksHandle);
    const CosmicTrackSeedCollection& trks = *trksHandle;
    bool virtualhit_good=false; 
    if (_virtualhit_good ){
      if (vdhits.isValid()){
	for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter){
	  const StepPointMC& hit = *iter;    
	 //TODO - test tracker acceptance?
	  double hitTimeUnfolded = _toff.timeWithOffsetsApplied(hit);
	  double hitTime         = fmod(hitTimeUnfolded,_mbtime);

	  if (hit.simParticle()->pdgId()==_VDPID && hitTime>(double) _T0MIN && hitTime<(double) _T0MAX ){
	    virtualhit_good=true;
	    break;
	  }
	}
      }
      if (! virtualhit_good ) return false;
    }

    bool trk_good=false; 
    if (_trk_good ){
      for ( size_t itrk=0; itrk< trks.size(); ++itrk ){
	CosmicTrackSeed const* trk = trks.at(itrk).get();
	CLHEP::Hep3Vector mom = (1,1,1);//trk->momentum(0);
	if ( //TODO trk feature cuts
	     (mom.mag() > _PMIN   )){
	  trk_good=true;
	}
      }
      if (_diagLevel > 0) std::cout << "Tracks Number: "  <<  trks.size() << " good?" << trk_good << endl;
      if ( ! trk_good ) return false;
    }

    
    return true;
      
  }
  void CosmicTriggerPreselect::endJob(){
    cout << "CosmicTriggerPreselect filter end job:" << _nProcessed << " events processed" << endl;
  }
  
}

using mu2e::CosmicTriggerPreselect;
DEFINE_ART_MODULE(CosmicTriggerPreselect);*/
