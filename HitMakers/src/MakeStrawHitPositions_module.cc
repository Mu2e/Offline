//
// A module to create strawhitpositions out of the strawhits
//
// $Id: MakeStrawHitPositions_module.cc,v 1.3 2014/04/28 13:34:43 brownd Exp $
// $Author: brownd $
// $Date: 2014/04/28 13:34:43 $
// 
//  Original Author: David Brown, LBNL
//  changes from G. Pezzullo, P. Murat

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"

// art includes.
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// root
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TMarker.h"
#include "TList.h"
#include "TLegend.h"

// C++ includes.
#include <iostream>
#include <float.h>

//#include "CalPatRec/inc/THackData.hh"

using namespace std;

namespace mu2e {


  void printShP(StrawHitPosition shp){
    double x,y,z;
    x = shp.pos().x();
    y = shp.pos().y();
    z = shp.pos().z();
    
    printf("%5.3f  %5.3f  %5.3f\n", x,y,z );
  }

  class MakeStrawHitPositions : public art::EDProducer {

  public:
    explicit MakeStrawHitPositions(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.

    void produce( art::Event& e);
    void beginJob();
    void beginRun(art::Run& Run);
    void printHits(const StrawHit& Hit, StrawHitPosition& Pos, int &banner);
  private:

    // Diagnostics level.
    int _debugLevel;
    int _printHits;
    bool _misscluster; // allow missing cluster for time division
    // Name of the StrawHit collection
    string _shLabel;
    //    THackData* fHackData;
  };

  MakeStrawHitPositions::MakeStrawHitPositions(fhicl::ParameterSet const& pset) :

    // Parameters
    _debugLevel(pset.get<int>("debugLevel",0)),
    _printHits(pset.get<int>("printHits",0)),
    _misscluster(pset.get<bool>("AllowMissingClusters",true)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH"))
 {
    // Tell the framework what we make.
    produces<StrawHitPositionCollection>();
  }

   void MakeStrawHitPositions::beginJob(){
    // create diagnostics if requested
    if(_debugLevel > 0){
      art::ServiceHandle<art::TFileService> tfs;
    }

  }

  void MakeStrawHitPositions::beginRun(art::Run& Run) {
  }

  void MakeStrawHitPositions::printHits(const StrawHit& Hit, StrawHitPosition& Pos, int &banner) {
    double x,y,z;
    Tracker const& tracker = getTrackerOrThrow();

    Straw const& straw = tracker.getStraw(Hit.strawIndex());
    x = Pos.pos().x();
    y = Pos.pos().y();
    z = Pos.pos().z();

    if(banner==0){
      printf("-----------------------------------------------------------------------------------");
      printf("------------------------------\n");
      printf("   x       y     z   SHID    Station Panel Layer Straw     Flags      Time          dt       eDep \n");
      printf("-----------------------------------------------------------------------------------");
      printf("------------------------------\n");
      banner = false;
    }

    printf("%5.3f %5.3f %5.3f %5i  %5i  %5i   %5i   %5i   0x%s %8.3f   %8.3f   %9.6f\n",
	   x,y,z,
	   Hit.strawIndex().asInt(),
	   straw.id().getPlane(),
	   straw.id().getPanel(),
	   straw.id().getLayer(),
	   straw.id().getStraw(),
	   Pos.flag().hex().data(),
	   Hit.time(),
	   Hit.dt(),
	   Hit.energyDep() );

  }

  void MakeStrawHitPositions::produce(art::Event& event) {
    Tracker const& tracker = getTrackerOrThrow();
    ConditionsHandle<StrawResponse> srep = ConditionsHandle<StrawResponse>("ignored");

    if ( _debugLevel > 0 ) cout << "MakeStrawHitPositions: produce() begin; event " << event.id().event() << endl;

    //    art::Handle<mu2e::StrawHitCollection> strawhitsH; 
    auto strawhitsH = event.getValidHandle<mu2e::StrawHitCollection>(_shLabel);
    const StrawHitCollection* strawhits(0);
    if (strawhitsH.product() != 0)
      //    if(event.getByLabel(_shLabel,strawhitsH))
      strawhits = strawhitsH.product();
    if(strawhits == 0){
      throw cet::exception("RECO")<<"mu2e::MakeStrawHitPositions: No StrawHit collection found for label " <<  _shLabel << endl;
    }
    // create a collection of StrawHitPosition, and intialize them using the time division
    size_t nsh = strawhits->size();
    unique_ptr<StrawHitPositionCollection> shpcol(new StrawHitPositionCollection);
    shpcol->reserve(nsh);

 //01 - 13 - 2014 gianipez added some printout
    int banner(0);
    static const double invsqrt12 = 1.0/sqrt(12.0);

    for(size_t ish=0;ish<nsh;++ish){
      StrawHit const& hit = strawhits->at(ish);
      Straw const& straw = tracker.getStraw(hit.strawIndex());
// get distance along wire from the straw center and it's estimated error
      float dw, dwerr;
      bool td = srep->wireDistance(hit,straw.getHalfLength(),dw,dwerr);
// create and fill the position struct
      StrawHitPosition shp;
      shp._pos = straw.getMidPoint()+dw*straw.getDirection();
      shp._wdir = straw.getDirection();
      shp._wdist = dw;
      shp._wres = dwerr;
      // crude initial estimate of the 
      shp._tres = straw.getRadius()*invsqrt12;
// if time division worked, flag the position accordingly
      if(td || _misscluster)shp._flag.merge(StrawHitFlag::tdiv); 
      if (_printHits>0) {
	printHits(hit,shp, banner);
	banner=1;
      }

      shpcol->push_back(shp);
    }

    event.put(std::move(shpcol));
  } // end MakeStrawHitPositions::produce.
} // end namespace mu2e

using mu2e::MakeStrawHitPositions;
DEFINE_ART_MODULE(MakeStrawHitPositions)

