//
// Reconstruction-level PID determination from a track.  The current implementation
// just uses the calorimeter matching information (TrkCaloHit), eventually it will
// also process dE/dx, etc.
//
// Original author: Dave Brown (LBNL)
//

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
// utilities
#include "Mu2eUtilities/inc/MVATools.hh"
#include "TrkDiag/inc/InfoStructHelper.hh"
// data
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkCaloHitPID.hh"
#include "RecoDataProducts/inc/RecoQual.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <float.h>
#include <vector>
using namespace std;

namespace mu2e {
  class TrackPID : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<int> debug{ Name("debugLevel"),
	Comment("Debug Level"), 0};
      fhicl::Atom<float> MaxDE{ Name("MaxDE"),
	Comment("Maximum difference between calorimeter cluster EDep energy and the track energy (assuming electron mass)")};
      fhicl::Atom<float> DT{ Name("DeltaTOffset"),
	Comment("Track - Calorimeter time offset")}; // this should be a condition FIXME!
      fhicl::Atom<art::InputTag> KSC { Name("KalSeedCollection"),
	Comment("KalSeedCollection producer")};
      fhicl::Table<typename MVATools::Config> MVAConfig{fhicl::Name("MVAConfig"),
	fhicl::Comment("MVA Configuration") };
    };

   using Parameters = art::EDProducer::Table<Config>;
   explicit TrackPID(const Parameters& conf);

  private:
   void produce(art::Event& event) override;
    bool _debug;
    float _maxde, _dtoffset;
    art::InputTag _kalSeedTag;
    MVATools* _tchmva;
  };

  TrackPID::TrackPID(const Parameters& config ) :
    art::EDProducer(config),
    _debug(config().debug()),
    _maxde(config().MaxDE()),
    _dtoffset(config().DT()),
    _kalSeedTag(config().KSC()),
    _tchmva(new MVATools(config().MVAConfig()))
  {
    produces<TrkCaloHitPIDCollection>();
    produces<RecoQualCollection>();
    _tchmva->initMVA();
    if(_debug> 0)_tchmva->showMVA();
  }

  void TrackPID::produce(art::Event& event ) {
    mu2e::GeomHandle<mu2e::Calorimeter> calo;
    // create output
    unique_ptr<TrkCaloHitPIDCollection> tchpcol(new TrkCaloHitPIDCollection());
    unique_ptr<RecoQualCollection> rqcol(new RecoQualCollection());
    // get the KalSeeds
    art::Handle<KalSeedCollection> kalSeedHandle;
    event.getByLabel(_kalSeedTag, kalSeedHandle);
    const auto& kalSeeds = *kalSeedHandle;

    for (const auto& kseed : kalSeeds) {
      TrkCaloHitPID tchpid;
      tchpid.setMVAStatus(MVAStatus::unset);
      tchpid.setMVAValue(-1.0);

      static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
      if (kseed.status().hasAllProperties(goodfit)){
	if(kseed.hasCaloCluster() &&
	    kseed.caloHit().flag().hasAllProperties(StrawHitFlag::active)){
	  auto const& tchs = kseed.caloHit();
	  auto const& cc = tchs.caloCluster();
	  XYZVec trkmom;
	  auto ikseg = kseed.nearestSegment(tchs.trkLen());
	  if(ikseg != kseed.segments().end())
	    ikseg->mom(ikseg->localFlt(tchs.trkLen()),trkmom);
	  else
	    throw cet::exception("RECO")<<"mu2e::TrackPID: KalSeed segment missing" << endl;
	  // compute the energy difference, assuming an electron mass
	  tchpid[TrkCaloHitPID::DeltaE] = cc->energyDep() - sqrt(trkmom.Mag2());
	  tchpid[TrkCaloHitPID::ClusterLen] = tchs.hitLen();
	  // move into detector coordinates.  Yikes!!
	  XYZVec cpos = Geom::toXYZVec(calo->geomUtil().mu2eToTracker(calo->geomUtil().diskFFToMu2e( cc->diskId(), cc->cog3Vector())));
	  tchpid[TrkCaloHitPID::RPOCA] = sqrt(cpos.Perp2());
	  // compute transverse direction WRT position
	  cpos.SetZ(0.0);
	  trkmom.SetZ(0.0);
	  tchpid[TrkCaloHitPID::TrkDir] = cpos.Dot(trkmom)/sqrt(cpos.Mag2()*trkmom.Mag2());
	  // the following includes the (Calibrated) light-propagation time delay.  It should eventually be put in the reconstruction FIXME!
	  // This velocity should come from conditions FIXME!
	  tchpid[TrkCaloHitPID::DeltaT] = tchs.t0().t0()-tchs.time()- std::min((float)200.0,std::max((float)0.0,tchs.hitLen()))*0.005 - _dtoffset;
	  // hard cut on the energy difference.  This rejects cosmic rays which hit the calo and produce an upstream-going track that is then
	  // reconstructed as a downstream particle associated to this cluster
	  if(tchpid[TrkCaloHitPID::DeltaE] < _maxde){
	    // evaluate the MVA
	    tchpid.setMVAValue(_tchmva->evalMVA(tchpid.values()));
	    tchpid.setMVAStatus(MVAStatus::calculated);
	  }
	}
      }
      tchpcol->push_back(tchpid);
      rqcol->push_back(RecoQual(tchpid.status(),tchpid.MVAValue()));
    }
    // put the output products into the event
    event.put(move(tchpcol));
    event.put(move(rqcol));
  }
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackPID);
