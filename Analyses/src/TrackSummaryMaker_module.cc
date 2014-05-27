// Andrei Gaponenko, 2014

#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include "cetlib/exception.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "KalmanTests/inc/KalRepPtrCollection.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalHit.hh"
#include "TrkBase/TrkHelixUtils.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"

#include "RecoDataProducts/inc/TrackSummary.hh"
#include "Mu2eUtilities/inc/toHepPoint.hh"

namespace mu2e {

  class TrackSummaryMaker : public art::EDProducer {
  public:
    explicit TrackSummaryMaker(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt) override;
  private:
    art::InputTag trackInput_;
  };

  //================================================================
  TrackSummaryMaker::TrackSummaryMaker(const fhicl::ParameterSet& pset)
    : trackInput_(pset.get<std::string>("trackInput"))
  {
    produces<TrackSummaryCollection>();
  }

  //================================================================
  void TrackSummaryMaker::produce(art::Event& event) {
    GeomHandle<DetectorSystem> det;
    GeomHandle<VirtualDetector> vdg;

    std::unique_ptr<TrackSummaryCollection> output(new TrackSummaryCollection());

    auto ih = event.getValidHandle<KalRepPtrCollection>(trackInput_);
    for(const auto& krep : *ih) {
      if(krep->fitCurrent()){
        TrackSummary sum(krep->fitStatus().success(),
                         krep->charge(), krep->nActive(),
                         krep->nDof(), krep->chisq(),
                         krep->t0().t0(), krep->t0().t0Err(),
                         krep->flt0());

        // The following code is based on KalFitMC.cc
        CLHEP::Hep3Vector entpos = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
        double zent = entpos.z();
        double firsthitfltlen = krep->firstHit()->kalHit()->hitOnTrack()->fltLen();
        double lasthitfltlen = krep->lastHit()->kalHit()->hitOnTrack()->fltLen();
        double entlen = std::min(firsthitfltlen,lasthitfltlen);
        if(!TrkHelixUtils::findZFltlen(krep->traj(),zent,entlen,0.1)) {
          throw cet::exception("RUNTIME")<<"Error from findZFltlen()\n";
        }

        double loclen(0.0);
        TrackSummary::HelixParams helix(*krep->localTrajectory(entlen,loclen));
        TrackSummary::TrackStateAtPoint st(helix,
                                           fromHepPoint(krep->position(entlen)),
                                           krep->momentum(entlen),
                                           krep->momentumErr(entlen).covMatrix(),
                                           krep->arrivalTime(entlen),
                                           entlen
                                           );

        sum.addState(st);
        output->emplace_back(sum);
      }
      else {
        throw cet::exception("BADINPUT")<<"TrackSummaryMaker: do not know what to do with a fitCurrent==0 track\n";
      }
    }

    event.put(std::move(output));
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TrackSummaryMaker);
