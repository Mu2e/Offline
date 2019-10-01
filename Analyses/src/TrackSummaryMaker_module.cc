// Andrei Gaponenko, 2014

#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "Mu2eUtilities/inc/toHepPoint.hh"
#include "RecoDataProducts/inc/TrackSummary.hh"
#include "RecoDataProducts/inc/TrackSummaryRecoMap.hh"

namespace mu2e {

  class TrackSummaryMaker : public art::EDProducer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> trackInput{ Name("trackInput"), Comment("The input collection.") };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit TrackSummaryMaker(const Parameters& conf);

    void produce(art::Event& evt) override;
  private:
    art::InputTag trackInput_;
  };

  //================================================================
  TrackSummaryMaker::TrackSummaryMaker(const Parameters& conf)
    : art::EDProducer{conf}
    , trackInput_(conf().trackInput())
  {
    produces<TrackSummaryCollection>();
    produces<TrackSummaryRecoMap>();
  }

  //================================================================
  void TrackSummaryMaker::produce(art::Event& event) {
    GeomHandle<DetectorSystem> det;
    GeomHandle<VirtualDetector> vdg;

    std::unique_ptr<TrackSummaryCollection> output(new TrackSummaryCollection());
    std::unique_ptr<TrackSummaryRecoMap> recomap(new TrackSummaryRecoMap());

    const art::ProductID trackSummaryPID = event.getProductID<TrackSummaryCollection>();
    const art::EDProductGetter *trackSummaryGetter = event.productGetter(trackSummaryPID);

    auto ih = event.getValidHandle<KalRepPtrCollection>(trackInput_);
    for(unsigned itrack=0; itrack<ih->size(); ++itrack) {
      const auto& krep = (*ih)[itrack];
      if(krep->fitCurrent()){
        TrackSummary sum(krep->fitStatus().success(),
                         krep->charge(), krep->nActive(),
                         krep->nDof(), krep->chisq(),
                         krep->t0().t0(), krep->t0().t0Err(),
                         krep->flt0());

        // The following code is based on KalFitMC.cc
        CLHEP::Hep3Vector entpos = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
        double zent = entpos.z();
        double firsthitfltlen = krep->firstHit()->kalHit()->hit()->fltLen();
        double lasthitfltlen = krep->lastHit()->kalHit()->hit()->fltLen();
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

        recomap->addSingle(art::Ptr<KalRepPtr>(ih, itrack),
                           art::Ptr<TrackSummary>(trackSummaryPID, output->size(), trackSummaryGetter));

        output->emplace_back(sum);
      }
      else {
        throw cet::exception("BADINPUT")<<"TrackSummaryMaker: do not know what to do with a fitCurrent==0 track\n";
      }
    }

    event.put(std::move(output));
    event.put(std::move(recomap));
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TrackSummaryMaker);
