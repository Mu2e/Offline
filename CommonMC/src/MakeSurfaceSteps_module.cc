//
//  This module creates SurfaceStep objects from G4 StepPointMCs
//
//  Original author: David Brown (LBNL) 7/2024
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Offline/GlobalConstantsService/inc/ParticleData.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

namespace mu2e {

  class MakeSurfaceSteps : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int> debug{ Name("debugLevel"), Comment("Debug Level"), 0};
        fhicl::Atom<art::InputTag> vdstepmcs { Name("VDStepPointMCs"), Comment("Virtual Detector StepPointMC collection")};
        fhicl::Atom<art::InputTag> absstepmcs { Name("AbsorberStepPointMCs"), Comment("Absorber StepPointMC collection")};
        fhicl::Atom<art::InputTag> ststepmcs { Name("TargetStepPointMCs"), Comment("Stopping target StepPointMC collection")};
        fhicl::Atom<double> maxdgap{ Name("MaxDistGap"), Comment("Maximum dstance gap between aggregated StepPointMCs")};
        fhicl::Atom<double> maxtgap{ Name("MaxTimeGap"), Comment("Maximum time gap between aggregated StepPointMCs")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit MakeSurfaceSteps(const Parameters& conf);
    private:
      void produce(art::Event& e) override;
      int debug_;
      double maxdgap_;
      double maxtgap_;
      art::ProductToken<StepPointMCCollection> vdstepmcs_, absstepmcs_, ststepmcs_;
      std::map<VirtualDetectorId,SurfaceId> vdmap_; // map of VDIds to surfaceIds
  };

  MakeSurfaceSteps::MakeSurfaceSteps(const Parameters& config )  :
    art::EDProducer{config},
    debug_(config().debug()),
    maxdgap_(config().maxdgap()),
    maxtgap_(config().maxtgap()),
    vdstepmcs_{ consumes<StepPointMCCollection>(config().vdstepmcs())},
    absstepmcs_{ consumes<StepPointMCCollection>(config().absstepmcs())},
    ststepmcs_{ consumes<StepPointMCCollection>(config().ststepmcs())}
    {
      produces <SurfaceStepCollection>();
      // build the VDId -> SId map by hand. This should come from a service TODO
      vdmap_[VirtualDetectorId(VirtualDetectorId::TT_FrontHollow)] = SurfaceId("TT_Front");
      vdmap_[VirtualDetectorId(VirtualDetectorId::TT_Mid)] = SurfaceId("TT_Mid");
      vdmap_[VirtualDetectorId(VirtualDetectorId::TT_MidInner)] = SurfaceId("TT_Mid");
      vdmap_[VirtualDetectorId(VirtualDetectorId::TT_Back)] = SurfaceId("TT_Back");
      vdmap_[VirtualDetectorId(VirtualDetectorId::TT_OutSurf)] = SurfaceId("TT_Outer");
      vdmap_[VirtualDetectorId(VirtualDetectorId::TT_InSurf)] = SurfaceId("TT_Inner");
    }

  void MakeSurfaceSteps::produce(art::Event& event) {
    GeomHandle<DetectorSystem> det;
    // create output
    std::unique_ptr<SurfaceStepCollection> ssc(new SurfaceStepCollection);
    // start with virtual detectors; these are copied directly
    auto const& vdspmccol_h = event.getValidHandle<StepPointMCCollection>(vdstepmcs_);
    auto const& vdspmccol = *vdspmccol_h;
    for(auto const& vdspmc : vdspmccol) {
      // only some VDs are kept
      auto isid = vdmap_.find(vdspmc.virtualDetectorId());
      if(isid != vdmap_.end())ssc->emplace_back(isid->second,vdspmc,det); // no aggregation of VD hits
    }
    auto nvdsteps = ssc->size();
    if(debug_ > 0)std::cout << "Added " << nvdsteps << " VD steps " << std::endl;
    // now absorbers: IPA, OPA, and TSDA
    // colate adjacent steppointMCs from the same particle. The following assumes the StepPointMCs from the same SimParticle are adjacent and (roughly) time-ordered. This is observationally true currently, but if G4 ever evolves so that its not, this code will need reworking
    auto const& absspmccol_h = event.getValidHandle<StepPointMCCollection>(absstepmcs_);
    auto const& absspmccol = *absspmccol_h;
    if(debug_ > 0)std::cout << "Absorber StepPointMC collection has " << absspmccol.size() << " Entries" << std::endl;
    // match steps by their sim particle. There is only 1 volume for IPA
    SurfaceStep absorberstep;
    for(auto const& spmc : absspmccol) {
      bool added(false);
      // decide if this step is contiguous to existing steps already aggregated
      // first, test SimParticle (surface is defined by the collection)
      if(absorberstep.simParticle() == spmc.simParticle()){
        if(debug_ > 2) std::cout << "Same SimParticle" << std::endl;
        // then time;
        auto tgap = fabs(absorberstep.time()-spmc.time());
        auto dgap = (absorberstep.endPosition() - XYZVectorF(det->toDetector(spmc.position()))).R();
        if(spmc.time() < absorberstep.time())throw cet::exception("Simulation") << " StepPointMC times out-of-order" << std::endl;
        if(debug_ > 2) std::cout << "Time gap " << tgap << " Distance gap " << dgap << std::endl;
        if(tgap < maxtgap_ && dgap < maxdgap_){
          // accumulate this step into the existing surface step
          if(debug_ > 1)std::cout <<"Added step" << std::endl;
          absorberstep.addStep(spmc,det);
          added = true;
        }
      }
      if(!added){
        // if the existing step is valid, save it
        if(absorberstep.surfaceId() != SurfaceIdDetail::unknown && absorberstep.simParticle().isNonnull())ssc->push_back(absorberstep);
        // start a new SurfaceStep for this step
        // hack around the fact that the IPA, OPA, and TSDA are all in the same collection ('absorber')
        SurfaceId sid;
        if(spmc.position().Z < -5960) // TSDA
          sid = SurfaceId(SurfaceIdDetail::IPA);
        else if(spmc.position().rho() > 310) // OPA
          sid = SurfaceId(SurfaceIdDetail::OPA);
        else
          sid = SurfaceId(SurfaceIdDetail::IPA);
        absorberstep = SurfaceStep(sid,spmc,det);
        if(debug_ > 1)std::cout <<"New step" << std::endl;
      }
    }
    // save last step
    if(absorberstep.surfaceId() != SurfaceIdDetail::unknown && absorberstep.simParticle().isNonnull())ssc->push_back(absorberstep);
    auto nabsorbersteps = ssc->size() - nvdsteps;
    if(debug_ > 0)std::cout << "Added " << nabsorbersteps << " IPA Steps "<< std::endl;
    // same for stopping target. Here, the foil number (= volumeid) matters
    auto const& stspmccol_h = event.getValidHandle<StepPointMCCollection>(ststepmcs_);
    auto const& stspmccol = *stspmccol_h;
    if(debug_ > 0)std::cout << "Target StepPointMC collection has " << stspmccol.size() << " Entries" << std::endl;
    // match steps by their sim particle. There is only 1 volume for IPA
    SurfaceStep ststep;
    for(auto const& spmc : stspmccol) {
      bool added(false);
      // decide if this step is contiguous to existing steps already aggregated
      // first, test SimParticle (surface is defined by the collection)
      if(ststep.simParticle() == spmc.simParticle() && ststep.surfaceId().index()== (int)spmc.volumeId()){
        if(debug_ > 2) std::cout << "Same SimParticle" << std::endl;
        // then time;
        auto tgap = fabs(ststep.time()-spmc.time());
        auto dgap = (ststep.endPosition() - XYZVectorF(det->toDetector(spmc.position()))).R();
        if(spmc.time() < ststep.time())throw cet::exception("Simulation") << " StepPointMC times out-of-order" << std::endl;
        if(debug_ > 2) std::cout << "Time gap " << tgap << " Distance gap " << dgap << std::endl;
        if(tgap < maxtgap_ && dgap < maxdgap_){
          // accumulate this step into the existing surface step
          if(debug_ > 1)std::cout <<"Added step" << std::endl;
          ststep.addStep(spmc,det);
          added = true;
        }
      }
      if(!added){
        // if the existing step is valid, save it
        if(ststep.surfaceId() != SurfaceIdDetail::unknown && ststep.simParticle().isNonnull())ssc->push_back(ststep);
        // start a new SurfaceStep for this step
        // hack around the problem that the ST wires are mixed in with the ST foils
        SurfaceId stid(SurfaceIdDetail::ST_Foils,spmc.volumeId());
        if(det->toDetector(spmc.position()).rho() > 75.001)stid = SurfaceId(SurfaceIdDetail::ST_Wires,spmc.volumeId());
        ststep = SurfaceStep(stid,spmc,det);
        if(debug_ > 1)std::cout <<"New step" << std::endl;
      }
    }
    if(ststep.surfaceId() != SurfaceIdDetail::unknown && ststep.simParticle().isNonnull())ssc->push_back(ststep);
    auto nststeps = ssc->size() - nabsorbersteps - nvdsteps;
    if(debug_ > 0)std::cout << "Added " << nststeps << " stopping target Steps "<< std::endl;
    // finish
    event.put(move(ssc));
  }

}
DEFINE_ART_MODULE(mu2e::MakeSurfaceSteps)
