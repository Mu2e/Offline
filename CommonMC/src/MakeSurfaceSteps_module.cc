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
#include <utility>
#include <cmath>

namespace mu2e {

  class MakeSurfaceSteps : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int> debug{ Name("debugLevel"), Comment("Debug Level"), 0};
        fhicl::Atom<art::InputTag> vdstepmcs { Name("VDStepPointMCs"), Comment("Virtual Detector StepPointMC collection")};
        fhicl::Atom<art::InputTag> ipastepmcs { Name("IPAStepPointMCs"), Comment("IPAStepPointMC collection")};
        fhicl::Atom<art::InputTag> ststepmcs { Name("TargetStepPointMCs"), Comment("Stopping target StepPointMC collection")};
        fhicl::Atom<double> maxdgap{ Name("MaxDistGap"), Comment("Maximum dstance gap between aggregated StepPointMCs")};
        fhicl::Atom<double> maxtgap{ Name("MaxTimeGap"), Comment("Maximum time gap between aggregated StepPointMCs")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit MakeSurfaceSteps(const Parameters& conf);

    private:
      void produce(art::Event& e) override;
      using SSPair = std::pair<SurfaceId,cet::map_vector_key>; // key for pair of surfaceId, SimParticle
      using SPMCP = art::Ptr<StepPointMC>;
      using SPMCPV = std::vector<SPMCP>;
      using SPSMap = std::map< SSPair , SPMCPV >; // steps by surface, SimParticle
      using SPMCCH = art::Handle<StepPointMCCollection>;
      using SPMCCHV = std::vector< SPMCCH >;
      void fillMap(SPMCCH const& spmcch, SPSMap& spsmap);
      void fillStep(SPMCPV const& spmcptrs,
          ParticleData const& pdata, cet::map_vector_key pid, SurfaceStep& sgs);
      int debug_;
      double maxdgap_;
      double maxtgap_;
      art::InputTag vdstepmcs_, ipastepmcs_, ststepmcs_;
      std::map<VirtualDetectorId,SurfaceId> vdmap_; // map of VDIds to surfaceIds
  };

  MakeSurfaceSteps::MakeSurfaceSteps(const Parameters& config )  :
    art::EDProducer{config},
    debug_(config().debug()),
    maxdgap_(config().maxdgap()),
    maxtgap_(config().maxtgap()),
    vdstepmcs_(config().vdstepmcs()),
    ipastepmcs_(config().ipastepmcs()),
    ststepmcs_(config().ststepmcs())
    {
      consumesMany<StepPointMCCollection>();
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
    // create output
    std::unique_ptr<SurfaceStepCollection> ssc(new SurfaceStepCollection);
    // start with virtual detectors; these are copied directly
    auto const& vdspmccol_h = event.getValidHandle<StepPointMCCollection>(vdstepmcs_);
    auto const& vdspmccol = *vdspmccol_h;
    for(auto const& vdspmc : vdspmccol) {
      // only some VDs are kept
      auto isid = vdmap_.find(vdspmc.virtualDetectorId());
      if(isid != vdmap_.end())ssc->emplace_back(isid->second,vdspmc); // no aggregation of VD hits
    }
    auto nvdsteps = ssc->size();
    if(debug_ > 0)std::cout << "Added " << nvdsteps << " VD steps " << std::endl;
    // now IPA
    // colate adjacent IPA steppointMCs from the same particle. The following assumes the StepPointMCs from the same SimParticle are adjacent and (roughly) time-ordered. This is observationally true currently, but if G4 ever evolves so that its not, this code will need reworking
    auto const& ipaspmccol_h = event.getValidHandle<StepPointMCCollection>(ipastepmcs_);
    auto const& ipaspmccol = *ipaspmccol_h;
    if(debug_ > 0)std::cout << "IPA StepPointMC collection has " << ipaspmccol.size() << " Entries" << std::endl;
    // match steps by their sim particle. There is only 1 volume for IPA
    SurfaceStep ipastep;
    for(auto const& spmc : ipaspmccol) {
      bool added(false);
      // decide if this step is contiguous to existing steps already aggregated
      // first, test SimParticle (surface is defined by the collection)
      if(ipastep.simParticle() == spmc.simParticle()){
        if(debug_ > 2) std::cout << "Same SimParticle" << std::endl;
        // then time;
        auto tgap = fabs(ipastep.time()-spmc.time());
        auto dgap = (ipastep.endPosition() - XYZVectorF(spmc.position())).R();
        if(spmc.time() < ipastep.time())throw cet::exception("Simulation") << " StepPointMC times out-of-order" << std::endl;
        if(debug_ > 2) std::cout << "Time gap " << tgap << " Distance gap " << dgap << std::endl;
        if(tgap < maxtgap_ && dgap < maxdgap_){
          // accumulate this step into the existing surface step
          if(debug_ > 1)std::cout <<"Added step" << std::endl;
          ipastep.addStep(spmc,maxdgap_,maxtgap_);
          added = true;
        }
      }
      if(!added){
        // if the existing step is valid, save it
        if(ipastep.surfaceId() != SurfaceIdDetail::unknown && ipastep.simParticle().isNonnull())ssc->push_back(ipastep);
        // start a new SurfaceStep for this step
        ipastep = SurfaceStep(SurfaceId(SurfaceIdDetail::IPA),spmc);
        if(debug_ > 1)std::cout <<"New step" << std::endl;
      }
    }
    // save last step
    if(ipastep.surfaceId() != SurfaceIdDetail::unknown && ipastep.simParticle().isNonnull())ssc->push_back(ipastep);
    auto nipasteps = ssc->size() - nvdsteps;
    if(debug_ > 0)std::cout << "Added " << nipasteps << " IPA Steps "<< std::endl;
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
      if(ststep.simParticle() == spmc.simParticle()){
        if(debug_ > 2) std::cout << "Same SimParticle" << std::endl;
        // then time;
        auto tgap = fabs(ststep.time()-spmc.time());
        auto dgap = (ststep.endPosition() - XYZVectorF(spmc.position())).R();
        if(spmc.time() < ststep.time())throw cet::exception("Simulation") << " StepPointMC times out-of-order" << std::endl;
        if(debug_ > 2) std::cout << "Time gap " << tgap << " Distance gap " << dgap << std::endl;
        if(tgap < maxtgap_ && dgap < maxdgap_){
          // accumulate this step into the existing surface step
          if(debug_ > 1)std::cout <<"Added step" << std::endl;
          ststep.addStep(spmc,maxdgap_,maxtgap_);
          added = true;
        }
      }
      if(!added){
        // if the existing step is valid, save it
        if(ststep.surfaceId() != SurfaceIdDetail::unknown && ststep.simParticle().isNonnull())ssc->push_back(ststep);
        // start a new SurfaceStep for this step
        ststep = SurfaceStep(SurfaceId(SurfaceIdDetail::ST_Foils,spmc.volumeId()),spmc);
        if(debug_ > 1)std::cout <<"New step" << std::endl;
      }
    }
    if(ststep.surfaceId() != SurfaceIdDetail::unknown && ststep.simParticle().isNonnull())ssc->push_back(ststep);
    auto nststeps = ssc->size() - nipasteps - nvdsteps;
    if(debug_ > 0)std::cout << "Added " << nststeps << " stopping target Steps "<< std::endl;
    // finish
    event.put(move(ssc));
  }

}
DEFINE_ART_MODULE(mu2e::MakeSurfaceSteps)
