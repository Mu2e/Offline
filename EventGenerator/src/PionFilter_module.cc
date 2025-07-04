// Purpose: Event filter for RPC simulations
// author: S Middleton  2024
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/OptionalAtom.h"

// Mu2e includes.
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/SumOfWeights.hh"
#include "Offline/Mu2eUtilities/inc/SimParticleGetTau.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include <iostream>
#include <string>
#include "TH1F.h"
#include "TTree.h"
using namespace std;
namespace mu2e {

  class PionFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("Diagnostic print level"), 0};
        fhicl::OptionalAtom<double> tmin{Name("tmin"), Comment("Selected pion minimum end time")};
        fhicl::OptionalAtom<double> tmax{Name("tmax"), Comment("Selected pion maximum end time")};
        fhicl::OptionalAtom<int> maxPions{Name("maxPions"), Comment("Maximum number of pion stops")};
        fhicl::Atom<int> processCode{Name("processCode"), Comment("Pion end process code to select")};
        fhicl::Atom<bool> isNull{Name("isNull"), Comment("Skip filtering is turned on"), false};
      };
      explicit PionFilter(const art::EDFilter::Table<Config>& config);
      virtual bool filter(art::Event& event) override;
      virtual bool beginSubRun(art::SubRun& sr) override;
      virtual bool endSubRun(art::SubRun& sr) override;
      virtual void beginJob() override;
      virtual void endJob() override;

    private:
      const SimParticleCollection* SimCol_;
      int diagLevel_;
      double tmin_;
      double tmax_;
      int maxPions_;
      int processCode_;
      bool isNull_;
      SumOfWeights total_;
      SumOfWeights selected_;
  };

  PionFilter::PionFilter(const art::EDFilter::Table<Config>& config) :
     EDFilter{config}
    , diagLevel_{config().diagLevel()}
    , processCode_{config().processCode()}
    , isNull_{config().isNull()}
  {
    if(!config().tmin(tmin_)) tmin_ = -1.e10;
    if(!config().tmax(tmax_)) tmax_ =  1.e10;
    if(!config().maxPions(maxPions_)) maxPions_ = -1;
    produces<SumOfWeights, art::InSubRun>("total");
    produces<SumOfWeights, art::InSubRun>("selected");
  }

  void PionFilter::beginJob(){
      art::ServiceHandle<art::TFileService> tfs;
  }

  bool PionFilter::filter(art::Event& evt) {
      if(isNull_) return true;
      bool passed = false;
      std::vector<art::Handle<SimParticleCollection>> vah = evt.getMany<SimParticleCollection>();
      const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>();
      const std::vector<int> decayOffCodes = {PDGCode::pi_plus, PDGCode::pi_minus};
      int npions(0);
      for (auto const& ah : vah) { //always one collection
        for(const auto& aParticle : *ah){
          art::Ptr<SimParticle> pp(ah, aParticle.first.asUint());

          // check if this is a pion of interest
          if( pp->stoppingCode() == processCode_ and std::abs(pp->pdgId()) == PDGCode::pi_plus){
            const float globalTime = pp->endGlobalTime();
            const float tau = SimParticleGetTau::calculate(pp, decayOffCodes, gc);
            const float weight = std::exp(-tau);

            // count found pions
            total_.add(weight);
            ++npions;

            // check additional filters
            if(globalTime > tmin_ and globalTime < tmax_ ){
              passed = true;
              selected_.add(weight);
            }
          }
        }
      }

      // check global filters
      passed &= maxPions_ < 0 || npions <= maxPions_;

      // return the result
      return passed;
  }

  bool PionFilter::beginSubRun(art::SubRun&) {
    total_   .reset();
    selected_.reset();
    return true;
  }

  bool PionFilter::endSubRun(art::SubRun& sr) {
    sr.put(std::unique_ptr<SumOfWeights>(new SumOfWeights(total_   .sum(), total_   .count())), "total"   , art::fullSubRun());
    sr.put(std::unique_ptr<SumOfWeights>(new SumOfWeights(selected_.sum(), selected_.count())), "selected", art::fullSubRun());
    return true;
  }

  void PionFilter::endJob(){
     if(diagLevel_ > 0 ){
       std::cout<<"Total weight for all stops "<<total_.sum()<<std::endl;
       std::cout<<"Total stops "<<total_.count()<<std::endl;
       std::cout<<"Selected weight for chosen stops "<<selected_.sum()<<std::endl;
       std::cout<<"Selected stops "<<selected_.count()<<std::endl;
    }
  }
}

using mu2e::PionFilter;
DEFINE_ART_MODULE(PionFilter)
