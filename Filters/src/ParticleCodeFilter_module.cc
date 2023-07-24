//
// Simple module to select events with SimParticles with specific codes
//
//  Original author: Dave Brown (LBNL)
//
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/types/Tuple.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include <string>

namespace mu2e {
  class ParticleCodeFilter : public art::EDFilter {
    public:
      using ParticleCodeConfig = fhicl::Sequence<fhicl::Tuple<int,std::string,std::string>>;

      struct ModuleConfig {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment ("Printout Level"), 0 };
        fhicl::Atom<art::InputTag> simParticles { Name("SimParticles"), Comment("SimParticle collection") };
        ParticleCodeConfig codeConfig { Name("ParticleCodes"), Comment("particle code to select: PDG, creation, and termination code \n"
            "(select 'uninitialized' to disable testing creation and/or termination codes)") };
      };

      struct ParticleCodeSelector {
        PDGCode::type pdgCode_;
        ProcessCode::enum_type creationCode_, terminationCode_;
        bool select(SimParticle const& part) const {
          return part.pdgId() == pdgCode_
            && ( creationCode_ == ProcessCode::uninitialized || part.creationCode() == creationCode_)
            && ( terminationCode_ == ProcessCode::uninitialized || part.stoppingCode() == terminationCode_);
        }
      };

      using Parameters = art::EDFilter::Table<ModuleConfig>;

      explicit ParticleCodeFilter(const Parameters& conf);
      bool filter(art::Event& event) override;
    private:
      int printLevel_;
      art::ProductToken<SimParticleCollection> const simParticles_;
      std::vector<ParticleCodeSelector> selectors_;
  };

  ParticleCodeFilter::ParticleCodeFilter(const Parameters& conf) : art::EDFilter{conf}
  , printLevel_(conf().printLevel())
    , simParticles_(consumes<SimParticleCollection>(conf().simParticles())) {
      for(auto const& pconfig : conf().codeConfig()) {
        ParticleCodeSelector pselector;
        pselector.pdgCode_ = static_cast<PDGCode::type>(std::get<0>(pconfig));
        pselector.creationCode_ = ProcessCode::findByName(std::get<1>(pconfig));
        pselector.terminationCode_ = ProcessCode::findByName(std::get<2>(pconfig));
        if(printLevel_ > 0) std::cout << "Creating selector of PDGcode " << pselector.pdgCode_
          << " creation code " << pselector.creationCode_
            << " termination code " << pselector.terminationCode_ << std::endl;
        selectors_.push_back(pselector);
      }
    }

  bool ParticleCodeFilter::filter(art::Event& event) {
    bool retval(false);
    auto simH = event.getValidHandle(simParticles_);
    auto const& simps = *simH;
    for (auto const& simp : simps) {
      if(printLevel_ > 1) std::cout << "Testing particle PDGcode " << simp.second.pdgId()
        << " creationcode =" << simp.second.creationCode()
          << " termination code =" << simp.second.stoppingCode() << std::endl;
      for( auto const& selector : selectors_ ){
        if(selector.select(simp.second)){
          if(printLevel_ > 0) std::cout << "Found matching particle " << std::endl;
          retval = true;
          break;
        }
      }
    }
    return retval;
  }

}
DEFINE_ART_MODULE(mu2e::ParticleCodeFilter)
