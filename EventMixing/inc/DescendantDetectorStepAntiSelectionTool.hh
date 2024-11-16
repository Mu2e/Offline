// Ed Callaghan
// Select detector steps which are not downstream of a certain process
// November 2024

#ifndef EventMixing_DescendantDetectorStepAntiSelectionTool_hh
#define EventMixing_DescendantDetectorStepAntiSelectionTool_hh

// stl
#include <memory>

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

// fhiclcpp
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/EventMixing/inc/DetectorStepSelectionTool.hh"
#include "Offline/EventMixing/inc/SimParticleSelectionTool.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"

namespace mu2e{
  class DescendantDetectorStepAntiSelectionTool: public DetectorStepSelectionTool{
    public:
      struct Config{
        fhicl::DelegatedParameter selection{
          fhicl::Name("particle_selection"),
          fhicl::Comment("Tool configuration to reject SimParticles")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      DescendantDetectorStepAntiSelectionTool(const Parameters&);
     ~DescendantDetectorStepAntiSelectionTool() = default;

      virtual bool Select(const CaloShowerStep&) override final;
      virtual bool Select(const CrvStep&)        override final;
      virtual bool Select(const StrawGasStep&)   override final;

    protected:
      std::unique_ptr<SimParticleSelectionTool> _selection;

      // templated to reduce code surface
      template<typename T>
      bool select(const T&);

      bool particle_match(const SimParticle&);
      bool heritage_match(const SimParticle&);

    private:
      /**/
  };

  // actively select steps from particles which do _not_ match configuration
  template<typename T>
  bool DescendantDetectorStepAntiSelectionTool::select(const T& step){
    const auto& particle = step.simParticle();
    bool matched = this->heritage_match(*particle);
    bool rv = !matched;
    return rv;
  }
} // namespace mu2e

#endif
