// Ed Callaghan
// Select detector steps which are not downstream of a certain process
// November 2024

#ifndef EventMixing_ProcessVolumeDetectorStepAntiSelectionTool_hh
#define EventMixing_ProcessVolumeDetectorStepAntiSelectionTool_hh

// stl
#include <unordered_set>

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"

// mu2e
#include "Offline/EventMixing/inc/DetectorStepSelectionTool.hh"
#include "Offline/EventMixing/inc/PseudoCylindricalVolumeLookupTool.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"

namespace mu2e{
  class ProcessVolumeDetectorStepAntiSelectionTool: public DetectorStepSelectionTool{
    public:
      struct Config{
        fhicl::Sequence<std::string> processes{
          fhicl::Name("processes"),
          fhicl::Comment("Steps of particles descendent from these processes are removed")
        };
        fhicl::Sequence<std::string> volumes{
          fhicl::Name("volumes"),
          fhicl::Comment("Steps of particles descendent from processes in these volumes are removed")
        };
        fhicl::Atom<double> momentum_threshold{
          fhicl::Name("momentum_threshold"),
          fhicl::Comment("Steps of particles descendent from a particle of below this momentum are removed")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      ProcessVolumeDetectorStepAntiSelectionTool(const Parameters&);
     ~ProcessVolumeDetectorStepAntiSelectionTool() = default;

      virtual bool Select(const CaloShowerStep&) override final;
      virtual bool Select(const CrvStep&)        override final;
      virtual bool Select(const StrawGasStep&)   override final;

    protected:
      std::unordered_set<ProcessCode::enum_type> _processCodes;
      std::unordered_set<std::string> _volumes;
      double _momentum_threshold;
      // this will be obviated by direct queries via PhysicalVolumeMultiHelper
      std::unique_ptr<PseudoCylindricalVolumeLookupTool> _lookup;

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
  bool ProcessVolumeDetectorStepAntiSelectionTool::select(const T& step){
    const auto& particle = step.simParticle();
    bool matched = this->heritage_match(*particle);
    bool rv = !matched;
    return rv;
  }
} // namespace mu2e

#endif
