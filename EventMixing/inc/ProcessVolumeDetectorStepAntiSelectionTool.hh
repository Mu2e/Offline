// Ed Callaghan
// Select detector steps which are not downstream of a certain process
// November 2024

#ifndef EventMixing_ProcessVolumeDetectorStepAntiSelectionTool_hh
#define EventMixing_ProcessVolumeDetectorStepAntiSelectionTool_hh

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/EventMixing/inc/DetectorStepSelectionTool.hh"
#include "Offline/EventMixing/inc/InnerProtonAbsorberPseudoVolumeLookupTool.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"

namespace mu2e{
  class ProcessVolumeDetectorStepAntiSelectionTool: public DetectorStepSelectionTool{
    public:
      struct Config{
        fhicl::Atom<std::string> process{
          fhicl::Name("process"),
          fhicl::Comment("Steps of particles descendent from this process are removed")
        };
        fhicl::Atom<std::string> volume{
          fhicl::Name("volume"),
          fhicl::Comment("Steps of particles descendent from process in this volume are removed")
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
      ProcessCode _processCode;
      std::string _volume;
      double _momentum_threshold;
      // this will be obviated by direct queries via PhysicalVolumeMultiHelper
      std::unique_ptr<InnerProtonAbsorberPseudoVolumeLookupTool> _lookup;

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
