// Ed Callaghan
// Tool to select SimParticles created by a specific G4 physics process
// November 2024

#ifndef EventMixing_ProcessVolumeSimParticleSelectionTool_hh
#define EventMixing_ProcessVolumeSimParticleSelectionTool_hh

// stl
#include <memory>

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/EventMixing/inc/InnerProtonAbsorberPseudoVolumeLookupTool.hh"
#include "Offline/EventMixing/inc/SimParticleSelectionTool.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

namespace mu2e{
class ProcessVolumeSimParticleSelectionTool: public SimParticleSelectionTool{
    public:
      /**/
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
          fhicl::Comment("Steps of particles descendent from a particle below this momentum are removed")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      ProcessVolumeSimParticleSelectionTool(const Parameters&);
     ~ProcessVolumeSimParticleSelectionTool() = default;

      virtual bool Select(const SimParticle&) override final;

    protected:
      ProcessCode _processCode;
      std::string _volume;
      double _momentum_threshold;
      // this will be obviated by direct queries via PhysicalVolumeMultiHelper
      std::unique_ptr<InnerProtonAbsorberPseudoVolumeLookupTool> _lookup;

    private:
      /**/
  };
} // namespace mu2e

#endif
