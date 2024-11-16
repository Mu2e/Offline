// Ed Callaghan
// Tool to select SimParticles downstread by a specific Mu2e generator
// November 2024

#ifndef EventMixing_GenIdSimParticleSelectionTool_hh
#define EventMixing_GenIdSimParticleSelectionTool_hh

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/EventMixing/inc/SimParticleSelectionTool.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

namespace mu2e{
class GenIdSimParticleSelectionTool: public SimParticleSelectionTool{
    public:
      /**/
      struct Config{
        fhicl::Atom<std::string> genId{
          fhicl::Name("GenId"),
          fhicl::Comment("Steps of particles descendent from this generator")
        };
        fhicl::Atom<double> momentum_threshold{
          fhicl::Name("momentum_threshold"),
          fhicl::Comment("Steps of particles descendent from a particle below this momentum are removed")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      GenIdSimParticleSelectionTool(const Parameters&);
     ~GenIdSimParticleSelectionTool() = default;

      virtual bool Select(const SimParticle&) override final;

    protected:
      GenId _genId;
      double _momentum_threshold;

    private:
      /**/
  };
} // namespace mu2e

#endif
