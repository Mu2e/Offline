// Sample elements from a material based on the relative volumes occupied
// in the material
// Ed Callaghan, 2024

#ifndef Mu2eG4_AtomicVolumeSamplerTool_hh
#define Mu2eG4_AtomicVolumeSamplerTool_hh

// stl
#include <string>

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/Mu2eG4/inc/WeightedElementSamplerTool.hh"

namespace mu2e{
  class AtomicVolumeSamplerTool:public WeightedElementSamplerTool{

    public:
      struct Config{
        fhicl::Atom<std::string> name{
          fhicl::Name("material"),
          fhicl::Comment("Material name")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      AtomicVolumeSamplerTool(const Parameters&);
     ~AtomicVolumeSamplerTool(){ /**/ };
    protected:
      void finish_initialize();
    private:
      /**/
  };
}

#endif
