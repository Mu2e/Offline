// art tool to sample different elements, according to some rule
// Ed Callaghan, 2024

#ifndef EventGenerator_ElementSamplerTool_hh
#define EventGenerator_ElementSamplerTool_hh

// stl
#include <iostream>
#include <string>
#include <vector>

// art
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// cetlib_except
#include "cetlib_except/exception.h"

// geant4
#include "Geant4/G4Material.hh"

namespace mu2e{
  class ElementSamplerTool{
    public:
      ElementSamplerTool(std::string name);
      virtual ~ElementSamplerTool();

      virtual void UseRandomEngine(art::RandomNumberGenerator::base_engine_t&) = 0;
      virtual std::string sample_element() = 0;
      std::string Sample();
      std::vector<std::string>& Elements();
    protected:
      std::string _name;
      G4Material* _material;
      std::vector<std::string> _elements;
      bool _initialized;

      virtual void finish_initialize() = 0;
      void initialize();

    private:
      /**/
  };
}; // namespace mu2e

#endif
