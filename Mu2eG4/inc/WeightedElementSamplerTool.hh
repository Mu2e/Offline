// art tool to sample different elements, according to a vector of weights
// Ed Callaghan, 2024

#ifndef Mu2eG4_WeightedElementSamplerTool_hh
#define Mu2eG4_WeightedElementSamplerTool_hh

// stl
#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

// art
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// clhep
#include "CLHEP/Random/RandFlat.h"

// mu2e
#include "Offline/Mu2eG4/inc/ElementSamplerTool.hh"

namespace mu2e{
  class WeightedElementSamplerTool: public ElementSamplerTool{
    public:
      WeightedElementSamplerTool(std::string);
      virtual void UseRandomEngine(art::RandomNumberGenerator::base_engine_t&);

      double GetWeight(std::string);
    protected:
      std::map<std::string, double> _weights;
      std::vector<double> _cumulative_weights;
      std::unique_ptr<CLHEP::RandFlat> _uniform;

      void normalize_weights();
      std::string sample_element();
    private:
      /**/
  };
}

#endif
