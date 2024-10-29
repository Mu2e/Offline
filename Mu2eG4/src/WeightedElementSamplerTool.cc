// art tool to sample different elements, according to a vector of weights
// Ed Callaghan, 2024

#include "Offline/Mu2eG4/inc/WeightedElementSamplerTool.hh"

namespace mu2e{
  WeightedElementSamplerTool::WeightedElementSamplerTool(std::string name):
      ElementSamplerTool(name){
    /**/
  }

  // propagate random engine to CLHEP distribution
  void WeightedElementSamplerTool::UseRandomEngine(art::RandomNumberGenerator::base_engine_t& engine){
    _uniform = std::make_unique<CLHEP::RandFlat>(engine);
  }

  double WeightedElementSamplerTool::GetWeight(std::string element){
    double rv = _weights.at(element);
    return rv;
  }

  std::string WeightedElementSamplerTool::sample_element(){
    // lazy initialization, to defer until after G4 has constructed materials
    if (!_initialized){
      this->initialize();
    }

    // sample a uniform on (0,1)
    double uniform = _uniform->fire();

    // and backtrack to which element this represents
    auto itr = std::upper_bound(_cumulative_weights.begin(),
                                _cumulative_weights.end(),
                                uniform);

    // sanity check
    if (itr == _cumulative_weights.end()){
      std::string msg = "random sample inconsistent with cumulative weights";
      throw cet::exception("WeightedElementSamplerTool") << msg << std::endl;
    }

    // look up this element and return
    size_t idx = std::distance(_cumulative_weights.begin(), itr);
    auto rv = _elements[idx];
    return rv;
  };

  // normalize elemental weights to sum to 1.0
  void WeightedElementSamplerTool::normalize_weights(){
    // first, compute the total weight...
    double sum = 0.0;
    for (const auto& pair: _weights){
      sum += pair.second;
    }

    // ...then, scale everything accordingly
    for (const auto& pair: _weights){
      auto key = pair.first;
      _weights[key] = pair.second / sum;
    }

    // compute cumulative weights, ordered as list of element names
    _cumulative_weights.resize(_elements.size());
    _cumulative_weights[0] = _weights[_elements[0]];
    for (size_t i = 1 ; i < _elements.size() ; i++){
      double last = _cumulative_weights[i-1];
      double curr = _weights[_elements[i]];
      _cumulative_weights[i] = last + curr;
    }
  }
}
