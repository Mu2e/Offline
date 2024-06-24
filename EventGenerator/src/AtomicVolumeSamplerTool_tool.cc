// Sample elements from a material based on the relative volumes occupied
// in the material

// stl
#include <algorithm>
#include <iostream>
#include <string>

// art
#include <art/Utilities/ToolConfigTable.h>
#include <art/Utilities/ToolMacros.h>

// cetlib_except
#include <cetlib_except/exception.h>

// fhiclcpp
#include <fhiclcpp/types/Atom.h>
#include <fhiclcpp/types/Comment.h>
#include <fhiclcpp/types/Name.h>

// mu2e
#include <Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh>
#include <Offline/EventGenerator/inc/ElementSamplerTool.hh>
#include <Offline/GlobalConstantsService/inc/PhysicsParams.hh>

namespace mu2e{
  class AtomicVolumeSamplerTool:public ElementSamplerTool{
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
      std::map<std::string, double> _weights;
      std::vector<double> _cumulative_weights;

      void normalize_weights();
      void finish_initialize();
      std::string sample_element();
    private:
      /**/
  };

  AtomicVolumeSamplerTool::AtomicVolumeSamplerTool(const Parameters& config):
      ElementSamplerTool(config().name()){
    /**/
  }

  void AtomicVolumeSamplerTool::finish_initialize(){
    // fill weight vector with volume fractions
    // G4Material::GetFractionVector returns mass fractions, rho^(m)_i
    // then rho^(v)_i = rho^(m)_i * (atomic volume / atomic mass)
    // then normalize rho^(v)_i -> rho^(v)_i / sum_i rho^(v)_i
    auto handle = GlobalConstantsHandle<PhysicsParams>();
    auto material = G4Material::GetMaterial(_name, false);
    auto elements = *material->GetElementVector();
    auto mass_fractions = material->GetFractionVector();
    for (size_t i = 0 ; i < elements.size() ; i++){
      auto element = elements[i];
      auto mass_fraction = *(mass_fractions + i);
      std::string name = static_cast<std::string>(element->GetName());
      double atomic_radius = handle->getAtomicRadius(name);
      double atomic_volume = pow(atomic_radius, 3); // proportional to volume
      double atomic_mass = handle->getAtomicMass(name);
      double volume_fraction = mass_fraction * (atomic_volume / atomic_mass);
      _weights[name] = volume_fraction;
    }
    this->normalize_weights();
    _initialized = true;
  }

  std::string AtomicVolumeSamplerTool::sample_element(){
    // lazy initialization, to defer until after G4 has constructed materials
    if (!_initialized){
      this->initialize();
    }

    // sample a uniform on (0,1)
    double uniform = G4UniformRand();

    // and backtrack to which element this represents
    auto itr = std::upper_bound(_cumulative_weights.begin(),
                                _cumulative_weights.end(),
                                uniform);

    // sanity check
    if (itr == _cumulative_weights.end()){
      std::string msg = "random sample inconsistent with cumulative weights";
      throw cet::exception("AtomicVolumeSamplerTool") << msg << std::endl;
    }

    // look up this element and return
    size_t idx = std::distance(_cumulative_weights.begin(), itr);
    auto rv = _elements[idx];
    return rv;
  };

  // normalize elemental weights to sum to 1.0
  void AtomicVolumeSamplerTool::normalize_weights(){
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
}; // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::AtomicVolumeSamplerTool)
