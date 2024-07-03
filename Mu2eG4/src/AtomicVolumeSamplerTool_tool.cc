// Sample elements from a material based on the relative volumes occupied
// in the material
// Ed Callaghan, 2024

// stl
#include <iostream>
#include <string>

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// cetlib_except
#include "cetlib_except/exception.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eG4/inc/AtomicVolumeSamplerTool.hh"

namespace mu2e{
  AtomicVolumeSamplerTool::AtomicVolumeSamplerTool(const Parameters& config):
      WeightedElementSamplerTool(config().name()){
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
}; // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::AtomicVolumeSamplerTool)
