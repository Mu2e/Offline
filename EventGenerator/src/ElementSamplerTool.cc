// art tool to sample different elements, according to some rule
// Ed Callaghan, 2024

// mu2e
#include "Offline/EventGenerator/inc/ElementSamplerTool.hh"

namespace mu2e{
  ElementSamplerTool::ElementSamplerTool(std::string name){
    _name = name;
    _initialized = false;
  }

  ElementSamplerTool::~ElementSamplerTool(){
    /**/
  }

  std::string ElementSamplerTool::Sample(){
    auto rv = this->sample_element();
    return rv;
  }

  std::vector<std::string>& ElementSamplerTool::Elements(){
    auto& rv = _elements;
    return rv;
  }

  void ElementSamplerTool::initialize(){
    _material = G4Material::GetMaterial(_name, false);

    // check that material is defined
    if (_material == NULL){
      std::string msg = "Undefined G4Material: " + _name;
      throw cet::exception("ElementSamplerTool") << msg << std::endl;
    }

    // populate own container of element names
    auto ptr = _material->GetElementVector();
    auto wrapped = *ptr;
    for (const auto& element: wrapped){
      std::string name = element->GetName();
      _elements.push_back(name);
    }

    // defer to subclass for any extra actions
    this->finish_initialize();
  }
}
