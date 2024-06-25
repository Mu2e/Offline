// art tool to sample different elements, according to some rule
// Ed Callaghan, 2024

#ifndef EventGenerator_ElementSamplerTool_hh
#define EventGenerator_ElementSamplerTool_hh

// stl
#include <iostream>
#include <string>
#include <vector>

// cetlib_except
#include <cetlib_except/exception.h>

// geant4
#include <Geant4/G4Material.hh>
#include <Geant4/Randomize.hh>

namespace mu2e{
  class ElementSamplerTool{
    public:
      ElementSamplerTool(std::string name){
        _name = name;
        _initialized = false;
      }
      virtual ~ElementSamplerTool(){ /**/ };

      virtual std::string sample_element() = 0;
      std::string Sample(){
        auto rv = this->sample_element();
        return rv;
      }
      std::vector<std::string>& Elements(){
        auto& rv = _elements;
        return rv;
      }
    protected:
      std::string _name;
      G4Material* _material;
      std::vector<std::string> _elements;
      bool _initialized;

      virtual void finish_initialize() = 0;
      void initialize(){
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
      };

    private:
      /**/
  };
}; // namespace mu2e

#endif
