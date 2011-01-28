#ifndef ConstructMaterials_H
#define ConstructMaterials_H 1
//
// Construct materials requested by the run-time configuration system.
//
// $Id: ConstructMaterials.hh,v 1.2 2011/01/28 23:51:58 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2011/01/28 23:51:58 $
//
// Original author Rob Kutschke
//

#include <string>
#include <memory>

// Forward references in global namespace.
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Element;

#include "G4String.hh"

namespace mu2e {

  // Return type of the isNeeded() method.
  class CheckedG4String{
  public:
    CheckedG4String ():
      doit(false),
      name(){
    }
    CheckedG4String ( bool b, G4String const& n):
      doit(b),
      name(n){
    }
    bool doit;
    G4String name;
  };

  // Forward references within mu2e namespace.
  class SimpleConfig;

  class ConstructMaterials{
  public:
    
    ConstructMaterials();
    ~ConstructMaterials();

    // Construct all of the materials.
    void construct();
    
  private:


    // Construct predefined G4 materials.
    void constructG4Materials( SimpleConfig const& config );

    // Construct some additional materials specific to Mu2e.
    void constructMu2eMaterials( SimpleConfig const& config );

    // Wrapper around FindOrBuildElement.
    G4Element* getElementOrThrow( G4String const& name);

    // Do we need to build this material?
    CheckedG4String isNeeded(  std::vector<std::string> const& V, std::string const& name);
    
    // Check that a material name is not already in use.
    void uniqueMaterialOrThrow( G4String const& name);
    
    // Check to see if a string appear within a vector of strings.
    bool isRequested( std::vector<std::string> const& V, std::string const& name);
    
  };

} // end namespace mu2e
#endif

