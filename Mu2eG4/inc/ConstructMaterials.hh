#ifndef Mu2eG4_ConstructMaterials_hh
#define Mu2eG4_ConstructMaterials_hh
//
// Construct materials requested by the run-time configuration system.
//
// $Id: ConstructMaterials.hh,v 1.6 2014/09/04 14:22:45 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/04 14:22:45 $
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

namespace fhicl { class ParameterSet; }

namespace mu2e {

  // Return type of the isUnique() method.
  class CheckedG4String{
  public:
    CheckedG4String ()                 : name() {}
    CheckedG4String (G4String const& n): name(n){}
    G4String name;
  };

  // Forward references within mu2e namespace.
  class SimpleConfig;

  class ConstructMaterials{
  public:

    ConstructMaterials();
    explicit ConstructMaterials(const fhicl::ParameterSet& pset);

    ~ConstructMaterials();

    // Construct all of the materials.
    void construct();

  private:
    bool mu2eStandardDetector_;
    bool printElements_;
    bool printMaterials_;

    // Construct materials specific to Mu2e.
    void constructMu2eMaterials();  // Don't add any more materials here.
    void constructMu2eMaterials2(); // Because of the number of variables
    // created to hold the many materials used, we move some of the
    // material creation into a new function.  Add new materials to this function.
    // Wrapper around FindOrBuildElement.
    G4Element* getElementOrThrow( G4String const& name);

    // Check that a material name is not already in use.
    CheckedG4String uniqueMaterialOrThrow( G4String const& name);

  };

} // end namespace mu2e
#endif /* Mu2eG4_ConstructMaterials_hh */
