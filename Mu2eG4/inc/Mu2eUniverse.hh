#ifndef Mu2eG4_Mu2eUniverse_hh
#define Mu2eG4_Mu2eUniverse_hh
//
// (Pure virtual) Umbrela for the the Mu2e G4 world classes 
//
// $Id: Mu2eUniverse.hh,v 1.2 2012/11/19 23:03:24 genser Exp $
// $Author: genser $
// $Date: 2012/11/19 23:03:24 $
//
// Original author K. Genser to generalize Mu2eWorld
//
// Notes
//


// C++ includes
#include <vector>


// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
#include "GeometryService/inc/GeometryService.hh"

// G4 includes
#include "G4Types.hh"

//G4 forward reference
class G4VPhysicalVolume;

namespace mu2e {

  // Forward references within mu2e namespace.
  class SimpleConfig;

  class Mu2eUniverse {
  public:

    Mu2eUniverse();
    virtual ~Mu2eUniverse();

    // Construct everything.
    // The non-const return type is eventually required 
    // by G4VUserDetectorConstruction::Construct();
    virtual G4VPhysicalVolume * construct() = 0;
      
      virtual void constructSDandField() = 0;

  protected:

    // Utility functions.
    static void setUnits( std::vector<double>& V, G4double unit );

    // A helper function for debugging.  Print a subset of the physical volume store
    static void printPhys();

    // geometry service
    GeometryService const & _geom;

    // Stash a pointer to the config object so that all methods can get at it easily.
    SimpleConfig const & _config; // make it ref?? (some functions need to change before it...

    // Access to the G4HelperService.
    G4Helper * _helper;

    int  _verbosityLevel;
    int  _g4VerbosityLevel; // for non geometry related printouts

  };

} // end namespace mu2e
#endif /* Mu2eG4_Mu2eUniverse_hh */
