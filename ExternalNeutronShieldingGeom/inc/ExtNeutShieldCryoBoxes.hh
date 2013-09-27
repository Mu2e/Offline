#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldCryoBoxes_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldCryoBoxes_hh

//
//
// Original author David Norvil Brown
//

// The PS External Shield is extruded, upside-down "u" shape, made of 
// concrete.  

#include <vector>
#include <ostream>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {

  class ExtNeutShieldCryoBoxesMaker;

  class ExtNeutShieldCryoBoxes : virtual public Detector {

  public:

    // Use a vector of Hep2Vectors for the corners of the shape to be extruded
    const std::vector<std::vector<double>>& getDimensions() const { return _dimensions; }
    const std::vector<std::string>& materialNames() const { return _materialNames; }
    const std::vector<CLHEP::Hep3Vector> centersOfBoxes() const { return _centerPositions; }

  private:

    friend class ExtNeutShieldCryoBoxesMaker;

    // Private ctr: the class should only be constructed via ExtNeutShieldCryoBoxes::ExtNeutShieldCryoBoxesMaker.
    ExtNeutShieldCryoBoxes(const std::vector<std::vector<double>>& dims, const std::vector<std::string>& mats, const std::vector<CLHEP::Hep3Vector> sites )
      : _dimensions(dims),
        _materialNames(mats),
	_centerPositions(sites)
    { }

    // Or read back from persistent storage
    ExtNeutShieldCryoBoxes();
    template<class T> friend class art::Wrapper;

    // Current description based on Geometry 13 from G4Beamline, adapted by
    // D. Norvil Brown, 


    std::vector<std::vector<double>> _dimensions;
    std::vector<std::string> _materialNames;
    std::vector<CLHEP::Hep3Vector> _centerPositions;

  };

  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldCryoBoxes& enscb);

}

#endif/*ExternalNeutronShieldingGeom_ExtNeutShieldCryoBoxes_hh*/
