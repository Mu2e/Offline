#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldUpstream1b_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldUpstream1b_hh

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

  class ExtNeutShieldUpstream1bMaker;

  class ExtNeutShieldUpstream1b : virtual public Detector {

  public:

    // Use a vector of Hep2Vectors for the corners of the shape to be extruded
    const std::vector<CLHEP::Hep2Vector>& externalShieldOutline() const { return _externalShieldOutline; }
    const double& getLength() const { return _length; }
    const std::string materialName() const { return _materialName; }
    const CLHEP::Hep3Vector centerOfShield() const { return _centerPosition; }
    const double& getRotation() const { return _rotation; }

  private:

    friend class ExtNeutShieldUpstream1bMaker;

    // Private ctr: the class should only be constructed via ExtNeutShieldUpstream1b::ExtNeutShieldUpstream1bMaker.
    ExtNeutShieldUpstream1b(const std::vector<CLHEP::Hep2Vector>& shape, const double& leng, const std::string mat, const CLHEP::Hep3Vector site, const double& rot)
      : _externalShieldOutline(shape),
	_length(leng), _materialName(mat),
	_centerPosition(site),
	_rotation(rot)
    { }

    // Or read back from persistent storage
    ExtNeutShieldUpstream1b();
    template<class T> friend class art::Wrapper;

    // Current description based on Geometry 13 from G4Beamline, adapted by
    // D. Norvil Brown, 


    std::vector<CLHEP::Hep2Vector> _externalShieldOutline;
    double _length;
    std::string _materialName;
    CLHEP::Hep3Vector _centerPosition;
    double _rotation;

  };

  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldUpstream1b& ensu1b);

}

#endif/*ExternalNeutronShieldingGeom_ExtNeutShieldUpstream1b_hh*/
