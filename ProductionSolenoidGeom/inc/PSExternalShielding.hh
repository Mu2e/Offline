#ifndef ProductionSolenoidGeom_PSExternalShielding_hh
#define ProductionSolenoidGeom_PSExternalShielding_hh

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

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class PSExternalShieldingMaker;

  class PSExternalShielding : virtual public Detector {

  public:

    // Use a vector of Hep2Vectors for the corners of the shape to be extruded
    const std::vector<CLHEP::Hep2Vector>& externalShieldOutline() const { return _externalShieldOutline; }
    const double& getLength() const { return _length; }
    const std::string materialName() const { return _materialName; }
    const CLHEP::Hep3Vector centerOfShield() const { return _centerPosition; }

  private:

    friend class PSExternalShieldingMaker;

    // Private ctr: the class should only be constructed via PSExternalShielding::PSExternalShieldingMaker.
    PSExternalShielding(const std::vector<CLHEP::Hep2Vector>& shape, const double& leng, const std::string mat, const CLHEP::Hep3Vector site)
      : _externalShieldOutline(shape),
	_length(leng), _materialName(mat),
	_centerPosition(site)
    { }

    // Or read back from persistent storage
    PSExternalShielding();
    template<class T> friend class art::Wrapper;

    // Current description based on private email Rick Coleman -> 
    // D. Norvil Brown, describing the implementation of the external shield
    // in G4BeamLine.  

    std::vector<CLHEP::Hep2Vector> _externalShieldOutline;
    double _length;
    std::string _materialName;
    CLHEP::Hep3Vector _centerPosition;

  };

  std::ostream& operator<<(std::ostream& os, const PSExternalShielding& pse);

}

#endif/*ProductionSolenoidGeom_PSExternalShielding_hh*/
