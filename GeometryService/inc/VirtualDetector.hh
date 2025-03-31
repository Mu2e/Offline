#ifndef VirtualDetectorGeom_VirtualDetector_hh
#define VirtualDetectorGeom_VirtualDetector_hh

//
// Class to represent the virtual detectors
//
//

#include <map>
#include <memory>
#include <string>

// Includes from Mu2e
#include "Offline/Mu2eInterfaces/inc/Detector.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class VirtualDetector : virtual public Detector {

  friend class VirtualDetectorMaker;

  public:
    VirtualDetector();

    double getHalfLength() const { return _halfLength; }

    unsigned int nDet() const { return _local.size(); }

    bool exist(int i) const { return _local.find(i) != _local.end(); }

    // none of the following functions are protected agains illegal values: FIXME!!!

    // Get position in the parent frame
    CLHEP::Hep3Vector  const& getLocal(int i) const { return _local.find(i)->second; }

    // Get position in the global Mu2e frame
    CLHEP::Hep3Vector  const& getGlobal(int i) const { return _global.find(i)->second; }

    const CLHEP::HepRotation *  getRotation(int i) const { return _rotation.find(i)->second; }

    // Used to identify Geant volumes
    static std::string volumeName(int i);

    void addVirtualDetector(int id,
                            const CLHEP::Hep3Vector& parentCenterInMu2e,
                            const CLHEP::HepRotation* parentRotationInMu2e,
                            const CLHEP::Hep3Vector& vdCenterInParent);

  protected:
    static const std::string _baseName;

    double _halfLength;

    std::map<int,CLHEP::Hep3Vector> _local;

    std::map<int,CLHEP::Hep3Vector> _global;

    std::map<int,const CLHEP::HepRotation*> _rotation;

    std::map<int,std::string> _name;
};

}
#endif /* VirtualDetectorGeom_VirtualDetector_hh */
