#ifndef VirtualDetectorGeom_VirtualDetector_hh
#define VirtualDetectorGeom_VirtualDetector_hh

//
// Class to represent the virtual detectors
//
// $Id: VirtualDetector.hh,v 1.7 2011/12/14 00:30:01 gandr Exp $
// $Author: gandr $
//

#include <map>
#include <memory>
#include <string>

// Includes from Mu2e
#include "GeometryService/inc/Detector.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class VirtualDetector : public Detector {

  friend class VirtualDetectorMaker;

  public:
    VirtualDetector();
    ~VirtualDetector(){;}

    virtual std::string name() const { return "VirtualDetector";}

    double getHalfLength() const { return _halfLength; }

    unsigned int nDet() const { return _local.size(); }

    bool exist(int i) const { return _local.find(i) != _local.end(); }

    // none of the following functions are protected agains illegal values: FIXME!!!

    // Get position in the parent frame
    CLHEP::Hep3Vector  const& getLocal(int i) const { return _local.find(i)->second; }

    // Get position in the global Mu2e frame
    CLHEP::Hep3Vector  const& getGlobal(int i) const { return _global.find(i)->second; }

    CLHEP::HepRotation *  getRotation(int i) const { return _rotation.find(i)->second; }

    std::string const& name(int i) const { return _name.find(i)->second;}

    void addVirtualDetector(int id,
			    const std::string& name,
			    const CLHEP::Hep3Vector& parentCenterInMu2e,
			    CLHEP::HepRotation* parentRotationInMu2e,
			    const CLHEP::Hep3Vector& vdCenterInParent);

  protected:

    double _halfLength;

    std::map<int,CLHEP::Hep3Vector> _local;

    std::map<int,CLHEP::Hep3Vector> _global;

    std::map<int,CLHEP::HepRotation*> _rotation;

    std::map<int,std::string> _name;
};

}
#endif /* VirtualDetectorGeom_VirtualDetector_hh */
