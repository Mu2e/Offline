#ifndef VirtualDetector_HH
#define VirtualDetector_HH

//
// Class to represent the virtual detectors
//
#include <map>
#include <string>
#include <memory>

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
    ~VirtualDetector(){;};

    virtual std::string name() const { return "VirtualDetector";};
    
    double getHalfLength() const { return _halfLength; };

    unsigned int nDet() const { return _local.size(); }; 

    bool exist(int i) const { return _local.find(i) != _local.end(); };

    // Get position in the parent frame
    CLHEP::Hep3Vector  const& getLocal(int i) const { return _local.find(i)->second; };

    // Get position in the global Mu2e frame
    CLHEP::Hep3Vector  const& getGlobal(int i) const { return _global.find(i)->second; };

    CLHEP::HepRotation *  getRotation(int i) const { return _rotation.find(i)->second; };

    void addVirtualDetector(int, std::string,CLHEP::Hep3Vector,CLHEP::HepRotation*,CLHEP::Hep3Vector);

  protected:

    double _halfLength;

    std::map<int,CLHEP::Hep3Vector> _local;

    std::map<int,CLHEP::Hep3Vector> _global;

    std::map<int,CLHEP::HepRotation*> _rotation;

    std::map<int,std::string> _name;
};

}
#endif
