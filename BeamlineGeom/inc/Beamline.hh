#ifndef Beamline_HH
#define Beamline_HH

//
// Class to represent the transport solenoid
//
#include <memory>

// Includes from Mu2e
#include "GeometryService/inc/Detector.hh"
#include "BeamlineGeom/inc/TransportSolenoid.hh"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class Beamline : public Detector {

  friend class BeamlineMaker;

  public:
    Beamline();
    ~Beamline(){;};

    virtual std::string name() const { return "Beamline";}
    
    double solenoidOffset() const { return _solenoidOffset; };

    TransportSolenoid const& getTS() const { return _ts; };

  protected:

    // All dimensions in mm.

    double _solenoidOffset;

    // Beamline subsystems
    TransportSolenoid _ts;

};

}
#endif
