#ifndef BeamlineGeom_Beamline_hh
#define BeamlineGeom_Beamline_hh

//
// Class to represent the transport solenoid
//
// Includes from Mu2e
#include "Mu2eInterfaces/inc/Detector.hh"
#include "BeamlineGeom/inc/TransportSolenoid.hh"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class Beamline : virtual public Detector {

  friend class BeamlineMaker;

  public:
    Beamline():
      Detector(),
      _solenoidOffset(0.0),
      _ts(){
    }

    // use compiler-generated copy c'tor, copy assignment, and d'tor

    double solenoidOffset() const { return _solenoidOffset; };

    TransportSolenoid const& getTS() const { return _ts; };

  protected:

    // All dimensions in mm.

    double _solenoidOffset;

    // Beamline subsystems
    TransportSolenoid _ts;

};

}
#endif /* BeamlineGeom_Beamline_hh */
