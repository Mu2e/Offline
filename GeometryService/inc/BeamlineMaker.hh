#ifndef BeamlineGeom_BeamlineMaker_hh
#define BeamlineGeom_BeamlineMaker_hh
//
// Construct and return an Beamline.
//
//
//
// Original author Peter Shanahan
//

#include <memory>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

  class Beamline;
  class SimpleConfig;
  class TransportSolenoid;

  class BeamlineMaker {

  public:
    static std::unique_ptr<Beamline> make(const SimpleConfig& config);

  private:
    static void BuildBeamline  (const SimpleConfig&, Beamline*);
    static void BuildTSCryostat(const SimpleConfig&, Beamline*);
    static void BuildTSCoils   (const SimpleConfig&, Beamline*);
    static void BuildTSCollimators(const SimpleConfig&, TransportSolenoid* );
    static void BuildTSVacua      (const SimpleConfig&, TransportSolenoid* );
    static void BuildTSCAs        (const SimpleConfig&, Beamline& );
    static void BuildTSPolyLining (const SimpleConfig&, TransportSolenoid* );
    static void BuildPbarWindow   (const SimpleConfig&, TransportSolenoid* );
  };

}  //namespace mu2e

#endif /* BeamlineGeom_BeamlineMaker_hh */
