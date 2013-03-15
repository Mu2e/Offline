#ifndef BeamlineGeom_BeamlineMaker_hh
#define BeamlineGeom_BeamlineMaker_hh
//
// Construct and return an Beamline.
//
//
// $Id: BeamlineMaker.hh,v 1.6 2013/03/15 15:52:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:03 $
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
    static void BuildBeamline(const SimpleConfig&, Beamline*);
    static void BuildTS(const SimpleConfig&, Beamline*);
  };

}  //namespace mu2e

#endif /* BeamlineGeom_BeamlineMaker_hh */
