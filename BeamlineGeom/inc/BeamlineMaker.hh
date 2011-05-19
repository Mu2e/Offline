#ifndef BeamlineGeom_BeamlineMaker_hh
#define BeamlineGeom_BeamlineMaker_hh
//
// Construct and return an Beamline.
//
//
// $Id: BeamlineMaker.hh,v 1.4 2011/05/19 18:47:25 wb Exp $
// $Author: wb $
// $Date: 2011/05/19 18:47:25 $
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

  BeamlineMaker( SimpleConfig const& config );

  // use compiler-generated copy c'tor, copy assignment, and d'tor

  // This is deprecated and will go away soon.
  // Still needed for root graphics version.
  const Beamline& getBeamline() const { return *_beamline;}

  // This is the accessor that will remain.
  std::auto_ptr<Beamline> getBeamlinePtr() { return _beamline; }

private:

  void BuildBeamline(SimpleConfig const&);
  void BuildTS(SimpleConfig const&);

  //SimpleConfig const& _config;

  // pointer to the Mu2E Geometry Beamline being made
  std::auto_ptr<Beamline> _beamline;

  // Read these variables from config file, data base, etc.

};

}  //namespace mu2e

#endif /* BeamlineGeom_BeamlineMaker_hh */
