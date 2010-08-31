#ifndef BEAMLINEMAKER_HH
#define BEAMLINEMAKER_HH
//
// Construct and return an Beamline.
//
//
// $Id: BeamlineMaker.hh,v 1.1 2010/08/31 00:24:51 logash Exp $
// $Author: logash $ 
// $Date: 2010/08/31 00:24:51 $
//
// Original author Peter Shanahan
//

#include <vector>
#include <string>
#include <memory>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

class Beamline;
class TransportSolenoid;
class SimpleConfig;

class BeamlineMaker {

public:

  BeamlineMaker( SimpleConfig const& config );  

  ~BeamlineMaker ();

  // This is depracted and will go away soon.  
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

#endif //BEAMLINEMAKER_HH
