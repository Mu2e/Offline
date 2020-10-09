#ifndef GeometryService_StoppingTargetMaker_hh
#define GeometryService_StoppingTargetMaker_hh
//
// Construct and return an Target.
//
// Original author Peter Shanahan
//

#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

class StoppingTarget;
class SimpleConfig;

class StoppingTargetMaker{

public:

  StoppingTargetMaker(const CLHEP::Hep3Vector& detSysOrigin, SimpleConfig const& config );

  // Use compiler-generated copy c'tor, copy assignment, and d'tor

  // This is the accessor that will remain.
  std::unique_ptr<StoppingTarget> getTargetPtr() { return std::move(_targ); }

private:

  void BuildIt();
  void PrintConfig();// prints the target configuration as TargetMaker
                     // understands it...

  //  SimpleConfig const& _config;

  // pointer to the Mu2E Geometry Target being made
  std::unique_ptr<StoppingTarget> _targ;

  // variables needed to build the Target.  Read in from config file,
  // data base, etc.  These (and the TargetMaker object itself) only need
  // to persist long enough to make the Target.  After that, the definition
  // resides entirely in the Target object.

  CLHEP::Hep3Vector _detSysOrigin;

  double _z0InMu2e; // nominal center of foils.
  double _deltaZ ; // nominal spacing of foils.
  double _rIn; // inner radius of foils.  Currently 0.

  std::vector<double> _rOut; //outer Radii of foils
  std::vector<double> _halfThicknesses; //half thicknesses of foils
  std::vector<double> _xVars; //variation of x position from axis of foils
  std::vector<double> _yVars; //variation of y position from axis of foils
  std::vector<double> _zVars; //variation of z position from nominal of foils
  std::vector<double> _xCos; // x directional cosines of foils
  std::vector<double> _yCos; // y directional cosines of foils
  std::vector<std::string> _materials; //material of foils
  std::string _fillMaterial; // material of enclosing cylinder

  // stopping foil target supporting structure
  bool   _foilTarget_supportStructure;
  int    _foilTarget_supportStructure_number; // number of support wires per foil
  double _foilTarget_supportStructure_angleOffset; // angle of first support wire wrt vertical
  double _foilTarget_supportStructure_radius; // radius of the support wires
//  double _foilTarget_supportStructure_length; // length of the support wires
  std::string _foilTarget_supportStructure_fillMaterial;
  double _foilTarget_supportStructure_rOut; //radius the structure goes out to
};

}  //namespace mu2e

#endif /* GeometryService_StoppingTargetMaker_hh */
