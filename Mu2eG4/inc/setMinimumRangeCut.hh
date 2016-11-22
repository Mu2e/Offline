#ifndef Mu2eG4_setMinimumRangeCut_hh
#define Mu2eG4_setMinimumRangeCut_hh
//
// Set the G4 minimum range cut as specified in the geometry file.
//
// $Id: setMinimumRangeCut.hh,v 1.1 2012/06/04 19:28:01 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/06/04 19:28:01 $
//
//-----------------------------------------------------------------------------

namespace fhicl { class ParameterSet; }

class G4VUserPhysicsList;

namespace mu2e{

  class SimpleConfig;

  void setMinimumRangeCut( SimpleConfig const& config, G4VUserPhysicsList* mPL );
  void setMinimumRangeCut(const fhicl::ParameterSet& pset, G4VUserPhysicsList* mPL);

}  // end namespace mu2e

#endif /* Mu2eG4_setMinimumRangeCut_hh */
