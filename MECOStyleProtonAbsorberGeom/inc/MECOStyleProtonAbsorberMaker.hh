#ifndef MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh
#define MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh
//
// Construct and return an MECOStyleProtonAbsorber.
//
//
// $Id: MECOStyleProtonAbsorberMaker.hh,v 1.3 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
//
// Original author MyeongJae Lee
//

#include <memory>
#include <string>
#include <vector>

namespace mu2e {

class MECOStyleProtonAbsorber;
class SimpleConfig;

class MECOStyleProtonAbsorberMaker{

public:

  MECOStyleProtonAbsorberMaker( SimpleConfig const& _config );

  void BuildIt ( SimpleConfig const& _config);

  void PrintConfig () ;

  // Use compiler-generated copy c'tor, copy assignment, and d'tor


  // This is the accessor that will remain.
  std::unique_ptr<MECOStyleProtonAbsorber> getMECOStyleProtonAbsorberPtr() { return std::move(_pabs); }

private:

  void Build (SimpleConfig const & _config);

  // pointer to the Mu2E Geometry MECOStyleProtonAbsorber being made
  std::unique_ptr<MECOStyleProtonAbsorber> _pabs;

  
  
};

}  //namespace mu2e

#endif /* MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh */
