#ifndef MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh
#define MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh
//
// Construct and return an MECOStyleProtonAbsorber.
//
//
// $Id: MECOStyleProtonAbsorberMaker.hh,v 1.1 2012/05/10 23:40:59 mjlee Exp $
// $Author: mjlee $
// $Date: 2012/05/10 23:40:59 $
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
  std::auto_ptr<MECOStyleProtonAbsorber> getMECOStyleProtonAbsorberPtr() { return _pabs; }

private:

  // pointer to the Mu2E Geometry MECOStyleProtonAbsorber being made
  std::auto_ptr<MECOStyleProtonAbsorber> _pabs;

  
  
};

}  //namespace mu2e

#endif /* MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh */
