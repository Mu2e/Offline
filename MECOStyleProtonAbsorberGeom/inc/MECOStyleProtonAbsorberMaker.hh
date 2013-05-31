#ifndef MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh
#define MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh
//
// Construct and return an MECOStyleProtonAbsorber.
//
//
// $Id: MECOStyleProtonAbsorberMaker.hh,v 1.4 2013/05/31 18:07:18 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 18:07:18 $
//
// Original author MyeongJae Lee
//

#include <memory>
#include <string>
#include <vector>

namespace mu2e {

class MECOStyleProtonAbsorber;
class SimpleConfig;
class DetectorSolenoid;
class Target;

class MECOStyleProtonAbsorberMaker{

public:

  MECOStyleProtonAbsorberMaker(SimpleConfig const& _config, const DetectorSolenoid& ds, const Target& target);

  void BuildIt ( SimpleConfig const& _config);

  void PrintConfig () ;

  // Use compiler-generated copy c'tor, copy assignment, and d'tor


  // This is the accessor that will remain.
  std::unique_ptr<MECOStyleProtonAbsorber> getMECOStyleProtonAbsorberPtr() { return std::move(_pabs); }

private:

  void Build (SimpleConfig const & _config);

  // pointer to the Mu2E Geometry MECOStyleProtonAbsorber being made
  std::unique_ptr<MECOStyleProtonAbsorber> _pabs;

  const DetectorSolenoid *_ds;
  const Target *_target;

};

}  //namespace mu2e

#endif /* MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh */
