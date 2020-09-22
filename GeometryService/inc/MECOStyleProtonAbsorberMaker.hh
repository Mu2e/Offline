#ifndef MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh
#define MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh
//
// Construct and return an MECOStyleProtonAbsorber.
//
//
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
class StoppingTarget;

class MECOStyleProtonAbsorberMaker{

public:

  MECOStyleProtonAbsorberMaker(SimpleConfig const& _config, const DetectorSolenoid& ds, const StoppingTarget& target);

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
  const StoppingTarget *_target;
  int _IPAVersion;
};

}  //namespace mu2e

#endif /* MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberMaker_hh */
