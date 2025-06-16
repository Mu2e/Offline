#ifndef GeometryService_inc_ProtonAbsorberMaker_hh
#define GeometryService_inc_ProtonAbsorberMaker_hh
//
// Construct and return an ProtonAbsorber.
//
//
//
// Original author MyeongJae Lee
//

#include <memory>
#include <string>
#include <vector>

namespace mu2e {

class ProtonAbsorber;
class SimpleConfig;
class DetectorSolenoid;
class StoppingTarget;

class ProtonAbsorberMaker{

public:

  ProtonAbsorberMaker(SimpleConfig const& _config, const DetectorSolenoid& ds, const StoppingTarget& target);

  void BuildIt ( SimpleConfig const& _config);

  void PrintConfig () ;

  // Use compiler-generated copy c'tor, copy assignment, and d'tor


  // This is the accessor that will remain.
  std::unique_ptr<ProtonAbsorber> getProtonAbsorberPtr() { return std::move(_pabs); }

private:

  void Build (SimpleConfig const & _config);

  // pointer to the Mu2E Geometry ProtonAbsorber being made
  std::unique_ptr<ProtonAbsorber> _pabs;

  const DetectorSolenoid *_ds;
  const StoppingTarget *_target;
  int _IPAVersion;
};

}  //namespace mu2e

#endif /* GeometryService_inc_ProtonAbsorberMaker_hh */
