//
// This class transfers user inputs from a SimpleConfig text file to a
// BFieldConfig object.
//
// Andrei Gaponenko, 2012, cut-and-pasting code from BFieldManagerMaker.

#ifndef BFieldGeom_BFieldConfigMaker_hh
#define BFieldGeom_BFieldConfigMaker_hh

#include <memory>

// Forward declaration is not sufficient b/c of potential ~unique_ptr in a client code.
#include "BFieldGeom/inc/BFieldConfig.hh"

namespace mu2e {

  class SimpleConfig;
  class Beamline;

  class BFieldConfigMaker {

  public:

    explicit BFieldConfigMaker(const SimpleConfig& config, const Beamline& beamg);

    // Transfer ownership of the BFManager.
    std::unique_ptr<BFieldConfig> getBFieldConfig() { return std::move(bfconf_); }

  private:
    // Hold the object while we are creating it. The GeometryService will take ownership.
    std::unique_ptr<BFieldConfig> bfconf_;

    void addGMC(const SimpleConfig& config, BFieldConfig* bfconf, const std::string& configFNKey, const std::string& configDimKey);
  };

} // namespace mu2e

#endif /* BFieldGeom_BFieldConfigMaker_hh */
