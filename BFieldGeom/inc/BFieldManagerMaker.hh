#ifndef BFieldGeom_BFieldManagerMaker_hh
#define BFieldGeom_BFieldManagerMaker_hh
//
// Build a magnetic field manager.
//
// $Id: BFieldManagerMaker.hh,v 1.9 2012/02/16 04:59:38 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/16 04:59:38 $
//

// Includes from C++
#include <memory>

// Includes from CLHEP
#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eUtilities/inc/ConfigFileLookupPolicy.hh"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;
  class BFieldManager;
  class BFMap;

  class BFieldManagerMaker{

  public:

    friend class BFieldManagerMakerMaker;

    BFieldManagerMaker( const SimpleConfig& config );
    ~BFieldManagerMaker();

    // Transfer ownership of the BFManager.
    std::auto_ptr<BFieldManager> getBFieldManager() { return _bfmgr; }

  private:

    // Helper object to turn filenames into full paths using MU2E_SEARCH_PATH.
    ConfigFileLookupPolicy _resolveFullPath;

    // The runtime configuration for the geometry subsystem.
    const SimpleConfig& _config;

    // Hold the object while we are creating it. The GeometryService will take ownership.
    std::auto_ptr<BFieldManager> _bfmgr;

    // Create a new magnetic field map, get information from config file.
    void loadGMC( const std::string& key,
                  const std::string& fileKey,
                  const std::string& dimensionKey );

    // Create a new magnetic field map, get information from config file.
    void loadG4BL( const std::string& key,
                   const std::string& fileKey );

    // Read a MECO GMC format map.
    void readGMCMap( const std::string& filename,
                     BFMap& bfmap );

    // Read a G4BL text format map.
    void readG4BLMap( const std::string& filename,
                      BFMap& bfmap, CLHEP::Hep3Vector offset );

    // Read a G4BL map that was stored using writeG4BLBinary.
    void readG4BLBinary( const std::string& headerFilename,
                         BFMap& bfmap );

    // Write an existing BFMap in binary format.
    void writeG4BLBinary(const std::string& key);

    // Special case: when the DS has a uniform field.
    void loadUniformDS();

    // Compute the size of the array needed to hold the raw data of the field map.
    int computeArraySize( int fd, const std::string& filename );

  }; // end class BFieldManagerMaker

} // end namespace mu2e
#endif /* BFieldGeom_BFieldManagerMaker_hh */
