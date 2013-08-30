#ifndef BFieldGeom_BFieldManagerMaker_hh
#define BFieldGeom_BFieldManagerMaker_hh
//
// Build a magnetic field manager.
//
// $Id: BFieldManagerMaker.hh,v 1.16 2013/08/30 22:25:22 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/08/30 22:25:22 $
//

// Includes from C++
#include <memory>

// Includes from CLHEP
#include "CLHEP/Vector/ThreeVector.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "BFieldGeom/inc/BFieldConfig.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "BFieldGeom/inc/BFInterpolationStyle.hh"

namespace mu2e {

  // Forward reference.
  class BFMap;

  class BFieldManagerMaker {

  public:

    friend class BFieldManagerMakerMaker;

    explicit BFieldManagerMaker(const BFieldConfig& config);

    // Transfer ownership of the BFManager.
    std::unique_ptr<BFieldManager> getBFieldManager() { return std::move(_bfmgr); }

  private:

    // Helper object to turn filenames into full paths using MU2E_SEARCH_PATH.
    ConfigFileLookupPolicy _resolveFullPath;

    // Hold the object while we are creating it. The GeometryService will take ownership.
    std::unique_ptr<BFieldManager> _bfmgr;

    void loadG4BL(BFieldManager::MapContainerType *whichMap,
                  const BFieldConfig::FileSequenceType& files,
                  double scaleFactor,
                  BFInterpolationStyle interpStyle);

    // Create a new magnetic field map, get information from config file.
    void loadG4BL(BFieldManager::MapContainerType *whichMap,
                  const std::string& key,
                  const std::string& resolvedFileName,
                  double scaleFactor,
                  BFInterpolationStyle interpStyle);

    // Create and fill a new magnetic field map
    void readGMCMap(const std::string& mapKey,
                    const std::string& resolvedFileName,
                    const std::vector<int>& dim,
                    double scaleFactor,
                    BFInterpolationStyle interpStyle
                    );

    // Read a G4BL text format map.
    void readG4BLMap( const std::string& filename,
                      BFMap& bfmap, CLHEP::Hep3Vector offset );

    // Read a G4BL map that was stored using writeG4BLBinary.
    void readG4BLBinary( const std::string& headerFilename,
                         BFMap& bfmap );

    // Write an existing BFMap in binary format.
    void writeG4BLBinary(const BFMap& bf, const std::string& outputfile);

    // Compute the size of the array needed to hold the raw data of the field map.
    int computeArraySize( int fd, const std::string& filename );

    void flipMap(BFMap *bf);

  }; // end class BFieldManagerMaker

} // end namespace mu2e
#endif /* BFieldGeom_BFieldManagerMaker_hh */
