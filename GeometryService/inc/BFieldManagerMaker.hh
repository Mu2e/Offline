#ifndef BFieldGeom_BFieldManagerMaker_hh
#define BFieldGeom_BFieldManagerMaker_hh
//
// Build a magnetic field manager.
//
//

// Includes from C++
#include <memory>

// Includes from CLHEP
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "Offline/BFieldGeom/inc/BFGridMap.hh"
#include "Offline/BFieldGeom/inc/BFInterpolationStyle.hh"
#include "Offline/BFieldGeom/inc/BFMap.hh"
#include "Offline/BFieldGeom/inc/BFParamMap.hh"
#include "Offline/BFieldGeom/inc/BFieldConfig.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"

namespace mu2e {

    // Forward reference.
    class BFMap;

    class BFieldManagerMaker {
       public:

        typedef std::vector<std::shared_ptr<BFMap>> MapContainerType;

        explicit BFieldManagerMaker(const BFieldConfig& config);

        // Transfer ownership of the BFManager.
        std::unique_ptr<BFieldManager> getBFieldManager() { return std::move(_bfmgr); }

       private:
        // Helper object to turn filenames into full paths using MU2E_SEARCH_PATH.
        ConfigFileLookupPolicy _resolveFullPath;

        // Hold the object while we are creating it. The GeometryService will take ownership.
        std::unique_ptr<BFieldManager> _bfmgr;

        // Hold the types of the inner and outer maps (if they differ)
        std::vector<BFMapType> _innerTypes;
        std::vector<BFMapType> _outerTypes;

        // Load a series of parametric magnetic field maps.
        void loadParam(MapContainerType& whichMap,
                       const BFieldConfig::FileSequenceType& files,
                       std::vector<BFMapType> mapTypeList,
                       BFInterpolationStyle interpStyle,
                       double scaleFactor,
                       bool flipBFieldMap);

        // Create a new parametric magnetic field map, get information from config file.
        void loadParam(MapContainerType& whichMap,
                       const std::string& key,
                       const std::string& resolvedFileName,
                       double scaleFactor);

        // Load a series of G4BL magnetic field maps.
        void loadG4BL(MapContainerType& whichMap,
                      const BFieldConfig::FileSequenceType& files,
                      double scaleFactor,
                      BFInterpolationStyle interpStyle,
                      bool flipBFieldMap);

        // Create a new magnetic field map, get information from config file.
        void loadG4BL(MapContainerType& whichMap,
                      const std::string& key,
                      const std::string& resolvedFileName,
                      double scaleFactor,
                      BFInterpolationStyle interpStyle,
                      bool flipBFieldMap);

        // Read a G4BL text format map.
        void readG4BLMap(const std::string& filename, BFGridMap& bfmap, CLHEP::Hep3Vector offset);

        // Read a G4BL map that was stored using writeG4BLBinary.
        void readG4BLBinary(const std::string& headerFilename, BFGridMap& bfmap);

        // Read a CSV with values for parametric map.
        void readParamFile(const std::string& filename, BFParamMap& bfmap);

        // Write an existing BFMap in binary format.
        void writeG4BLBinary(const BFGridMap& bf, const std::string& outputfile);
        void flipMap(BFGridMap& bf);

    };  // end class BFieldManagerMaker

}  // end namespace mu2e
#endif /* BFieldGeom_BFieldManagerMaker_hh */
