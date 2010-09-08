#ifndef BFieldManagerMaker_HH
#define BFieldManagerMaker_HH
//
// Build a magnetic field manager.
//
// $Id: BFieldManagerMaker.hh,v 1.2 2010/09/08 00:07:27 logash Exp $
// $Author: logash $ 
// $Date: 2010/09/08 00:07:27 $
//

// Includes from C++
#include <memory>
// Includes from CLHEP
#include "CLHEP/Vector/ThreeVector.h"

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

    const SimpleConfig& _config;

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

    // Read a MECO G4BL format map.
    void readG4BLMap( const std::string& filename,
		      BFMap& bfmap, CLHEP::Hep3Vector offset );

    // Special case: when the DS has a uniform field.
    void loadUniformDS();

    // Compute the size of the array needed to hold the raw data of the field map.
    int computeArraySize( int fd, const std::string& filename );

  }; // end class BFieldManagerMaker

} // end namespace mu2e
#endif
