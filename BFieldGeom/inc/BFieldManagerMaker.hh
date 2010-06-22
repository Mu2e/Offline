#ifndef BFieldManagerMaker_HH
#define BFieldManagerMaker_HH
//
// Build a magnetic field manager.
//
// $Id: BFieldManagerMaker.hh,v 1.1 2010/06/22 16:44:25 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/06/22 16:44:25 $
//

// Includes from C++
#include <memory>

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

    // Read a MECO GMC format map.
    void readGMCMap( const std::string& filename,
                     BFMap& bfmap );

    // Special case: when the DS has a uniform field.
    void loadUniformDS();

    // Compute the size of the array needed to hold the raw data of the field map.
    int computeArraySize( int fd, const std::string& filename );

  }; // end class BFieldManagerMaker

} // end namespace mu2e
#endif
