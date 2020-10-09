#ifndef Mu2eG4_constructDS_hh
#define Mu2eG4_constructDS_hh
//
// Free function to create the Detector Solenoid
//
//
// Original author KLG
//

class G4Colour;

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;
  class TubsParams;

  void constructDS(const VolumeInfo& parent,
                   const SimpleConfig& config
                   );

  // limited utility function
  void placeTubeCore ( const std::string & name,
                       double radiusFract,
                       double radiusDFract,
                       double dPhiFraction,
                       const std::string & material,
                       const G4Colour & color,
                       const VolumeInfo& parent,
                       const TubsParams & parentParams,
                       const std::string & lookupToken,
                       const SimpleConfig & config, // to be removed?
		       const int zNotPhi = 0
                       );

  TubsParams calculateTubeCoreParams (const TubsParams& parentParams,
                                      double radiusFract,
                                      double radiusDFract,
                                      double dPhiFraction,
                                      int verbosityLevel=0,
				      const int zNotPhi = 0
				      );

}

#endif /* Mu2eG4_constructDS_hh */
