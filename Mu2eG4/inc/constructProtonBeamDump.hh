// Andrei Gaponenko, 2011

#ifndef CONSTRUCTPROTONBEAMDUMP_HH
#define CONSTRUCTPROTONBEAMDUMP_HH

namespace mu2e {
  
  class VolumeInfo;
  class SimpleConfig;
  class SensitiveDetectorHelper;
  
  void constructProtonBeamDump(const VolumeInfo& parent,
			       const SimpleConfig& config,
                               const SensitiveDetectorHelper& sdHelper
			       );
}

#endif /* CONSTRUCTPROTONBEAMDUMP_HH */
