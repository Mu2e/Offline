#ifndef Mu2eG4_constructDS_hh
#define Mu2eG4_constructDS_hh
//
// Free function to create the Detector Solenoid 
//
// $Id: constructDS.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructDS(VolumeInfo   const & parent, 
                   SimpleConfig const * const _config
                   );

}

#endif /* Mu2eG4_constructDS_hh */
