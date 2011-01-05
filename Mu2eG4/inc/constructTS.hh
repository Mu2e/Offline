#ifndef constructTS_HH
#define constructTS_HH
//
// Free function to create  Transport Solenoid
//
// $Id: constructTS.hh,v 1.1 2011/01/05 21:04:31 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:31 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructTS(VolumeInfo   const & parent, 
                   SimpleConfig const * const _config
                   );

}

#endif
