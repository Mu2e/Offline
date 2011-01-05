#ifndef constructSteel_HH
#define constructSteel_HH
//
// Free function to create Hall Steel
//
// $Id: constructSteel.hh,v 1.1 2011/01/05 21:04:31 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:31 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructSteel(VolumeInfo   const & parent, 
                      SimpleConfig const * const _config
                      );

}

#endif
