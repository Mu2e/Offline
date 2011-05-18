#ifndef Mu2eG4_constructSteel_hh
#define Mu2eG4_constructSteel_hh
//
// Free function to create Hall Steel
//
// $Id: constructSteel.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
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

#endif /* Mu2eG4_constructSteel_hh */
