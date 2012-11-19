#ifndef Mu2eG4_constructSteel_hh
#define Mu2eG4_constructSteel_hh
//
// Free function to create Hall Steel
//
// $Id: constructSteel.hh,v 1.4 2012/11/19 23:03:24 genser Exp $
// $Author: genser $
// $Date: 2012/11/19 23:03:24 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructSteel(VolumeInfo   const & parent,
                      SimpleConfig const & _config
                      );

}

#endif /* Mu2eG4_constructSteel_hh */
