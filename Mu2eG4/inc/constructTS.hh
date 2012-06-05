#ifndef Mu2eG4_constructTS_hh
#define Mu2eG4_constructTS_hh
//
// Free function to create  Transport Solenoid
//
// $Id: constructTS.hh,v 1.4 2012/06/05 16:19:33 genser Exp $
// $Author: genser $
// $Date: 2012/06/05 16:19:33 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructTS(VolumeInfo   const & parent,
                   SimpleConfig const & _config
                   );

}

#endif /* Mu2eG4_constructTS_hh */
