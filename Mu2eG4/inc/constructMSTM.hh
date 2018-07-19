#ifndef Mu2eG4_constructMSTM_hh
#define Mu2eG4_constructMSTM_hh
//
// Free function to create the Detector Solenoid
//
// $Id: constructMSTM.hh,v 1.1 2014/06/05 21:05:05 genser Exp $
// $Author: genser $
// $Date: 2014/06/05 21:05:05 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructMSTM(const VolumeInfo& parent,
                     const SimpleConfig& _config
                     );

}

#endif /* Mu2eG4_constructMSTM_hh */
