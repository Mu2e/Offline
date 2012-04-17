#ifndef Mu2eG4_constructDirt_hh
#define Mu2eG4_constructDirt_hh
//
// Free function to create the earthen overburden.
//
// $Id: constructDirt.hh,v 1.4 2012/04/17 19:56:56 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/17 19:56:56 $
//
// Original author KLG
//

namespace mu2e {

  class SimpleConfig;
  class VolumeInfo;

  void constructDirt(const VolumeInfo& parent, const SimpleConfig& config);
}

#endif /* Mu2eG4_constructDirt_hh */
