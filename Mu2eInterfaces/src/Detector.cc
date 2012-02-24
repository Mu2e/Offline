//
// A base class for detector components.
//
// $Id: Detector.cc,v 1.1 2012/02/24 16:36:36 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 16:36:36 $
//
// Original author Rob Kutschke
//

#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e
{
  Detector::~Detector() { }
  void Detector::update() { }
  std::string Detector::name() const { return "NONAME"; }
}
