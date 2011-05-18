//
// A base class for detector components.
//
// $Id: Detector.cc,v 1.2 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
//
// Original author Rob Kutschke
//

#include "GeometryService/inc/Detector.hh"

namespace mu2e
{
  Detector::~Detector() { }
  void Detector::update() { }
  std::string Detector::name() const { return "NONAME"; }
}
