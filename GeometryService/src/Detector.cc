//
// A base class for detector components.
//
// $Id: Detector.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
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
