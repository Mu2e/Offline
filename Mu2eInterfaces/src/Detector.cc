//
// A base class for detector components.
//
//
// Original author Rob Kutschke
//

#include "Offline/Mu2eInterfaces/inc/Detector.hh"

namespace mu2e
{
  Detector::~Detector() = default;
  Detector::Detector(const Detector&) = default;
  Detector& Detector::operator=(const Detector &) = default;
}
