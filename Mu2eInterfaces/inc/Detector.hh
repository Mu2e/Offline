#ifndef GeometryService_Detector_hh
#define GeometryService_Detector_hh

//
// A base class for detector components.
//
// $Id: Detector.hh,v 1.3 2012/02/24 20:55:48 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 20:55:48 $
//
// Original author Rob Kutschke
//


#include <string>

namespace mu2e
{
  class Detector
  {
  public:
    virtual ~Detector();
  };
}

#endif /* GeometryService_Detector_hh */
