#ifndef GeometryService_Detector_hh
#define GeometryService_Detector_hh

//
// A base class for detector components.
//
// $Id: Detector.hh,v 1.2 2012/02/24 16:37:09 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 16:37:09 $
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
    //virtual std::string name() const = 0;
  };
}

#endif /* GeometryService_Detector_hh */
