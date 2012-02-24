#ifndef GeometryService_Detector_hh
#define GeometryService_Detector_hh

//
// A base class for detector components.
//
// $Id: Detector.hh,v 1.1 2012/02/24 16:36:36 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 16:36:36 $
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
    virtual void update();
    virtual std::string name() const;
  };
}

#endif /* GeometryService_Detector_hh */
