#ifndef GeometryService_Detector_hh
#define GeometryService_Detector_hh

//
// A base class for detector components.
//
// $Id: Detector.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:35 $
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
