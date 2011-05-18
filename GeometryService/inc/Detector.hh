#ifndef GeometryService_Detector_hh
#define GeometryService_Detector_hh

//
// A base class for detector components.
//
// $Id: Detector.hh,v 1.3 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
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
