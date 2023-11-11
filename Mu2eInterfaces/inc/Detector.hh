#ifndef GeometryService_Detector_hh
#define GeometryService_Detector_hh

//
// A base class for detector components.
//
//
// Original author Rob Kutschke
//


#include <string>

namespace mu2e
{
  class Detector
  {
  public:
    Detector() {}
    virtual ~Detector();
    Detector(const Detector&);
    Detector& operator=(const Detector &);
  };
}

#endif /* GeometryService_Detector_hh */
