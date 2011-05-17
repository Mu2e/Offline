#ifndef GeometryService_GeomHandle_hh
#define GeometryService_GeomHandle_hh

//
// A safe pointer to the geometry information for a 
// detector component.
//
// $Id: GeomHandle.hh,v 1.2 2011/05/17 15:36:00 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:36:00 $
//
// Original author Rob Kutschke
//

#include "GeometryService/inc/GeometryService.hh"

namespace mu2e {
  template <typename DET>
  class GeomHandle
  {
  public:
    GeomHandle()
    {
      art::ServiceHandle<GeometryService> sg;
      _detector = sg->getElement<DET>();
    }
    ~GeomHandle() { }
    
    DET const * operator->() const { return _detector;}
    DET const & operator*()  const { return *_detector;}
    DET const * operator->() { return _detector;}
    DET const & operator*()  { return *_detector;}
    
  private:
    GeomHandle(const GeomHandle&);
    GeomHandle& operator=(const GeomHandle&);
    
    // unnecessary
    DET* operator&();
    
    DET* _detector;
  };
}

#endif
