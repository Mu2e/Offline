#ifndef GeometryService_GeomHandle_hh
#define GeometryService_GeomHandle_hh

//
// A safe pointer to the geometry information for a 
// detector component.
//
// $Id: GeomHandle.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
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
      edm::Service<GeometryService> sg;
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
