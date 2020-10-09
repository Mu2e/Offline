#ifndef GeometryService_GeomHandle_hh
#define GeometryService_GeomHandle_hh

//
// A safe pointer to the geometry information for a
// detector component.
//
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
    DET const * get()        const { return _detector;}
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

#endif /* GeometryService_GeomHandle_hh */
