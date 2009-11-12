#ifndef ConditionsService_ConditionsHandle_hh
#define ConditionsService_ConditionsHandle_hh

//
// A safe pointer to a ConditionsEntity.
//
// $Id: ConditionsHandle.hh,v 1.1 2009/11/12 00:51:08 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/12 00:51:08 $
//
// Original author Rob Kutschke
//

#include "ConditionsService/inc/ConditionsService.hh"

namespace mu2e {
  template <typename ENTITY>
  class ConditionsHandle
  {
  public:
    ConditionsHandle()
    {
      edm::Service<ConditionsService> sg;
      _entity = sg->getElement<ENTITY>();
    }
    ~ConditionsHandle() { }
    
    ENTITY const * operator->() const { return _entity;}
    ENTITY const & operator*()  const { return *_entity;}
    ENTITY const * operator->() { return _entity;}
    ENTITY const & operator*()  { return *_entity;}
    
  private:
    ConditionsHandle(const ConditionsHandle&);
    ConditionsHandle& operator=(const ConditionsHandle&);
    
    // unnecessary
    ENTITY* operator&();
    
    ENTITY* _entity;
  };
}

#endif
