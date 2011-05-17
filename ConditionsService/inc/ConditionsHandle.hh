#ifndef ConditionsService_ConditionsHandle_hh
#define ConditionsService_ConditionsHandle_hh

//
// A safe pointer to a ConditionsEntity.
//
// $Id: ConditionsHandle.hh,v 1.5 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
//
// Original author Rob Kutschke
//

#include <string>

#include "ConditionsService/inc/ConditionsService.hh"

namespace mu2e {
  template <typename ENTITY>
  class ConditionsHandle
  {
  public:
    ConditionsHandle( std::string const& key,
                      std::string const& version="current")
    {
      art::ServiceHandle<ConditionsService> sg;
      _entity = sg->getElement<ENTITY>(key,version);
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

#endif /* ConditionsService_ConditionsHandle_hh */
