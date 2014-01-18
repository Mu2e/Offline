#ifndef ConditionsService_ConditionsHandle_hh
#define ConditionsService_ConditionsHandle_hh

//
// A safe pointer to a ConditionsEntity.
//
// $Id: ConditionsHandle.hh,v 1.7 2014/01/18 17:31:59 brownd Exp $
// $Author: brownd $
// $Date: 2014/01/18 17:31:59 $
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

    ConditionsHandle() : _entity(0) {}
    ConditionsHandle(const ConditionsHandle& other ) : _entity(other._entity) {}
    ConditionsHandle& operator=(const ConditionsHandle& other) {
      if(this != &other){
	_entity = other._entity;
      }
      return *this;
    }

  private:
    // unnecessary
    ENTITY* operator&();

    ENTITY* _entity;
  };
}

#endif /* ConditionsService_ConditionsHandle_hh */
