#ifndef ConditionsService_ConditionsHandle_hh
#define ConditionsService_ConditionsHandle_hh

//
// A safe pointer to a ConditionsEntity.
//
//
// Original author Rob Kutschke
//

#include <string>

#include "Offline/ConditionsService/inc/ConditionsService.hh"

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
