#ifndef ConditionsService_ConditionsEntity_hh
#define ConditionsService_ConditionsEntity_hh
//
// A base class for objects held by the conditions data system.
//
//
// Original author Rob Kutschke
//


#include <string>

namespace mu2e
{
  class ConditionsEntity
  {
  public:
    ConditionsEntity() {}
    virtual ~ConditionsEntity() = default;
  };
}

#endif /* ConditionsService_ConditionsEntity_hh */
