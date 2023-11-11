//
// A base class for objects held by the conditions data system.
//
//
//
// Original author Rob Kutschke
//

#include "Offline/Mu2eInterfaces/inc/ConditionsEntity.hh"

namespace mu2e
{
  ConditionsEntity::~ConditionsEntity() = default;
  ConditionsEntity::ConditionsEntity(const ConditionsEntity&) = default;
  ConditionsEntity& ConditionsEntity::operator=(const ConditionsEntity &) = default;
}
