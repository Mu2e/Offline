//
// A base class for objects held by the conditions data system.
//
//
// $Id: ConditionsEntity.cc,v 1.1 2009/11/12 00:51:08 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/12 00:51:08 $
//
// Original author Rob Kutschke
//

#include "ConditionsService/inc/ConditionsEntity.hh"

namespace mu2e
{
  ConditionsEntity::~ConditionsEntity() { }
  void ConditionsEntity::update() { }
  std::string ConditionsEntity::name() const { return "NONAME"; }
}
