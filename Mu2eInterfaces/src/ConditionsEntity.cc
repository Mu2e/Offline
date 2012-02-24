//
// A base class for objects held by the conditions data system.
//
//
// $Id: ConditionsEntity.cc,v 1.1 2012/02/24 16:36:36 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 16:36:36 $
//
// Original author Rob Kutschke
//

#include "Mu2eInterfaces/inc/ConditionsEntity.hh"

namespace mu2e
{
  ConditionsEntity::~ConditionsEntity() { }
  void ConditionsEntity::update() { }
  std::string ConditionsEntity::name() const { return "NONAME"; }
}
