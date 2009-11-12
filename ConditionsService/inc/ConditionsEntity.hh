#ifndef ConditionsEntity_H
#define ConditionsEntity_H
//
// A base class for objects held by the conditions data system.
//
// $Id: ConditionsEntity.hh,v 1.1 2009/11/12 00:51:08 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/12 00:51:08 $
//
// Original author Rob Kutschke
//


#include <string>

namespace mu2e
{
  class ConditionsEntity
  {
  public:
    virtual ~ConditionsEntity();
    virtual void update();
    virtual std::string name() const;
  };
}

#endif
