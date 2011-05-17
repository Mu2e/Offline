#ifndef ConditionsService_ConditionsEntity_hh
#define ConditionsService_ConditionsEntity_hh
//
// A base class for objects held by the conditions data system.
//
// $Id: ConditionsEntity.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
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

#endif /* ConditionsService_ConditionsEntity_hh */
