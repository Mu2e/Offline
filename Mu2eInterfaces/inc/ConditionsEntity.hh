#ifndef ConditionsService_ConditionsEntity_hh
#define ConditionsService_ConditionsEntity_hh
//
// A base class for objects held by the conditions data system.
//
// $Id: ConditionsEntity.hh,v 1.1 2012/02/24 16:36:36 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 16:36:36 $
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
