#ifndef ConditionsService_ConditionsEntity_hh
#define ConditionsService_ConditionsEntity_hh
//
// A base class for objects held by the conditions data system.
//
// $Id: ConditionsEntity.hh,v 1.2 2012/02/24 16:37:09 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 16:37:09 $
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
    //virtual std::string name() const = 0;
  };
}

#endif /* ConditionsService_ConditionsEntity_hh */
