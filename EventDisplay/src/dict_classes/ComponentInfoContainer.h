//
// Container class for component info objects.
//
// $Id: ComponentInfoContainer.h,v 1.2 2012/09/14 17:17:34 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2012/09/14 17:17:34 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_dict_classes_ComponentInfoContainer_h
#define EventDisplay_src_dict_classes_ComponentInfoContainer_h


#include "EventDisplay/src/dict_classes/ComponentInfo.h"
#ifndef __CINT__
#include "boost/shared_ptr.hpp"
#endif

namespace mu2e_eventdisplay
{

  class ComponentInfoContainer
  {
#ifndef __CINT__
    boost::shared_ptr<ComponentInfo> _componentInfo;
#endif

    ComponentInfoContainer();
    ComponentInfoContainer(const ComponentInfoContainer &);
    ComponentInfoContainer& operator=(const ComponentInfoContainer &);

    public:
#ifndef __CINT__
    ComponentInfoContainer(const boost::shared_ptr<ComponentInfo> c) : _componentInfo(c) {}
#endif

    virtual ~ComponentInfoContainer() {}

#ifndef __CINT__
    const boost::shared_ptr<ComponentInfo> getComponentInfo() const {return _componentInfo;}
#endif

    ClassDef(ComponentInfoContainer,1);
  };

}
#endif /* EventDisplay_src_dict_classes_ComponentInfoContainer_h */

