//
// Container class for component info objects.
//
// $Id: ComponentInfoContainer.h,v 1.1 2011/09/08 03:54:45 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2011/09/08 03:54:45 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_dict_classes_ComponentInfoContainer_h
#define EventDisplay_src_dict_classes_ComponentInfoContainer_h


#include "ComponentInfo.h"
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

    ClassDef(ComponentInfoContainer,0);
  };

}
#endif /* EventDisplay_src_dict_classes_ComponentInfoContainer_h */
