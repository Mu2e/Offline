#ifndef Mu2eG4_InitEnvToolBase_hh
#define Mu2eG4_InitEnvToolBase_hh

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  class InitEnvToolBase {
  protected:
    std::string _name;
  public:

    InitEnvToolBase() noexcept = default ;
    virtual ~InitEnvToolBase()  noexcept = default ;

    virtual const std::string& name() { return _name; };

    virtual int construct(VolumeInfo const & ParentVInfo, SimpleConfig const& Config) ;

  };
}

#endif

