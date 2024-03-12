#ifndef Mu2eG4_InitEnvToolBase_hh
#define Mu2eG4_InitEnvToolBase_hh

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  class InitEnvToolBase {
  public:

    InitEnvToolBase() noexcept = default ;

    virtual ~InitEnvToolBase()  noexcept = default ;

    virtual int construct(VolumeInfo const & ParentVInfo, SimpleConfig const& Config) = 0;

  };
}

#endif
