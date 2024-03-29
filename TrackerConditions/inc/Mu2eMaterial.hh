#ifndef TrackerConditions_Mu2eMaterial_hh
#define TrackerConditions_Mu2eMaterial_hh
//
// This holds pointers to BTrk objects that
// describe materails that BTrk will use.
// BTrk has singletons and to reflect that, there
// is only ever one instance of this service entity.
// Initialized with Mu2eMaterialMaker
//

#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/Mu2eBTrk/inc/FileFinder.hh"
#include "Offline/Mu2eBTrk/inc/ParticleInfo.hh"
#include "Offline/Mu2eBTrk/inc/DetStrawType.hh"
#include "BTrk/MatEnv/MatDBInfo.hh"
#include <memory>

namespace mu2e {

  class Mu2eMaterial : public ProditionsEntity {

    public:

      typedef std::shared_ptr<Mu2eMaterial> ptr_t;
      typedef std::shared_ptr<const Mu2eMaterial> cptr_t;
      friend class Mu2eMaterialMaker;
      constexpr static const char* cxname = {"Mu2eMaterial"};

      Mu2eMaterial():ProditionsEntity(cxname)  {}
      virtual ~Mu2eMaterial() {}

      DetStrawType const* strawType() const { return _strawtype.get(); }

      void print( std::ostream& ) const;

    private:

      // types for straw elements
      std::unique_ptr<DetStrawType> _strawtype; // straw materials description
      // materials
      std::string _gasmatname, _wallmatname, _wirematname;

      // this hold pointers to material descriptions inside of
      // BTrk material singleton
      MatDBInfo _mat;

      std::unique_ptr<FileFinder> _fileFinder;
      std::unique_ptr<ParticleInfo> _particleInfo;

  };

} // namespace mu2e

#endif /* TrackerConditions_Mu2eMaterial_hh */
