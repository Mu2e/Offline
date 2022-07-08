#ifndef TrackerConditions_Mu2eDetector_hh
#define TrackerConditions_Mu2eDetector_hh
//
// This holds the straw geometry as a set of
// DetStrawElem, which is the format needed by BTrk.
// DetStrawElem has a pointer to a DetStrawType,
// which describes the material in a format BTrk needs, and
// a pointer to a Straw, including aligned wire geometry.
// Initialized with Mu2eDetectorMaker
//

#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/Mu2eBTrk/inc/DetStrawElem.hh"
#include <memory>
#include <map>

namespace mu2e {

  class Mu2eDetector : public ProditionsEntity {

    public:

      typedef std::shared_ptr<Mu2eDetector> ptr_t;
      typedef std::shared_ptr<const Mu2eDetector> cptr_t;
      friend class Mu2eDetectorMaker;
      constexpr static const char* cxname = {"Mu2eDetector"};

      Mu2eDetector(): ProditionsEntity(cxname) {}
      virtual ~Mu2eDetector();

      const DetStrawElem* strawElem(Straw const& straw) const {
        return strawElem(straw.id());
      }
      const DetStrawElem* strawElem(StrawId const& strawid) const;

      void print( std::ostream& ) const;

    private:

      // map between straw index and detector elements
      std::map<StrawId,DetStrawElem*> _strawmap;

  };

} // namespace mu2e

#endif /* TrackerConditions_Mu2eDetector_hh */
