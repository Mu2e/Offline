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

#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "Mu2eBTrk/inc/DetStrawElem.hh"
#include <memory>
#include <map>

namespace mu2e {

  class Mu2eDetector : public ProditionsEntity {

  public:

    typedef std::shared_ptr<Mu2eDetector> ptr_t;
    typedef std::shared_ptr<const Mu2eDetector> cptr_t;
    friend class Mu2eDetectorMaker;

    Mu2eDetector(): _name("Mu2eDetector") {}
    virtual ~Mu2eDetector();

    const DetStrawElem* strawElem(Straw const& straw) const {
      return strawElem(straw.id());
    }
    const DetStrawElem* strawElem(StrawId const& strawid) const;

    std::string const& name() const { return _name; }
    void print( std::ostream& ) const;

  private:
    std::string _name;

    // map between straw index and detector elements 
    std::map<StrawId,DetStrawElem*> _strawmap;

  };

} // namespace mu2e

#endif /* TrackerConditions_Mu2eDetector_hh */
