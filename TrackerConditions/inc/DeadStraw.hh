#ifndef TrackerConditions_DeadStraw_hh
#define TrackerConditions_DeadStraw_hh
//
// A list of dead straws.
// Initialized with DeadStrawMaker
//

// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"

// C++ includes
#include <set>
#include <memory>

namespace mu2e {

// struct to define a range of deadness in a straw.  This is assumed to be equidistant
// from the straw center, to model radiation damage due to dose (which is concentrated
// at the lowest radius = straw center)
//

  class DeadStraw : public ProditionsEntity {

  public:

    struct Range {
      StrawId _strawId; // which straw
      double _range; // how much is dead
      // constructor.  Default range is longer than the longest straw, meaning the
      // straw is 100% dead
      Range(StrawId const& id,double range=1000.0) : 
	                            _strawId(id), _range(range) {}
      Range(Straw const& straw,double range=1000.0) : Range(straw.id(),range) {}
      // test if a given position (relative to the straw center) is dead
      bool isDead(double hitpos) const { return fabs(hitpos)<_range; }
      // define comparator to ignore range: this allows finding an 
      // element based on strawid quickly
      bool operator < (Range const& other) const { 
	return _strawId < other._strawId; }
    };

    typedef std::shared_ptr<DeadStraw> ptr_t;
    typedef std::shared_ptr<const DeadStraw> cptr_t;
    typedef std::set<Range> set_t;

    DeadStraw(): _name("DeadStraw") {}
    DeadStraw(set_t const& deadstraws):
      _name("DeadStraw"),_deadstraws(deadstraws) {}

    virtual ~DeadStraw() {}

    // Accessors; test if a straw is dead.  The position argument refers to 
    // distance from the straw center in mm
    bool isAlive( StrawId id, double hitpos=0.0 ) const { return !(isDead(id,hitpos)); }
    bool isDead ( StrawId id, double hitpos=0.0 ) const;

    std::string const& name() const { return _name; }
    void print( std::ostream& ) const;

  private:
    std::string _name;
    set_t _deadstraws; // sparse list of dead straws.  Ordering allows rapid lookup.

  };

} // namespace mu2e

#endif /* TrackerConditions_DeadStraw_hh */
