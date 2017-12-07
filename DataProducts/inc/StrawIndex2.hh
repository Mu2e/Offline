#ifndef TrackerGeom_StrawIndex2_hh
#define TrackerGeom_StrawIndex2_hh

//
// Non contiguous uint16_t identifier of one straw used to access 
// std::array<Straw const*,TTracker::_maxRedirect> _allStraws2_p;
//
// Original author Krzysztof Genser based on StrawIndex by Rob Kutschke
//
// This class exists to allow compile time checks that
// the code is using the correct sort of index when indexing
// into a container that requires a StrawIndex2 as an index.
//

#include <ostream>

namespace mu2e {


  class StrawIndex2{

  public: 
    enum NoStraw {NO_STRAW = -999999}; // do we still need it and what it should be???

  public:

    // No default c'tor by design.

    // No automatic conversion of int to StrawIndex2.
    explicit StrawIndex2(uint16_t idx):
      _idx(idx){
    }

    // Compiler generated versions are OK for:
    // copy c'tor, destructor, operator=

    // Return the value as an int or as an unsigned in
    // Do not want automatic conversion to an int.
    int      asInt() const { return static_cast<int>(_idx);}
    unsigned asUint16() const { return _idx;}

    bool operator==( StrawIndex2 const& rhs) const{
      return (_idx == rhs._idx);
    }

    bool operator<( StrawIndex2 const& rhs) const{
      return ( _idx < rhs._idx);
    }
    bool operator>( StrawIndex2 const& rhs) const{
      return ( _idx > rhs._idx);
    }
  private:

    uint16_t _idx;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StrawIndex2 const& i){
    ost << i.asInt();
    return ost;
  }

  inline bool operator!=( StrawIndex2 const& lhs,
                          StrawIndex2 const& rhs) {
    return !( lhs == rhs);
  }

}
#endif /* TrackerGeom_StrawIndex2_hh */
