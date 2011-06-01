#ifndef TrackerGeom_StrawIndex_hh
#define TrackerGeom_StrawIndex_hh

//
// Dense integer identifier of one straw.
// Has values 0...(N-1), where N is the number
// of straws in the system.  This works for both the LTracker
// and the TTracker.
//
// $Id: StrawIndex.hh,v 1.9 2011/06/01 16:02:58 mf Exp $
// $Author: mf $
// $Date: 2011/06/01 16:02:58 $
//
// Original author Rob Kutschke
//
// This class exists to allow compile time checks that
// the code is using the correct sort of index when indexing
// into a container that requires a StrawIndex as an index.
// For now it is a bit of an experiment.
//
// Most users should never need to access a StrawIndex as an
// integer or set one from an integer.  If this goes according
// to plan, only the code that creates the LTracker will need
// to set one from an integer.
//
// Added functionallity allowing this class to be used effeciently
// when traversing collectoins of straws in a panel -- the NO_STRAW
// indicator.

#include <ostream>

namespace mu2e {


  class StrawIndex{

  public: 
    enum NoStraw {NO_STRAW = -999999};

  public:

    // No default c'tor by design.

    // No automatic conversion of int to StrawIndex.
    explicit StrawIndex(int idx):
      _idx(idx){
    }

    // Compiler generated versions are OK for:
    // copy c'tor, destructor, operator=

    // Return the value as an int or as an unsigned in
    // Do not want automatic conversion to an int.
    int      asInt() const { return _idx;}
    unsigned asUint() const { return static_cast<unsigned>(_idx);}

    bool operator==( StrawIndex const& rhs) const{
      return (_idx == rhs._idx);
    }

    bool operator<( StrawIndex const& rhs) const{
      return ( _idx < rhs._idx);
    }
    bool operator>( StrawIndex const& rhs) const{
      return ( _idx > rhs._idx);
    }
  private:

    int _idx;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StrawIndex const& i){
    ost << i.asInt();
    return ost;
  }

  inline bool operator!=( StrawIndex const& lhs,
                          StrawIndex const& rhs) {
    return !( lhs == rhs);
  }

}
#endif /* TrackerGeom_StrawIndex_hh */
