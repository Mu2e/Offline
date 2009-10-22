#ifndef STRAWINDEX_HH
#define STRAWINDEX_HH

//
// Dense integer identifier of one straw.
// Has values 0...(N-1), where N is the number
// of straws in the system.
//
// $Id: StrawIndex.hh,v 1.2 2009/10/22 16:21:05 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/10/22 16:21:05 $
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

#include <ostream>

namespace mu2e {


  class StrawIndex{

  public:

    // No default c'tor by design.
    
    // Compiler generated versions are OK for:
    // copy c'tor, destructor, operator=

    // Return the value as an int.
    int asInt() const { return _idx;}

    // Set value from an int.
    static StrawIndex fromInt(int i) { return StrawIndex(i); }

    bool operator==( StrawIndex const& rhs) const{
      return (_idx == rhs._idx);
    }

    bool operator!=( StrawIndex const& rhs) const{
      return !( *this == rhs);
    }

    bool operator<( StrawIndex const& rhs) const{
      return ( _idx < rhs._idx);
    }

  private:

    // 
    StrawIndex(int idx):
      _idx(idx){
    }

    int _idx;
  };

  inline std::ostream& operator<<( std::ostream& ost,
				   StrawIndex const& i){
    ost << i.asInt();
    return ost;
  }
  

}

#endif
