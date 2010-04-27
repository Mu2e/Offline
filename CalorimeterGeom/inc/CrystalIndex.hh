#ifndef CRYSTALINDEX_HH
#define CRYSTALINDEX_HH

//
// Dense integer identifier of one crystal
// Has values 0...(N-1), where N is the number
// of crystals in the system.  This works for both the LTracker
// and the TTracker having been copied over.  
//
// $Id: CrystalIndex.hh,v 1.3 2010/04/27 18:47:06 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/27 18:47:06 $
//
// Original author Rob Kutschke
//
// This class exists to allow compile time checks that
// the code is using the correct sort of index when indexing
// into a container that requires a CrystalIndex as an index.
// For now it is a bit of an experiment.
//
// Most users should never need to access a CrystalIndex as an
// integer or set one from an integer.  If this goes according
// to plan, only the code that creates the Calorimeter will need
// to set one from an integer.
//

#include <ostream>

namespace mu2e {
  namespace calorimeter{

  class CrystalIndex{

  public:

    // No default c'tor by design.

    // No automatic conversion of int to CrystalIndex.
    explicit CrystalIndex(int idx):
      _idx(idx){
    }
    
    // Compiler generated versions are OK for:
    // copy c'tor, destructor, operator=

    // Return the value as an int.
    // Do not want automatic conversion to an int.
    int asInt() const { return _idx;}

    bool operator==( CrystalIndex const& rhs) const{
      return (_idx == rhs._idx);
    }

    bool operator<( CrystalIndex const& rhs) const{
      return ( _idx < rhs._idx);
    }

  private:

    int _idx;
  };

  inline std::ostream& operator<<( std::ostream& ost,
				   CrystalIndex const& i){
    ost << i.asInt();
    return ost;
  }

  inline bool operator!=( CrystalIndex const& lhs, 
			  CrystalIndex const& rhs) {
      return !( lhs == rhs);
  }

  
  } // namespace calorimeter
} //namespace mu2e

#endif
