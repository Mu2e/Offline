#ifndef STRAWINDEX_HH
#define STRAWINDEX_HH

//
// Helper struct for cluster finding.  
// Holds a StrawIndex and an index into a container of hits.
// 
// $Id: Candidate.hh,v 1.2 2010/05/18 20:28:07 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 20:28:07 $
//
// Original author Rob Kutschke
//

#include <ostream>

namespace mu2e {


  struct Candidate{

    // Identifier of the hit straw ( StrawIndex or VolumeId ).
    int id;

    // Index in the original hit container.
    int hitIndex;

    explicit Candidate( int id_=-1, int hitIndex_=-1):
      id(id_),
      hitIndex(hitIndex_){
    }

    // Compiler generated versions are OK for:
    // copy c'tor, destructor, operator=

    // The required behaviour of the comparison operators is that they
    // compare only the values of the id.
    bool operator==( Candidate const& rhs) const{
      return (id == rhs.id);
    }

    bool operator!=( Candidate const& rhs) const{
      return !( *this == rhs );
    }

    bool operator<( Candidate const& rhs) const{
      return ( id < rhs.id);
    }

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   Candidate const& c){
    ost << "( " 
        << c.id << "," 
        << c.hitIndex <<")";
    return ost;
  }  

}

#endif
