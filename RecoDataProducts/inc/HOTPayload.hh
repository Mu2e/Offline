#ifndef RecoDataProducts_HOTPayload_hh
#define RecoDataProducts_HOTPayload_hh
//
// Peristent data for one hit-on-track object.
//
// $Id: HOTPayload.hh,v 1.1 2012/07/03 03:27:24 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/03 03:27:24 $
//
// Contact author Rob Kutschke
//

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes
//#include "DataProducts/inc/StrawIndex.hh"

namespace mu2e {

  class HOTPayload{
  public:

    double t0() const { return t0_; }

    HOTPayload(){}
    // Accept compiler written d'tor, copy c'tor and copy assignment.

    HOTPayload( double at0):
      t0_(at0){
    }

    // Accessors

    // Print contents of the object.
    void print( std::ostream& ost, bool doEndl = true ) const;

  private:

    double t0_;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   HOTPayload const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* RecoDataProducts_HOTPayload_hh */
