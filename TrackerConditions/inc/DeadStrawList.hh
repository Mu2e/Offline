#ifndef HitMakers_DeadStrawList_hh
#define HitMakers_DeadStrawList_hh
//
// Maintain a list of dead straws.  This should be moved to the conditions system.
//
// $Id: DeadStrawList.hh,v 1.1 2014/05/30 19:15:32 brownd Exp $
// $Author: brownd $
// $Date: 2014/05/30 19:15:32 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "DataProducts/inc/StrawIndex.hh"

// art and its tool chain
#include "fhiclcpp/ParameterSet.h"

// C++ includes
#include <vector>
#include <iosfwd>

namespace mu2e {

  class DeadStrawList{

  public:
    DeadStrawList( fhicl::ParameterSet const& pset );
    // Accept compiler supplied d'tor, copy assignment and copy c'tor.

    // Reset the state when conditions information changes.
    void reset( fhicl::ParameterSet const& pset );

    // Accessors
    bool isAlive( StrawIndex i ) const { return  _alive.at(i.asUint()); }
    bool isDead ( StrawIndex i ) const { return !_alive.at(i.asUint()); }

    void print( std::ostream& ) const;

  private:

    // Control printout
    int _verbosity;

    // True if the straw is alive; indexed by the StrawIndex.
    std::vector<bool> _alive;

  };

} // namespace mu2e

#endif /* HitMakers_DeadStrawList_hh */
