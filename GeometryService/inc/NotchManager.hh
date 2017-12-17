#ifndef GeomPrimitives_NotchManager_HH
#define GeomPrimitives_NotchManager_HH
// NotchManager is a class that keeps track of notches.  Specifically,
// it is intended to keep track of notches in the walls/floors/etc of
// the Mu2e building.  This allows us to add this possibility with
// a fairly small footprint in memory and code.  Current code should not
// change its functionality or form (except where the notches are specifically
// added), so this solution should be minimally invasive, while allowing
// us extra realism in our building geometry.
// There is no reason this NotchManager couldn't also be used to track notches
// used in other geometry elements, such as the external shielding.

#include "GeomPrimitives/inc/Notch.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include <vector>
#include <string>
#include <unordered_map>

namespace mu2e{

  class NotchManager {
  public:
    NotchManager():hasLoaded_(false){}
    ~NotchManager(){}

    void loadNotches( const SimpleConfig& config );

    const bool hasNotches( const std::string& part ) const;

    const std::vector<Notch>& getNotchVector( const std::string& part ) const;

  private:

    std::unordered_map<std::string,std::vector<Notch> > theMap_;
    std::vector<Notch> emptyVec_;
    bool hasLoaded_;
  }; // end NotchManager class definition

} // end namespace mu2e

#endif //  GeomPrimitives_NotchManager_HH
