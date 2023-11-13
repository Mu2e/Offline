#ifndef GeomPrimitives_NotchHoleManager_HH
#define GeomPrimitives_NotchHoleManager_HH
// NotchHoleManager is a class that keeps track of notches.  Specifically,
// it is intended to keep track of notches in the walls/floors/etc of
// the Mu2e building.  This allows us to add this possibility with
// a fairly small footprint in memory and code.  Current code should not
// change its functionality or form (except where the notches are specifically
// added), so this solution should be minimally invasive, while allowing
// us extra realism in our building geometry.
// There is no reason this NotchHoleManager couldn't also be used to track notches
// used in other geometry elements, such as the external shielding.

// changed for holes too :Sridhar Tripathy


#include "Offline/GeomPrimitives/inc/Hole.hh"
#include "Offline/GeomPrimitives/inc/Notch.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include <vector>
#include <string>
#include <unordered_map>

namespace mu2e{

  class NotchHoleManager {
  public:
    NotchHoleManager():hasLoaded_(false){}
    ~NotchHoleManager(){}

    void loadNotches( const SimpleConfig& config );
    void loadHoles( const SimpleConfig& config );

    const bool hasNotches( const std::string& part ) const;
    const bool hasHoles( const std::string& part ) const;

    const std::vector<Notch>& getNotchVector( const std::string& part ) const;
    const std::vector<Hole>& getHoleVector( const std::string& part ) const;

  private:

    std::unordered_map<std::string,std::vector<Notch> > theMap_;
    std::unordered_map<std::string,std::vector<Hole> > theMapH_;
    std::vector<Notch> emptyVec_;
    std::vector<Hole> emptyVecH_;
    bool hasLoaded_;
    bool hasLoadedH_;
  }; // end NotchHoleManager class definition

} // end namespace mu2e

#endif //  GeomPrimitives_NotchHoleManager_HH
