#ifndef G4GEOMETRY_OPTIONS
#define G4GEOMETRY_OPTIONS

#include <map>
#include <string>
#include <vector>

namespace mu2e {

  class SimpleConfig;

  class G4GeometryOptions {
  public:

    typedef std::vector<std::string> VS;
    typedef std::pair<bool,VS> Ordering;
    typedef std::vector<Ordering> OrderingList;

    G4GeometryOptions( const SimpleConfig& config );

    // Disable copy c'tor and copy assignment
    G4GeometryOptions           (const G4GeometryOptions&) = delete;
    G4GeometryOptions& operator=(const G4GeometryOptions&) = delete;

    void loadOrderingStrings( OrderingList& orderingList,
                              const SimpleConfig& config,
                              const std::string& varString );
                              

    void loadEntry    ( const SimpleConfig& config,
                        const std::string& volName,
                        const std::string& prefix );

    bool isSolid            ( const std::string& volName ) const;
    bool isVisible          ( const std::string& volName ) const;
    bool doSurfaceCheck     ( const std::string& volName ) const;
    bool forceAuxEdgeVisible( const std::string& volName ) const;
    bool placePV            ( const std::string& volName ) const;

  private:
    
    static void mapInserter( std::map<std::string,bool>& map,
                             const SimpleConfig& c,
                             const std::string& volName,
                             const std::string& prefix,
                             const bool default_value );
    
    static bool queryMap( const std::map<std::string,bool>& map,
                          const OrderingList& orderingList,
                          const std::string& volName,
                          const bool default_value );

    static std::pair<bool,bool> flagOverridden( const std::string& volName,
                                                const OrderingList& orderingList );

    bool defaultIsSolid_;
    bool defaultIsVisible_;
    bool defaultDoSurfaceCheck_;
    bool defaultForceAuxEdgeVisible_;
    bool defaultPlacePV_;

    std::map<std::string,bool> isSolid_;
    std::map<std::string,bool> isVisible_;
    std::map<std::string,bool> doSurfaceCheck_;
    std::map<std::string,bool> forceAuxEdgeVisible_;
    std::map<std::string,bool> placePV_;

    OrderingList orderingIsSolid_;
    OrderingList orderingIsVisible_;
    OrderingList orderingDoSurfaceCheck_;
    OrderingList orderingForceAuxEdgeVisible_;
    OrderingList orderingPlacePV_;

  };

} /* end mu2e namespace*/

#endif /*G4GEOMETRY_OPTIONS*/
