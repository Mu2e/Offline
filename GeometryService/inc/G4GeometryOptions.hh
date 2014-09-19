#ifndef G4GEOMETRY_OPTIONS
#define G4GEOMETRY_OPTIONS
//
// G4 geometry options look-up facility, to be used in conjunction
// with SimpleConfig.
//
// $Id: G4GeometryOptions.hh,v 1.3 2014/09/19 20:06:25 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/19 20:06:25 $
//
// Original author: Kyle Knoepfel
//
// This method is used for setting and overriding various flags that
// are specified when creating volumes in G4.  Ideally, it would go in
// the G4Helper service, but it is tied to GeometryService because of
// SimpleConfig and linkage loops.
//
// The idiom of this helper is the following:
//
//   (1) A SimpleConfig file can specify the following assignments:
//  
//         bool <var_prefix>.isVisible           = [ true or false ];
//         bool <var_prefix>.isSolid             = [ true or false ];
//         bool <var_prefix>.forceAuxEdgeVisible = [ true or false ]; --> DISCOURAGED
//         bool <var_prefix>.placePV             = [ true or false ]; --> DISCOURAGED
//         bool <var_prefix>.doSurfaceCheck      = [ true or false ]; --> DISCOURAGED
//
//  (2) The various flags are loaded into the option maps by the
//      following syntax within a .cc file:
//
//         G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
//         geomOption->loadEntry( configFile, "MATCHING_TOKEN", <var_prefix> );
//
//      where the "MATCHING_TOKEN" is specified by the User in terms
//      of what you want the querying functions to look for.  Normally
//      the value of "MATCHING_TOKEN" is the physical volume name, but
//      it does not have to be (see next point).  If "loadEntry" is
//      not included for a given volume, then the 5 flags above
//      default to global values.
//
//  (3) To access the flags, the following can be done:
// 
//         const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
//         geomOptions->isVisible( "MATCHING_TOKEN" );
//         geomOptions->isSolid  ( "MATCHING_TOKEN" );
//         etc.
//
//      If one were to do the following (the following is pseudo-code):
//
//         vector<VolumeParams> volumes;  // A vector with a lot of volume parameters
//
//         for ( const auto& volParams ; volumes ) {
//
//            finishNesting(  volParams,
//                            ...  
//                            geomOptions->isVisible( volParams.volumeName );
//                            ... );
//         }
//
//      such a query could take a long time.  For that reason, the
//      "MATCHING_TOKEN" value does not need to match that of the
//      volume name to be created.  The following can be much faster:
//
//         vector<VolumeParams> volumes;   // A vector with a lot of volume parameters
//         bool isVisible = geomOptions->isVisible( "Straw" ); // look-up once.
//         for ( const auto& volParams ; volumes ) {
//
//            finishNesting(  volParams,
//                            ...  
//                            isVisible
//                            ... );
//         }      
//
//      Note that an individual volume (e.g. straw) can be viewed by
//      specifying an override (see point 5).
//
//  (4) The (e.g.) isVisible() facility will search through the
//      corresponding map for a match.  If no match is found---i.e. an
//      entry corresponding to the requested "MATCHING_TOKEN" does not
//      exist---the default isVisible value is returned.  
//
//  (5) The value returned from step 4 can be overridden by specifying
//      override commands in Mu2eG4/geom/g4_userOptions.txt (e.g.):
//      
//         bool g4.doSurfaceCheck = true;
//         vector<string> g4.doSurfaceCheck.drop  = {"*"};
//         vector<string> g4.doSurfaceCheck.keep  = {"PSShield*"};
//         vector<string> g4.doSurfaceCheck.order = { "g4.doSurfaceCheck.drop", 
//                                                    "g4.doSurfaceCheck.keep" };
//       
//      In this case, the default "doSurfaceCheck" value is true, but
//      the doSurfaceCheck's for all volumes are disabled by the drop
//      "*" command, since "*" matches to all volumes.  All volumes
//      that match "PSShield*" then have their surface checks enabled.
//      Note that the commands in "drop" and "keep" always override
//      the default g4.doSurfaceCheck value.
//
//      The actual drop/keep commands are not implemented unless they
//      are specified in the *.order vector in the order desired. 
//
//      Additional drop/keep commands can be added.  The only
//      requirement is that their suffixes must of the form *.keep* or
//      *.drop*.

// C++ includes
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


    void loadEntry( const SimpleConfig& config,
		    const std::string& volName,
		    const std::string& prefix );

    bool isSolid            ( const std::string& volName ) const;
    bool isVisible          ( const std::string& volName ) const;
    bool doSurfaceCheck     ( const std::string& volName ) const;
    bool forceAuxEdgeVisible( const std::string& volName ) const;
    bool placePV            ( const std::string& volName ) const;

  private:
    
    static void loadOrderingStrings( OrderingList& orderingList,
				     const SimpleConfig& config,
				     const std::string& varString );

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

    // I don't like having 5 copies of the same structure.  Might be
    // better to do an enum-specified map, like:
    //
    //    enum enum_type { ISSOLID, ISVISIBLE, SURFACECHECK, FORCEAUX, PLACEPV };
    //    std::map<enum_type,OrderingList> map_of_flags_;
    //    
    // but for an OrderingList map, that would resolve to:
    //
    //    std::map<enum_type,std::vector<std::pair<bool,std::vector<std::string>>>>;
    //
    // kind of ugly.  Should probably make a class that houses the
    // global defaults, local-default maps, and OrderingLists.  Maybe
    // in the future...

    // Global defaults
    bool defaultIsSolid_;
    bool defaultIsVisible_;
    bool defaultDoSurfaceCheck_;
    bool defaultForceAuxEdgeVisible_;
    bool defaultPlacePV_;

    // Local defaults
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
