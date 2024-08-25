#ifndef G4GEOMETRY_OPTIONS
#define G4GEOMETRY_OPTIONS
//
// G4 geometry options look-up facility, to be used in conjunction
// with SimpleConfig.
//
//
// Original author: Kyle Knoepfel
//
// This method is used for setting and overriding various flags that
// are specified when creating volumes in G4.  Ideally, it would go in
// the Mu2eG4Helper service, but it is tied to GeometryService because of
// SimpleConfig and linkage loops.
//
// The idiom of this helper is the following:
//
//   (1) A SimpleConfig file can specify the following assignments:
//
//         bool <var_prefix>.isVisible           = [ true or false ];
//         bool <var_prefix>.isSolid             = [ true or false ];
//         bool <var_prefix>.forceAuxEdgeVisible = [ true or false ];
//         bool <var_prefix>.placePV             = [ true or false ];
//         bool <var_prefix>.doSurfaceCheck      = [ true or false ];
//
//  (2) The various flags are loaded into the option maps by the
//      following syntax within a .cc file:
//
//         G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
//         geomOption->loadEntry( configFile, "MATCHING_TOKEN", <var_prefix> );
//
//      where the "MATCHING_TOKEN" is specified by the User in terms
//      of what you want the querying functions to look for.  Normally
//      the value of "MATCHING_TOKEN" applies to several volumes, but
//      it could be chosen for each volume. If "loadEntry" is
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
//  (4) The (e.g.) visible() facility will first search through the
//      corresponding map for a match.  If no match is found---i.e. an
//      entry corresponding to the requested "MATCHING_TOKEN" does not
//      exist---the default visible value is returned.
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

#include <map>
#include <string>
#include <vector>
#include <regex>


namespace mu2e {

  class SimpleConfig;



  class G4GeometryOptData {

     public:

       typedef std::vector<std::string>   VS;
       typedef std::pair<bool,std::regex> Ordering;
       typedef std::vector<Ordering>      OrderingList;


       G4GeometryOptData( bool defaultValue, const std::string& name );

       void loadOrderingStrings( const SimpleConfig& config, const std::string& varString );
       void mapInserter        ( const std::string& volName, bool value );
       bool queryMap           ( const std::string& volName             ) const;
       bool default_value() const {return default_;}


     private:

       std::pair<bool,bool> flagOverridden( const std::string& volName ) const;

       std::string                 name_;
       std::map<std::string, bool> map_;
       OrderingList                ordering_;
       bool                        default_;
  };




  class G4GeometryOptions {

    public:

      G4GeometryOptions( const SimpleConfig& config );
      ~G4GeometryOptions() = default;

      // Disable copy c'tor and copy assignment
      G4GeometryOptions           (const G4GeometryOptions&) = delete;
      G4GeometryOptions& operator=(const G4GeometryOptions&) = delete;
      G4GeometryOptions           (G4GeometryOptions&&)      = delete;
      G4GeometryOptions& operator=(G4GeometryOptions&&)      = delete;

      void loadEntry( const SimpleConfig& config, const std::string& volName, const std::string& prefix );

      bool isSolid            ( const std::string& volName ) const;
      bool isVisible          ( const std::string& volName ) const;
      bool doSurfaceCheck     ( const std::string& volName ) const;
      bool forceAuxEdgeVisible( const std::string& volName ) const;
      bool placePV            ( const std::string& volName ) const;


    private:

      G4GeometryOptData dataSurfaceCheck_;
      G4GeometryOptData dataIsVisible_;
      G4GeometryOptData dataIsSolid_;
      G4GeometryOptData dataForceAuxEdge_;
      G4GeometryOptData dataPlacePV_;

  };

}

#endif /*G4GEOMETRY_OPTIONS*/
