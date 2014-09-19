// C++ includes
#include <regex>

// Mu2e includes
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"

// Framework include files
#include "cetlib/exception.h"

// Boost includes
#include "boost/algorithm/string/replace.hpp"

namespace mu2e {


  //======================================================================
  G4GeometryOptions::G4GeometryOptions( const SimpleConfig& config ) 
    : defaultIsSolid_            ( config.getBool("g4.isSolid"       ) )
    , defaultIsVisible_          ( config.getBool("g4.isVisible"     ) )
    , defaultDoSurfaceCheck_     ( config.getBool("g4.doSurfaceCheck") )
    , defaultForceAuxEdgeVisible_( config.getBool("g4.forceAuxEdgeVisible") )
    , defaultPlacePV_            ( config.getBool("g4.placePV") )
  {

    loadOrderingStrings( orderingIsSolid_,             config, "g4.isSolid.order"            );
    loadOrderingStrings( orderingIsVisible_,           config, "g4.isVisible.order"          );
    loadOrderingStrings( orderingDoSurfaceCheck_,      config, "g4.doSurfaceCheck.order"     );
    loadOrderingStrings( orderingForceAuxEdgeVisible_, config, "g4.forceAuxEdgeVisible.order");
    loadOrderingStrings( orderingPlacePV_,             config, "g4.placePV.order"            );

  }

  //======================================================================
  void G4GeometryOptions::loadVolume( const SimpleConfig& config, 
                                      const std::string& volName, 
                                      const std::string& prefix ) {
    
    mapInserter( doSurfaceCheck_     , config, volName, prefix+".doSurfaceCheck"     , defaultDoSurfaceCheck_      );
    mapInserter( isSolid_            , config, volName, prefix+".isSolid"            , defaultIsSolid_             );
    mapInserter( isVisible_          , config, volName, prefix+".isVisible"          , defaultIsVisible_           );
    mapInserter( forceAuxEdgeVisible_, config, volName, prefix+".forceAuxEdgeVisible", defaultForceAuxEdgeVisible_ );
    mapInserter( placePV_            , config, volName, prefix+".placePV"            , defaultPlacePV_             );

  }
  
  //======================================================================
  bool G4GeometryOptions::isSolid( const std::string& volName ) const {
    return queryMap( isSolid_, 
                     orderingIsSolid_, 
                     volName, 
                     defaultIsSolid_ );
  }
  
  bool G4GeometryOptions::isVisible( const std::string& volName ) const {
    return queryMap( isVisible_, 
                     orderingIsVisible_, 
                     volName, 
                     defaultIsVisible_ );
  }
  
  bool G4GeometryOptions::doSurfaceCheck( const std::string& volName ) const {
    return queryMap( doSurfaceCheck_, 
                     orderingDoSurfaceCheck_, 
                     volName, 
                     defaultDoSurfaceCheck_ );
  }
  bool G4GeometryOptions::forceAuxEdgeVisible( const std::string& volName ) const {
    return queryMap( forceAuxEdgeVisible_, 
                     orderingForceAuxEdgeVisible_, 
                     volName, 
                     defaultForceAuxEdgeVisible_ );
  }
  bool G4GeometryOptions::placePV( const std::string& volName ) const {
    return queryMap( placePV_, 
                     orderingPlacePV_, 
                     volName, 
                     defaultPlacePV_ );
  }

  //======================================================================
  void G4GeometryOptions::loadOrderingStrings( OrderingList& orderingList,
                                               const SimpleConfig& config, 
                                               const std::string& var ) {
    
    VS vs;
    config.getVectorString( var, vs, VS() );

    if ( vs.empty() ) return;

    for ( const std::string& orderingVar : vs ) {

      const bool keepFlag = std::regex_match( orderingVar, std::regex("(.*)\\.(keep)(.*)") );
      const bool dropFlag = std::regex_match( orderingVar, std::regex("(.*)\\.(drop)(.*)") );

      if ( !keepFlag && !dropFlag )
        throw cet::exception("G4GeometryOptions") << "Ordering lists must be prefixed with \"*.keep\" or \"*.drop\"";
      
      VS volTokens;
      config.getVectorString( orderingVar, volTokens );
      
      if ( volTokens.empty() ) continue;            
      
      // Replace wildcards with suitable regex characters
      std::for_each( volTokens.begin(), 
                     volTokens.end(), 
                     [](std::string& token){
                       boost::replace_all( token, "*", "(.*)" );
                     } );
      
      orderingList.emplace_back( keepFlag, volTokens );
      
    }
    
  }

  //======================================================================
  void G4GeometryOptions::mapInserter( std::map<std::string,bool>& map,
                                       const SimpleConfig& config,
                                       const std::string& volName,
                                       const std::string& var,
                                       const bool default_value ){
    
    auto test = map.emplace( volName,
                              config.getBool( var, default_value )
                              );
    if ( !test.second ) 
      cet::exception("GEOM") << volName << " not successfully inserted into map using var " << var << "\n" ;
    
  }

  //======================================================================
  bool G4GeometryOptions::queryMap( const std::map<std::string,bool>& map,                                 
                                    const OrderingList& orderingList,
                                    const std::string& volName,
                                    const bool default_value ) {

    auto findEntry           = map.find( volName );

    const bool initial_value = ( findEntry != map.end() ) ? findEntry->second : default_value;

    auto overrideFlag        = flagOverridden( volName, orderingList );

    return overrideFlag.first ? overrideFlag.second : initial_value;
  }

  //======================================================================
  std::pair<bool,bool> G4GeometryOptions::flagOverridden( const std::string& volName,
                                                          const OrderingList& orderingList ) {
    
    bool overridden(false), value(false);

    for ( const Ordering& ordering : orderingList ) {
      for ( const std::string& volToken : ordering.second ) {
        if ( std::regex_match( volName, std::regex( volToken ) ) ) {
          overridden = true;
          value      = ordering.first; 
        }
      }
    }

    return std::make_pair( overridden, value );

  }

}
