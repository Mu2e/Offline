// G4 geometry options look-up facility, to be used in conjunction with SimpleConfig.
//
//
// Original author: Kyle Knoepfel, modified by B. Echenard

#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "cetlib_except/exception.h"

#include "boost/algorithm/string/replace.hpp"
#include <regex>



namespace mu2e {


  //======================================================================
  G4GeometryOptData::G4GeometryOptData(bool defaultValue, const std::string& name) 
    : name_(name),map_(),ordering_(),default_(defaultValue)
  {}

  //======================================================================
  void G4GeometryOptData::loadOrderingStrings(const SimpleConfig& config, const std::string& varString) 
  {
      VS vs;
      config.getVectorString( varString, vs, VS() );
      if ( vs.empty() ) return;

      for (const auto& orderingVar : vs )
      {
          const bool keepFlag = std::regex_match( orderingVar, std::regex("(.*)\\.(keep)(.*)") );
          const bool dropFlag = std::regex_match( orderingVar, std::regex("(.*)\\.(drop)(.*)") );

          if ( !keepFlag && !dropFlag )
            throw cet::exception("G4GeometryOptions") << "Ordering lists must be prefixed with \"*.keep\" or \"*.drop\"";

          VS volTokens;
          config.getVectorString( orderingVar, volTokens );

          if ( volTokens.empty() ) continue;            

          // Replace wildcards with suitable regex characters and store regex object
          for (auto& volToken : volTokens)
          {
              boost::replace_all( volToken, "*", "(.*)" );
              ordering_.emplace_back(keepFlag, std::regex(volToken));            
          }
      }     
  }

  
  //======================================================================
  void G4GeometryOptData::mapInserter(const std::string& volName, bool value)
  {    
     if (map_.find(volName)==map_.end()) map_.emplace(volName, value);
  }


  //======================================================================
  bool G4GeometryOptData::queryMap(const std::string& volName) const
  {
     auto findEntry = map_.find(volName);
     //if (findEntry == map_.end())  std::cout<<"Warning on G4GeometryOptData, "<<volName<<" not defined in "<<name_<<std::endl;

     const bool initial_value = ( findEntry != map_.end() ) ? findEntry->second : default_;

     auto overrideFlag = flagOverridden(volName);     
     return overrideFlag.first ? overrideFlag.second : initial_value;
  }


  //======================================================================
  std::pair<bool,bool> G4GeometryOptData::flagOverridden(const std::string& volName) const
  {        
     bool overridden(false), value(false);

     for (const Ordering& ordering : ordering_)
     {
         if (std::regex_match(volName,ordering.second))
         {
            overridden = true;
            value      = ordering.first; 
         }           
     }

     return std::make_pair(overridden,value);
  }











  //======================================================================
  G4GeometryOptions::G4GeometryOptions( const SimpleConfig& config ) 
    : dataSurfaceCheck_ ( config.getBool("g4.doSurfaceCheck"      ), "doSurfaceCheck" )
    , dataIsVisible_    ( config.getBool("g4.visible"           ),   "visible" )
    , dataIsSolid_      ( config.getBool("g4.solid"             ),   "solid" )
    , dataForceAuxEdge_ ( config.getBool("g4.forceAuxEdgeVisible" ), "forceAuxEdgeVisible" )
    , dataPlacePV_      ( config.getBool("g4.placePV"             ), "placePV" )
    
  {    
    dataSurfaceCheck_.loadOrderingStrings( config,"g4.doSurfaceCheck.order"      );
    dataIsVisible_.loadOrderingStrings   ( config,"g4.visible.order"           );
    dataIsSolid_.loadOrderingStrings     ( config,"g4.solid.order"             );
    dataForceAuxEdge_.loadOrderingStrings( config,"g4.forceAuxEdgeVisible.order" );
    dataPlacePV_.loadOrderingStrings     ( config,"g4.placePV.order" );
  }

  //======================================================================
  void G4GeometryOptions::loadEntry( const SimpleConfig& config, const std::string& volName, const std::string& prefix )                                      
  {
    dataSurfaceCheck_.mapInserter( volName, config.getBool(prefix+".doSurfaceCheck",     dataSurfaceCheck_.default_value()) ); 
    dataIsVisible_.mapInserter   ( volName, config.getBool(prefix+".visible",            dataIsVisible_.default_value())    ); 
    dataIsSolid_.mapInserter     ( volName, config.getBool(prefix+".solid",              dataIsSolid_.default_value())      ); 
    dataForceAuxEdge_.mapInserter( volName, config.getBool(prefix+".forceAuxEdgeVisible",dataForceAuxEdge_.default_value()) ); 
    dataPlacePV_.mapInserter     ( volName, config.getBool(prefix+".placePV",            dataPlacePV_.default_value()) ); 
  }
  
  //======================================================================
  bool G4GeometryOptions::doSurfaceCheck( const std::string& volName ) const 
  {
    return dataSurfaceCheck_.queryMap(volName);
  }

  bool G4GeometryOptions::isVisible( const std::string& volName ) const 
  {
    return dataIsVisible_.queryMap(volName);
  }

  bool G4GeometryOptions::isSolid( const std::string& volName ) const 
  {
    return dataIsSolid_.queryMap(volName);
  }
  
  bool G4GeometryOptions::forceAuxEdgeVisible( const std::string& volName ) const 
  {
    return dataForceAuxEdge_.queryMap(volName);
  }
  bool G4GeometryOptions::placePV( const std::string& volName ) const 
  {
    return dataPlacePV_.queryMap(volName);
  }
  

}

