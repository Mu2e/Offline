//
// Construct Mu2e building
//
//
// Original author: Andrei Gaponenko

// Mu2e include files
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/GeometryService/inc/Mu2eEnvelope.hh"
#include "Offline/GeometryService/inc/Mu2eHallMaker.hh"
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"

// C++ include files
#include <sstream>
#include <iostream>

// Framework include files
#include "cetlib_except/exception.h"

// CLHEP include files
#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {

  std::unique_ptr<Mu2eHall> Mu2eHallMaker::makeBuilding( G4GeometryOptions& geomOptions,
                                                         const SimpleConfig& c )
  {
    std::unique_ptr<Mu2eHall> b (new Mu2eHall());
    loadSolids( b->bldgSolids_, geomOptions, c, "bldg.prefix.list" );
    return b;
  }

  //==================================================================
  void Mu2eHallMaker::makeDirt( Mu2eHall& b,
                                G4GeometryOptions& geomOptions,
                                const SimpleConfig& c,
                                const Mu2eEnvelope& mu2eEnv ) {
    loadSolids           ( b.dirtSolids_, geomOptions, c, "dirt.prefix.list" );
    replaceBoundaryValues( b.dirtSolids_, c, "dirt.prefix.list", "y", mu2eEnv.xmin(), mu2eEnv.xmax() );
    replaceBoundaryValues( b.dirtSolids_, c, "dirt.prefix.list", "x", mu2eEnv.zmin(), mu2eEnv.zmax() );
  }
 //==================================================================
  void Mu2eHallMaker::makeRotated( Mu2eHall& b,
                                G4GeometryOptions& geomOptions,
                                const SimpleConfig& c,
                                const Mu2eEnvelope& mu2eEnv ) {
    loadRotSolids           ( b.rotatedSolids_, geomOptions, c, "rotated.prefix.list" );
  }

  //==================================================================
  void Mu2eHallMaker::makeTrapDirt( Mu2eHall& b,
                                    G4GeometryOptions& geomOptions,
                                    const SimpleConfig& c,
                                    const Mu2eEnvelope& mu2eEnv ) {
    loadTrapSolids      ( b.dirtTrapSolids_, geomOptions, c, "dirt.trap.prefix.list" );
    replaceBoundaryValues( b.dirtTrapSolids_, c, "dirt.trap.prefix.list", "y", mu2eEnv.xmin(), mu2eEnv.xmax() );
    replaceBoundaryValues( b.dirtTrapSolids_, c, "dirt.trap.prefix.list", "x", mu2eEnv.zmin(), mu2eEnv.zmax() );
  }

  //==================================================================
  void Mu2eHallMaker::loadSolids( std::map<std::string,ExtrudedSolid>& solidMap,
                                  G4GeometryOptions& geomOptions,
                                  const SimpleConfig& c,
                                  const std::string& varPrefixStr )
  {
    std::vector<std::string> varNames;
    c.getVectorString( varPrefixStr, varNames );

    for ( const auto& prefix : varNames ) {

      CLHEP::Hep3Vector offset
        (
         c.getDouble( prefix+".offsetFromMu2eOrigin.x" ),
         c.getDouble( prefix+".offsetFromFloorSurface.y" )+c.getDouble( "yOfFloorSurface.below.mu2eOrigin"),
         c.getDouble( prefix+".offsetFromMu2eOrigin.z" )
         );

      std::vector<double> x,y;
      c.getVectorDouble( prefix+".xPositions", x );
      c.getVectorDouble( prefix+".yPositions", y );

      const std::string volName = c.getString( prefix+".name" );

      std::string loadPrefix = prefix;
      std::string dot = ".";
      std::size_t place1 = prefix.find(dot);
      std::size_t place2 = std::string::npos;
      if ( place1 != std::string::npos ) place2 = prefix.find(dot,place1+1);
      if ( place2 != std::string::npos ) loadPrefix = prefix.substr(0,place2);

      solidMap[volName] = ExtrudedSolid( volName,
                                         c.getString( prefix+".material"),
                                         offset,
                                         c.getDouble( prefix+".yHalfThickness" ),
                                         getPairedVector(x,y) );

      geomOptions.loadEntry( c, volName, loadPrefix );

    }

  }


  //==================================================================
  std::vector<CLHEP::Hep2Vector>
  Mu2eHallMaker::getPairedVector( const std::vector<double>& x,
                                  const std::vector<double>& y ){

    assert ( x.size() == y.size() );

    std::vector<CLHEP::Hep2Vector> vCLHEP;

    for ( std::size_t i(0) ; i < x.size() ; ++i ) {
      vCLHEP.emplace_back( x.at(i), y.at(i) );
    }

    return vCLHEP;

  }

  //==================================================================
  void Mu2eHallMaker::loadRotSolids( std::map<std::string,RotExtrudedSolid>& solidMap,
                                      G4GeometryOptions& geomOptions,
                                      const SimpleConfig& c,
                                      const std::string& varPrefixStr )
  {
    std::vector<std::string> varNames;
    c.getVectorString( varPrefixStr, varNames, varNames ); //default is empty list

    for ( const auto& prefix : varNames ) {
      CLHEP::Hep3Vector offset
        (
         c.getDouble( prefix+".offsetFromMu2eOrigin.x" ),
         c.getDouble( prefix+".offsetFromFloorSurface.y" )+c.getDouble( "yOfFloorSurface.below.mu2eOrigin"),
         c.getDouble( prefix+".offsetFromMu2eOrigin.z" )
         );
      std::vector<double> x,y, angles;
      c.getVectorDouble( prefix+".xPositions", x );
      c.getVectorDouble( prefix+".yPositions", y );
      c.getVectorDouble( prefix+".angles", angles, {0.,0.,0.}, 3); //if given rotation angles, must have all 3


      CLHEP::HepRotation rot(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
      rot.rotateX(angles[0]);
      rot.rotateY(angles[1]);
      rot.rotateZ(angles[2]);
      const std::string volName = c.getString( prefix+".name" );

      std::string loadPrefix = prefix;
      std::string dot = ".";
      std::size_t place1 = prefix.find(dot);
      std::size_t place2 = std::string::npos;
      if ( place1 != std::string::npos ) place2 = prefix.find(dot,place1+1);
      if ( place2 != std::string::npos ) loadPrefix = prefix.substr(0,place2);

      solidMap[volName] = RotExtrudedSolid( volName,
                                       c.getString( prefix+".material"),
                                       offset,
                                       c.getDouble( prefix+".yHalfThickness" ),
                                       getPairedVector(x,y),
                                       rot);


      geomOptions.loadEntry( c, volName, loadPrefix );

    }

  }

  //==================================================================
  void Mu2eHallMaker::loadTrapSolids( std::map<std::string,GenericTrap>& solidMap,
                                      G4GeometryOptions& geomOptions,
                                      const SimpleConfig& c,
                                      const std::string& varPrefixStr )
  {
    std::vector<std::string> varNames;
    c.getVectorString( varPrefixStr, varNames, varNames ); //default is empty list

    for ( const auto& prefix : varNames ) {
      CLHEP::Hep3Vector offset
        (
         c.getDouble( prefix+".offsetFromMu2eOrigin.x" ),
         c.getDouble( prefix+".offsetFromFloorSurface.y" )+c.getDouble( "yOfFloorSurface.below.mu2eOrigin"),
         c.getDouble( prefix+".offsetFromMu2eOrigin.z" )
         );
      std::vector<double> x,y, angles;
      c.getVectorDouble( prefix+".xPositions", x );
      c.getVectorDouble( prefix+".yPositions", y );
      c.getVectorDouble( prefix+".angles", angles, {0.,0.,0.}, 3); //if given rotation angles, must have all 3


      CLHEP::HepRotation rot(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
      rot.rotateX(angles[0]);
      rot.rotateY(angles[1]);
      rot.rotateZ(angles[2]);
      const std::string volName = c.getString( prefix+".name" );

      std::string loadPrefix = prefix;
      std::string dot = ".";
      std::size_t place1 = prefix.find(dot);
      std::size_t place2 = std::string::npos;
      if ( place1 != std::string::npos ) place2 = prefix.find(dot,place1+1);
      if ( place2 != std::string::npos ) loadPrefix = prefix.substr(0,place2);

      solidMap[volName] = GenericTrap( volName,
                                       c.getString( prefix+".material"),
                                       offset,
                                       c.getDouble( prefix+".yHalfThickness" ),
                                       getPairedVector(x,y),
                                       rot);


      geomOptions.loadEntry( c, volName, loadPrefix );

    }

  }

  //==================================================================
  void Mu2eHallMaker::replaceBoundaryValues(  std::map<std::string,ExtrudedSolid>& dirtMap,
                                              const SimpleConfig& c,
                                              const std::string& varPrefixStr,
                                              const std::string& dim,
                                              const double min, const double max ) {

    std::vector<std::string> varNames;
    c.getVectorString( varPrefixStr, varNames );

    for ( const auto& prefix : varNames ) {

      std::vector<int> vr;
      c.getVectorInt( prefix+"."+dim+"replace", vr, std::vector<int>() );
      if ( vr.empty() ) continue;

      const std::string volName = c.getString( prefix+".name" );
      auto vol = dirtMap.find( volName );

      if ( vol == dirtMap.end() ) throw cet::exception("GEOM") << "Dirt volume << " <<  volName << " >> not found!\n";

      for ( const int rindex : vr ) {

        CLHEP::Hep2Vector& vertex = vol->second.modifyVertex(rindex);

        // The offsets need to be removed because they were included
        // in the definiton of the min/max values for Mu2eEnvelope.
        //
        // IMPORTANT NOTE!
        //
        //    - This replacement works assuming that each vol. offset is equal to that
        //      which determined the mu2eEnvelope.  This is a reasonable
        //      assumption as long as no one did anything sinister by using
        //      different offsets from Mu2e origin for each building volume!

        const double zOffset =  vol->second.getOffsetFromMu2eOrigin().z();
        if ( dim=="x" ) vertex.setX( vertex.x() < 0 ? min-zOffset : max-zOffset );

        const double xOffset =  vol->second.getOffsetFromMu2eOrigin().x();
        if ( dim=="y" ) vertex.setY( vertex.y() < 0 ? min-xOffset : max-xOffset );
      }

    }

  }

  //==================================================================
  void Mu2eHallMaker::replaceBoundaryValues(  std::map<std::string,GenericTrap>& dirtMap,
                                              const SimpleConfig& c,
                                              const std::string& varPrefixStr,
                                              const std::string& dim,
                                              const double min, const double max ) {

    std::vector<std::string> varNames;
    c.getVectorString( varPrefixStr, varNames,varNames ); //defult is empty list

    for ( const auto& prefix : varNames ) {

      std::vector<int> vr;
      c.getVectorInt( prefix+"."+dim+"replace", vr, std::vector<int>() );
      if ( vr.empty() ) continue;

      const std::string volName = c.getString( prefix+".name" );
      auto vol = dirtMap.find( volName );

      if ( vol == dirtMap.end() ) throw cet::exception("GEOM") << "Dirt volume << " <<  volName << " >> not found!\n";

      for ( const int rindex : vr ) {

        CLHEP::Hep2Vector& vertex = vol->second.modifyVertex(rindex);

        // The offsets need to be removed because they were included
        // in the definiton of the min/max values for Mu2eEnvelope.
        //
        // IMPORTANT NOTE!
        //
        //    - This replacement works assuming that each vol. offset is equal to that
        //      which determined the mu2eEnvelope.  This is a reasonable
        //      assumption as long as no one did anything sinister by using
        //      different offsets from Mu2e origin for each building volume!

        const double zOffset =  vol->second.getOffsetFromMu2eOrigin().z();
        if ( dim=="x" ) vertex.setX( vertex.x() < 0 ? min-zOffset : max-zOffset );

        const double xOffset =  vol->second.getOffsetFromMu2eOrigin().x();
        if ( dim=="y" ) vertex.setY( vertex.y() < 0 ? min-xOffset : max-xOffset );
      }

    }

  }
}
