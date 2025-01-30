#ifndef GeometryService_Mu2eHallMaker_hh_HH
#define GeometryService_Mu2eHallMaker_hh_HH

#include "CLHEP/Vector/TwoVector.h"

#include <memory>

namespace mu2e {

  class ExtrudedSolid;
  class RotExtrudedSolid;
  class GenericTrap;
  class SimpleConfig;
  class Mu2eEnvelope;
  class G4GeometryOptions;

  class Mu2eHallMaker {
  public:

    static std::unique_ptr<Mu2eHall>
    makeBuilding(G4GeometryOptions& geomOptions,const SimpleConfig& config);

    static void
    makeDirt( Mu2eHall& mh,
              G4GeometryOptions& geomOptions,
              const SimpleConfig& config,
              const Mu2eEnvelope& mu2eEnv );
    static void
    makeRotated( Mu2eHall& mh,
              G4GeometryOptions& geomOptions,
              const SimpleConfig& config,
              const Mu2eEnvelope& mu2eEnv );
    static void
    makeTrapDirt( Mu2eHall& mh,
                  G4GeometryOptions& geomOptions,
                  const SimpleConfig& config,
                  const Mu2eEnvelope& mu2eEnv );

    static void
    loadSolids( std::map<std::string,ExtrudedSolid>& solidMap,
                G4GeometryOptions& geomOptions,
                const SimpleConfig& config,
                const std::string& varPrefixStr );
    static void
    loadRotSolids( std::map<std::string,RotExtrudedSolid>& solidMap,
                G4GeometryOptions& geomOptions,
                const SimpleConfig& config,
                const std::string& varPrefixStr );

    static void
    loadTrapSolids( std::map<std::string,GenericTrap>& solidMap,
                    G4GeometryOptions& geomOptions,
                    const SimpleConfig& config,
                    const std::string& varPrefixStr );

    static std::vector<CLHEP::Hep2Vector>
    getPairedVector( const std::vector<double>& v1,
                     const std::vector<double>& v2 );

    static void
    replaceBoundaryValues( std::map<std::string,ExtrudedSolid>& solidMap,
                           const SimpleConfig& config,
                           const std::string& varPrefixStr,
                           const std::string& dim,
                           const double min, const double max );
    static void
    replaceBoundaryValues( std::map<std::string,GenericTrap>& solidMap,
                           const SimpleConfig& config,
                           const std::string& varPrefixStr,
                           const std::string& dim,
                           const double min, const double max );
  };
}

#endif/* GeometryService_Mu2eHallMaker_hh_HH */
