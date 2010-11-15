#ifndef Mu2eWorld_H
#define Mu2eWorld_H 1
//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eWorld.hh,v 1.22 2010/11/15 23:27:53 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/11/15 23:27:53 $
//
// Original author Rob Kutschke
//
// Notes
// 1) The variable _volumeInfoList, holds some pointers to information
//    about volumes.  It also holds some position information.  The
//    data member centerInWorld is not always meaningful.  If you follow
//    the trail from a given volume back to the world volume, if you
//    do not find any rotations, then this data member is meaningful.
//    If you do find a rotation, then the data member is not meaningful.
//    This data member will only be filled for some of the upper level
//    volumes.   It won't be filled for straws and crystals.  It's purpose
//    is to make it easier to break up one giant method into many smaller
//    ones, by allowing code to look up its mother volume by name.
//    The bottom line is that this is not a fully general facility and
//    must be used with care.

#include <string>
#include <memory>
#include <vector>
#include <map>

// Forward references.
class G4Material;
class G4Mag_UsualEqRhs;
class G4UserLimits;

// Mu2e includes
#include "Mu2eG4/inc/WorldInfo.hh"
#include "Mu2eG4/inc/VolumeInfo.hh"
#include "Mu2eG4/inc/FieldMgr.hh"
#include "TrackerGeom/inc/TubsParams.hh"

//G4 includes 
#include "G4String.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

namespace mu2e {

  // Forward references within mu2e namespace.
  class SimpleConfig;

  class Mu2eWorld {
  public:
    
    Mu2eWorld();
    ~Mu2eWorld();

    // Construct everything.
    WorldInfo const* construct();

    WorldInfo const& getWorldInfo() const { return _info; }

    G4ThreeVector const& getCosmicReferencePoint() const{
      return _cosmicReferencePoint;
    }

    G4ThreeVector const& getMu2eOrigin() const{
      return _mu2eOrigin;
    }

    G4ThreeVector const& getMu2eDetectorOrigin() const{
      return _mu2eDetectorOrigin;
    }

    G4ThreeVector const& getPrimaryProtonGunOrigin() const{
      return _primaryProtonGunOrigin;
    }
 
    G4RotationMatrix const& getPrimaryProtonGunRotation() const{
      return _primaryProtonGunRotation;
    }

  private:

    // Do all of the work.
    void constructWorld();

    // Break the big task into many smaller ones.
    void defineMu2eOrigin();
    VolumeInfo constructDirt();
    VolumeInfo constructHall( const VolumeInfo& parent );
    void constructDS( const VolumeInfo& parent );
    void constructTS( const VolumeInfo& parent );
    void constructPS( const VolumeInfo& parent );
    void constructVD( );
    VolumeInfo constructTracker();
    VolumeInfo constructTarget();
    void constructProtonAbs();
    void constructCal();
    void constructMagnetYoke();
    void constructCRV();
    void constructSteel( const VolumeInfo& parent );
    void constructBFieldAndManagers();  // MECO field maps
    void constructBFieldAndManagers2(); // G4BL field maps
    void constructStepLimiters();

    // The world coordinates of the center of the cosmic ray reference plane.
    G4ThreeVector _cosmicReferencePoint;

    // The world coordinates of the origin of the Mu2e coordinate system.
    G4ThreeVector _mu2eOrigin;

    // World coorindates of the reference point for building the detectors:
    //  - on axis on the DS with z=12,000. in the Mu2e coordinate system.
    G4ThreeVector _mu2eDetectorOrigin;

    // Origin of the hall air volume in the Mu2e coordinate system.
    G4ThreeVector _hallOriginInMu2e;

    // Information about the world that can be passed to others.
    WorldInfo _info;

    // Utility functions.
    void setUnits( std::vector<double>& V, G4double unit );
    VolumeInfo nestBox ( std::string const& name,
                         double const halfDim[3],
                         G4Material* material,
                         G4RotationMatrix* rot,
                         G4ThreeVector const& offset,
                         const VolumeInfo& parent,
                         int copyNo,
                         bool isVisible = false,
                         G4Colour color = G4Colour::Black(),
                         bool forceSolid = false
                         );

    // Alternate argument list, using a vector for the half dimensions.
    VolumeInfo nestBox ( std::string const& name,
                         std::vector<double> const&  halfDim,
                         G4Material* material,
                         G4RotationMatrix* rot,
                         G4ThreeVector const& offset,
                         const VolumeInfo& parent,
                         int copyNo,
                         bool isVisible = false,
                         G4Colour color = G4Colour::Black(),
                         bool forceSolid = false
                         ){
      return nestBox( name, 
                      &halfDim[0],
                      material,
                      rot,
                      offset,
                      parent,
                      copyNo,
                      isVisible,
                      color,
                      forceSolid
                      );
    }

    VolumeInfo nestTubs2 ( std::string const& name,
                           double params[5],
                           G4Material* material,
                           G4RotationMatrix* rot,
                           const G4ThreeVector& offset,
                           const VolumeInfo& parent,
                           int copyNo,
                           bool isVisible = false,
                           G4Colour color = G4Colour::Black(),
                           bool forceSolid = false
                           );
  


    // Alternate argument list, using a vector for the parameters.
    inline VolumeInfo nestTubs2 ( std::string const& name,
                                  std::vector<double>&  params,
                                  G4Material* material,
                                  G4RotationMatrix* rot,
                                  const G4ThreeVector& offset,
                                  const VolumeInfo& parent,
                                  int copyNo,
                                  bool isVisible = false,
                                  G4Colour color = G4Colour::Black(),
                                  bool forceSolid = false
                                  ){
      return nestTubs2( name, 
                        &params[0],
                        material,
                        rot,
                        offset,
                        parent,
                        copyNo,
                        isVisible,
                        color,
                        forceSolid
                        );
    }

    // Alternate argument list, using a TubsParams object for the parameters.
    inline VolumeInfo nestTubs2 ( std::string const& name,
                                  TubsParams& params,
                                  G4Material* material,
                                  G4RotationMatrix* rot,
                                  const G4ThreeVector& offset,
                                  const VolumeInfo& parent,
                                  int copyNo,
                                  bool isVisible = false,
                                  G4Colour color = G4Colour::Black(),
                                  bool forceSolid = false
                                  ){
      return nestTubs2( name, 
                        &params.innerRadius,
                        material,
                        rot,
                        offset,
                        parent,
                        copyNo,
                        isVisible,
                        color,
                        forceSolid
                        );
    }

    VolumeInfo nestCons2 ( std::string const& name,
                           double params[7],
                           G4Material* material,
                           G4RotationMatrix* rot,
                           const G4ThreeVector& offset,
                           const VolumeInfo& parent,
                           int copyNo,
                           bool isVisible = false,
                           G4Colour color = G4Colour::Black(),
                           bool forceSolid = false
                           );
  
    VolumeInfo nestExtrudedSolid2 ( std::string const& name,
				    double hz,
				    std::vector<double> &x,
				    std::vector<double> &y,
				    G4Material* material,
				    G4RotationMatrix* rot,
				    const G4ThreeVector& offset,
				    const VolumeInfo& parent,
				    int copyNo,
				    bool isVisible = false,
				    G4Colour color = G4Colour::Black(),
				    bool forceSolid = false
				    );
  
    VolumeInfo nestTorus2 ( std::string const& name,
                            double halfDim[5],
                            G4Material* material,
                            G4RotationMatrix* rot,
                            G4ThreeVector const& offset,
                            const VolumeInfo& parent,
                            int copyNo,
                            bool isVisible = false,
                            G4Colour color = G4Colour::Black(),
                            bool forceSolid = false
                            );
  


    // Alternate argument list, using a vector for the half dimensions.
    inline VolumeInfo nestTorus2 ( std::string const& name,
                                   std::vector<double>&  halfDim,
                                   G4Material* material,
                                   G4RotationMatrix* rot,
                                   G4ThreeVector& offset,
                                   const VolumeInfo& parent,
                                   int copyNo,
                                   bool isVisible = false,
                                   G4Colour color = G4Colour::Black(),
                                   bool forceSolid = false
                                   ){
      return nestTorus2( name, 
                         &halfDim[0],
                         material,
                         rot,
                         offset,
                         parent,
                         copyNo,
                         isVisible,
                         color,
                         forceSolid
                         );
    }

    // Versions of the map [] operator that check for errors.
    VolumeInfo& locateVolInfo( const std::string key);
    void addVolInfo( const VolumeInfo& info );
    
    // Stash a pointer to the config object so that all methods can get at it easily.
    SimpleConfig const* _config;

    // Location in G4 world coordinates of the reference point for the Primary Proton Gun
    G4ThreeVector _primaryProtonGunOrigin;
    G4RotationMatrix _primaryProtonGunRotation;

    // Models of the DS magnetic field:
    // 0 - whole DS uses the field map.
    // 1 - upstream uses the full field map; downstream uses a uniform field.
    // 2 - whole DS uses a uniform field.
    enum DSFieldModel { dsModelFull, dsModelSplit, dsModelUniform};

    // Field managers for the different regions of magnetic field.
    // These have a lifetime equal to that of the G4 geometry.
    std::auto_ptr<FieldMgr> _dsFull;
    std::auto_ptr<FieldMgr> _dsUniform;
    std::auto_ptr<FieldMgr> _psFull;
    std::auto_ptr<FieldMgr> _tsFull;
    std::auto_ptr<FieldMgr> _hallFull;

    // Allow access to the volume information by volume name.  See note 1.
    std::map<std::string,VolumeInfo> _volumeInfoList;

  };

} // end namespace mu2e
#endif
