//
// Free function to create a geant4 test environment geometry
//
//
// Original author KLG 
//
// Notes:
//
// one can nest volume inside other volumes if needed
// see other construct... functions for examples
//

#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.

#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/InitEnvToolBase.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  class ConstructEnvCalo01: public InitEnvToolBase {
  public:
    ConstructEnvCalo01(const fhicl::ParameterSet& PSet);
    ~ConstructEnvCalo01();

    int construct(VolumeInfo const& ParentVInfo, SimpleConfig const& Config);
  private:
    
    void constructModule(std::string         const&  mNamePrefix,
			 VolumeInfo          const&  parentVInfo,
			 std::vector<double> const&  tHL,
			 std::vector<double> const&  tCInParent,
			 vector<double>      const&  mHL,
			 vector<string>      const&  mMat,
			 G4int                       mNumberOfLayers,
			 G4int                       verbosityLevel,
			 G4double                    mLCenterInParent,
			 G4int                       passiveVolumeStartingCopyNumber,
			 G4int                       activeVolumeStartingCopyNumber,
			 bool                const   isVisible,
			 G4Colour            const&  passiveVolumeColour,
			 G4Colour            const&  activeVolumeColour,
			 bool                const   forceSolid,
			 bool                const   forceAuxEdgeVisible,
			 bool                const   placePV,
			 bool                const   doSurfaceCheck
			 ); 
  };

//-----------------------------------------------------------------------------
  ConstructEnvCalo01::ConstructEnvCalo01(const fhicl::ParameterSet& PSet) {
    _name = "Calo01";
  }

//-----------------------------------------------------------------------------
  ConstructEnvCalo01::~ConstructEnvCalo01() {
    _name = "Calo01";
  }

//-----------------------------------------------------------------------------
  int ConstructEnvCalo01::construct(VolumeInfo const& parentVInfo, SimpleConfig const& config) {

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    // Extract calorimeter information from the config file and construct it

    // it will be a set of thin plates (G4Box'es)

    G4int  verbosityLevel  = config.getInt("calo.verbosityLevel",-1);
    G4bool isVisible     = config.getBool("calo.visible",true);
    G4bool forceSolid      = config.getBool("calo.solid",true);

    vector<double> tHL;
    config.getVectorDouble( "calo.transverseHalfLengths", tHL);

    vector<double> tCInParent;
    config.getVectorDouble( "calo.transverseCenterInWorld", tCInParent);

    // front Scintillator (just one layer, we could have used simple variables)

    vector<double> fSHL;
    config.getVectorDouble( "calo.frontScintLayerHalfLengths",fSHL);

    vector<string> fSMat;
    config.getVectorString( "calo.frontScintLayerMaterials", fSMat);

    G4Material* fSMaterial = findMaterialOrThrow(fSMat[0]);

    G4int activeVolumeCopyNumber = config.getInt("calo.activeVolumeStartingCopyNumber");

    G4double const fSParams[] = {tHL[0], tHL[1], fSHL[0]};

    G4double fSLCenterInParent = config.getDouble("calo.frontScintStartingLongitPosition");

    G4ThreeVector fSCenterInParent(tCInParent[0],tCInParent[1],fSLCenterInParent+fSHL[0]);

    //    G4Colour  orange  (.75, .55, .0);

    string vName("fS");

    if (verbosityLevel > 0 ) {
      cout << __func__ << " constructing: " << vName 
           << " " << "at " << fSCenterInParent
           << endl;
    }

    VolumeInfo fsVInfo(nestBox( vName,
                                fSParams,
                                fSMaterial,
                                0x0, // no rotation
                                fSCenterInParent,
                                parentVInfo,
                                activeVolumeCopyNumber++,
                                // non 0 copy nuber for volume tracking purposes
                                isVisible,
                                G4Colour::Yellow(),
                                forceSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck
                                ));

    // front module

    vector<double> fMHL;
    config.getVectorDouble( "calo.frontModuleLayerHalfLengths",fMHL);

    vector<string> fMMat;
    config.getVectorString( "calo.frontModuleLayerMaterials", fMMat);

    G4int passiveVolumeCopyNumber = config.getInt("calo.passiveVolumeStartingCopyNumber");
    G4double    fMLCenterInParent = config.getDouble("calo.frontModuleStartingLongitPosition");
    G4int        fMNumberOfLayers = config.getInt("calo.frontModuleNumberOfLayers");


    constructModule("fm",
		    parentVInfo,
		    tHL,
		    tCInParent,
		    fMHL,
		    fMMat,
		    fMNumberOfLayers,
		    verbosityLevel,
		    fMLCenterInParent,
		    passiveVolumeCopyNumber,
		    activeVolumeCopyNumber,
		    isVisible,
		    G4Colour::Gray(),
		    G4Colour::Yellow(),
		    forceSolid,
		    forceAuxEdgeVisible,
		    placePV,
		    doSurfaceCheck
		    );

    // rear module

    vector<double> rMHL;
    config.getVectorDouble( "calo.rearModuleLayerHalfLengths",rMHL);

    vector<string> rMMat;
    config.getVectorString( "calo.rearModuleLayerMaterials", rMMat);

    // we need to "calculate" the starting copy numbers

    passiveVolumeCopyNumber += fMNumberOfLayers;
    activeVolumeCopyNumber += fMNumberOfLayers;

    G4double    rMLCenterInParent = config.getDouble("calo.rearModuleStartingLongitPosition");
    G4int        rMNumberOfLayers = config.getInt("calo.rearModuleNumberOfLayers");

    constructModule("rm",
		    parentVInfo,
		    tHL,
		    tCInParent,
		    rMHL,
		    rMMat,
		    rMNumberOfLayers,
		    verbosityLevel,
		    rMLCenterInParent,
		    passiveVolumeCopyNumber,
		    activeVolumeCopyNumber,
		    isVisible,
		    G4Colour::Brown(),
		    G4Colour::Yellow(),
		    forceSolid,
		    forceAuxEdgeVisible,
		    placePV,
		    doSurfaceCheck
		    );
    return 0;
  }

  void ConstructEnvCalo01::constructModule(std::string const & mNamePrefix,
					   VolumeInfo const & parentVInfo,
					   std::vector<double> const & tHL,
					   std::vector<double> const & tCInParent,
					   vector<double> const & mHL,
					   vector<string> const & mMat,
					   G4int    mNumberOfLayers,
					   G4int    verbosityLevel,
					   G4double mLCenterInParent,
					   G4int passiveVolumeStartingCopyNumber,
					   G4int  activeVolumeStartingCopyNumber,
					   bool const isVisible,
					   G4Colour const & passiveVolumeColour,
					   G4Colour const &  activeVolumeColour,
					   bool const forceSolid,
					   bool const forceAuxEdgeVisible,
					   bool const placePV,
					   bool const doSurfaceCheck
					   ) { 
    //
    //  construct a calorimetric module (could be made to a n-layerd if needed)
    //  ( P - Passive and A - Active layers)

    G4int passiveVolumeCopyNumber = passiveVolumeStartingCopyNumber;
    G4int activeVolumeCopyNumber  = activeVolumeStartingCopyNumber;

    G4Material* mPMaterial = findMaterialOrThrow(mMat[0]); // passive
    G4Material* mAMaterial = findMaterialOrThrow(mMat[1]); // active

    G4double const mPParams[] = {tHL[0], tHL[1], mHL[0]}; // passive
    G4double const mAParams[] = {tHL[0], tHL[1], mHL[1]}; // active

    // we place all the volumes directly "in" the parent

    ostringstream vsPNumber("");
    vsPNumber.width(3);
    vsPNumber.fill('0');
    string vPName;

    ostringstream vsANumber("");
    vsANumber.width(3);
    vsANumber.fill('0');
    string vAName;

    G4double      mLayerStep  = 2.*(mHL[0]+mHL[1]);

    // we initialize the Z position to be one "step earlier" where it
    // should eventually be and then add the step to position of the
    // individual layers

    G4double      mPZCenterInParent = mLCenterInParent + mHL[0] - mLayerStep;
    G4double      mAZCenterInParent = mLCenterInParent - mHL[1];

    G4ThreeVector mPCenterInParent(tCInParent[0],tCInParent[1],mPZCenterInParent);
    G4ThreeVector mACenterInParent(tCInParent[0],tCInParent[1],mAZCenterInParent);

    for( G4int nl = 0; nl<mNumberOfLayers; ++nl ) {

      vsPNumber.str("");
      vsPNumber << passiveVolumeCopyNumber;
      vPName = mNamePrefix + "PL_" + vsPNumber.str();

      mPZCenterInParent += mLayerStep;
      mPCenterInParent.setZ(mPZCenterInParent);

      vsANumber.str("");
      vsANumber << activeVolumeCopyNumber;
      vAName = mNamePrefix + "AL_" + vsANumber.str();

      mAZCenterInParent += mLayerStep;
      mACenterInParent.setZ(mAZCenterInParent);

      if (verbosityLevel > 0 ) {
        cout << __func__ << " constructing: " << vPName 
             << " " << "at " << mPCenterInParent 
             << endl;
        cout << __func__ << " constructing: " << vAName 
             << " " << "at " << mACenterInParent 
             << endl;
      }

      // passive medium
      VolumeInfo mPVInfo(nestBox( vPName,
                                   mPParams,
                                   mPMaterial,
                                   0x0, // no rotation
                                   mPCenterInParent,
                                   parentVInfo,
                                   passiveVolumeCopyNumber++,
                                   // non 0 for volume tracking purposes
                                   isVisible,
                                   passiveVolumeColour,
                                   forceSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));

      // active medium
      VolumeInfo mAVInfo(nestBox( vAName,
                                   mAParams,
                                   mAMaterial,
                                   0x0, // no rotation
                                   mACenterInParent,
                                   parentVInfo,
                                   activeVolumeCopyNumber++,
                                   // non 0 for volume tracking purposes
                                   isVisible,
                                   activeVolumeColour,
                                   forceSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));

    }

  } // constructModule
}

DEFINE_ART_CLASS_TOOL(mu2e::ConstructEnvCalo01)
