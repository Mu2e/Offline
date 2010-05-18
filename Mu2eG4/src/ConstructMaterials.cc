//
// Construct materials requested by the run-time configuration system.
//
// $Id: ConstructMaterials.cc,v 1.3 2010/05/18 21:16:11 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/05/18 21:16:11 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) This code new's many objects of type G4Material.  The lifeime
//    of these objects is controlled within G4.  We must never
//    delete them.
//

// C++ includes
#include <iostream>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/ConstructMaterials.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// G4 includes
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"

using namespace std;

namespace mu2e {

  // Value used to request that all Mu2e specific materials be made.
  static const std::string DoAllValue = "DoAll";

  ConstructMaterials::ConstructMaterials(){
  }
  
  ConstructMaterials::~ConstructMaterials(){
  }

  // This is the main function.
  void ConstructMaterials::construct(){

    // Get access to the master geometry system and its run time config.
    edm::Service<GeometryService> geom;
    SimpleConfig const& config = geom->config();

    // Construct the requested materials.
    constructG4Materials( config );
    constructMu2eMaterials( config );
    
    // Print element table, if requested.
    if ( config.getBool("g4.printElements",false) ){
      edm::LogInfo  log("GEOM");
      log << *G4Element::GetElementTable();
    }
    
    // Print material table, if requested.
    if ( config.getBool("g4.printMaterials",false) ){
      edm::LogInfo  log("GEOM");
      log << *G4Material::GetMaterialTable();
    }
    
    //cout << "Big test: " << endl;
    //cout << *::G4Element::GetElementTable();
  }
  
  void ConstructMaterials::constructG4Materials( SimpleConfig const& config){
    
    // Get list of G4 predefined elements we need to load.
    vector<string> elementsToLoad;
    if ( config.hasName("g4.elements") ){
      config.getVectorString("g4.elements",elementsToLoad);
    }
    
    // G4 manager for elements and materials.
    G4NistManager* nistMan = G4NistManager::Instance();
    
    // Load the elements.
    for ( vector<string>::const_iterator i=elementsToLoad.begin();
          i!=elementsToLoad.end(); ++i ){
      getElementOrThrow(*i);
    }

    // Get list of G4 predefined materials we need to load.
    vector<string> materialsToLoad;
    if ( config.hasName("g4.materials") ){
      config.getVectorString("g4.materials",materialsToLoad);
    }

    // Load the materials.
    for ( vector<string>::const_iterator i=materialsToLoad.begin();
          i!=materialsToLoad.end(); ++i ){
      findMaterialOrThrow(*i);
    }
    
  }

  // Build the requested Mu2e specific materials.
  // Notes:
  // 1) In the methods: G4Material::AddElement( G4Element* elem, ... )
  //    elem is actually an intrusive reference counted pointer to an element.
  //    The element keeps track of how many times it is used - so the first
  //    arugment cannot be const*.
  void ConstructMaterials::constructMu2eMaterials(SimpleConfig const& config){

    G4NistManager* nistMan = G4NistManager::Instance();

    // List of requested Mu2e specific materials from the config file.
    vector<string> materialsToLoad;
    config.getVectorString("mu2e.materials",materialsToLoad);

    CheckedG4String mat = isNeeded(materialsToLoad, "ShieldingConcrete");
    if ( mat.doit ) {
      //
      // Concrete is 2.00 for fraction, but shielding concrete has reinforcing Iron bars:
      // www-esh.fnal.gov/TM1834_PDF_FILES/TM_1834_Revision_9.pdf
      //
      G4Material* ShieldingConcrete = 
        new G4Material( mat.name, 2.5*g/cm3, 6);
      G4Element* eO  = getElementOrThrow("O");
      G4Element* eSi = getElementOrThrow("Si");
      G4Element* eCa = getElementOrThrow("Ca");
      G4Element* eNa = getElementOrThrow("Na");
      G4Element* eFe = getElementOrThrow("Fe");
      G4Element* eAl = getElementOrThrow("Al");
      ShieldingConcrete->AddElement( eO ,  0.5200);
      ShieldingConcrete->AddElement( eSi,  0.3250);
      ShieldingConcrete->AddElement( eCa,  0.0600);
      ShieldingConcrete->AddElement( eNa,  0.0150);
      ShieldingConcrete->AddElement( eFe,  0.0400); 
      ShieldingConcrete->AddElement( eAl,  0.0400);
    }
    
    mat = isNeeded(materialsToLoad, "Polyethylene");
    if ( mat.doit ){
      G4Material* Polyethylene = 
        new G4Material( mat.name, 0.956*g/cm3, 2);
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eH  = getElementOrThrow("H");
      Polyethylene->AddElement( eC, 1);
      Polyethylene->AddElement( eH, 2);
    }

    mat = isNeeded(materialsToLoad, "IsoButane");
    if ( mat.doit ){
      G4Material* IsoButane = 
        new G4Material( mat.name, 0.00265*g/cm3, 2);
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eH  = getElementOrThrow("H");
      IsoButane->AddElement( eC, 4);
      IsoButane->AddElement( eH, 10);
    }

    mat = isNeeded(materialsToLoad, "StrawGas");
    if ( mat.doit ) {
      G4Material* StrawGas = 
        new G4Material(mat.name, 0.0028561*g/cm3, 3);
      G4Element* eAr = getElementOrThrow("Ar");
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eF  = getElementOrThrow("F");
      StrawGas->AddElement( eAr, 1);
      StrawGas->AddElement( eC,  1);
      StrawGas->AddElement( eF,  4);
    }
    
    mat = isNeeded(materialsToLoad, "Kapton");
    if ( mat.doit ){
      //
      // Kapton: from NIST: physics.nist.gov/cgi-bin/Star/compos.pl?matno=179
      //
      G4Material* Kapton = 
        new G4Material(mat.name, 1.42*g/cm3, 4);
      G4Element* eH  = getElementOrThrow("H");
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eN  = getElementOrThrow("N");
      G4Element* eO  = getElementOrThrow("O");
      Kapton->AddElement( eH, 0.026362);
      Kapton->AddElement( eC, 0.691133);
      Kapton->AddElement( eN, 0.073270);
      Kapton->AddElement( eO, 0.209235);
    }

    mat = isNeeded(materialsToLoad, "Scintillator");
    if ( mat.doit ){
      //  
      // Scintillator.  
      // We probably want several flavors of scintillator so that we can change the
      // detector by just changing a name in the config file.
      //
      G4Material* Sci = 
        new G4Material( mat.name, 1.032*g/cm3, 2);
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eH  = getElementOrThrow("H");
      Sci->AddElement( eC, 9);
      Sci->AddElement( eH, 10);
    }

    mat = isNeeded(materialsToLoad, "WAGVacuum");
    if ( mat.doit ){
      //
      // WAG at properties of Vacuum.  
      // May need several different levels of vacuum in different parts
      // of Mu2e.
      //
      G4double density     = universe_mean_density;
      G4double pressure    = 3.e-18*pascal;
      G4double temperature = 2.73*kelvin;
      G4Material* Vacuum = 
        new G4Material( mat.name, 1., 1.01 *g/mole,
                        density, kStateGas, temperature, pressure);
    }

    mat = isNeeded(materialsToLoad, "MBOverburden");
    if ( mat.doit ){
      //
      // MiniBoone model of the earthen overburden.  See Mu2e-doc
      //  
      G4Material* mbOverburden = new G4Material( mat.name, 2.15*g/cm3, 3);
      G4Element* eO  = getElementOrThrow("O");
      G4Element* eSi = getElementOrThrow("Si");
      G4Element* eAl = getElementOrThrow("Al");
      mbOverburden->AddElement( eO,  65);
      mbOverburden->AddElement( eSi, 20);
      mbOverburden->AddElement( eAl, 15);
    }

    mat = isNeeded(materialsToLoad, "ITGasMix");
    if ( mat.doit ){
      //He/C4H10-gas-mixture

      G4double a, z;
      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.0001786*g/cm3;
      G4double densityIsoB = 0.00267  *g/cm3;
      G4double fractionHe  = 90.0*perCent;

      density = fractionHe*densityHe + (1.0-fractionHe)*densityIsoB;

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature= 293.15*kelvin, pressure= 1*atmosphere);

      //      G4Element* He = new G4Element("He"       , "He", z=2.0, a= 4.002602 *g/mole);
      //      G4Element* C  = new G4Element("Carbonium", "C" , z=6.0, a= 12.0107  *g/mole);
      //      G4Element* H  = new G4Element("Hydrogen" , "H" , z=1.0, a= 1.00794  *g/mole);
      G4Element* He = getElementOrThrow("He");
      G4Element* C  = getElementOrThrow("C");
      G4Element* H  = getElementOrThrow("H");
      GasMix->AddElement(He, 0.9   );
      GasMix->AddElement(H , 0.0173);
      GasMix->AddElement(C , 0.0827);
    }

    mat = isNeeded(materialsToLoad, "CarbonFiber");
    if ( mat.doit ){
      G4double density;
      G4int nel;
      G4Material* CarbonFiber =
        new G4Material(mat.name, density = 2.265*g/cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CarbonFiber->AddElement(C, 100.0*perCent );
    }

    mat = isNeeded(materialsToLoad, "PolypropyleneFoam");
    if ( mat.doit ){
      //Polypropylene (CH3)
      G4double density;
      G4int nel;
      G4Material *Polypropylene = new G4Material(mat.name, density = 0.04*g/cm3, nel=2);
      G4Element* H  = getElementOrThrow("H");
      G4Element* C  = getElementOrThrow("C");
      Polypropylene->AddElement(H, 3 );
      Polypropylene->AddElement(C, 1 );
    }


    // Completed constructing Mu2e specific materials.

    // Check that all requested materials are present. 
    for ( vector<string>::const_iterator i=materialsToLoad.begin();
          i!=materialsToLoad.end(); ++i ){
      if ( *i != DoAllValue ){
        findMaterialOrThrow(*i);
      }
    }
  }

  // Decide if we need to build this material.
  // If additional tests are required, call them from within this method.
  CheckedG4String ConstructMaterials::isNeeded( vector<string> const& V, string const& s){

    // Default return value is not to build it.
    CheckedG4String val(false,s);

    // Throw if the material already exists.
    uniqueMaterialOrThrow(val.name);

    // Is this material requested explicitly?
    val.doit = isRequested( V, s );

    // Is this material requested implicitly?
    if ( !val.doit ) val.doit = isRequested(V, DoAllValue);

    return val;

  }

  // Return true if the requested string is present in the container. 
  // The match must be exact.
  bool ConstructMaterials::isRequested( vector<string> const& V, string const& s){
    for ( vector<string>::const_iterator i=V.begin(), e=V.end();
          i != e; ++i ){
      if ( *i == s ) return true;
    }
    return false;
  }

  // Check to see if the named material already exists.
  void ConstructMaterials::uniqueMaterialOrThrow( G4String const& name){
    if ( G4Material::GetMaterial(name,false) != 0 ){
      throw cms::Exception("GEOM")
        << "mu2e::ConstructMaterials::constructMu2eMaterials(): "
        << "The requested material is already defined: "
        << name
        << "\n";
    }
  }

  // Wrapper around FindOrBuildElement.
  // This is protection against spelling mistakes in the name.
  G4Element* ConstructMaterials::getElementOrThrow( G4String const& name){

    // G4 manager for elements and materials.
    G4NistManager* nistMan = G4NistManager::Instance();
    
    G4Element* answer = nistMan->FindOrBuildElement(name,true);
      
    // Throw if we could not find a requested element.
    if ( !answer ){
      throw cms::Exception("GEOM")
        << "mu2e::ConstructMaterials::constructMaterials(): "
        << "Could not load predefined G4 element named: "
        << name
        << "\n";
    }

    return answer;
  }

  
} // end namespace mu2e
