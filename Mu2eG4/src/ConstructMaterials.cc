//
// Construct materials requested by the run-time configuration system.
//
// $Id: ConstructMaterials.cc,v 1.23 2011/10/11 20:29:51 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/10/11 20:29:51 $
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
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

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
    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const& config = geom->config();

    // Construct the requested materials.
    constructMu2eMaterials( config );

    // Print element table, if requested.
    if ( config.getBool("g4.printElements",false) ){
      mf::LogInfo  log("GEOM");
      log << *G4Element::GetElementTable();
    }

    // Print material table, if requested.
    if ( config.getBool("g4.printMaterials",false) ){
      mf::LogInfo  log("GEOM");
      log << *G4Material::GetMaterialTable();
    }

  }

  // Build the requested Mu2e specific materials.
  // Notes:
  // 1) In the methods: G4Material::AddElement( G4Element* elem, ... )
  //    Each element object keeps track of how many time it is used in a material.
  //    Therefore the first argument cannot be a const pointer.
  void ConstructMaterials::constructMu2eMaterials(SimpleConfig const& config){

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

    // polyethylene data below as in John Caunt Scientific, also see shieldwerx
    mat = isNeeded(materialsToLoad, "Polyethylene092");
    if ( mat.doit ){
      G4Material* Polyethylene092 =
        new G4Material( mat.name, 0.92*g/cm3, 2);
      G4Material* mC  = findMaterialOrThrow("G4_C");
      G4Material* mH  = findMaterialOrThrow("G4_H");
      Polyethylene092->AddMaterial( mH, 0.1428);
      Polyethylene092->AddMaterial( mC, 0.8572);
    }

    mat = isNeeded(materialsToLoad, "Polyethylene0956");
    if ( mat.doit ){
      G4Material* Polyethylene0956 =
        new G4Material( mat.name, 0.956*g/cm3, 2);
      G4Material* mC  = findMaterialOrThrow("G4_C");
      G4Material* mH  = findMaterialOrThrow("G4_H");
      Polyethylene0956->AddMaterial( mH, 0.143711);
      Polyethylene0956->AddMaterial( mC, 0.856289);
    }

    //   note that G4 has:
    //   AddMaterial("G4_POLYETHYLENE", 0.94, 0, 57.4, 2);
    //   AddElementByWeightFraction( 1, 0.143711);
    //   AddElementByWeightFraction( 6, 0.856289);
    //   AddChemicalFormula("G4_POLYETHYLENE","(C_2H_4)_N-Polyethylene");

    // borated polyethylene data as in John Caunt Scientific
    mat = isNeeded(materialsToLoad, "Polyethylene092B050d095");
    if ( mat.doit ){
      G4Material* Polyethylene092B050d095 =
        new G4Material( mat.name, 0.95*g/cm3, 2);
      // we will use the Polyethylene092 and add B as a material
      G4Material* Polyethylene092 = findMaterialOrThrow("Polyethylene092");
      G4Material* mB              = findMaterialOrThrow("G4_B");
      double BPercentage = 5.0;
      Polyethylene092B050d095->AddMaterial( Polyethylene092, (100.-BPercentage)*perCent);
      Polyethylene092B050d095->AddMaterial( mB, BPercentage*perCent);
    }

    mat = isNeeded(materialsToLoad, "Polyethylene092B300d119");
    if ( mat.doit ){
      G4Material* Polyethylene092B300d119 =
        new G4Material( mat.name, 1.19*g/cm3, 2);
      // we will use the Polyethylene092 and add B as a material
      G4Material* Polyethylene092 = findMaterialOrThrow("Polyethylene092");
      G4Material* mB              = findMaterialOrThrow("G4_B");
      double BPercentage = 30.0;
      Polyethylene092B300d119->AddMaterial( Polyethylene092, (100.-BPercentage)*perCent);
      Polyethylene092B300d119->AddMaterial( mB, BPercentage*perCent);
    }

    mat = isNeeded(materialsToLoad, "Polyethylene092Li075d106");
    if ( mat.doit ){
      G4Material* Polyethylene092Li075d106 =
        new G4Material( mat.name, 1.06*g/cm3, 2);
      // we will use the Polyethylene092 and add Li as a material
      G4Material* Polyethylene092 = findMaterialOrThrow("Polyethylene092");
      G4Material* mLi              = findMaterialOrThrow("G4_Li");
      double LiPercentage = 7.5;
      Polyethylene092Li075d106->AddMaterial( Polyethylene092, (100.-LiPercentage)*perCent);
      Polyethylene092Li075d106->AddMaterial( mLi, LiPercentage*perCent);
    }

    // Stainless Steel (Medical Physics, Vol 25, No 10, Oct 1998) based on brachytherapy example
    mat = isNeeded(materialsToLoad, "StainlessSteel");
    if ( mat.doit ){
      G4Material* StainlessSteel =
        new G4Material( mat.name, 8.02*g/cm3, 5);
      StainlessSteel->AddMaterial(findMaterialOrThrow("G4_Mn"), 0.02);
      StainlessSteel->AddMaterial(findMaterialOrThrow("G4_Si"), 0.01);
      StainlessSteel->AddMaterial(findMaterialOrThrow("G4_Cr"), 0.19);
      StainlessSteel->AddMaterial(findMaterialOrThrow("G4_Ni"), 0.10);
      StainlessSteel->AddMaterial(findMaterialOrThrow("G4_Fe"), 0.68);
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
        new G4Material(mat.name, 0.0028561*g/cm3, 3); // it is OK not to use kStateGas
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

      // G4 takes ownership of this object and manages its lifetime.
      new G4Material( mat.name, 1., 1.01 *g/mole,
                      density, kStateGas, temperature, pressure);
    }


    //new added. Check if it is correct
    mat = isNeeded(materialsToLoad, "DSVacuum");
    if ( mat.doit ){

      G4Element* N = getElementOrThrow("N");

      G4double refDensity  = 1.1652e-3*g/cm3; // Nitrogen pressure at 0 celsius 1 bar
      G4double refPress    = 100e3*pascal; // 1 bar
      G4double refTemp     = 273.15*kelvin; // 0 celsius
      G4double pressure    = 133.322e-4*pascal; // 10e-4 Torr
      G4double temperature = 300.00*kelvin; // Temperature of the DS

      G4double density = refDensity*pressure*refTemp/(refPress*temperature);
      
      G4Material* DSVacuum =
	new G4Material(mat.name, density, 1, kStateGas, temperature, pressure);

      G4int nAtoms;
      DSVacuum->AddElement(N, nAtoms=2);
    }

    mat = isNeeded(materialsToLoad, "MBOverburden");
    if ( mat.doit ){
      //
      // MiniBoone model of the earthen overburden.  See Mu2e-doc-570.
      //
      G4Material* mbOverburden = new G4Material( mat.name, 2.15*g/cm3, 3);
      G4Element* eO  = getElementOrThrow("O");
      G4Element* eSi = getElementOrThrow("Si");
      G4Element* eAl = getElementOrThrow("Al");
      mbOverburden->AddElement( eO,  65);
      mbOverburden->AddElement( eSi, 20);
      mbOverburden->AddElement( eAl, 15);
    }

    mat = isNeeded(materialsToLoad, "ITGasHe_90Isob_10");
    if ( mat.doit ){

      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.000166 *g/cm3;
      G4double densityIsoB = 0.00249  *g/cm3;
      G4double fractionHe  = 90.0*perCent;

      density = fractionHe*densityHe + (1.0-fractionHe)*densityIsoB;

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature= 293.15*kelvin, pressure= 1*atmosphere);

      G4Element* He = getElementOrThrow("He");
      G4Element* C  = getElementOrThrow("C");
      G4Element* H  = getElementOrThrow("H");

      G4double atomicWeight_He =  4.002602 *g/mole;
      G4double atomicWeight_C  = 12.0107   *g/mole;
      G4double atomicWeight_H  =  1.00794  *g/mole;
      G4double pwHe = fractionHe*atomicWeight_He;
      G4double pwC  = (1.0-fractionHe) *  4.0*atomicWeight_C;
      G4double pwH  = (1.0-fractionHe) * 10.0*atomicWeight_H;
      G4double atomicWeightMix = pwHe + pwC + pwH ;

      pwHe/=atomicWeightMix;
      pwH/=atomicWeightMix;
      GasMix->AddElement(He, pwHe );
      GasMix->AddElement(H , pwH  );
      GasMix->AddElement(C , 1.0-pwHe-pwH  );
    }

    mat = isNeeded(materialsToLoad, "ITGasHe_90CF4_10");
    if ( mat.doit ){

      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.000166 *g/cm3;
      G4double densityCF4  = 0.003780 *g/cm3;
      G4double fractionHe  = 90.0*perCent;

      density = fractionHe*densityHe + (1.0-fractionHe)*densityCF4;

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature= 293.15*kelvin, pressure= 1*atmosphere);

      G4Element* He = getElementOrThrow("He");
      G4Element* C  = getElementOrThrow("C");
      G4Element* F  = getElementOrThrow("F");

      G4double atomicWeight_He =  4.002602  *g/mole;
      G4double atomicWeight_C  = 12.0107    *g/mole;
      G4double atomicWeight_F  = 18.9984032 *g/mole;
      G4double pwHe = fractionHe*atomicWeight_He;
      G4double pwC  = (1.0-fractionHe) *  1.0*atomicWeight_C;
      G4double pwF  = (1.0-fractionHe) *  4.0*atomicWeight_F;
      G4double atomicWeightMix = pwHe + pwC + pwF ;

      pwHe/=atomicWeightMix;
      pwF/=atomicWeightMix;
      GasMix->AddElement(He, pwHe );
      GasMix->AddElement(F , pwF  );
      GasMix->AddElement(C , 1.0-pwHe-pwF  );
    }

    mat = isNeeded(materialsToLoad, "ITGasMix");
    if ( mat.doit ){
      //He/C4H10-gas-mixture

      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.0001786*g/cm3;
      G4double densityIsoB = 0.00267  *g/cm3;
      G4double fractionHe  = 90.0*perCent;

      density = fractionHe*densityHe + (1.0-fractionHe)*densityIsoB;

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature= 293.15*kelvin, pressure= 1*atmosphere);

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
        new G4Material(mat.name, density = 1.384*g/cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CarbonFiber->AddElement(C, 100.0*perCent );
    }
    
    mat = isNeeded(materialsToLoad, "Lyso_01");  /// Alessandra
    if ( mat.doit ){
      G4double density;
      G4int nel;
      G4Material* Lyso_00 =
        new G4Material(mat.name, density = 7.4*g/cm3, nel=4);
      G4Element* Lu  = getElementOrThrow("Lu");
      G4Element* Si  = getElementOrThrow("Si");
      G4Element* O  = getElementOrThrow("O");
      G4Element* Y  = getElementOrThrow("Y");
      G4Element* Ce  = getElementOrThrow("Ce");
      Lyso_00->AddElement( Lu, 71.0*perCent );
      Lyso_00->AddElement( Si, 7.0*perCent );      
      Lyso_00->AddElement( O, 18.0*perCent );      
      Lyso_00->AddElement( Y, 4.0*perCent );
      G4Material* Lyso_01 =
        new G4Material(mat.name, density = 7.4*g/cm3, nel=2);
      Lyso_01->AddMaterial( Lyso_00, 99.85*perCent ); 
      Lyso_01->AddElement( Ce, 0.15*perCent );    
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

    mat = isNeeded(materialsToLoad, "CFoam_480");
    if ( mat.doit ){
      G4double density;
      G4int nel;
      G4Material *CFoam = new G4Material(mat.name, density = 0.480*g/cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CFoam->AddElement(C, 100.0*perCent );
    }

    mat = isNeeded(materialsToLoad, "CFoam_100");
    if ( mat.doit ){
      G4double density;
      G4int nel;
      G4Material *CFoam = new G4Material(mat.name, density = 0.100*g/cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CFoam->AddElement(C, 100.0*perCent );
    }

    mat = isNeeded(materialsToLoad, "CFoam_080");
    if ( mat.doit ){
      G4double density;
      G4int nel;
      G4Material *CFoam = new G4Material(mat.name, density = 0.080*g/cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CFoam->AddElement(C, 100.0*perCent );
    }

    mat = isNeeded(materialsToLoad, "CFoam");
    if ( mat.doit ){
      G4double density;
      G4int nel;
      G4Material *CFoam = new G4Material(mat.name, density = 0.030*g/cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CFoam->AddElement(C, 100.0*perCent );
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
      throw cet::exception("GEOM")
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
      throw cet::exception("GEOM")
        << "mu2e::ConstructMaterials::constructMaterials(): "
        << "Could not load predefined G4 element named: "
        << name
        << "\n";
    }

    return answer;
  }

} // end namespace mu2e
