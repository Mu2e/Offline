//
// Construct materials requested by the run-time configuration system.
//

// Original author Rob Kutschke
//
// Notes:
// 1) This code new's many objects of type G4Material.  The lifeime
//    of these objects is controlled within G4.  We must never
//    delete them.
// 2) (DNB Louisville).  Because of the large number of materials, we
//    have split the function for building Mu2e material into two
//    functions, which I have unimaginatively named
//    constructMu2eMaterials and constructMu2eMaterials2.  New materials
//    should be added to the latter.
//

// C++ includes
#include <iostream>
#include <iomanip>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Offline/Mu2eG4/inc/ConstructMaterials.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/setBirksConstant.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

// CLHEP includes
#include "CLHEP/Units/PhysicalConstants.h"

// G4 includes
#include "Geant4/G4GeometryManager.hh"
//#include "Geant4/G4PhysicalVolumeStore.hh"
//#include "Geant4/G4LogicalVolumeStore.hh"

#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialTable.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/globals.hh"
#include "Geant4/G4NistManager.hh"
#include "Geant4/G4VisAttributes.hh"

using namespace std;

namespace mu2e {

  ConstructMaterials::ConstructMaterials(const Mu2eG4Config::Top& config)
    : config_(config)
    , printElements_(config.debug().printElements())
    , printMaterials_(config.debug().printMaterials())
  {
    art::ServiceHandle<GeometryService> geom;
    mu2eStandardDetector_ = geom->isStandardMu2eDetector();
  }

  ConstructMaterials::~ConstructMaterials(){
  }

  // This is the main function.
  void ConstructMaterials::construct(){

    if (!mu2eStandardDetector_) {
      G4cout << __func__ << " Non standard mu2e configuration, Will NOT construct mu2e materials " << G4endl;
      // one could move this down to constructMu2eMaterials before using GeomHandle's
      return;
    }

    // Construct the requested materials.
    constructMu2eMaterials();
    constructMu2eMaterials2();

    setBirksConstant(config_);

    // Turning on density effect correction in ionization loss calculations
    // for selected materials deemed conductors
    useDensityEffectInIonizationLossCorrectionIfRequested();

    // Print element table, if requested.
    if (printElements_){
      mf::LogInfo  log("GEOM");
      log << *G4Element::GetElementTable();
    }

    // Print material table, if requested.
    if (printMaterials_){
      mf::LogInfo  log("GEOM");
      log << *G4Material::GetMaterialTable();
    }

  }

  // Build the requested Mu2e specific materials.
  // Notes:
  //
  // 1) In the methods: G4Material::AddElement( G4Element* elem, ... )
  //    Each element object keeps track of how many time it is used in a material.
  //    Therefore the first argument cannot be a const pointer.
  //
  // 2) After each mat assignment (mat = ... ), the material construction takes place within
  //    in a code body of local scope--i.e. within { ... }.  This is to allow for
  //    multiple C++ definitions of (e.g.) density, temperature, etc., which are objects
  //    that have G4-specific types (G4double).  Removing the braces would require removing
  //    the G4double typename in front of the variable strings for assignments after the first
  //    one, possibly making things a little confusing for people needing to define materials.
  //  **************** Added Note *********************
  //  Because of warnings from the compiler about too many variables in the
  //  construct function, we have split the constructMu2eMaterials function into
  //  two functions.  New materials should be added to the second function,
  //  named constructMu2eMaterials2.
  //  *************************************************
  void ConstructMaterials::constructMu2eMaterials(){


    CheckedG4String mat = uniqueMaterialOrThrow( "CONCRETE_MARS" );
    {
      G4Material* ConcreteMars = new G4Material(mat.name, 2.35*CLHEP::g/CLHEP::cm3, 9 );
      ConcreteMars->AddElement( getElementOrThrow("H") , 0.006); //Hydrogen
      ConcreteMars->AddElement( getElementOrThrow("C") , 0.030); //Carbon
      ConcreteMars->AddElement( getElementOrThrow("O") , 0.500); //Oxygen
      ConcreteMars->AddElement( getElementOrThrow("Na"), 0.010); //Sodium
      ConcreteMars->AddElement( getElementOrThrow("Al"), 0.030); //Aluminum
      ConcreteMars->AddElement( getElementOrThrow("Si"), 0.200); //Silicon
      ConcreteMars->AddElement( getElementOrThrow("K") , 0.010); //Potassium
      ConcreteMars->AddElement( getElementOrThrow("Ca"), 0.200); //Calcium
      ConcreteMars->AddElement( getElementOrThrow("Fe"), 0.014); //Iron
    }

    mat = uniqueMaterialOrThrow( "CONCRETE_MASONRY" );
    {
      G4Material* ConcreteMasonry = new G4Material(mat.name, 1.16*CLHEP::g/CLHEP::cm3, 9 );
      ConcreteMasonry->AddElement( getElementOrThrow("H") , 0.006); //Hydrogen
      ConcreteMasonry->AddElement( getElementOrThrow("C") , 0.030); //Carbon
      ConcreteMasonry->AddElement( getElementOrThrow("O") , 0.500); //Oxygen
      ConcreteMasonry->AddElement( getElementOrThrow("Na"), 0.010); //Sodium
      ConcreteMasonry->AddElement( getElementOrThrow("Al"), 0.030); //Aluminum
      ConcreteMasonry->AddElement( getElementOrThrow("Si"), 0.200); //Silicon
      ConcreteMasonry->AddElement( getElementOrThrow("K") , 0.010); //Potassium
      ConcreteMasonry->AddElement( getElementOrThrow("Ca"), 0.200); //Calcium
      ConcreteMasonry->AddElement( getElementOrThrow("Fe"), 0.014); //Iron
    }

    mat = uniqueMaterialOrThrow( "CONCRETE_CB4_07P" );
    {
      G4Material* ConcreteCB4_07P = new G4Material(mat.name, 2.35*CLHEP::g/CLHEP::cm3, 2 );
      ConcreteCB4_07P->AddMaterial( findMaterialOrThrow("CONCRETE_MARS") , 0.993); // MARS Concrete
      ConcreteCB4_07P->AddMaterial( findMaterialOrThrow("G4_BORON_CARBIDE") , 0.007); // Boron carbide
    }

    mat = uniqueMaterialOrThrow( "CONCRETE_CB4_1P" );
    {
      G4Material* ConcreteCB4_1P = new G4Material(mat.name, 2.35*CLHEP::g/CLHEP::cm3, 2 );
      ConcreteCB4_1P->AddMaterial( findMaterialOrThrow("CONCRETE_MARS") , 0.99); // MARS Concrete
      ConcreteCB4_1P->AddMaterial( findMaterialOrThrow("G4_BORON_CARBIDE") , 0.01); // Boron carbide
    }

    mat = uniqueMaterialOrThrow( "CONCRETE_CB4_3P" );
    {
      G4Material* ConcreteCB4_3P = new G4Material(mat.name, 2.35*CLHEP::g/CLHEP::cm3, 2 );
      ConcreteCB4_3P->AddMaterial( findMaterialOrThrow("CONCRETE_MARS") , 0.97); // MARS Concrete
      ConcreteCB4_3P->AddMaterial( findMaterialOrThrow("G4_BORON_CARBIDE") , 0.03); // Boron carbide
    }

    mat = uniqueMaterialOrThrow( "CONCRETE_CB4_5P" );
    {
      G4Material* ConcreteCB4_5P = new G4Material(mat.name, 2.35*CLHEP::g/CLHEP::cm3, 2 );
      ConcreteCB4_5P->AddMaterial( findMaterialOrThrow("CONCRETE_MARS") , 0.95); // MARS Concrete
      ConcreteCB4_5P->AddMaterial( findMaterialOrThrow("G4_BORON_CARBIDE") , 0.05); // Boron carbide
    }

    mat = uniqueMaterialOrThrow( "BARITE" );
    {
      G4Material* Barite = new G4Material(mat.name, 3.5*CLHEP::g/CLHEP::cm3, 9 );
      Barite->AddElement( getElementOrThrow("H") , 0.0069); //Hydrogen
      Barite->AddElement( getElementOrThrow("O") , 0.3386); //Oxygen
      Barite->AddElement( getElementOrThrow("Mg"), 0.0011); //Magnesium
      Barite->AddElement( getElementOrThrow("Al"), 0.0039); //Aluminum
      Barite->AddElement( getElementOrThrow("Si"), 0.0100); //Silicon
      Barite->AddElement( getElementOrThrow("S") , 0.1040); //Sulfur
      Barite->AddElement( getElementOrThrow("Ca"), 0.0478); //Calcium
      Barite->AddElement( getElementOrThrow("Fe"), 0.0457); //Iron
      Barite->AddElement( getElementOrThrow("Ba"), 0.4420); //Barium - should be 0.4431, but
                                                            //         G4 complains if fractional
                                                            //         weights don't add up to 1.
    }

    mat = uniqueMaterialOrThrow( "BARITE_CB4_05P" );
    {
      G4Material* BariteCB4_05P = new G4Material(mat.name, 3.5*CLHEP::g/CLHEP::cm3, 2 );
      BariteCB4_05P->AddMaterial( findMaterialOrThrow("BARITE") , 0.995); // Barite concrete
      BariteCB4_05P->AddMaterial( findMaterialOrThrow("G4_BORON_CARBIDE") , 0.005); // Boron carbide
    }

    mat = uniqueMaterialOrThrow( "BARITE_CB4_1P" );
    {
      G4Material* BariteCB4_1P = new G4Material(mat.name, 3.5*CLHEP::g/CLHEP::cm3, 2 );
      BariteCB4_1P->AddMaterial( findMaterialOrThrow("BARITE") , 0.99); // Barite concrete
      BariteCB4_1P->AddMaterial( findMaterialOrThrow("G4_BORON_CARBIDE") , 0.01); // Boron carbide
    }

    mat = uniqueMaterialOrThrow("HeavyConcrete");
    {
      G4Material* HeavyConcrete = new G4Material(mat.name, 3.295*CLHEP::g/CLHEP::cm3, 17);
      HeavyConcrete->AddElement( getElementOrThrow("H") , 0.01048482); //Hydrogen
      HeavyConcrete->AddElement( getElementOrThrow("B") , 0.00943758); //Boron
      HeavyConcrete->AddElement( getElementOrThrow("C") , 0.0129742);  //Carbon
      HeavyConcrete->AddElement( getElementOrThrow("O") , 0.27953541); //Oxygen
      HeavyConcrete->AddElement( getElementOrThrow("F") , 1.5175E-4);  //Fluorine
      HeavyConcrete->AddElement( getElementOrThrow("Na"), 3.7014E-4);  //Sodium
      HeavyConcrete->AddElement( getElementOrThrow("Mg"), 0.08298213); //Magnesium
      HeavyConcrete->AddElement( getElementOrThrow("Al"), 0.02769028); //Aluminum
      HeavyConcrete->AddElement( getElementOrThrow("Si"), 0.06317253); //Silicon
      HeavyConcrete->AddElement( getElementOrThrow("P") , 0.00176963); //Phosphorus
      HeavyConcrete->AddElement( getElementOrThrow("S") , 5.8275E-4);  //Sulfur
      HeavyConcrete->AddElement( getElementOrThrow("K") , 4.2024E-4);  //Potassium
      HeavyConcrete->AddElement( getElementOrThrow("Ca"), 0.03227609); //Calcium
      HeavyConcrete->AddElement( getElementOrThrow("Ti"), 5.457E-5);   //Titanium
      HeavyConcrete->AddElement( getElementOrThrow("Mn"), 0.00321757); //Manganese
      HeavyConcrete->AddElement( getElementOrThrow("Fe"), 0.47423935); //Iron
      HeavyConcrete->AddElement( getElementOrThrow("Sr"), 6.4097E-4);  //Strontium
    }

    mat = uniqueMaterialOrThrow( "ShieldingConcrete");
    {
      //
      // Concrete is 2.00 for fraction, but shielding concrete has reinforcing Iron bars:
      // www-esh.fnal.gov/TM1834_PDF_FILES/TM_1834_Revision_9.pdf
      //
      G4Material* ShieldingConcrete = new G4Material( mat.name, 2.5*CLHEP::g/CLHEP::cm3, 6);
      ShieldingConcrete->AddElement( getElementOrThrow("O") ,  0.5200);
      ShieldingConcrete->AddElement( getElementOrThrow("Si"),  0.3250);
      ShieldingConcrete->AddElement( getElementOrThrow("Ca"),  0.0600);
      ShieldingConcrete->AddElement( getElementOrThrow("Na"),  0.0150);
      ShieldingConcrete->AddElement( getElementOrThrow("Fe"),  0.0400);
      ShieldingConcrete->AddElement( getElementOrThrow("Al"),  0.0400);
    }

    mat = uniqueMaterialOrThrow( "Polyethylene");
    {
      G4Material* Polyethylene = new G4Material( mat.name, 0.956*CLHEP::g/CLHEP::cm3, 2);
      Polyethylene->AddElement( getElementOrThrow("C"), 1);
      Polyethylene->AddElement( getElementOrThrow("H"), 2);
    }


    mat = uniqueMaterialOrThrow( "Half_Poly" );
    {
      G4Material* HalfPoly = new G4Material( mat.name, 0.465*CLHEP::g/CLHEP::cm3, 1);
      HalfPoly->AddMaterial( findMaterialOrThrow("Polyethylene"), 1. );
    }

    // polyethylene data below as in John Caunt Scientific, also see shieldwerx
    mat = uniqueMaterialOrThrow( "Polyethylene092");
    {
      G4Material* Polyethylene092 = new G4Material( mat.name, 0.92*CLHEP::g/CLHEP::cm3, 2);
      Polyethylene092->AddMaterial( findMaterialOrThrow("G4_H"), 0.1428);
      Polyethylene092->AddMaterial( findMaterialOrThrow("G4_C"), 0.8572);
    }

    mat = uniqueMaterialOrThrow( "Polyethylene0956");
    {
      G4Material* Polyethylene0956 = new G4Material( mat.name, 0.956*CLHEP::g/CLHEP::cm3, 2);
      Polyethylene0956->AddMaterial( findMaterialOrThrow("G4_H"), 0.143711);
      Polyethylene0956->AddMaterial( findMaterialOrThrow("G4_C"), 0.856289);
    }

    mat = uniqueMaterialOrThrow( "IPAPolyethylene");
    {
      G4Material* IPAPolyethylene = new G4Material( mat.name, 0.954*CLHEP::g/CLHEP::cm3, 2);
      IPAPolyethylene->AddMaterial( findMaterialOrThrow("G4_H"), 0.11);
      IPAPolyethylene->AddMaterial( findMaterialOrThrow("G4_C"), 0.89); // Carbon doped Polytehylene, additional carbon 2-5% from MDS (DeWal DW 402B),  density measured by S. Krave 6/22/2021
    }

    mat = uniqueMaterialOrThrow( "Polyethylene096");
    {
      G4Material* Polyethylene096 = new G4Material( mat.name, 0.96*CLHEP::g/CLHEP::cm3, 2);
      Polyethylene096->AddMaterial( findMaterialOrThrow("G4_H"), 0.14);
      Polyethylene096->AddMaterial( findMaterialOrThrow("G4_C"), 0.86);
    }

    mat = uniqueMaterialOrThrow( "Polyethylene094");
    {
      G4Material* Polyethylene094 = new G4Material( mat.name, 0.94*CLHEP::g/CLHEP::cm3, 2);
      Polyethylene094->AddMaterial( findMaterialOrThrow("G4_H"), 0.14);
      Polyethylene094->AddMaterial( findMaterialOrThrow("G4_C"), 0.86);
    }

    mat = uniqueMaterialOrThrow( "Polyethylene0935");
    {
      G4Material* Polyethylene0935 = new G4Material( mat.name, 0.935*CLHEP::g/CLHEP::cm3, 2);
      Polyethylene0935->AddMaterial( findMaterialOrThrow("G4_H"), 0.14);
      Polyethylene0935->AddMaterial( findMaterialOrThrow("G4_C"), 0.86);
    }

    mat = uniqueMaterialOrThrow( "Polyethylene090");
    {
      G4Material* Polyethylene090 = new G4Material( mat.name, 0.90*CLHEP::g/CLHEP::cm3, 2);
      Polyethylene090->AddMaterial( findMaterialOrThrow("G4_H"), 0.14);
      Polyethylene090->AddMaterial( findMaterialOrThrow("G4_C"), 0.86);
    }

    // Not real, very thin Polyethylene
    mat = uniqueMaterialOrThrow( "Polyethylene0010");
    {
      G4Material* Polyethylene0010 = new G4Material( mat.name, 0.010*CLHEP::g/CLHEP::cm3, 2);
      Polyethylene0010->AddMaterial( findMaterialOrThrow("G4_H"), 0.143711);
      Polyethylene0010->AddMaterial( findMaterialOrThrow("G4_C"), 0.856289);
    }

    // Not real, very thin Polyethylene
    mat = uniqueMaterialOrThrow( "Polyethylene0020");
    {
      G4Material* Polyethylene0020 = new G4Material( mat.name, 0.020*CLHEP::g/CLHEP::cm3, 2);
      Polyethylene0020->AddMaterial( findMaterialOrThrow("G4_H"), 0.143711);
      Polyethylene0020->AddMaterial( findMaterialOrThrow("G4_C"), 0.856289);
    }

    //   note that G4 has:
    //   AddMaterial("G4_POLYETHYLENE", 0.94, 0, 57.4, 2);
    //   AddElementByWeightFraction( 1, 0.143711);
    //   AddElementByWeightFraction( 6, 0.856289);
    //   AddChemicalFormula("G4_POLYETHYLENE","(C_2H_4)_N-Polyethylene");

    // borated polyethylene data as in John Caunt Scientific
    mat = uniqueMaterialOrThrow( "Polyethylene092B050d095");
    {
      G4Material* Polyethylene092B050d095 = new G4Material( mat.name, 0.95*CLHEP::g/CLHEP::cm3, 2);
      // we will use the Polyethylene092 and add B as a material
      const double BPercentage = 5.0;
      Polyethylene092B050d095->AddMaterial(findMaterialOrThrow("Polyethylene092"), (100.-BPercentage)*CLHEP::perCent);
      Polyethylene092B050d095->AddMaterial(findMaterialOrThrow("G4_B")           , BPercentage*CLHEP::perCent);
    }

    mat = uniqueMaterialOrThrow( "Polyethylene092B300d119");
    {
      G4Material* Polyethylene092B300d119 = new G4Material( mat.name, 1.19*CLHEP::g/CLHEP::cm3, 2);
      // we will use the Polyethylene092 and add B as a material
      const double BPercentage = 30.0;
      Polyethylene092B300d119->AddMaterial(findMaterialOrThrow("Polyethylene092"), (100.-BPercentage)*CLHEP::perCent);
      Polyethylene092B300d119->AddMaterial(findMaterialOrThrow("G4_B")           , BPercentage*CLHEP::perCent);
    }

    mat = uniqueMaterialOrThrow( "Polyethylene092Li075d106");
    {
      G4Material* Polyethylene092Li075d106 =
        new G4Material( mat.name, 1.06*CLHEP::g/CLHEP::cm3, 2);
      // we will use the Polyethylene092 and add Li as a material
      const double LiPercentage = 7.5;
      Polyethylene092Li075d106->AddMaterial(findMaterialOrThrow("Polyethylene092"), (100.-LiPercentage)*CLHEP::perCent);
      Polyethylene092Li075d106->AddMaterial(findMaterialOrThrow("G4_Li")          , LiPercentage*CLHEP::perCent);
    }

    mat = uniqueMaterialOrThrow( "Polyetheretherketone");
    {
      G4Material* Polyetheretherketone = new G4Material(mat.name, 1.32*CLHEP::g/CLHEP::cm3, 3);
      Polyetheretherketone->AddMaterial(findMaterialOrThrow("G4_C"), 0.513514);
      Polyetheretherketone->AddMaterial(findMaterialOrThrow("G4_H"), 0.405405);
      Polyetheretherketone->AddMaterial(findMaterialOrThrow("G4_O"), 0.081081);
    }

    //  FEP, AKA Teflon FEP, AKA flourinated ethylene propylene.
    //  Used, among other places, in signal cable insulation.
    //  Here, specifically treating as an equal mix of
    //  tetraflouroethylene and hexaflouropropylene (the manufacturers
    //  tend to protect the specifics).  So that's C2F4 + C3F6 = C5F10

    mat = uniqueMaterialOrThrow( "TeflonFEP" );
    {
      G4Material* TeflonFEP = new G4Material(mat.name, 2.14*CLHEP::g/CLHEP::cm3, 2);
      TeflonFEP->AddElement( getElementOrThrow("C"), 5 );
      TeflonFEP->AddElement( getElementOrThrow("F"), 10);
    }


    // Stainless Steel (Medical Physics, Vol 25, No 10, Oct 1998) based on brachytherapy example
    // FIXME is there a better reference?
    mat = uniqueMaterialOrThrow( "StainlessSteel");
    {
      G4Material* StainlessSteel = new G4Material( mat.name, 8.02*CLHEP::g/CLHEP::cm3, 5);
      StainlessSteel->AddMaterial(findMaterialOrThrow("G4_Mn"), 0.02);
      StainlessSteel->AddMaterial(findMaterialOrThrow("G4_Si"), 0.01);
      StainlessSteel->AddMaterial(findMaterialOrThrow("G4_Cr"), 0.19);
      StainlessSteel->AddMaterial(findMaterialOrThrow("G4_Ni"), 0.10);
      StainlessSteel->AddMaterial(findMaterialOrThrow("G4_Fe"), 0.68);
    }


    // Stainless Steel 316 http://en.wikipedia.org/wiki/Marine_grade_stainless
    mat = uniqueMaterialOrThrow( "StainlessSteel316");
    {
      G4Material* StainlessSteel316 = new G4Material( mat.name, 8.00*CLHEP::g/CLHEP::cm3, 10);
      StainlessSteel316->AddMaterial(findMaterialOrThrow("G4_Cr"), 0.17    );
      StainlessSteel316->AddMaterial(findMaterialOrThrow("G4_Ni"), 0.12    );
      StainlessSteel316->AddMaterial(findMaterialOrThrow("G4_C"),  0.0008  );
      StainlessSteel316->AddMaterial(findMaterialOrThrow("G4_Mn"), 0.02    );
      StainlessSteel316->AddMaterial(findMaterialOrThrow("G4_Si"), 0.0075  );
      StainlessSteel316->AddMaterial(findMaterialOrThrow("G4_P"),  0.00045 );
      StainlessSteel316->AddMaterial(findMaterialOrThrow("G4_S"),  0.0003  );
      StainlessSteel316->AddMaterial(findMaterialOrThrow("G4_N"),  0.001   );
      StainlessSteel316->AddMaterial(findMaterialOrThrow("G4_Mo"), 0.025   );
      StainlessSteel316->AddMaterial(findMaterialOrThrow("G4_Fe"), 0.65495 );
    }

    // Stainless Steel 316L http://en.wikipedia.org/wiki/Marine_grade_stainless
    mat = uniqueMaterialOrThrow( "StainlessSteel316L");
    {
      G4Material* StainlessSteel316L = new G4Material( mat.name, 8.00*CLHEP::g/CLHEP::cm3, 10);
      StainlessSteel316L->AddMaterial(findMaterialOrThrow("G4_Cr"), 0.17    );
      StainlessSteel316L->AddMaterial(findMaterialOrThrow("G4_Ni"), 0.12    );
      StainlessSteel316L->AddMaterial(findMaterialOrThrow("G4_C"),  0.0003  );
      StainlessSteel316L->AddMaterial(findMaterialOrThrow("G4_Mn"), 0.02    );
      StainlessSteel316L->AddMaterial(findMaterialOrThrow("G4_Si"), 0.0075   );
      StainlessSteel316L->AddMaterial(findMaterialOrThrow("G4_P"),  0.00045 );
      StainlessSteel316L->AddMaterial(findMaterialOrThrow("G4_S"),  0.0003  );
      StainlessSteel316L->AddMaterial(findMaterialOrThrow("G4_N"),  0.001   );
      StainlessSteel316L->AddMaterial(findMaterialOrThrow("G4_Mo"), 0.025   );
      StainlessSteel316L->AddMaterial(findMaterialOrThrow("G4_Fe"), 0.65545 );
    }

    // A standard carbon-steel used for racks
    mat = uniqueMaterialOrThrow( "RackSteel" );
    {
      G4Material* RackSteel = new G4Material( mat.name, 8.05*CLHEP::g/CLHEP::cm3, 2);
      RackSteel->AddMaterial(findMaterialOrThrow("G4_Al"),0.985);
      RackSteel->AddMaterial(findMaterialOrThrow("G4_C"), 0.015);
    }

    // Inconel 718 alloy.  Taken from http://www.matweb.com/
    // Only building from elements that comprise more than 0.1% of it composition
    mat = uniqueMaterialOrThrow( "Inconel718" );
    {
      G4Material* Inconel718 = new G4Material( mat.name, 8.19*CLHEP::g/CLHEP::cm3, 11);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Al"),0.005);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Cr"),0.190);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Co"),0.010);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Cu"),0.003);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Fe"),0.170);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Mn"),0.003);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Mo"),0.030);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Ni"),0.527);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Nb"),0.052);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Si"),0.003);
      Inconel718->AddMaterial(findMaterialOrThrow("G4_Ti"),0.007);
    }

    // Bronze used in the HRS.  Formally, Bronze C63200.
    mat = uniqueMaterialOrThrow( "HRSBronze" );
    {
      G4Material* HRSBronze = new G4Material( mat.name, 7.64*CLHEP::g/CLHEP::cm3, 4);
      HRSBronze->AddMaterial(findMaterialOrThrow("G4_Cu"),0.820);
      HRSBronze->AddMaterial(findMaterialOrThrow("G4_Al"),0.090);
      HRSBronze->AddMaterial(findMaterialOrThrow("G4_Fe"),0.040);
      HRSBronze->AddMaterial(findMaterialOrThrow("G4_Ni"),0.050);
    }

    // Bronze C93800  from www.matweb.com
    mat = uniqueMaterialOrThrow( "BronzeC943" );
    {
      G4Material* BronzeC943 = new G4Material( mat.name, 9.29*CLHEP::g/CLHEP::cm3, 6);
      BronzeC943->AddMaterial(findMaterialOrThrow("G4_Cu"),0.700);
      BronzeC943->AddMaterial(findMaterialOrThrow("G4_Pb"),0.240);
      BronzeC943->AddMaterial(findMaterialOrThrow("G4_Sn"),0.0475);
      BronzeC943->AddMaterial(findMaterialOrThrow("G4_Sb"),0.005);
      BronzeC943->AddMaterial(findMaterialOrThrow("G4_Ni"),0.005);
      BronzeC943->AddMaterial(findMaterialOrThrow("G4_Zn"),0.0025);
    }

    // Bronze C60800 somewat based on www.matweb.com
    mat = uniqueMaterialOrThrow( "BronzeC608" );
    {
      G4Material* BronzeC608 = new G4Material( mat.name, 8.17*CLHEP::g/CLHEP::cm3, 5);
      BronzeC608->AddMaterial(findMaterialOrThrow("G4_Cu"),0.9310);
      BronzeC608->AddMaterial(findMaterialOrThrow("G4_Al"),0.0572);
      BronzeC608->AddMaterial(findMaterialOrThrow("G4_Fe"),0.0050);
      BronzeC608->AddMaterial(findMaterialOrThrow("G4_Pb"),0.0050);
      BronzeC608->AddMaterial(findMaterialOrThrow("G4_As"),0.0018);
    }

    // Bronze C94500  from https://alloys.copper.org/alloy/C94500
    mat = uniqueMaterialOrThrow( "BronzeC945" );
    {
      G4Material* BronzeC945 = new G4Material( mat.name, 9.40*CLHEP::g/CLHEP::cm3, 3);
      BronzeC945->AddMaterial(findMaterialOrThrow("G4_Cu"),0.78);
      BronzeC945->AddMaterial(findMaterialOrThrow("G4_Pb"),0.16);
      BronzeC945->AddMaterial(findMaterialOrThrow("G4_Sn"),0.06);
    }

    // Bronze C93800  from https://alloys.copper.org/alloy/C93800
    mat = uniqueMaterialOrThrow( "BronzeC938" );
    {
      G4Material* BronzeC938 = new G4Material( mat.name, 9.40*CLHEP::g/CLHEP::cm3, 6);
      BronzeC938->AddMaterial(findMaterialOrThrow("G4_Cu"),0.76);
      BronzeC938->AddMaterial(findMaterialOrThrow("G4_Pb"),0.145);
      BronzeC938->AddMaterial(findMaterialOrThrow("G4_Sn"),0.069);
      BronzeC938->AddMaterial(findMaterialOrThrow("G4_Ni"),0.01);
      BronzeC938->AddMaterial(findMaterialOrThrow("G4_Zn"),0.008);
      BronzeC938->AddMaterial(findMaterialOrThrow("G4_Sb"),0.008);
    }


    // C360 brass
    mat = uniqueMaterialOrThrow( "BrassC360" );
    {
      G4Material* BrassC360 = new G4Material(mat.name, 8.50*CLHEP::g/CLHEP::cm3,3);
      BrassC360->AddMaterial(findMaterialOrThrow("G4_Cu"),0.615);
      BrassC360->AddMaterial(findMaterialOrThrow("G4_Zn"),0.354);
      BrassC360->AddMaterial(findMaterialOrThrow("G4_Pb"),0.031);
    }

    // C64200 from https://alloys.copper.org/alloy/C64200
    mat = uniqueMaterialOrThrow( "BronzeC642" );
    {
      G4Material* BronzeC642 = new G4Material( mat.name, 7.70*CLHEP::g/CLHEP::cm3, 3);
      BronzeC642->AddMaterial(findMaterialOrThrow("G4_Cu"),0.922);
      BronzeC642->AddMaterial(findMaterialOrThrow("G4_Al"),0.063);
      BronzeC642->AddMaterial(findMaterialOrThrow("G4_Si"),0.015);
    }

    // A mix made to represent the MBS spherical support
    mat = uniqueMaterialOrThrow( "MBSSupportMix" );
    {
      G4Material* MBSSupportMix = new G4Material(mat.name, 1.26*CLHEP::g/CLHEP::cm3, 2);
      MBSSupportMix->AddMaterial(findMaterialOrThrow("StainlessSteel"),0.8);
      MBSSupportMix->AddMaterial(findMaterialOrThrow("HRSBronze"),0.2);
    }

    // Construction Aluminum
    //http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA5083O
    //http://ppd-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=1112;filename=MD-ENG-109.pdf;version=1
    mat = uniqueMaterialOrThrow( "A95083");
    {
      G4Material* A95083 = new G4Material( mat.name, 2.66*CLHEP::g/CLHEP::cm3, 9);
      A95083->AddMaterial(findMaterialOrThrow("G4_Al"), 0.9400);
      A95083->AddMaterial(findMaterialOrThrow("G4_Mg"), 0.0450);
      A95083->AddMaterial(findMaterialOrThrow("G4_Mn"), 0.0070);
      A95083->AddMaterial(findMaterialOrThrow("G4_Fe"), 0.0020);
      A95083->AddMaterial(findMaterialOrThrow("G4_Si"), 0.0020);
      A95083->AddMaterial(findMaterialOrThrow("G4_Cr"), 0.0015);
      A95083->AddMaterial(findMaterialOrThrow("G4_Zn"), 0.0013);
      A95083->AddMaterial(findMaterialOrThrow("G4_Ti"), 0.0007);
      A95083->AddMaterial(findMaterialOrThrow("G4_Cu"), 0.0005);
    }

    // 1100 Aluminum
    // https://en.wikipedia.org/wiki/1100_aluminium_alloy
    // http://www.matweb.com/search/DataSheet.aspx?MatGUID=2ca5a0592e4147848bdbd40d1ff1a056&ckck=1
    mat = uniqueMaterialOrThrow( "A1100");
    {
      G4Material* A1100 = new G4Material( mat.name, 2.71*CLHEP::g/CLHEP::cm3, 6);
      A1100->AddMaterial(findMaterialOrThrow("G4_Al"), 0.99275);
      A1100->AddMaterial(findMaterialOrThrow("G4_Fe"), 0.00275);
      A1100->AddMaterial(findMaterialOrThrow("G4_Si"), 0.00275);
      A1100->AddMaterial(findMaterialOrThrow("G4_Cu"), 0.00100);
      A1100->AddMaterial(findMaterialOrThrow("G4_Mn"), 0.00025);
      A1100->AddMaterial(findMaterialOrThrow("G4_Zn"), 0.0005);
    }

    // 6105 Aluminum, "small" extrusion density
    // 8020 T-Slotted Profile 3030
    // https://en.wikipedia.org/wiki/6105_aluminium_alloy
    // https://www.matweb.com/search/datasheet.aspx?matguid=9d1c81ac4e2b4e5590e5781f842b4446&ckck=1
    // This has reduced density to make up for the fact that the Al extrusion is not a solid block
    mat = uniqueMaterialOrThrow("A6105SmallExtrusion");
    {
      G4Material* A6105SmallExtrusion = new G4Material(mat.name, 2.69*0.51*CLHEP::g/CLHEP::cm3, 9);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Al"), 0.973);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Cr"), 0.001);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Cu"), 0.001);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Fe"), 0.0035);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Mg"), 0.008);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Mn"), 0.0015);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Si"), 0.01);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Ti"), 0.001);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Zn"), 0.001);
    }

    // 6105 Aluminum, "large" extrusion density
    // 8020 T-Slotted Profile 3030
    // https://en.wikipedia.org/wiki/6105_aluminium_alloy
    // https://www.matweb.com/search/datasheet.aspx?matguid=9d1c81ac4e2b4e5590e5781f842b4446&ckck=1
    // This has reduced density to make up for the fact that the Al extrusion is not a solid block
    mat = uniqueMaterialOrThrow("A6105LargeExtrusion");
    {
      G4Material* A6105SmallExtrusion = new G4Material(mat.name, 2.69*0.36*CLHEP::g/CLHEP::cm3, 9);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Al"), 0.973);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Cr"), 0.001);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Cu"), 0.001);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Fe"), 0.0035);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Mg"), 0.008);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Mn"), 0.0015);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Si"), 0.01);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Ti"), 0.001);
      A6105SmallExtrusion->AddMaterial(findMaterialOrThrow("G4_Zn"), 0.001);
    }

    // NbTi
    mat = uniqueMaterialOrThrow( "NbTi"); // FIXME verify it
    {
      G4Material* NbTi = new G4Material( mat.name, 6.5*CLHEP::g/CLHEP::cm3, 2);
      NbTi->AddMaterial(findMaterialOrThrow("G4_Nb"), 0.65);
      NbTi->AddMaterial(findMaterialOrThrow("G4_Ti"), 0.35);
    }

    // NbTiCu
    mat = uniqueMaterialOrThrow( "NbTiCu"); // FIXME verify it
    {
      G4Material* NbTiCu = new G4Material( mat.name, 7.69*CLHEP::g/CLHEP::cm3, 2);
      NbTiCu->AddMaterial(findMaterialOrThrow("NbTi"),  0.45);
      NbTiCu->AddMaterial(findMaterialOrThrow("G4_Cu"), 0.55);
    }

    // DS1CoilMix - DS Coil material mix, representing NbTi and Al stiffener.
    // From the dimensions in TDR figure 6.76, the NbTi is 17% by volume
    // while Al makes up 83% by volume.  Percentages by mass are different.
    mat = uniqueMaterialOrThrow( "DS1CoilMix");
    {
      G4Material* DS1CoilMix = new G4Material( mat.name, 3.35*CLHEP::g/CLHEP::cm3,2);
      DS1CoilMix->AddMaterial(findMaterialOrThrow("G4_Al"), 0.67); //proportion
      DS1CoilMix->AddMaterial(findMaterialOrThrow("NbTi"),  0.33); //by mass
    }

    // DS2CoilMix - DS Coil material mix, representing NbTi and Al stiffener.
    // From dimensions in TDR figure 6.77, the NbTi is 8.7% by volume and Al is
    // 91.3% by volume.
    mat = uniqueMaterialOrThrow( "DS2CoilMix");
    {
      G4Material* DS2CoilMix = new G4Material( mat.name, 3.031*CLHEP::g/CLHEP::cm3,2);
      DS2CoilMix->AddMaterial(findMaterialOrThrow("G4_Al"), 0.813);//proportion
      DS2CoilMix->AddMaterial(findMaterialOrThrow("NbTi"),  0.187); //by mass
    }

    // AL999Ni001 by volume ?
    mat = uniqueMaterialOrThrow( "AL999Ni001"); // FIXME verify it
    {
      G4Material* AL999Ni001 = new G4Material( mat.name, 2.706*CLHEP::g/CLHEP::cm3, 2);
      AL999Ni001->AddMaterial(findMaterialOrThrow("G4_Al"), 0.9967);
      AL999Ni001->AddMaterial(findMaterialOrThrow("G4_Ni"), 0.0033);
    }

    // http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/Epotek-301-1.html
    // C_19_H_20_O_4

    mat = uniqueMaterialOrThrow( "C_19_H_20_O_4");
    {
      G4Material* C_19_H_20_O_4 = new G4Material( mat.name, 1.16*CLHEP::g/CLHEP::cm3, 3);
      C_19_H_20_O_4->AddElement( getElementOrThrow("C"), 19);
      C_19_H_20_O_4->AddElement( getElementOrThrow("H"), 20);
      C_19_H_20_O_4->AddElement( getElementOrThrow("O"),  4);
    }

    // C_10_H_18_O_4

    mat = uniqueMaterialOrThrow( "C_10_H_18_O_4");
    {
      G4Material* C_10_H_18_O_4 = new G4Material( mat.name, 1.10*CLHEP::g/CLHEP::cm3, 3);
      C_10_H_18_O_4->AddElement( getElementOrThrow("C"), 10);
      C_10_H_18_O_4->AddElement( getElementOrThrow("H"), 18);
      C_10_H_18_O_4->AddElement( getElementOrThrow("O"),  4);
    }

    // C_9_H_22_N_2

    mat = uniqueMaterialOrThrow( "C_9_H_22_N_2");
    {
      G4Material* C_9_H_22_N_2 = new G4Material( mat.name, 0.865*CLHEP::g/CLHEP::cm3, 3);
      C_9_H_22_N_2->AddElement( getElementOrThrow("C"),  9);
      C_9_H_22_N_2->AddElement( getElementOrThrow("H"), 22);
      C_9_H_22_N_2->AddElement( getElementOrThrow("N"),  2);
    }

    // http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/Epotek-301-1.html
    mat = uniqueMaterialOrThrow( "Epotek301");
    {
      G4Material* Epotek301 = new G4Material( mat.name, 1.19*CLHEP::g/CLHEP::cm3, 3);
      Epotek301->AddMaterial(findMaterialOrThrow("C_19_H_20_O_4"), 0.56);
      Epotek301->AddMaterial(findMaterialOrThrow("C_10_H_18_O_4"), 0.24);
      Epotek301->AddMaterial(findMaterialOrThrow("C_9_H_22_N_2"),  0.20);
    }

    // http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/E_glass.html
    mat = uniqueMaterialOrThrow( "EGlass");
    {
      G4Material* EGlass = new G4Material (mat.name, 2.61*CLHEP::g/CLHEP::cm3, 10);
      EGlass->AddMaterial(findMaterialOrThrow("G4_SILICON_DIOXIDE"), 0.54);
      EGlass->AddMaterial(findMaterialOrThrow("G4_CALCIUM_OXIDE"), 0.19 );
      EGlass->AddMaterial(findMaterialOrThrow("G4_ALUMINUM_OXIDE"), 0.13 );
      EGlass->AddMaterial(findMaterialOrThrow("G4_MAGNESIUM_OXIDE"), 0.025 );
      EGlass->AddMaterial(findMaterialOrThrow("G4_BORON_OXIDE"), 0.075 );
      EGlass->AddMaterial(findMaterialOrThrow("G4_TITANIUM_DIOXIDE"), 0.008 );
      EGlass->AddMaterial(findMaterialOrThrow("G4_SODIUM_MONOXIDE"), 0.01 );
      EGlass->AddMaterial(findMaterialOrThrow("G4_POTASSIUM_OXIDE"), 0.01 );
      EGlass->AddMaterial(findMaterialOrThrow("G4_FERRIC_OXIDE"), 0.005 );
      EGlass->AddMaterial(findMaterialOrThrow("G4_F"), 0.007 );
    }

    // G10 http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/SDD_G10FR4.html
    // http://pdg.lbl.gov/2002/atomicrpp.pdf
    mat = uniqueMaterialOrThrow( "G10");
    {
      G4Material* G10 = new G4Material( mat.name, 1.7*CLHEP::g/CLHEP::cm3, 2);
      G10->AddMaterial(findMaterialOrThrow("G4_SILICON_DIOXIDE"), 0.60);//FIXME do e-glass etc...
      G10->AddMaterial(findMaterialOrThrow("Epotek301"), 0.40);
    }

    mat = uniqueMaterialOrThrow( "G10Lite" );
    {
      // Just a low-density version of G10 to represent space half-filled
      // with G10
      G4Material* G10Lite = new G4Material( mat.name, 0.85*CLHEP::g/CLHEP::cm3,1);
      G10Lite->AddMaterial(findMaterialOrThrow("G10"),1.0);
    }


    // DC 704 diffusion pump oil.  The DC stands for Dow-Corning, though
    // many brands produce the same "704" pump oil.  Chemical
    // name is tetraphenyl tetramethyl trisiloxane, C28H32O2Si3
    // Information from http://www.sigmaaldrich.com/catalog/product/aldrich/445975?lang=en&region=US  on 24 June 2017.

    mat = uniqueMaterialOrThrow( "DC704" );
    {
      G4Material* DC704 = new G4Material( mat.name, 1.07*CLHEP::g/CLHEP::cm3,4);
      DC704->AddElement(getElementOrThrow("C"),28);
      DC704->AddElement(getElementOrThrow("H"),32);
      DC704->AddElement(getElementOrThrow("O"),2);
      DC704->AddElement(getElementOrThrow("Si"),3);
    }

    mat = uniqueMaterialOrThrow( "Electronics" );
    {
      // This material based on measurements at Argonne of some Mu2e
      // electronics by Gary Drake, conveyed by Vitaly Pronskikh
      G4Material* Electronics = new G4Material( mat.name, 0.58*CLHEP::g/CLHEP::cm3, 5);
      Electronics->AddElement( getElementOrThrow("Si"),             0.103);
      Electronics->AddElement( getElementOrThrow("Al"),             0.034);
      Electronics->AddElement( getElementOrThrow("Au"),             0.034);
      Electronics->AddElement( getElementOrThrow("Cu"),             0.691);
      Electronics->AddMaterial( findMaterialOrThrow("Polyethylene"), 0.138);
    }

    // These are needed (i) for MARS to have unique materials, and (ii)
    // in case electronics board materials vary.

    mat = uniqueMaterialOrThrow( "ElectronicsFEB" );
    {
      G4Material* ElectronicsFEB = new G4Material(mat.name, 0.58*CLHEP::g/CLHEP::cm3, 1);
      ElectronicsFEB->AddMaterial( findMaterialOrThrow("Electronics"), 1.0);
    }

    mat = uniqueMaterialOrThrow( "ElectronicsCMB" );
    {
      G4Material* ElectronicsCMB = new G4Material(mat.name, 0.58*CLHEP::g/CLHEP::cm3, 1);
      ElectronicsCMB->AddMaterial( findMaterialOrThrow("Electronics"), 1.0);
    }

    mat = uniqueMaterialOrThrow( "RackElectronics" );
    {
      // This material represents a typical instrumented rack and is currently
      // a placeholder, awaiting measurements.
      // This version is based on PCB Electronics material (above),
      // Aluminum, air, and steel.
      G4Material* RackElectronics = new G4Material( mat.name, 0.19*CLHEP::g/CLHEP::cm3, 7);
      RackElectronics->AddMaterial( findMaterialOrThrow("Electronics"),0.1464);
      RackElectronics->AddMaterial( findMaterialOrThrow( "G4_AIR" ),   0.0045);
      RackElectronics->AddMaterial( findMaterialOrThrow( "G4_Al" ),    0.1170);
      RackElectronics->AddMaterial( findMaterialOrThrow( "RackSteel" ),0.6644);
      RackElectronics->AddMaterial( findMaterialOrThrow( "StainlessSteel" ), 0.0003);
      RackElectronics->AddMaterial( findMaterialOrThrow( "Polyethylene" ), 0.0162);
      RackElectronics->AddMaterial( findMaterialOrThrow( "G4_Cu" ), 0.0512 );
    }


    // Superconducting Cable Insulation
    mat = uniqueMaterialOrThrow( "SCCableInsulation");
    {
      G4Material* SCCableInsulation = new G4Material( mat.name, 1.54*CLHEP::g/CLHEP::cm3, 3);
      SCCableInsulation->AddMaterial(findMaterialOrThrow("G4_KAPTON"), 0.18);
      SCCableInsulation->AddMaterial(findMaterialOrThrow("Epotek301"), 0.16);
      SCCableInsulation->AddMaterial(findMaterialOrThrow("G10"), 0.66);
    }

    //http://panel7.xor.aps.anl.gov/~dufresne/UofM/techinfo/kapton.html
    //The chemical formula of Kapton is C22H10N205,  its density is 1.43
    // also see below, do we need more than one?

    // Superconducting Cable
    mat = uniqueMaterialOrThrow( "SCCable"); // FIXME verify it
    {
      G4Material* SCCable = new G4Material( mat.name, 3.95*CLHEP::g/CLHEP::cm3, 3);
      SCCable->AddMaterial(findMaterialOrThrow("SCCableInsulation"), 0.04);
      SCCable->AddMaterial(findMaterialOrThrow("AL999Ni001"),        0.43);
      SCCable->AddMaterial(findMaterialOrThrow("NbTiCu"),            0.53);
    }


    // TS Collimator Cu mix (includes fiberglass and poly wrap/fill,
    // assumed about 2%)
    mat = uniqueMaterialOrThrow( "CollCu" );
    {
      G4Material* CollCu = new G4Material( mat.name, 8.815*CLHEP::g/CLHEP::cm3, 2);
      CollCu->AddMaterial( findMaterialOrThrow("G4_Cu"),0.98 );
      CollCu->AddMaterial( findMaterialOrThrow("G10"),  0.02 );
    }



    mat = uniqueMaterialOrThrow( "IsoButane");
    {
      G4Material* IsoButane = new G4Material( mat.name, 0.00265*CLHEP::g/CLHEP::cm3, 2);
      IsoButane->AddElement( getElementOrThrow("C"), 4);
      IsoButane->AddElement( getElementOrThrow("H"), 10);
    }

    mat = uniqueMaterialOrThrow( "StrawGasArCF4");
    {
      G4Material* StrawGasArCF4 = new G4Material(mat.name, 0.0028561*CLHEP::g/CLHEP::cm3, 3); // it is OK not to use kStateGas
      StrawGasArCF4->AddElement( getElementOrThrow("Ar"), 1);
      StrawGasArCF4->AddElement( getElementOrThrow("C") , 1);
      StrawGasArCF4->AddElement( getElementOrThrow("F") , 4);
    }

    mat = uniqueMaterialOrThrow( "TrackerManifold"); // materials and proportions defined in doc888v7
    {
      G4Material* TrackerManifold = new G4Material( mat.name, 1.95*CLHEP::g/CLHEP::cm3, 2);
      TrackerManifold->AddMaterial(findMaterialOrThrow("G4_POLYVINYL_CHLORIDE"), 0.355);
      TrackerManifold->AddMaterial(findMaterialOrThrow("G4_Al"), 0.645);
    }

    mat = uniqueMaterialOrThrow( "StrawWall"); // materials and proportions defined in doc888v7
    {
      G4Material* StrawWall = new G4Material( mat.name, 1.43*CLHEP::g/CLHEP::cm3, 3);
      StrawWall->AddMaterial(findMaterialOrThrow("G4_MYLAR"), 0.97);
      StrawWall->AddMaterial(findMaterialOrThrow("G4_Au"), 0.018);
      StrawWall->AddMaterial(findMaterialOrThrow("G4_Al"), 0.012);
    }

    mat = uniqueMaterialOrThrow( "StrawGas");
    {
      G4double density;
      G4double temperature = 293.15*CLHEP::kelvin;
      G4double pressure = 1.*CLHEP::atmosphere;
      G4int nel;

      G4double densityAr   = 0.00166 *CLHEP::g/CLHEP::cm3; //from PDG
      G4double densityCO2  = 0.00184 *CLHEP::g/CLHEP::cm3; //from PDG
      G4double fractionAr  = 80.0*CLHEP::perCent;

      density = fractionAr*densityAr + (1.0-fractionAr)*densityCO2;

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature, pressure);

      G4Element* Ar = getElementOrThrow("Ar");
      G4Element* C  = getElementOrThrow("C");
      G4Element* O  = getElementOrThrow("O");

      G4double atomicWeight_Ar =  39.948 *CLHEP::g/CLHEP::mole;
      G4double atomicWeight_C  = 12.0107 *CLHEP::g/CLHEP::mole;
      G4double atomicWeight_O  = 15.9994 *CLHEP::g/CLHEP::mole;
      G4double pwAr = fractionAr*atomicWeight_Ar;
      G4double pwC  = (1.0-fractionAr) *  1.0*atomicWeight_C;
      G4double pwO  = (1.0-fractionAr) *  2.0*atomicWeight_O;
      G4double atomicWeightMix = pwAr + pwC + pwO ;

      pwAr/=atomicWeightMix;
      pwO/=atomicWeightMix;
      GasMix->AddElement(Ar, pwAr );
      GasMix->AddElement(O , pwO  );
      GasMix->AddElement(C , 1.0-pwAr-pwO  );
    }

    mat = uniqueMaterialOrThrow( "Kapton");
    {
      //
      // Kapton: from NIST: physics.nist.gov/cgi-bin/Star/compos.pl?matno=179
      //
      G4Material* Kapton = new G4Material(mat.name, 1.42*CLHEP::g/CLHEP::cm3, 4);
      Kapton->AddElement( getElementOrThrow("H"), 0.026362);
      Kapton->AddElement( getElementOrThrow("C"), 0.691133);
      Kapton->AddElement( getElementOrThrow("N"), 0.073270);
      Kapton->AddElement( getElementOrThrow("O"), 0.209235);
    }

    mat = uniqueMaterialOrThrow( "Scintillator");
    {
      //
      // Scintillator.
      // We probably want several flavors of scintillator so that we can change the
      // detector by just changing a name in the config file.
      //
      G4Material* Sci = new G4Material( mat.name, 1.032*CLHEP::g/CLHEP::cm3, 2);
      Sci->AddElement( getElementOrThrow("C"), 9);
      Sci->AddElement( getElementOrThrow("H"), 10);
    }

    // These are materials used for backfill around Mu2e Building
    mat = uniqueMaterialOrThrow( "Calcite" );
    {
      G4Material* Calcite = new G4Material( mat.name, 2.71*CLHEP::g/CLHEP::cm3, 3 );
      Calcite->AddElement( getElementOrThrow("Ca"), 1);
      Calcite->AddElement( getElementOrThrow("C"), 1);
      Calcite->AddElement( getElementOrThrow("O"), 3);
    }

    mat = uniqueMaterialOrThrow( "Dolomite" );
    {
      G4Material* Dolomite = new G4Material( mat.name, 2.84*CLHEP::g/CLHEP::cm3, 4 );
      Dolomite->AddElement( getElementOrThrow("Ca"), 1);
      Dolomite->AddElement( getElementOrThrow("Mg"), 1);
      Dolomite->AddElement( getElementOrThrow("C"), 2);
      Dolomite->AddElement( getElementOrThrow("O"), 6);
    }

    // The following are some "engineered" fill materials, with information
    // provided by Tom Hamernik of Fermilab.  The materials are limestone,
    // which is typically a mixture of Calcite and Dolomite, just defined.

    mat = uniqueMaterialOrThrow( "CA7Backfill" );
    {
      G4Material* CA7Backfill = new G4Material( mat.name, 1.64*CLHEP::g/CLHEP::cm3, 2 );
      CA7Backfill->AddMaterial( findMaterialOrThrow("Calcite"),0.60 );
      CA7Backfill->AddMaterial( findMaterialOrThrow("Dolomite"),0.40);
    }

    mat = uniqueMaterialOrThrow( "RockBackfill" );
    {
      G4Material* RockBackfill = new G4Material( mat.name, 2.0*CLHEP::g/CLHEP::cm3, 2 );
      RockBackfill->AddMaterial( findMaterialOrThrow("Calcite"),0.60 );
      RockBackfill->AddMaterial( findMaterialOrThrow("Dolomite"),0.40);
    }

    // Vacuum(s)
    mat = uniqueMaterialOrThrow( "WAGVacuum");
    {
      //
      // This is the lowest density vacuum allowed by G4.
      G4double density     = CLHEP::universe_mean_density;
      G4double pressure    = 3.e-18*CLHEP::pascal;
      G4double temperature = 2.73*CLHEP::kelvin;

      // G4 takes ownership of this object and manages its lifetime.
      new G4Material( mat.name, 1., 1.01 *CLHEP::g/CLHEP::mole,
                      density, kStateGas, temperature, pressure);
    }


    // Presume that the residual gas in the DS will be leakage from the straws.
    mat = uniqueMaterialOrThrow( "DSVacuum");
    {

      const double oneTorr = CLHEP::atmosphere/760.;
      GeomHandle<DetectorSolenoid> ds;

      G4Material* StrawLeak = findMaterialOrThrow("StrawGas");

      G4double temperature = 300.00*CLHEP::kelvin; // Temperature of the DS
      G4double pressure    = oneTorr * ds->vac_pressure();
      G4double refTemp     = StrawLeak->GetTemperature();
      G4double refPress    = StrawLeak->GetPressure();

      G4double density = StrawLeak->GetDensity()*pressure*refTemp/(refPress*temperature);

      G4Material* DSVacuum =
        new G4Material(mat.name, density, StrawLeak->GetNumberOfElements(),
                       kStateGas, temperature, pressure);

      for (size_t i = 0 ; i < StrawLeak->GetNumberOfElements(); ++i) {
        DSVacuum->AddElement(const_cast<G4Element*>(StrawLeak->GetElementVector()->at(i)), StrawLeak->GetFractionVector()[i]);
      }

    }

    // Presume that the residual gas in the PS will be purged with Argon
    mat = uniqueMaterialOrThrow( "PSVacuum");
    {

      const double oneTorr = CLHEP::atmosphere/760.;
      GeomHandle<PSVacuum> ps;

      G4Material* gas = findMaterialOrThrow(ps->vacuumG4Material());

      G4double temperature = 300.00*CLHEP::kelvin;
      G4double pressure    = oneTorr * ps->vacuumPressure();
      G4double refTemp     = gas->GetTemperature();
      G4double refPress    = gas->GetPressure();

      G4double density = gas->GetDensity()*pressure*refTemp/(refPress*temperature);

      G4Material* PSVacuum =
        new G4Material(mat.name, density, gas->GetNumberOfElements(),
                       kStateGas, temperature, pressure);

      for (size_t i = 0 ; i < gas->GetNumberOfElements(); ++i) {
        PSVacuum->AddElement(const_cast<G4Element*>(gas->GetElementVector()->at(i)), gas->GetFractionVector()[i]);
      }

    }



    mat = uniqueMaterialOrThrow( "MBOverburden");
    {
      //
      // MiniBoone model of the earthen overburden.  See Mu2e-doc-570.
      //
      G4Material* mbOverburden = new G4Material( mat.name, 2.15*CLHEP::g/CLHEP::cm3, 3);
      G4Element* eO  = getElementOrThrow("O");
      G4Element* eSi = getElementOrThrow("Si");
      G4Element* eAl = getElementOrThrow("Al");
      mbOverburden->AddElement( eO,  65);
      mbOverburden->AddElement( eSi, 20);
      mbOverburden->AddElement( eAl, 15);
    }

    mat = uniqueMaterialOrThrow( "ITGasHe_95Isob_5");
    {

      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.000166 *CLHEP::g/CLHEP::cm3;
      G4double densityIsoB = 0.00249  *CLHEP::g/CLHEP::cm3;
      G4double fractionHe  = 95.0*CLHEP::perCent;

      density = fractionHe*densityHe + (1.0-fractionHe)*densityIsoB;

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature= 293.15*CLHEP::kelvin,
                                           pressure= 1.*CLHEP::atmosphere);

      G4Element* He = getElementOrThrow("He");
      G4Element* C  = getElementOrThrow("C");
      G4Element* H  = getElementOrThrow("H");

      G4double atomicWeight_He =  4.002602 *CLHEP::g/CLHEP::mole;
      G4double atomicWeight_C  = 12.0107   *CLHEP::g/CLHEP::mole;
      G4double atomicWeight_H  =  1.00794  *CLHEP::g/CLHEP::mole;
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

    mat = uniqueMaterialOrThrow( "ITGasHe_90Isob_10");
    {

      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.000166 *CLHEP::g/CLHEP::cm3;
      G4double densityIsoB = 0.00249  *CLHEP::g/CLHEP::cm3;
      G4double fractionHe  = 90.0*CLHEP::perCent;

      density = fractionHe*densityHe + (1.0-fractionHe)*densityIsoB;

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature= 293.15*CLHEP::kelvin,
                                           pressure= 1.*CLHEP::atmosphere);

      G4Element* He = getElementOrThrow("He");
      G4Element* C  = getElementOrThrow("C");
      G4Element* H  = getElementOrThrow("H");

      G4double atomicWeight_He =  4.002602 *CLHEP::g/CLHEP::mole;
      G4double atomicWeight_C  = 12.0107   *CLHEP::g/CLHEP::mole;
      G4double atomicWeight_H  =  1.00794  *CLHEP::g/CLHEP::mole;
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

    mat = uniqueMaterialOrThrow( "ITGasHe_75Isob_25_400mbar");
    {

      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.000166 *CLHEP::g/CLHEP::cm3;
      G4double densityIsoB = 0.00249  *CLHEP::g/CLHEP::cm3;
      G4double fractionHe  = 75.0*CLHEP::perCent;

      density = fractionHe*densityHe + (1.0-fractionHe)*densityIsoB;
      pressure = 0.4*CLHEP::bar;
      density *= pressure/(1.0*CLHEP::atmosphere);

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature= 293.15*CLHEP::kelvin, pressure);

      G4Element* He = getElementOrThrow("He");
      G4Element* C  = getElementOrThrow("C");
      G4Element* H  = getElementOrThrow("H");

      G4double atomicWeight_He =  4.002602 *CLHEP::g/CLHEP::mole;
      G4double atomicWeight_C  = 12.0107   *CLHEP::g/CLHEP::mole;
      G4double atomicWeight_H  =  1.00794  *CLHEP::g/CLHEP::mole;
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

    mat = uniqueMaterialOrThrow( "ITGasHe_90CF4_10");
    {

      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.000166 *CLHEP::g/CLHEP::cm3;
      G4double densityCF4  = 0.003780 *CLHEP::g/CLHEP::cm3;
      G4double fractionHe  = 90.0*CLHEP::perCent;

      density = fractionHe*densityHe + (1.0-fractionHe)*densityCF4;

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature= 293.15*CLHEP::kelvin,
                                           pressure= 1.*CLHEP::atmosphere);

      G4Element* He = getElementOrThrow("He");
      G4Element* C  = getElementOrThrow("C");
      G4Element* F  = getElementOrThrow("F");

      G4double atomicWeight_He =  4.002602  *CLHEP::g/CLHEP::mole;
      G4double atomicWeight_C  = 12.0107    *CLHEP::g/CLHEP::mole;
      G4double atomicWeight_F  = 18.9984032 *CLHEP::g/CLHEP::mole;
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

    mat = uniqueMaterialOrThrow( "ITGasMix");
    {
      //He/C4H10-gas-mixture

      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.0001786*CLHEP::g/CLHEP::cm3;
      G4double densityIsoB = 0.00267  *CLHEP::g/CLHEP::cm3;
      G4double fractionHe  = 90.0*CLHEP::perCent;

      density = fractionHe*densityHe + (1.0-fractionHe)*densityIsoB;

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature= 293.15*CLHEP::kelvin,
                                           pressure= 1.*CLHEP::atmosphere);

      G4Element* He = getElementOrThrow("He");
      G4Element* C  = getElementOrThrow("C");
      G4Element* H  = getElementOrThrow("H");

      GasMix->AddElement(He, 0.9   );
      GasMix->AddElement(H , 0.0173);
      GasMix->AddElement(C , 0.0827);
    }

    mat = uniqueMaterialOrThrow( "ITGasVacuum");
    {
      //
      // This is the lowest density vacuum allowed by G4.
      G4double density     = CLHEP::universe_mean_density;
      G4double pressure    = 3.e-18*CLHEP::pascal;
      G4double temperature = 2.73*CLHEP::kelvin;

      // G4 takes ownership of this object and manages its lifetime.
      new G4Material( mat.name, 1., 1.01 *CLHEP::g/CLHEP::mole,
                      density, kStateGas, temperature, pressure);
    }

    //Material for the MBS ring to protect the calorimeter
    mat = uniqueMaterialOrThrow( "MBSCalShieldRing" );
    {
      //Shield is 90% tungsten 10% copper and 18 g/cm^3
      G4Material* MBSCalShieldRing = new G4Material( mat.name, 18.0*CLHEP::g/CLHEP::cm3, 2);
      MBSCalShieldRing->AddElement( getElementOrThrow("W"),  90.0*CLHEP::perCent);
      MBSCalShieldRing->AddElement( getElementOrThrow("Cu"), 10.0*CLHEP::perCent);
    }

    // Completed constructing Mu2e specific materials, first function.
    // Add new materials et the end of the  second function for Mu2e specific materials.

  }


  void ConstructMaterials::constructMu2eMaterials2(){


    // An alias for the stopping target material
    CheckedG4String mat = uniqueMaterialOrThrow("StoppingTarget_"+GlobalConstantsHandle<PhysicsParams>()->getStoppingTargetMaterial());
    if ( true /* Always load the stopping target material */ ){
      G4Material* met = findMaterialOrThrow("G4_"+GlobalConstantsHandle<PhysicsParams>()->getStoppingTargetMaterial());
      G4Material* tgt = new G4Material(mat.name, met->GetDensity(), 1);
      tgt->AddMaterial(met, 1.);
    }

    mat = uniqueMaterialOrThrow( "CarbonFiber_resin");
    {
      G4double density;
      G4int nel;
      G4Material* CFresin =
        new G4Material(mat.name, density = 1.1*CLHEP::g/CLHEP::cm3, nel=3);
      G4int natoms;
      CFresin->AddElement(getElementOrThrow("H"),natoms=5);
      CFresin->AddElement(getElementOrThrow("C"),natoms=5);
      CFresin->AddElement(getElementOrThrow("O"),natoms=2);
    }

    mat = uniqueMaterialOrThrow( "CarbonFiber");
    {
      G4double density, fiberFrac=46.0*CLHEP::perCent;
      G4int nel, natoms;
      G4Material* CFresin = findMaterialOrThrow("CarbonFiber_resin");
      CheckedG4String mat0 = uniqueMaterialOrThrow("CFibers");
      G4Material* CFibers = new G4Material(mat0.name, density = 1.8*CLHEP::g/CLHEP::cm3,nel=1);
      CFibers->AddElement(getElementOrThrow("C"),natoms=1);

      density = fiberFrac*CFibers->GetDensity()+(1.0-fiberFrac)*CFresin->GetDensity();
      G4Material* CarbonFiber =
        new G4Material(mat.name, density /*= 1.384*CLHEP::g/CLHEP::cm3*/, nel=2);
      CarbonFiber->AddMaterial(CFibers, fiberFrac );
      CarbonFiber->AddMaterial(CFresin, (1.0-fiberFrac) );
    }

    mat = uniqueMaterialOrThrow("PET_P100");
    G4Material* PET_P100 = new G4Material( mat.name, 0.11*CLHEP::g/CLHEP::cm3, 3);
    PET_P100->AddElement( getElementOrThrow("C"), 10);
    PET_P100->AddElement( getElementOrThrow("H"),  8);
    PET_P100->AddElement( getElementOrThrow("O"),  4);

    mat = uniqueMaterialOrThrow( "Lyso_01");  /// Alessandra
    {
      G4double density;
      G4int nel;
      CheckedG4String mat0 = uniqueMaterialOrThrow( "Lyso_00");
      G4Material* Lyso_00 =
        new G4Material(mat0.name, density = 7.4*CLHEP::g/CLHEP::cm3, nel=4);
      G4Element* Lu  = getElementOrThrow("Lu");
      G4Element* Si  = getElementOrThrow("Si");
      G4Element* O  = getElementOrThrow("O");
      G4Element* Y  = getElementOrThrow("Y");
      G4Element* Ce  = getElementOrThrow("Ce");
      Lyso_00->AddElement( Lu, 71.0*CLHEP::perCent );
      Lyso_00->AddElement( Si, 7.0*CLHEP::perCent );
      Lyso_00->AddElement( O, 18.0*CLHEP::perCent );
      Lyso_00->AddElement( Y, 4.0*CLHEP::perCent );
      G4Material* Lyso_01 =
        new G4Material(mat.name, density = 7.4*CLHEP::g/CLHEP::cm3, nel=2);
      Lyso_01->AddMaterial( Lyso_00, 99.85*CLHEP::perCent );
      Lyso_01->AddElement( Ce, 0.15*CLHEP::perCent );
    }

    mat = uniqueMaterialOrThrow( "CuW1090");
    G4Material* CuW1090 = new G4Material(mat.name, 17.3*CLHEP::g/CLHEP::cm3, 2);
    CuW1090->AddMaterial( findMaterialOrThrow("G4_W"),0.90);
    CuW1090->AddMaterial( findMaterialOrThrow("G4_Cu"),0.10);


    mat = uniqueMaterialOrThrow("CorrugatedPolypropylene");
    {
      G4double standardPolypropyleneDensity = 0.946*CLHEP::g/CLHEP::cm3;
      G4double corrugatedEffectiveDensity = 0.2 * standardPolypropyleneDensity;
      G4Material* CorrugatedPolypropylene = new G4Material( mat.name, corrugatedEffectiveDensity, 2);
      CorrugatedPolypropylene->AddElement( getElementOrThrow("C"), 3);
      CorrugatedPolypropylene->AddElement( getElementOrThrow("H"), 6);
    }


    //G10-FR4 used for printed board of the I-Tracker
    // G10 http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/SDD_G10FR4.html
    // http://pdg.lbl.gov/2002/atomicrpp.pdf
    mat = uniqueMaterialOrThrow( "G10_FR4");
    {
      G4double density;

      G4Material* G10_FR4 =
        new G4Material(mat.name, density = 1.8*CLHEP::g/CLHEP::cm3, 2);
      G10_FR4->AddMaterial(findMaterialOrThrow("EGlass"), 0.60);
      G10_FR4->AddMaterial(findMaterialOrThrow("Epotek301"), 0.40);
    }


    mat = uniqueMaterialOrThrow( "Electronics2" );
    {
      // This material represents the passive part of the board in the Calorimeter crates
      G4double density;

      G4Material* Electronics2 =
        new G4Material(mat.name, density = 4.52*CLHEP::g/CLHEP::cm3, 2);
      Electronics2->AddMaterial(findMaterialOrThrow("G10_FR4"), 0.26);
      Electronics2->AddMaterial(findMaterialOrThrow("G4_Cu"), 0.74);
    }

    mat = uniqueMaterialOrThrow( "AluminumHoneycomb");
    {
      //Honeycomb used to fill the source-panel and the inner step margins of the calorimeter
      G4double density;
      G4int nel;
      G4Material *AluminumHoneycomb = new G4Material(mat.name, density = 0.03*CLHEP::g/CLHEP::cm3, nel=1);
      AluminumHoneycomb->AddMaterial(findMaterialOrThrow("G4_Al"), 100.0*CLHEP::perCent);
    }

    mat = uniqueMaterialOrThrow( "PolypropyleneFoam");
    {
      //Polypropylene (CH3)
      G4double density;
      G4int nel;
      G4Material *Polypropylene = new G4Material(mat.name, density = 0.04*CLHEP::g/CLHEP::cm3, nel=2);
      G4Element* H  = getElementOrThrow("H");
      G4Element* C  = getElementOrThrow("C");
      Polypropylene->AddElement(H, 3 );
      Polypropylene->AddElement(C, 1 );
    }

    mat = uniqueMaterialOrThrow( "BeFoam_018");
    {
      G4double density;
      G4int nel;
      G4Material *BeFoam_018 = new G4Material(mat.name, density = 0.018*CLHEP::g/CLHEP::cm3, nel=1);
      G4Element* Be  = getElementOrThrow("Be");
      BeFoam_018->AddElement(Be, 100.0*CLHEP::perCent );
    }

    mat = uniqueMaterialOrThrow( "CFoam_332");
    {
      G4double density;
      G4int nel;
      G4Material *CFoam = new G4Material(mat.name, density = 0.332*CLHEP::g/CLHEP::cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CFoam->AddElement(C, 100.0*CLHEP::perCent );
    }

    mat = uniqueMaterialOrThrow( "CFoam_166");
    {
      G4double density;
      G4int nel;
      G4Material *CFoam = new G4Material(mat.name, density = 0.166*CLHEP::g/CLHEP::cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CFoam->AddElement(C, 100.0*CLHEP::perCent );
    }

    mat = uniqueMaterialOrThrow( "CFoam_080");
    {
      G4double density;
      G4int nel;
      G4Material *CFoam = new G4Material(mat.name, density = 0.080*CLHEP::g/CLHEP::cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CFoam->AddElement(C, 100.0*CLHEP::perCent );
    }

    mat = uniqueMaterialOrThrow( "CFoam");
    {
      G4double density;
      G4int nel;
      G4Material *CFoam = new G4Material(mat.name, density = 0.030*CLHEP::g/CLHEP::cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CFoam->AddElement(C, 100.0*CLHEP::perCent );
    }

    mat = uniqueMaterialOrThrow( "KptFoam_030");
    {
      G4double density;
      G4int nel;
      G4Material *KptFoam = new G4Material(mat.name, density = 0.030*CLHEP::g/CLHEP::cm3, nel=1);
      KptFoam->AddMaterial(findMaterialOrThrow("G4_KAPTON"), 100.0*CLHEP::perCent );
    }

    mat = uniqueMaterialOrThrow( "ZirconiumHydridePolyethylene");
    {
      G4Material* ZirconiumHydridePolyethylene =
        new G4Material( mat.name, 3.67*CLHEP::g/CLHEP::cm3, 5);
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eH  = getElementOrThrow("H");
      G4Element* eB  = getElementOrThrow("B");
      G4Element* eO  = getElementOrThrow("O");
      G4Element* eZr  = getElementOrThrow("Zr");

      ZirconiumHydridePolyethylene->AddElement( eC, 8.9*CLHEP::perCent);
      ZirconiumHydridePolyethylene->AddElement( eH, 3.4*CLHEP::perCent);
      ZirconiumHydridePolyethylene->AddElement( eB, 0.5*CLHEP::perCent);
      ZirconiumHydridePolyethylene->AddElement( eO, 2.2*CLHEP::perCent);
      ZirconiumHydridePolyethylene->AddElement( eZr, 85.0*CLHEP::perCent);
    }

    mat = uniqueMaterialOrThrow( "StrawWallEq");
    {
      G4double density;
      G4int nel;
      G4Material* strwMl = findMaterialOrThrow("G4_MYLAR");
      G4Material* strwMet1 = findMaterialOrThrow("G4_Au");
      G4Material* strwMet2 = findMaterialOrThrow("G4_Al");

      density = 1.4325*CLHEP::g/CLHEP::cm3;
      G4Material* stWallEq =
        new G4Material(mat.name, density, nel=3);
      stWallEq->AddMaterial(strwMl, 96.95e-2 );
      stWallEq->AddMaterial(strwMet1, 1.80e-2 );
      stWallEq->AddMaterial(strwMet2, 1.25e-2 );
    }


    // scaled W for hayman studies
    mat = uniqueMaterialOrThrow("G4_W_Hayman");
    {  //220 mm with gaps; this models as lower density without gaps
      G4int nel;
      G4Material* tung = findMaterialOrThrow("G4_W");
      G4Material*  G4_W_Hayman = new G4Material(mat.name, (80./110.0)*tung->GetDensity(), nel = 1);
      G4Element* Tung  = getElementOrThrow("W");
      G4_W_Hayman->AddElement(Tung, 100.0*CLHEP::perCent );
    }

    // various densities of Al and Be to permit staging of pbar window studies without changing geometry
    // between stages
    mat = uniqueMaterialOrThrow( "G4_Be_Quarter");
    {
      G4int nel;
      G4Material* G4_Be_Quarter = new G4Material(mat.name, 0.25*1.85*CLHEP::g/CLHEP::cm3, nel = 1);
      G4Element* Be  = getElementOrThrow("Be");
      G4_Be_Quarter->AddElement(Be, 100.0*CLHEP::perCent );
    }
    mat = uniqueMaterialOrThrow( "G4_Be_Half");
    {
      G4int nel;
      G4Material*  G4_Be_Half = new G4Material(mat.name, 0.50*1.85*CLHEP::g/CLHEP::cm3, nel = 1);
      G4Element* Be  = getElementOrThrow("Be");
      G4_Be_Half->AddElement(Be, 100.0*CLHEP::perCent );
    }
    mat = uniqueMaterialOrThrow( "G4_Be_Standard");
    {
      G4int nel;
      G4Material*  G4_Be_Standard = new G4Material(mat.name, 1.0*1.85*CLHEP::g/CLHEP::cm3, nel = 1);
      G4Element* Be  = getElementOrThrow("Be");
      G4_Be_Standard->AddElement(Be, 100.0*CLHEP::perCent );
    }
    mat = uniqueMaterialOrThrow( "G4_Be_Double");
    {
      G4int nel;
      G4Material*  G4_Be_Double = new G4Material(mat.name, 2.0*1.85*CLHEP::g/CLHEP::cm3, nel = 1);
      G4Element* Be  = getElementOrThrow("Be");
      G4_Be_Double->AddElement(Be, 100.0*CLHEP::perCent );
    }
    mat = uniqueMaterialOrThrow( "G4_Be_Triple");
    {
      G4int nel;
      G4Material*  G4_Be_Triple = new G4Material(mat.name, 3.0*1.85*CLHEP::g/CLHEP::cm3, nel = 1);
      G4Element* Be  = getElementOrThrow("Be");
      G4_Be_Triple->AddElement(Be, 100.0*CLHEP::perCent );
    }
    mat = uniqueMaterialOrThrow( "G4_Al_Quarter");
    {
      G4int nel;
      G4Material*  G4_Al_Quarter = new G4Material(mat.name, 0.25*2.70*CLHEP::g/CLHEP::cm3, nel = 1);
      G4Element* Al  = getElementOrThrow("Al");
      G4_Al_Quarter->AddElement(Al, 100.0*CLHEP::perCent );
    }
    mat = uniqueMaterialOrThrow( "G4_Al_Half");
    {
      G4int nel;
      G4Material*  G4_Al_Half = new G4Material(mat.name, 0.50*2.70*CLHEP::g/CLHEP::cm3, nel = 1);
      G4Element* Al  = getElementOrThrow("Al");
      G4_Al_Half->AddElement(Al, 100.0*CLHEP::perCent );

    }
    mat = uniqueMaterialOrThrow( "G4_Al_Standard");
    {
      G4int nel;
      G4Material*  G4_Al_Standard = new G4Material("G4_Al_Standard", 1.0*2.70*CLHEP::g/CLHEP::cm3, nel = 1);
      G4Element* Al  = getElementOrThrow("Al");
      G4_Al_Standard->AddElement(Al, 100.0*CLHEP::perCent );
    }
    mat = uniqueMaterialOrThrow( "G4_Al_Double");
    {
      G4int nel;
      G4Material*  G4_Al_Double = new G4Material(mat.name, 2.0*2.70*CLHEP::g/CLHEP::cm3, nel = 1);
      G4Element* Al  = getElementOrThrow("Al");
      G4_Al_Double->AddElement(Al, 100.0*CLHEP::perCent );
    }
    mat = uniqueMaterialOrThrow( "G4_Al_Triple");
    {
      G4int nel;
      G4Material*  G4_Al_Triple = new G4Material(mat.name, 3.0*2.70*CLHEP::g/CLHEP::cm3, nel = 1);
      G4Element* Al  = getElementOrThrow("Al");
      G4_Al_Triple->AddElement(Al, 100.0*CLHEP::perCent );
    }


    //Material for the Calorimeter cable runs bulk
    mat = uniqueMaterialOrThrow( "CalCableRunOuter" );
    {

      G4Material* CalCableRunOuter =
        new G4Material( mat.name, 6.8440*CLHEP::g/CLHEP::cm3, 4);
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eCu = getElementOrThrow("Cu");
      G4Element* eAg = getElementOrThrow("Ag");
      G4Element* eF  = getElementOrThrow("F");

      //Add elements by mass fraction
      CalCableRunOuter->AddElement( eAg, 52.603*CLHEP::perCent);
      CalCableRunOuter->AddElement( eCu, 34.825*CLHEP::perCent);
      CalCableRunOuter->AddElement( eF,  9.553*CLHEP::perCent);
      CalCableRunOuter->AddElement( eC,  3.019*CLHEP::perCent);
    }

    //Material for the Calorimeter cable runs fiber optic cable
    mat = uniqueMaterialOrThrow( "CalCableRunFiber" );
    {

      G4Material* CalCableRunFiber =
        new G4Material( mat.name, 1.6575*CLHEP::g/CLHEP::cm3, 6);
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eSi = getElementOrThrow("Si");
      G4Element* eO  = getElementOrThrow("O");
      G4Element* eAl = getElementOrThrow("Al");
      G4Element* eF  = getElementOrThrow("F");
      G4Element* eH  = getElementOrThrow("H");

      //Add elements by mass fraction
      CalCableRunFiber->AddElement( eC,   34.093*CLHEP::perCent);
      CalCableRunFiber->AddElement( eO,   28.578*CLHEP::perCent);
      CalCableRunFiber->AddElement( eAl,  28.265*CLHEP::perCent);
      CalCableRunFiber->AddElement( eH,   5.676*CLHEP::perCent);
      CalCableRunFiber->AddElement( eSi,  3.172*CLHEP::perCent);
      CalCableRunFiber->AddElement( eF,   0.216*CLHEP::perCent);
    }

    //Material for the Calorimeter cable runs bulk
    mat = uniqueMaterialOrThrow( "TrkCableRunOuter" );
    {

      G4Material* TrkCableRunOuter =
        new G4Material( mat.name, 4.6254*CLHEP::g/CLHEP::cm3, 6);
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eCu = getElementOrThrow("Cu");
      G4Element* eH  = getElementOrThrow("H");
      G4Element* eN  = getElementOrThrow("N");
      G4Element* eSi = getElementOrThrow("Si");
      G4Element* eO  = getElementOrThrow("O");

      //Add elements by mass fraction
      TrkCableRunOuter->AddElement( eCu, 71.458*CLHEP::perCent);
      TrkCableRunOuter->AddElement( eC,  10.259*CLHEP::perCent);
      TrkCableRunOuter->AddElement( eH,  2.099*CLHEP::perCent);
      TrkCableRunOuter->AddElement( eN,  0.197*CLHEP::perCent);
      TrkCableRunOuter->AddElement( eO,  6.161*CLHEP::perCent);
      TrkCableRunOuter->AddElement( eSi, 9.826*CLHEP::perCent);
    }

    //Material for the Calorimeter cable runs fiber optic cable
    mat = uniqueMaterialOrThrow( "TrkCableRunFiber" );
    {

      G4Material* TrkCableRunFiber =
        new G4Material( mat.name, 1.6575*CLHEP::g/CLHEP::cm3, 6);
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eSi = getElementOrThrow("Si");
      G4Element* eO  = getElementOrThrow("O");
      G4Element* eAl = getElementOrThrow("Al");
      G4Element* eF  = getElementOrThrow("F");
      G4Element* eH  = getElementOrThrow("H");

      //Add elements by mass fraction
      TrkCableRunFiber->AddElement( eC,   34.093*CLHEP::perCent);
      TrkCableRunFiber->AddElement( eO,   28.578*CLHEP::perCent);
      TrkCableRunFiber->AddElement( eAl,  28.265*CLHEP::perCent);
      TrkCableRunFiber->AddElement( eH,   5.676*CLHEP::perCent);
      TrkCableRunFiber->AddElement( eSi,  3.172*CLHEP::perCent);
      TrkCableRunFiber->AddElement( eF,   0.216*CLHEP::perCent);
    }

    mat = uniqueMaterialOrThrow( "MLI"); // assuming 15 ~1.4g/cm^3 mylar layers with ~13 um thickness each becomes ~5 mm thick blanket
    {
      G4Material* mli = new G4Material( mat.name, 0.055*CLHEP::g/CLHEP::cm3, 1);
      mli->AddMaterial(findMaterialOrThrow("G4_MYLAR"), 1.0);
    }

    mat = uniqueMaterialOrThrow( "ST_Wires"); // assuming 6% gold 94% tungsten, from docdb-31260 5-7% gold plating expected
    {
      G4Material* wires = new G4Material( mat.name, 19.25*CLHEP::g/CLHEP::cm3, 2); //
      wires->AddMaterial(findMaterialOrThrow("G4_W"), 0.94);
      wires->AddMaterial(findMaterialOrThrow("G4_Au"), 0.06);
    }

    //Carbon steel
    mat = uniqueMaterialOrThrow( "MildSteel"); // DocDB-42993
    {
      G4Material* MildSteel = new G4Material(mat.name, 7.86*CLHEP::g/CLHEP::cm3, 4); //used an example density of steel
      MildSteel->AddMaterial(findMaterialOrThrow("G4_Mn"), 0.0100);
      MildSteel->AddMaterial(findMaterialOrThrow("G4_Si"), 0.0025);
      MildSteel->AddMaterial(findMaterialOrThrow("G4_C" ), 0.0025);
      MildSteel->AddMaterial(findMaterialOrThrow("G4_Fe"), 0.9850);
    }

    // Completed constructMu2eMaterials2(), second function for
    // building all Mu2e materials.



    mat = uniqueMaterialOrThrow("ProductionTargetTungstenLa2_O3");
    {

     G4Material* ProductionTargetTungstenLa2_O3 = new G4Material(mat.name
                                                                 ,18.75*CLHEP::g/CLHEP::cm3
                                                                 ,2);
     constexpr double wPercentage = 99.;

     // first define lanthanum oxide
     //     G4Element* lanthanum = new G4Element("Lanthanum","La",57.,138.905*CLHEP::g/CLHEP::mole);
     //G4Element* oxygen = new G4Element("Oxygen","O",8.,15.999*CLHEP::g/CLHEP::mole);
     G4Material* La2_O3 = new G4Material("La2_O3", 6.51*CLHEP::g/CLHEP::cm3  ,2);
     La2_O3->AddElement(getElementOrThrow("La"),2);
     La2_O3->AddElement(getElementOrThrow("O"),3);
     //
     // and now tungsten
     ProductionTargetTungstenLa2_O3->AddMaterial(La2_O3,(100. - wPercentage)*CLHEP::perCent);
     ProductionTargetTungstenLa2_O3->AddElement(getElementOrThrow("W"),wPercentage*CLHEP::perCent);
    }


    mat = uniqueMaterialOrThrow("LaBr3Ce");
    {
     G4Material* LaBr3Ce = new G4Material(mat.name, 5.08*CLHEP::g/CLHEP::cm3, 2);
     G4Material* LaBr3 = new G4Material("LaBr3", 5.06*CLHEP::g/CLHEP::cm3  ,2);
     LaBr3->AddElement(getElementOrThrow("La"),1);
     LaBr3->AddElement(getElementOrThrow("Br"),3);

     G4Element* elCe  = new G4Element("Cerium"  ,"Ce" , 58., 140.116*CLHEP::g/CLHEP::mole);

     LaBr3Ce -> AddElement(elCe, 1.9*CLHEP::perCent);
     LaBr3Ce -> AddMaterial(LaBr3, 98.1*CLHEP::perCent);
    }

    mat = uniqueMaterialOrThrow("BP");   //Borated polyethylene
    {
     G4Material* BP  = new G4Material(mat.name, 1.04*CLHEP::g/CLHEP::cm3, 2);

     G4Element* elB  = new G4Element("Boron"  ,"B" , 5., 10.81*CLHEP::g/CLHEP::mole);
     G4Material* Poly = findMaterialOrThrow("G4_POLYETHYLENE");

     BP -> AddElement(elB, 5*CLHEP::perCent);
     BP -> AddMaterial(Poly, 95*CLHEP::perCent);
    }

    mat = uniqueMaterialOrThrow( "ClosePackedExtMonSteelShot");
    {
      //7.85 g/cm3 is the density of steel and 63.5% is the densest packing of spheres
      constexpr double density = 7.85 * .635;
      G4Material* ClosePackedExtMonSteelShot = new G4Material(mat.name, density*CLHEP::g/CLHEP::cm3, 1);
      ClosePackedExtMonSteelShot -> AddMaterial(findMaterialOrThrow("MildSteel"), 100*CLHEP::perCent);
    }

    //Information from https://en.wikipedia.org/wiki/Difluoromethane
    mat = uniqueMaterialOrThrow("R32");
    {
     G4Material* R32 = new G4Material(mat.name, 1.1*CLHEP::g/CLHEP::cm3, 3);

     G4Element* eC = getElementOrThrow("C");
     G4Element* eH  = getElementOrThrow("H");
     G4Element* eF = getElementOrThrow("F");

     R32->AddElement( eC,   1);
     R32->AddElement( eH,   2);
     R32->AddElement( eF,   2);

    }

    //Information from https://en.wikipedia.org/wiki/Pentafluoroethane
    mat = uniqueMaterialOrThrow("R125");
    {
     G4Material* R125 = new G4Material(mat.name, 1.53*CLHEP::g/CLHEP::cm3, 3);

     G4Element* eC = getElementOrThrow("C");
     G4Element* eH  = getElementOrThrow("H");
     G4Element* eF = getElementOrThrow("F");

     R125->AddElement( eC,   2);
     R125->AddElement( eH,   1);
     R125->AddElement( eF,   5);

    }

    //Information from https://en.wikipedia.org/wiki/2,3,3,3-Tetrafluoropropene
    mat = uniqueMaterialOrThrow("R1234yf");
    {
     G4Material* R1234yf = new G4Material(mat.name, 1.1*CLHEP::g/CLHEP::cm3, 3);

     G4Element* eC = getElementOrThrow("C");
     G4Element* eH  = getElementOrThrow("H");
     G4Element* eF = getElementOrThrow("F");

     R1234yf->AddElement( eC,   3);
     R1234yf->AddElement( eH,   2);
     R1234yf->AddElement( eF,   4);

    }

    //Information from https://en.wikipedia.org/wiki/1,1,1,2-Tetrafluoroethane
    mat = uniqueMaterialOrThrow("R134a");
    {
     G4Material* R134a = new G4Material(mat.name, 1.206*CLHEP::g/CLHEP::cm3, 3);

     G4Element* eC = getElementOrThrow("C");
     G4Element* eH  = getElementOrThrow("H");
     G4Element* eF = getElementOrThrow("F");

     R134a->AddElement( eC,   2);
     R134a->AddElement( eH,   2);
     R134a->AddElement( eF,   4);

    }

    //Information from https://www.honeywell-refrigerants.com/europe/wp-content/uploads/2017/10/FPR-029-2017-09_Solstice_452A_A4_2892017.pdf
    mat = uniqueMaterialOrThrow("R452A");
    {
     G4Material* R452A = new G4Material(mat.name, 1.1488*CLHEP::g/CLHEP::cm3, 3);

     R452A->AddMaterial( findMaterialOrThrow("R1234yf"), 30.0*CLHEP::perCent);
     R452A->AddMaterial( findMaterialOrThrow("R32"),     11.0*CLHEP::perCent);
     R452A->AddMaterial( findMaterialOrThrow("R125"),    59.0*CLHEP::perCent);

    }

    //information from https://www.opteon.com/en/-/media/files/opteon/opteon-xp40-prodinfo.pdf?la=en&rev=f82b8f89deec4f19bdbc7c17a04fe314
    mat = uniqueMaterialOrThrow("R449A");
    {
     G4Material* R449A = new G4Material(mat.name, 1.1141*CLHEP::g/CLHEP::cm3, 4);

     R449A->AddMaterial( findMaterialOrThrow("R32"),     24.3*CLHEP::perCent);
     R449A->AddMaterial( findMaterialOrThrow("R125"),    24.7*CLHEP::perCent);
     R449A->AddMaterial( findMaterialOrThrow("R1234yf"), 25.3*CLHEP::perCent);
     R449A->AddMaterial( findMaterialOrThrow("R134a"),   25.7*CLHEP::perCent);
    }

    //Information from https://julabo.us/wp-content/uploads/2023/03/Julabo-USA-SDS-Thermal-C5-1-1.pdf
    // and https://en.wikipedia.org/wiki/Polydimethylsiloxane
    // Assuming n = 1 for molecular formula
    mat = uniqueMaterialOrThrow("C5Coolant");
    {
     G4Material* C5Coolant = new G4Material(mat.name, 0.965*CLHEP::g/CLHEP::cm3, 4);

     G4Element* eH  = getElementOrThrow("H");
     G4Element* eC  = getElementOrThrow("C");
     G4Element* eSi = getElementOrThrow("Si");
     G4Element* eO  = getElementOrThrow("O");

     C5Coolant->AddElement( eH,   18);
     C5Coolant->AddElement( eC,   6);
     C5Coolant->AddElement( eSi,  2);
     C5Coolant->AddElement( eO,   1);

    }
    // Add new materials before this line

  }


  // Check to see if the named material already exists.
  CheckedG4String ConstructMaterials::uniqueMaterialOrThrow( G4String const& name){
    if ( G4Material::GetMaterial(name,false) != 0 ){
      throw cet::exception("GEOM")
        << "mu2e::ConstructMaterials::constructMu2eMaterials(): "
        << "The requested material is already defined: "
        << name
        << "\n";
    }
    return name;
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

  // Turning on density effect correction in ionization loss calculations
  // for selected materials deemed conductors
  void ConstructMaterials::useDensityEffectInIonizationLossCorrectionIfRequested() {

    if (config_.physics().useDensityEffectInIonizationLossCalc()) {

      std::vector<std::string> conductors = config_.physics().conductingMaterials();

      // conductor is defined as having GetFreeElectronDensity() > 0.0
      // one can set that density with SetFreeElectronDensity(val)

      if (config_.debug().diagLevel() > 0) {
        G4cout << "ConstructMaterials::" << __func__
               << " Using density effect correction in ionization loss calculations"
               << " for selected materials deemed conductors "
               << G4endl;
      }
      G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
      for(size_t i=0; i!=theMaterialTable->size(); ++i) {
        G4Material* theMaterial = (*theMaterialTable)[i];
        if (config_.debug().diagLevel() > 0) {
          G4String cond( ( theMaterial->GetFreeElectronDensity() > 0. ) ?
                             " conductor" : " not conductor" );
          G4cout << "ConstructMaterials::" <<  __func__
                 << " "
                 << theMaterial->GetName()
                 << ", has free electron density of "
                 << theMaterial->GetFreeElectronDensity()
                 << cond
                 << G4endl;
        }
        if (std::find(conductors.begin(), conductors.end(), theMaterial->GetName())
            != conductors.end() ) {
          G4NistManager::Instance()->SetDensityEffectCalculatorFlag(theMaterial, true);
          if (config_.debug().diagLevel() > 0) {
            G4cout << "ConstructMaterials::" <<  __func__
                   << " Using correction in calculations for "
                   << theMaterial->GetName()
                   << ", its free electron density is "
                   << theMaterial->GetFreeElectronDensity()
                   << G4endl;

          }
        }
      }
    }
  }

} // end namespace mu2e
