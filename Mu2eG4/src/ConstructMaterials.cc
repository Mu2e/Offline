//
// Construct materials requested by the run-time configuration system.
//
// $Id: ConstructMaterials.cc,v 1.35 2012/09/17 17:00:50 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/09/17 17:00:50 $
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
#include <iomanip>

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

    CheckedG4String mat = isNeeded(materialsToLoad, "HeavyConcrete");
    if ( mat.doit ) {
      G4Material* HeavyConcrete = new G4Material(mat.name, 3.295*g/cm3, 17);
      G4Element* eH  = getElementOrThrow("H");
      G4Element* eB  = getElementOrThrow("B");
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eO  = getElementOrThrow("O");
      G4Element* eF  = getElementOrThrow("F");
      G4Element* eNa = getElementOrThrow("Na");
      G4Element* eMg = getElementOrThrow("Mg");
      G4Element* eAl = getElementOrThrow("Al");
      G4Element* eSi = getElementOrThrow("Si");
      G4Element* eP  = getElementOrThrow("P");
      G4Element* eS  = getElementOrThrow("S");
      G4Element* eK  = getElementOrThrow("K");
      G4Element* eCa = getElementOrThrow("Ca");
      G4Element* eTi = getElementOrThrow("Ti");
      G4Element* eMn = getElementOrThrow("Mn");
      G4Element* eFe = getElementOrThrow("Fe");
      G4Element* eSr = getElementOrThrow("Sr");
      HeavyConcrete->AddElement( eH , 0.01048482); //Hydrogen
      HeavyConcrete->AddElement( eB , 0.00943758); //Boron
      HeavyConcrete->AddElement( eC , 0.0129742);  //Carbon
      HeavyConcrete->AddElement( eO , 0.27953541); //Oxygen
      HeavyConcrete->AddElement( eF , 1.5175E-4);  //Fluorine
      HeavyConcrete->AddElement( eNa, 3.7014E-4);  //Sodium
      HeavyConcrete->AddElement( eMg, 0.08298213); //Magnesium
      HeavyConcrete->AddElement( eAl, 0.02769028); //Aluminum
      HeavyConcrete->AddElement( eSi, 0.06317253); //Silicon
      HeavyConcrete->AddElement( eP , 0.00176963); //Phosphorus
      HeavyConcrete->AddElement( eS , 5.8275E-4);  //Sulfur
      HeavyConcrete->AddElement( eK , 4.2024E-4);  //Potassium
      HeavyConcrete->AddElement( eCa, 0.03227609); //Calcium
      HeavyConcrete->AddElement( eTi, 5.457E-5);   //Titanium
      HeavyConcrete->AddElement( eMn, 0.00321757); //Manganese
      HeavyConcrete->AddElement( eFe, 0.47423935); //Iron
      HeavyConcrete->AddElement( eSr, 6.4097E-4);  //Strontium
    }

    mat = isNeeded(materialsToLoad, "ShieldingConcrete");
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
    // FIXME is there a better reference?
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

    // Construction Aluminum
    //http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA5083O
    //http://ppd-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=1112;filename=MD-ENG-109.pdf;version=1
    mat = isNeeded(materialsToLoad, "A95083");
    if ( mat.doit ){
      G4Material* A95083 =
        new G4Material( mat.name, 2.66*g/cm3, 9);
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


    // NbTi
    mat = isNeeded(materialsToLoad, "NbTi"); // FIXME verify it
    if ( mat.doit ){
      G4Material* NbTi =
        new G4Material( mat.name, 6.5*g/cm3, 2);
      NbTi->AddMaterial(findMaterialOrThrow("G4_Nb"), 0.65);
      NbTi->AddMaterial(findMaterialOrThrow("G4_Ti"), 0.35);
    }

    // NbTiCu
    mat = isNeeded(materialsToLoad, "NbTiCu"); // FIXME verify it
    if ( mat.doit ){
      G4Material* NbTiCu =
        new G4Material( mat.name, 7.69*g/cm3, 2);
      NbTiCu->AddMaterial(findMaterialOrThrow("NbTi"),  0.45);
      NbTiCu->AddMaterial(findMaterialOrThrow("G4_Cu"), 0.55);
    }

    // AL999Ni001 by volume ?
    mat = isNeeded(materialsToLoad, "AL999Ni001"); // FIXME verify it
    if ( mat.doit ){
      G4Material* AL999Ni001 =
        new G4Material( mat.name, 2.706*g/cm3, 3);
      AL999Ni001->AddMaterial(findMaterialOrThrow("G4_Al"), 0.9967);
      AL999Ni001->AddMaterial(findMaterialOrThrow("G4_Ni"), 0.0033);
    }

    // http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/Epotek-301-1.html
    // C_19_H_20_O_4

    mat = isNeeded(materialsToLoad, "C_19_H_20_O_4");
    if ( mat.doit ){
      G4Material* C_19_H_20_O_4 =
        new G4Material( mat.name, 1.16*g/cm3, 3);

      G4Element* eC  = getElementOrThrow("C");
      G4Element* eH  = getElementOrThrow("H");
      G4Element* eO  = getElementOrThrow("O");

      C_19_H_20_O_4->AddElement( eC, 19);
      C_19_H_20_O_4->AddElement( eH, 20);
      C_19_H_20_O_4->AddElement( eO,  4);

    }

    // C_10_H_18_O_4

    mat = isNeeded(materialsToLoad, "C_10_H_18_O_4");
    if ( mat.doit ){
      G4Material* C_10_H_18_O_4 =
        new G4Material( mat.name, 1.10*g/cm3, 3);

      G4Element* eC  = getElementOrThrow("C");
      G4Element* eH  = getElementOrThrow("H");
      G4Element* eO  = getElementOrThrow("O");

      C_10_H_18_O_4->AddElement( eC, 10);
      C_10_H_18_O_4->AddElement( eH, 18);
      C_10_H_18_O_4->AddElement( eO,  4);

    }

    // C_9_H_22_N_2

    mat = isNeeded(materialsToLoad, "C_9_H_22_N_2");
    if ( mat.doit ){
      G4Material* C_9_H_22_N_2 =
        new G4Material( mat.name, 0.865*g/cm3, 3);

      G4Element* eC  = getElementOrThrow("C");
      G4Element* eH  = getElementOrThrow("H");
      G4Element* eN  = getElementOrThrow("N");

      C_9_H_22_N_2->AddElement( eC,  9);
      C_9_H_22_N_2->AddElement( eH, 22);
      C_9_H_22_N_2->AddElement( eN,  2);

    }

    // http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/Epotek-301-1.html
    mat = isNeeded(materialsToLoad, "Epotek301");
    if ( mat.doit ){
      G4Material* Epotek301 =
        new G4Material( mat.name, 1.19*g/cm3, 3);
      Epotek301->AddMaterial(findMaterialOrThrow("C_19_H_20_O_4"), 0.56);
      Epotek301->AddMaterial(findMaterialOrThrow("C_10_H_18_O_4"), 0.24);
      Epotek301->AddMaterial(findMaterialOrThrow("C_9_H_22_N_2"),  0.20);
    }

   // http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/E_glass.html
    mat = isNeeded(materialsToLoad, "EGlass");
    if ( mat.doit ){
      G4Material* EGlass =
	new G4Material (mat.name, 2.61*g/cm3, 10);
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
    mat = isNeeded(materialsToLoad, "G10");
    if ( mat.doit ){
      G4Material* G10 =
        new G4Material( mat.name, 1.7*g/cm3, 2);
      G10->AddMaterial(findMaterialOrThrow("G4_SILICON_DIOXIDE"), 0.60);//FIXME do e-glass etc...
      G10->AddMaterial(findMaterialOrThrow("Epotek301"), 0.40);
    }

    // Superconducting Cable Insulation
    mat = isNeeded(materialsToLoad, "SCCableInsulation");
    if ( mat.doit ){
      G4Material* SCCableInsulation =
        new G4Material( mat.name, 1.54*g/cm3, 3);
      SCCableInsulation->AddMaterial(findMaterialOrThrow("G4_KAPTON"), 0.18);
      SCCableInsulation->AddMaterial(findMaterialOrThrow("Epotek301"), 0.16);
      SCCableInsulation->AddMaterial(findMaterialOrThrow("G10"), 0.66);
    }

    //http://sector7.xor.aps.anl.gov/~dufresne/UofM/techinfo/kapton.html
    //The chemical formula of Kapton is C22H10N205,  its density is 1.43
    // also see below, do we need more than one?

    // Superconducting Cable
    mat = isNeeded(materialsToLoad, "SCCable"); // FIXME verify it
    if ( mat.doit ){
      G4Material* SCCable =
        new G4Material( mat.name, 3.95*g/cm3, 3);
      SCCable->AddMaterial(findMaterialOrThrow("SCCableInsulation"), 0.04);
      SCCable->AddMaterial(findMaterialOrThrow("AL999Ni001"),        0.43);
      SCCable->AddMaterial(findMaterialOrThrow("NbTiCu"),            0.53);
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

    mat = isNeeded(materialsToLoad, "StrawGasArCF4");

    if ( mat.doit ) {
      
      G4Material* StrawGasArCF4 =
        new G4Material(mat.name, 0.0028561*g/cm3, 3); // it is OK not to use kStateGas
      G4Element* eAr = getElementOrThrow("Ar");
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eF  = getElementOrThrow("F");
      StrawGasArCF4->AddElement( eAr, 1);
      StrawGasArCF4->AddElement( eC,  1);
      StrawGasArCF4->AddElement( eF,  4);

    }
    
    mat = isNeeded(materialsToLoad, "StrawGas");
    if ( mat.doit ) {
     
      G4double density;
      G4double temperature = 293.15*kelvin;
      G4double pressure = 1*atmosphere;
      G4int nel;

      G4double densityAr   = 0.00166 *g/cm3; //from PDG
      G4double densityCO2  = 0.00184 *g/cm3; //from PDG
      G4double fractionAr  = 80.0*perCent;

      density = fractionAr*densityAr + (1.0-fractionAr)*densityCO2;

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature, pressure);

      G4Element* Ar = getElementOrThrow("Ar");
      G4Element* C  = getElementOrThrow("C");
      G4Element* O  = getElementOrThrow("O");

      G4double atomicWeight_Ar =  39.948  *g/mole;
      G4double atomicWeight_C  = 12.0107    *g/mole;
      G4double atomicWeight_O  = 15.9994 *g/mole;
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
      // This is the lowest density vacuum allowed by G4.
      G4double density     = universe_mean_density;
      G4double pressure    = 3.e-18*pascal;
      G4double temperature = 2.73*kelvin;

      // G4 takes ownership of this object and manages its lifetime.
      new G4Material( mat.name, 1., 1.01 *g/mole,
                      density, kStateGas, temperature, pressure);
    }


    // Presume that the residual gas in the DS will be leakage from the straws,
    // pumped down to 10^{-4} torr.
    mat = isNeeded(materialsToLoad, "DSVacuum");
    if ( mat.doit ){

      G4Material* StrawLeak = findMaterialOrThrow("StrawGas");

      G4double temperature = 300.00*kelvin; // Temperature of the DS
      G4double pressure    = 133.322e-4*pascal; // 10e-4 Torr
      G4double refTemp     = StrawLeak->GetTemperature();
      G4double refPress    = StrawLeak->GetPressure();

      G4double density = StrawLeak->GetDensity()*pressure*refTemp/(refPress*temperature);

      G4Material* DSVacuum =
	new G4Material(mat.name, density, StrawLeak->GetNumberOfElements(),
		       kStateGas, temperature, pressure);

      for (size_t i = 0 ; i < StrawLeak->GetNumberOfElements(); ++i) {
	DSVacuum->AddElement(StrawLeak->GetElementVector()->at(i), StrawLeak->GetFractionVector()[i]);
      }

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

    mat = isNeeded(materialsToLoad, "ITGasHe_95Isob_5");
    if ( mat.doit ){

      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.000166 *g/cm3;
      G4double densityIsoB = 0.00249  *g/cm3;
      G4double fractionHe  = 95.0*perCent;

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

    mat = isNeeded(materialsToLoad, "ITGasHe_75Isob_25_400mbar");
    if ( mat.doit ){

      G4double density, temperature, pressure;
      G4int nel;

      G4double densityHe   = 0.000166 *g/cm3;
      G4double densityIsoB = 0.00249  *g/cm3;
      G4double fractionHe  = 75.0*perCent;

      density = fractionHe*densityHe + (1.0-fractionHe)*densityIsoB;
      pressure = 0.4*bar;
      density *= pressure/(1.0*atmosphere);

      G4Material *GasMix = new G4Material( mat.name, density, nel=3,
                                           kStateGas, temperature= 293.15*kelvin, pressure);

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


    //G10-FR4 used for printed board of the I-Tracker
    // G10 http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/SDD_G10FR4.html
    // http://pdg.lbl.gov/2002/atomicrpp.pdf
    mat = isNeeded(materialsToLoad, "G10_FR4");
    if ( mat.doit ) {
      G4double density;

      G4Material* G10_FR4 =
        new G4Material(mat.name, density = 1.8*g/cm3, 2);
      G10_FR4->AddMaterial(findMaterialOrThrow("EGlass"), 0.60);
      G10_FR4->AddMaterial(findMaterialOrThrow("Epotek301"), 0.40);
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

    mat = isNeeded(materialsToLoad, "CFoam_332");
    if ( mat.doit ){
      G4double density;
      G4int nel;
      G4Material *CFoam = new G4Material(mat.name, density = 0.332*g/cm3, nel=1);
      G4Element* C  = getElementOrThrow("C");
      CFoam->AddElement(C, 100.0*perCent );
    }

    mat = isNeeded(materialsToLoad, "CFoam_166");
    if ( mat.doit ){
      G4double density;
      G4int nel;
      G4Material *CFoam = new G4Material(mat.name, density = 0.166*g/cm3, nel=1);
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

    mat = isNeeded(materialsToLoad, "ZirconiumHydridePolyethylene");
    if ( mat.doit ){
      G4Material* ZirconiumHydridePolyethylene =
        new G4Material( mat.name, 3.67*g/cm3, 5);
      G4Element* eC  = getElementOrThrow("C");
      G4Element* eH  = getElementOrThrow("H");
      G4Element* eB  = getElementOrThrow("B");
      G4Element* eO  = getElementOrThrow("O");
      G4Element* eZr  = getElementOrThrow("Zr");

      ZirconiumHydridePolyethylene->AddElement( eC, 8.9*perCent);
      ZirconiumHydridePolyethylene->AddElement( eH, 3.4*perCent);
      ZirconiumHydridePolyethylene->AddElement( eB, 0.5*perCent);
      ZirconiumHydridePolyethylene->AddElement( eO, 2.2*perCent);
      ZirconiumHydridePolyethylene->AddElement( eZr, 85.0*perCent);
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
