//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.cc,v 1.22 2013/10/25 21:37:29 genser Exp $
// $Author: genser $
// $Date: 2013/10/25 21:37:29 $
//
// Original author Peter Shanahan
//
// Notes:


// C++ includes
#include <iostream>
#include <string>

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/constructStoppingTarget.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

// G4 includes
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolumeStore.hh"

using namespace std;

namespace mu2e {
    VolumeInfo constructStoppingTarget( VolumeInfo   const& parent,
                                      SimpleConfig const& config ){
 //TODO - add this param to the ST txt
      string design_name = config.getString("stoppingTarget.designName", "foil");
      std::cout << "design Name : " << design_name << std::endl;
      if(design_name == "screen"){ 
        return constructStoppingTarget_screen(parent, config);
      }else if(design_name == "cylinder"){
        return constructStoppingTarget_cylinder(parent, config);
      }else if(design_name == "hexagon"){
        return constructStoppingTarget_hexagon(parent, config);
      }else{
        return constructStoppingTarget_foil(parent, config);
      }
    }
    /*
      Implementation of foil related stopping target designs, including
      1. foilhole (37-folder foil with 21mm inner radius hole).
    */
    VolumeInfo constructStoppingTarget_foil( VolumeInfo   const& parent,
                                      SimpleConfig const& config ){
   
    
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "stoppingTarget", "stoppingTarget");

    const bool stoppingTargetIsVisible = geomOptions->isVisible("stoppingTarget"); 
    const bool stoppingTargetIsSolid   = geomOptions->isSolid("stoppingTarget"); 
    const bool forceAuxEdgeVisible     = geomOptions->forceAuxEdgeVisible("stoppingTarget"); 
    const bool doSurfaceCheck          = geomOptions->doSurfaceCheck("stoppingTarget"); 
    const bool placePV                 = geomOptions->placePV("stoppingTarget"); 


    int verbosity(config.getInt("stoppingTarget.verbosity",0));

    if ( verbosity > 1 ) std::cout << "In constructStoppingTarget" << std::endl;
    // Master geometry for the Target assembly
    GeomHandle<StoppingTarget> target;

    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    TubsParams targetMotherParams(0., target->cylinderRadius(), target->cylinderLength()/2.);
    //if ( verbosity > 1 ) {
      std::cout<<"Foil ST with length "<<target->cylinderLength()<<" radius "<<target->cylinderRadius()<<std::endl;
      std::cout<<"Target Center is "<<target->centerInMu2e()<<std::endl;
    
    //}
    VolumeInfo targetInfo = nestTubs("StoppingTargetMother",
                                     targetMotherParams,
                                     findMaterialOrThrow(target->fillMaterial()),
                                     0,
                                     target->centerInMu2e() - parent.centerInMu2e(),
                                     parent,
                                     0,
                                     false/*visible*/,
                                     G4Colour::Black(),
                                     false/*solid*/,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    // now create the individual target foils

    G4VPhysicalVolume* pv;

    for (int itf=0; itf<target->nFoils(); ++itf) {
        
        TargetFoil foil=target->foil(itf);
        if(itf==0) cout<< foil.centerInMu2e()<<std::endl;
        if(itf==target->nFoils()-1) cout<< foil.centerInMu2e()<<std::endl;
        VolumeInfo foilInfo;
        G4Material* foilMaterial = findMaterialOrThrow(foil.material());

        std::ostringstream os;
        os << std::setfill('0') << std::setw(2) << itf;
        foilInfo.name = "Foil_" + os.str();

        if ( verbosity > 0 )  std::cout << __func__ << " " << foilInfo.name << std::endl;

        foilInfo.solid = new G4Tubs(foilInfo.name
                                    ,foil.rIn()
                                    ,foil.rOut()
                                    ,foil.halfThickness()
                                    ,0.
                                    ,CLHEP::twopi
                                    );

        foilInfo.logical = new G4LogicalVolume( foilInfo.solid
                                                , foilMaterial
                                                , foilInfo.name
                                                );
        

        // rotation matrix...
        G4RotationMatrix* rot = 0; //... will have to wait

        G4ThreeVector foilOffset(foil.centerInMu2e() - targetInfo.centerInMu2e());
        if ( verbosity > 1 ) std::cout << "foil " 
                                  << itf
                                  << " centerInMu2e=" 
                                  << foil.centerInMu2e()
                                  << ", offset="<< foilOffset<< std::endl;

        // G4 manages the lifetime of this object.
        pv = new G4PVPlacement( rot,
                                foilOffset,
                                foilInfo.logical,
                                "Target"+foilInfo.name,
                                targetInfo.logical,
                                0,
                                itf,
                                false);

        doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);


double totMass  = 0.0;
    double foil_rIn = config.getDouble("stoppingTarget.inner_radius", 0);
    
    for (int itf=0; itf<target->nFoils(); ++itf)
    {

        TargetFoil foil=target->foil(itf);


        totMass += M_PI * (foil.rOut()*foil.rOut() - foil_rIn*foil_rIn) * 2 * foil.halfThickness();
        
     }
      totMass *= 2.7 * 0.001;
      std::cout<<"Mass of target = "<<totMass<<std::endl;
        if (!stoppingTargetIsVisible) {
          foilInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
        } else {
          G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
          visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
          visAtt->SetForceSolid(stoppingTargetIsSolid);
          foilInfo.logical->SetVisAttributes(visAtt);
        }
      }// target foils

    for (int itf=0; itf<target->nSupportStructures(); ++itf) {
	TargetFoilSupportStructure supportStructure=target->supportStructure(itf);
      double stoppingTarget_z0InMu2e = config.getDouble("stoppingTarget.z0InMu2e", 5871.);
      cout << "-------- Design Details --------" << endl;
      cout << "target start : "<< stoppingTarget_z0InMu2e<< endl;
      cout << "target length : "<< target->cylinderLength() << endl;
      cout << "target center : "<<target->centerInMu2e() << endl;
      cout<< "parent center : " << parent.centerInMu2e() << endl;
        VolumeInfo supportStructureInfo;

        G4Material* supportStructureMaterial = findMaterialOrThrow(supportStructure.material());

        std::ostringstream os;
        os << std::setfill('0') << std::setw(2) << itf;
        supportStructureInfo.name = "FoilSupportStructure_" + os.str();

        if ( verbosity > 0 )  std::cout << __func__ << " " << supportStructureInfo.name 
                                          << std::endl;

        supportStructureInfo.solid = new G4Tubs(supportStructureInfo.name
                                    ,0
                                    ,supportStructure.radius()
                                    ,supportStructure.length()/2. // G4Tubs useses half-lengths as this parameter
                                    ,0.
                                    ,CLHEP::twopi
                                    );

        supportStructureInfo.logical = new G4LogicalVolume( supportStructureInfo.solid
                                                , supportStructureMaterial
                                                , supportStructureInfo.name
                                                );
        

	if ( verbosity > 1 ) std::cout << "supportStructure.support_id() = "
                                       << supportStructure.support_id() 
                                       << "    target->nSupportStructures() = " 
                                       << target->nSupportStructures() 
                                       << "     target->nFoils() = " 
                                       << target->nFoils() 
                                       << "     supportStructure.length() = " 
                                       << supportStructure.length() 
                                       << std::endl;

        // rotation matrices to rotate the orientation of the
        // supporting wires. First rotate into xy-plane by 90deg
        // rotation around y-axis, then rotate within xy-plane by
        // appropiate rotation around z-axis

        CLHEP::HepRotationY secRy(-M_PI/2.);
        CLHEP::HepRotationZ secRz( -supportStructure.support_id() * 360.*CLHEP::deg / 
                                   (target->nSupportStructures()/target->nFoils()) 
                                   - 90.*CLHEP::deg - supportStructure.angleOffset()*CLHEP::deg);
        G4RotationMatrix* supportStructure_rotMatrix = reg.add(G4RotationMatrix(secRy*secRz));

	if ( verbosity > 1 ) std::cout << "supportStructure_rotMatrix = " 
                                       << *supportStructure_rotMatrix << std::endl;

        // vector where to place to support tube
        // first find target center
        G4ThreeVector supportStructureOffset(supportStructure.centerInMu2e() - targetInfo.centerInMu2e()); 

        if ( verbosity > 1 ) std::cout << supportStructureInfo.name << " "
                                  << itf 
                                  << " centerInMu2e="
                                  << supportStructure.centerInMu2e()
                                  << ", offset=" 
                                  << supportStructureOffset
                                  << std::endl;

        if ( verbosity > 1 ) std::cout << __func__ << " "
                                       << supportStructureInfo.name
                                       <<  std::endl;

        G4ThreeVector 
          vector_supportStructure_Orientation( (supportStructure.length()/2.+supportStructure.foil_outer_radius()) * 
                                               std::cos(supportStructure.support_id() * 360.*CLHEP::deg / 
                                                        (target->nSupportStructures()/target->nFoils()) + 
                                                        90.*CLHEP::deg + 
                                                        supportStructure.angleOffset()*CLHEP::deg), 
                                               (supportStructure.length()/2.+supportStructure.foil_outer_radius()) * 
                                               std::sin(supportStructure.support_id() * 360.*CLHEP::deg / 
                                                        (target->nSupportStructures()/target->nFoils()) + 
                                                        90.*CLHEP::deg + supportStructure.angleOffset()*CLHEP::deg), 0);
        
	if ( verbosity > 1 ) std::cout << "vector_supportStructure_Orientation = " 
                                       << vector_supportStructure_Orientation << std::endl;

	supportStructureOffset += vector_supportStructure_Orientation; // second add vector to support wire tube center

	if ( verbosity > 1 ) std::cout << "supportStructureOffset += vector_supportStructure_Orientation = " 
                                       << supportStructureOffset << std::endl;

        pv = new G4PVPlacement( supportStructure_rotMatrix,
                                supportStructureOffset,
                                supportStructureInfo.logical,
                                supportStructureInfo.name,
                                targetInfo.logical,
                                0,
                                itf,
                                false);

        doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

        if (!stoppingTargetIsVisible) {
          supportStructureInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
        } else {
          G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Blue()));
          visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
          visAtt->SetForceSolid(stoppingTargetIsSolid);
          supportStructureInfo.logical->SetVisAttributes(visAtt);
        }
      }// target foils support structures

    return targetInfo;
  }
  
  /*
  Implementation of screen related stopping target designs, including
  1. screen_default: .0021'' diam of wire, .0029 width of opening
  2. screen_mesh: change the mesh of screen_default
  */
VolumeInfo constructStoppingTarget_screen( VolumeInfo   const& parent,
                                      SimpleConfig const& config ){
    std::cout << "Debug: Enter 1 " << std::endl;
    const bool forceAuxEdgeVisible  = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV              = true;

    int verbosity(config.getInt("target.verbosity",0));

    if ( verbosity > 1 ) std::cout<<"In constructStoppingTarget"<<std::endl;
    std::cout << "Debug: Enter 2 " << std::endl;
    // Master geometry for the Target assembly
    GeomHandle<StoppingTarget> target;

    G4VSensitiveDetector* stSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::StoppingTarget());

    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    //    TubsParams targetMotherParams(0., target->cylinderRadius(), target->cylinderLength()/2.); //MR20150616 replace by hardcoded tube for primary volume
    TubsParams targetMotherParams(0., 300., 1000.);//MR20150616 hardcoded tube for primary volume

    VolumeInfo targetInfo = nestTubs("StoppingTargetMother",
                                     targetMotherParams,
                                     findMaterialOrThrow(target->fillMaterial()),
                                     0,
                                     target->centerInMu2e() - parent.centerInMu2e(),
                                     parent,
                                     0,
                                     false/*visible*/,
                                     G4Colour::Black(),
                                     false/*solid*/,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );
    G4VPhysicalVolume* pv;

    ////////////////////////////////////////////////
    // START create screen
    ////////////////////////////////////////////////
    double outer_radius = config.getDouble("stoppingTarget.outer_radius");
    double inner_radius = config.getDouble("stoppingTarget.inner_radius");
    int number_screenlayers = config.getInt("stoppingTarget.number_screenlayers");
    double string_radius = config.getDouble("stoppingTarget.string_radius");
    double string_target_opening_size = config.getDouble("stoppingTarget.string_target_opening_size");
    double stoppingTarget_z0InMu2e = config.getDouble("stoppingTarget.z0InMu2e", 5871.);
    double target_length = config.getDouble("stoppingTarget.target_length", 800.);

//CBE All numbers should be configurable, or extracted from other configurable parameters
    double screenlayers_deltaz = target_length / (number_screenlayers - 1); // spacing of screen layers in the stopping target volume
    double string_spacing_dx = string_target_opening_size + 2 * string_radius;
    int Nstrings_per_layer = (int)round(outer_radius * 2 / string_spacing_dx);
    double total_target_mass(0);
    double total_target_mass_GEANT4(0);
    std::cout << "============================================================================" << std::endl;
    std::cout << "START Target info" << std::endl;
    std::cout << "target->nFoils() == " << target->nFoils() << std::endl;
    std::cout << "number_screenlayers == " << number_screenlayers << std::endl;
    std::cout << "screenlayers_deltaz == " << screenlayers_deltaz << std::endl;
    std::cout << "Nstrings_per_layer == " << Nstrings_per_layer << std::endl;
    std::cout << "stoppingTarget_z0InMu2e == " << stoppingTarget_z0InMu2e << std::endl;
    std::cout << "string_radius == " << string_radius << std::endl;
    std::cout << "string_target_opening_size == " << string_target_opening_size << std::endl;
    std::cout << "string_spacing_dx == " << string_spacing_dx << std::endl;
    std::cout << "STOP Target info" << std::endl;
    std::cout << "============================================================================" << std::endl;
    //MR20160824 START
    VolumeInfo stringInfoLeft, stringInfoRight;
    G4Material* stringMaterial = findMaterialOrThrow("StoppingTarget_Al");
    stringInfoLeft.name = "String";
    stringInfoRight.name = "String";
    // rotation matrix...
    // G4RotationMatrix* rot = 0;
    CLHEP::HepRotationZ secRz(-M_PI/2.); // 90 degree rotation around z-axis
    CLHEP::HepRotationX secRx(-M_PI/2.); // 90 degree rotation around x-axis

    G4RotationMatrix* rotation_matrix_strings_parallel_xaxis = new G4RotationMatrix(secRx * secRz);
    G4RotationMatrix* rotation_matrix_strings_parallel_yaxis = new G4RotationMatrix(secRx);

    // now create the individual strings per screen layer
    for (int counter_stringpair=0; counter_stringpair < Nstrings_per_layer/2; counter_stringpair++) {
        double screenlayer_mass=0;
        // compute halfLength of string as a function of the radius at this position
        double center_distance = string_spacing_dx / 2. + counter_stringpair * string_spacing_dx;
        double string_halfLength = sqrt( pow(outer_radius, 2) - pow(center_distance, 2) );
        if(center_distance < inner_radius)
            string_halfLength -= sqrt( pow(inner_radius, 2) - pow(center_distance, 2) );
        // create string wire with halfLength according to the position in the disc
        
//CBE Do you need to create two? Given the symmetry of the problem, a single volume should suffice, but it will be placed at different positions.
// I note the two volumes are identical!
	
	      stringInfoLeft.solid = new G4Tubs(stringInfoLeft.name
                ,0. // inner radius
                ,string_radius // outer radius
                ,string_halfLength/2.0 // half length of the string dependent on position
                ,0.
                ,CLHEP::twopi
                );
        
        stringInfoRight.solid = new G4Tubs(stringInfoRight.name
                ,0. // inner radius
                ,string_radius // outer radius
                ,string_halfLength/2.0 // half length of the string dependent on position
                ,0.
                ,CLHEP::twopi
                );
        
        // create logic volume for the string
        stringInfoLeft.logical = new G4LogicalVolume(stringInfoLeft.solid
                ,stringMaterial
                ,stringInfoLeft.name
                );
        if(stSD) stringInfoLeft.logical->SetSensitiveDetector(stSD);
        
        stringInfoRight.logical = new G4LogicalVolume(stringInfoRight.solid
                ,stringMaterial
                ,stringInfoRight.name
                );
        if(stSD) stringInfoRight.logical->SetSensitiveDetector(stSD);


        // loop over screen layers
        for (int counter_screenlayer=0; counter_screenlayer < number_screenlayers; counter_screenlayer++) {
            total_target_mass_GEANT4 += 4 * ( stringInfoLeft.logical->GetMass() + stringInfoRight.logical->GetMass() ) / CLHEP::Avogadro;
            // create strings pairwise at +y and -y (strings parallel to x-axis), and  +x and -x (strings parallel to y-axis)
            double delta_x = string_halfLength / 2.0;
            if(center_distance < inner_radius) delta_x += sqrt( pow(inner_radius, 2) - pow(center_distance, 2) );
            double delta_y = center_distance;
            double delta_z = ( (number_screenlayers-1)*screenlayers_deltaz/2 ) + screenlayers_deltaz + (int(counter_screenlayer)-number_screenlayers)*screenlayers_deltaz;

//CBE once again, no "magic numbers" should be there
            
            G4ThreeVector stringRight_x_n_offset(+delta_x, -delta_y, delta_z-string_radius);
            G4ThreeVector stringRight_x_p_offset(+delta_x, +delta_y, delta_z-string_radius);
            G4ThreeVector stringRight_y_n_offset(-delta_y, +delta_x, delta_z+string_radius);
            G4ThreeVector stringRight_y_p_offset(+delta_y, +delta_x, delta_z+string_radius);

            G4ThreeVector stringLeft_x_n_offset(-delta_x, -delta_y, delta_z-string_radius);
            G4ThreeVector stringLeft_x_p_offset(-delta_x, delta_y, delta_z-string_radius);
            G4ThreeVector stringLeft_y_n_offset(-delta_y, -delta_x, delta_z+string_radius);
            G4ThreeVector stringLeft_y_p_offset(+delta_y, -delta_x, delta_z+string_radius);


//CBE do we really need to use boost_lexical here? Just declare a bunch of strings and that should do it.	    
            pv = new G4PVPlacement( rotation_matrix_strings_parallel_xaxis,
                    stringRight_x_n_offset,
                    stringInfoRight.logical,
                    ("TargetString_stringpair_Right_" + boost::lexical_cast<std::string>(counter_stringpair) + "_negativeside_parallel_xaxis").c_str(),
                    targetInfo.logical,
                    0,
                    counter_stringpair,
                    false);
            doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

            pv = new G4PVPlacement( rotation_matrix_strings_parallel_xaxis,
                    stringLeft_x_n_offset,
                    stringInfoLeft.logical,
                    ("TargetString_stringpair_Left_" + boost::lexical_cast<std::string>(counter_stringpair) + "_negativeside_parallel_xaxis").c_str(),
                    targetInfo.logical,
                    0,
                    counter_stringpair,
                    false);
            doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

            // place positive side string at +x for layer parallel to xaxis (downstream)

            pv = new G4PVPlacement( rotation_matrix_strings_parallel_xaxis,
                    stringRight_x_p_offset,
                    stringInfoRight.logical,
                    ("TargetString_stringpair_Right_" + boost::lexical_cast<std::string>(counter_stringpair) + "_negativeside_parallel_xaxis").c_str(),
                    targetInfo.logical,
                    0,
                    counter_stringpair,
                    false);
            doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

            pv = new G4PVPlacement( rotation_matrix_strings_parallel_xaxis,
                    stringLeft_x_p_offset,
                    stringInfoLeft.logical,
                    ("TargetString_stringpair_Left_" + boost::lexical_cast<std::string>(counter_stringpair) + "_negativeside_parallel_xaxis").c_str(),
                    targetInfo.logical,
                    0,
                    counter_stringpair,
                    false);
            doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

            // place negative side string at -x for layer parallel to yaxis (upstream)

            pv = new G4PVPlacement( rotation_matrix_strings_parallel_yaxis,
                    stringRight_y_n_offset,
                    stringInfoRight.logical,
                    ("TargetString_stringpair_Right_" + boost::lexical_cast<std::string>(counter_stringpair) + "_negativeside_parallel_xaxis").c_str(),
                    targetInfo.logical,
                    0,
                    counter_stringpair,
                    false);
            doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);


            pv = new G4PVPlacement( rotation_matrix_strings_parallel_yaxis,
                    stringLeft_y_n_offset,
                    stringInfoLeft.logical,
                    ("TargetString_stringpair_Left_" + boost::lexical_cast<std::string>(counter_stringpair) + "_negativeside_parallel_xaxis").c_str(),
                    targetInfo.logical,
                    0,
                    counter_stringpair,
                    false);
            doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

            // place positive side string at +x for layer parallel to yaxis (upstream)

            pv = new G4PVPlacement( rotation_matrix_strings_parallel_yaxis,
                    stringRight_y_p_offset,
                    stringInfoRight.logical,
                    ("TargetString_stringpair_Right_" + boost::lexical_cast<std::string>(counter_stringpair) + "_negativeside_parallel_xaxis").c_str(),
                    targetInfo.logical,
                    0,
                    counter_stringpair,
                    false);
            doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);


            pv = new G4PVPlacement( rotation_matrix_strings_parallel_yaxis,
                    stringLeft_y_p_offset,
                    stringInfoLeft.logical,
                    ("TargetString_stringpair_Left_" + boost::lexical_cast<std::string>(counter_stringpair) + "_negativeside_parallel_xaxis").c_str(),
                    targetInfo.logical,
                    0,
                    counter_stringpair,
                    false);
            doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

            // compute mass of this screenlayer. The density of aluminum is 2.7 g/cm^3
            screenlayer_mass += (4 * M_PI * string_radius * string_radius * 2 * string_halfLength) * (2.7 / 1000.) ; // four strings places in this screen layer


            if (!config.getBool("target.visible",true)) {
                stringInfoLeft.logical->SetVisAttributes(G4VisAttributes::Invisible);
                stringInfoRight.logical->SetVisAttributes(G4VisAttributes::Invisible);
            } else {
                G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
                visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
                visAtt->SetForceSolid(config.getBool("target.solid",true));
                stringInfoLeft.logical->SetVisAttributes(visAtt);
                stringInfoRight.logical->SetVisAttributes(visAtt);
            }

            if (counter_stringpair==0) std::cout << "Created screen layer " << counter_screenlayer << std::endl;
        }
        std::cout << "stopping target screenlayer_mass == " << screenlayer_mass << " gram" << std::endl;
        total_target_mass+=screenlayer_mass;
    }
    std::cout << "stopping target total mass computed by GEANT4 == " << total_target_mass_GEANT4 << std::endl;
    std::cout << "stopping target total mass computed manually == " << total_target_mass << " gram" << std::endl;
       return targetInfo;
  }
  
  /*
  Implementation of cylinder related designs, including
  1. cylinder_default: 7 cylinders
  2. cylinder_mesh: change mesh and number of cylinders of cylinder_default
  */
  VolumeInfo constructStoppingTarget_cylinder( VolumeInfo const& parent,
                                      SimpleConfig const& config){

    const bool forceAuxEdgeVisible  = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV              = true;

    int verbosity(config.getInt("target.verbosity",10));

    if ( verbosity > 1 ) std::cout<<"In constructStoppingTarget"<<std::endl;
    // Master geometry for the Target assembly
    GeomHandle<StoppingTarget> target;

    G4VSensitiveDetector* stSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::StoppingTarget()); //TODO

    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    //TubsParams targetMotherParams(0., target->cylinderRadius(), target->cylinderLength()/2.); //SM2020 added back in//MR20150616 replace by hardcoded tube for primary volume
   //CBE Lost, no hardcoded numbers allowed! in the whole function!
   TubsParams targetMotherParams(0., 300., 1000.); //MR20150616 hardcoded tube for primary volume
   if ( verbosity > 1 ) {
      std::cout<<"Cylinder ST with length "<<target->cylinderLength()<<" radius "<<target->cylinderRadius()<<std::endl;
      std::cout<<"Target Center is "<<target->centerInMu2e()<<std::endl;
    
    }
    VolumeInfo targetInfo = nestTubs("StoppingTargetMother",
                                     targetMotherParams,
                                     findMaterialOrThrow(target->fillMaterial()),
                                     0,
                                     target->centerInMu2e() - parent.centerInMu2e(),
                                     parent,
                                     0,
                                     false/*visible*/,
                                     G4Colour::Black(),
                                     false/*solid*/,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    G4VPhysicalVolume* pv;

    ////////////////////////////////////////////////
    // START create target foils
    ////////////////////////////////////////////////

    /* Structure parameters */
    double wire_diameter = config.getDouble("stoppingTarget.string_radius", 0.02665) * 2.0;
    double width_of_Opening = config.getDouble("stoppingTarget.string_target_opening_size", 0.07366);
    double target_length = config.getDouble("stoppingTarget.target_length");

    double mesh_per_mm = 1.0 / (width_of_Opening + wire_diameter);
    int number_disks = (target_length - wire_diameter) / (wire_diameter + width_of_Opening);
    double stoppingTarget_z0InMu2e = config.getDouble("stoppingTarget.z0InMu2e", 5871.);
    cout << "-------- Design Details --------" << endl;
    cout << "Number of disks : " << number_disks << endl;
    cout << "mesh per mm : " << mesh_per_mm << endl;
    cout << "wire diameter : " <<wire_diameter <<endl;
    cout << "target start : "<< stoppingTarget_z0InMu2e <<endl;
    cout << "target center : "<<target_length/2 + stoppingTarget_z0InMu2e<<endl;
      cout << "target length : "<< target_length<< endl;
      cout << "target center : "<<target->centerInMu2e() << endl;
      cout<< "parent center : " << parent.centerInMu2e() << endl;
    //double dummyStoppingTarget_z0 = 5471.;

    /* Disk parameters */
    //double disk_rIn = 0.1;
    //double disk_rOut = 0.12;
    //vector<double> disk_rOut_list {190.5, 170.18, 149.86, 129.54, 109.22, 88.9};
    vector<double> disk_rOut_list;
    config.getVectorDouble("stoppingTarget.radius_list", disk_rOut_list);
    double disk_diff_rOut_rIn = sqrt(M_PI) * wire_diameter / 2.0;
    double disk_halfHeight = disk_diff_rOut_rIn / 2; /* assume cross area of disk is square to approximate circle */
    double deltaZBetweenTargets = width_of_Opening + disk_halfHeight * 2.0;
    cout << "disk_halfHeight : " << disk_halfHeight << endl;
    cout << "deltaZBetweenTargets : " << deltaZBetweenTargets << endl; 

    /* String parameters */
    double string_rIn = 0.0;
    double string_rOut = wire_diameter / 2;
    double string_halfHeight = (number_disks - 1) * deltaZBetweenTargets / 2
                                + disk_halfHeight;

    G4Material* Material = findMaterialOrThrow("StoppingTarget_Al");
    /* compute Mass of this design */
    double totMass = 0.0;
    // disk mass 
    double disk_totMass = 0.0;
    for (const auto& disk_rOut : disk_rOut_list) {
      double disk_rIn = disk_rOut - disk_diff_rOut_rIn;
      double single_disk_mass = M_PI * (disk_rOut*disk_rOut - disk_rIn*disk_rIn) * (2 * disk_halfHeight);
      disk_totMass += single_disk_mass;
      cout << "Disk Mass is : " << single_disk_mass * number_disks * 0.001 * 2.7 << endl; 
    }
    disk_totMass *= number_disks;
    // string mass
    double string_totMass = 0.0;
    double string_singleMass = M_PI * (wire_diameter/2)*(wire_diameter/2) * (2 * string_halfHeight);
    for (const auto& disk_rOut : disk_rOut_list) {
        int number_strings = 2 * M_PI * disk_rOut * mesh_per_mm;
        string_totMass += string_singleMass * number_strings;
        cout << "String Mass is : " << string_singleMass * number_strings * 0.001 * 2.7 << endl;
    }
    // total mass, density of Aluminum is 2.7 (g/cm3)
    totMass = (string_totMass + disk_totMass) * 0.001 * 2.7;
    cout << "Total Mass(gram) of design is : " << totMass << endl;




    /////////////////////////////////////////////
    // Loop all values in disk_rOut_list, to construct
    // many cylinders, each is made of meshes.
    /////////////////////////////////////////////
    for (const auto& disk_rOut : disk_rOut_list) {
      double disk_rIn = disk_rOut - disk_diff_rOut_rIn;
      int number_strings = 2 * M_PI * disk_rOut * mesh_per_mm;
      /////////////////////////////////////////////
      // Define a single cylinder made of meshes
      // Step 1: compute location and rotation.
      // Step 2: construction, place them in correct location and rotation.
      /////////////////////////////////////////////

      /* Step 1: location and rotation computation */
      /* disk location and rotation */
      vector<G4ThreeVector> diskInfo_location_list;
      for (int i = 0; i < number_disks; ++i) {
        double disk_px = 0.0;
        double disk_py = 0.0;
        double disk_pz = disk_halfHeight + i * deltaZBetweenTargets;
        diskInfo_location_list.push_back(G4ThreeVector(disk_px, disk_py, disk_pz));
       
      }
      /* string location and rotation */
      vector<G4ThreeVector> stringInfo_location_list;
      for (int i = 0; i < number_strings; ++i) {
        double theta = 2.0 * M_PI / number_strings;
        double string_px = ( disk_rOut + string_rOut + 0.0001 ) * cos(theta * i);
        double string_py = ( disk_rOut + string_rOut + 0.0001 ) * sin(theta * i);
        double string_pz = string_halfHeight;
        stringInfo_location_list.push_back(G4ThreeVector(string_px, string_py, string_pz));
      }


      // Step 2: construction
      // Define geometries and logical for disk, and
      // place them in correct location and rotation.
      /* disk construction */
      for (int i = 0; i < number_disks; ++i) {
        VolumeInfo diskInfo;
        string disk_no = to_string(disk_rOut);
        disk_no.push_back('_');
        disk_no += to_string(i);
        diskInfo.name = "Disk" + disk_no;
        diskInfo.solid = new G4Tubs(diskInfo.name,
                                      disk_rIn,
                                      disk_rOut,
                                      disk_halfHeight,
                                      0.,
                                      CLHEP::twopi);
        diskInfo.logical = new G4LogicalVolume(diskInfo.solid,
                                               Material,
                                               diskInfo.name);
        if(stSD) diskInfo.logical->SetSensitiveDetector(stSD);
        if(i==0) cout<<" Location of first disk "<<diskInfo_location_list[0]<<endl;
        if(i==number_disks-1) cout<<" Location of last disk "<<diskInfo_location_list[number_disks-1]<<endl;
        pv = new G4PVPlacement(0,
                               diskInfo_location_list[i],
                               diskInfo.logical,
                               diskInfo.name,
                               targetInfo.logical,
                               false,
                               0);

        doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);
        if (!config.getBool("stoppingTarget.visible",true)) {
            diskInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
          } else {
            G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
            visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
            visAtt->SetForceSolid(config.getBool("stoppingTarget.solid",true));
            diskInfo.logical->SetVisAttributes(visAtt);
          }
      }
     /* string construction */
     for (int i = 0; i < number_strings; ++i) {
       VolumeInfo stringInfo;
       string disk_no = to_string(disk_rOut);
       disk_no.push_back('_');
       disk_no += to_string(i);
       stringInfo.name = "String" + disk_no;
       stringInfo.solid = new G4Tubs(stringInfo.name,
                                     string_rIn,
                                     string_rOut,
                                     string_halfHeight,
                                     0.,
                                     CLHEP::twopi);
       stringInfo.logical = new G4LogicalVolume(stringInfo.solid,
                                              Material,
                                              stringInfo.name);
       if(stSD) stringInfo.logical->SetSensitiveDetector(stSD);
       
       pv = new G4PVPlacement(0,
                              stringInfo_location_list[i],
                              stringInfo.logical,
                              stringInfo.name,
                              targetInfo.logical,
                              false,
                              0);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);
       if (!config.getBool("stoppingTarget.visible",true)) {
           stringInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
         } else {
           G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
           visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
           visAtt->SetForceSolid(config.getBool("stoppingTarget.solid",true));
           stringInfo.logical->SetVisAttributes(visAtt);
         }
     }
  }
    return targetInfo;
  }


  /*
  Implementation of hexagon related stopping target designs, including
  1. hexagon_default: 19 cylinder combinations
  2. hexagon_mesh: change mesh of hexagon_default
  */
  VolumeInfo constructStoppingTarget_hexagon( VolumeInfo const& parent,
                                      SimpleConfig const& config){

    const bool forceAuxEdgeVisible  = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV              = true;

    int verbosity(config.getInt("target.verbosity",0));

    if ( verbosity > 1 ) std::cout<<"In constructStoppingTarget"<<std::endl;
    // Master geometry for the Target assembly
    GeomHandle<StoppingTarget> target;

    G4VSensitiveDetector* stSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::StoppingTarget());


    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    //    TubsParams targetMotherParams(0., target->cylinderRadius(), target->cylinderLength()/2.); //MR20150616 replace by hardcoded tube for primary volume
    TubsParams targetMotherParams(0., 300., 1000.);//MR20150616 hardcoded tube for primary volume

    VolumeInfo targetInfo = nestTubs("StoppingTargetMother",
                                     targetMotherParams,
                                     findMaterialOrThrow(target->fillMaterial()),
                                     0,
                                     target->centerInMu2e() - parent.centerInMu2e(),
                                     parent,
                                     0,
                                     false/*visible*/,
                                     G4Colour::Black(),
                                     false/*solid*/,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    G4VPhysicalVolume* pv;

    ////////////////////////////////////////////////
    // START create target foils
    ////////////////////////////////////////////////

    /* Structure parameters */
    double wire_diameter = config.getDouble("stoppingTarget.string_radius", 0.02665) * 2.0;
    double width_of_Opening = config.getDouble("stoppingTarget.string_target_opening_size", 0.07366);
    double length = config.getDouble("stoppingTarget.hexagon_length", 676.0);
    double outer_radius = config.getDouble("stoppingTarget.outer_radius", 21.0);
    double target_length = config.getDouble("stoppingTarget.target_length");

    double mesh_per_mm = 1.0 / (width_of_Opening + wire_diameter);
    double delta_pos = (target_length - length)/2.0;
    int number_disks = (length - wire_diameter) / (wire_diameter + width_of_Opening);
    cout << "Number of disks : " << number_disks << endl;
    //double dummyStoppingTarget_z0 = 5471.;


    /* Disk parameters */
    vector<double> disk_rOut_list;
    vector<double> x_list;
    vector<double> y_list;
    config.getVectorDouble("stoppingTarget.radius_list", disk_rOut_list);
    config.getVectorDouble("stoppingTarget.x_list", x_list);
    config.getVectorDouble("stoppingTarget.y_list", y_list);
    /*
    vector<double> disk_rOut_list(19, outer_radius);
    vector<double> x_list {0., 2., -2., 4., -4., 1., -1., 1., -1., 2., -2., 2., -2., 0., 0., 3., -3., 3., -3.};
    vector<double> y_list {0., 0., 0., 0., 0., sqrt(3), sqrt(3), -sqrt(3), -sqrt(3), 2.*sqrt(3), 2.*sqrt(3), 
         -2.*sqrt(3), -2.*sqrt(3), 2.*sqrt(3), -2.*sqrt(3), sqrt(3), sqrt(3), -sqrt(3), -sqrt(3)};
    */
    double disk_diff_rOut_rIn = sqrt(M_PI) * wire_diameter / 2.0;
    double disk_halfHeight = disk_diff_rOut_rIn / 2; /* assume cross area of disk is square to approximate circle */
    double deltaZBetweenTargets = width_of_Opening + disk_halfHeight * 2.0;
    int number_vertical_screen_strings = 2 * (outer_radius - wire_diameter) / (width_of_Opening + wire_diameter);
    CLHEP::HepRotationZ secRz(-M_PI/2.); // 90 degree rotation around z-axis
    CLHEP::HepRotationX secRx(-M_PI/2.); // 90 degree rotation around x-axis
    G4RotationMatrix* rotation_matrix_strings_parallel_xaxis = new G4RotationMatrix(secRx * secRz);
    G4RotationMatrix* rotation_matrix_strings_parallel_yaxis = new G4RotationMatrix(secRx);

    cout << "disk_halfHeight : " << disk_halfHeight << endl;
    cout << "deltaZBetweenTargets : " << deltaZBetweenTargets << endl; 

    /* String parameters */
    double string_rIn = 0.0;
    double string_rOut = wire_diameter / 2;
    double string_halfHeight = (number_disks - 1) * deltaZBetweenTargets / 2
                                + disk_halfHeight;

    G4Material* Material = findMaterialOrThrow("StoppingTarget_Al");
    /* compute Mass of this design */
    double totMass = 0.0;
    // disk mass 
    double disk_totMass = 0.0;
    for (const auto& disk_rOut : disk_rOut_list) {
      double disk_rIn = disk_rOut - disk_diff_rOut_rIn;
      double single_disk_mass = M_PI * (disk_rOut*disk_rOut - disk_rIn*disk_rIn) * (2 * disk_halfHeight);
      disk_totMass += single_disk_mass;
      cout << "Disk Mass is : " << single_disk_mass * number_disks * 0.001 * 2.7 << endl; 
    }
    disk_totMass *= number_disks;
    // string mass
    double string_totMass = 0.0;
    double string_singleMass = M_PI * (wire_diameter/2)*(wire_diameter/2) * (2 * string_halfHeight);
    for (const auto& disk_rOut : disk_rOut_list) {
        int number_strings = 2 * M_PI * disk_rOut * mesh_per_mm;
        string_totMass += string_singleMass * number_strings;
        cout << "String Mass is : " << string_singleMass * number_strings * 0.001 * 2.7 << endl;
    }
//CBE I would make a comment about magic numbers, but I guess you already know!
    // screen mass
    double screen_totMass = 0.0;
    for(int i=0; i<number_vertical_screen_strings; ++i) {
        double loc_ = 2. * outer_radius * (i+1) / number_vertical_screen_strings;
        double len_ = 2 * sqrt( pow(outer_radius, 2) - pow(outer_radius-loc_, 2) );
        double mass_ = M_PI * (wire_diameter/2)*(wire_diameter/2) * len_;
        screen_totMass += 4. * mass_;  // front/end; vertical/horizontal, thus 4 totally
    }
    //screen_totMass *= 20;
    screen_totMass *= 18; // There are 19 cylinders totally, but the centered one does not have screen, thus 18 screens for cylinders together.
    // total mass, density of Aluminum is 2.7 (g/cm3)
    totMass = (string_totMass + disk_totMass + screen_totMass) * 0.001 * 2.7;
    cout << "Total Mass(gram) of design is : " << totMass << endl;




    /////////////////////////////////////////////
    // Loop all values in disk_rOut_list, to construct
    // many cylinders, each is made of meshes.
    /////////////////////////////////////////////
    for (unsigned int circle_idx=0; circle_idx<disk_rOut_list.size(); ++circle_idx) {
      double disk_rOut = disk_rOut_list[circle_idx];
      double x_loc = x_list[circle_idx];
      double y_loc = y_list[circle_idx];
      double disk_rIn = disk_rOut - disk_diff_rOut_rIn;
      int number_strings = 2 * M_PI * disk_rOut * mesh_per_mm;
      /////////////////////////////////////////////
      // Define a single cylinder made of meshes
      // Step 1: compute location and rotation.
      // Step 2: construction, place them in correct location and rotation.
      /////////////////////////////////////////////

      /* Step 1: location and rotation computation */
      /* disk location and rotation */
      vector<G4ThreeVector> diskInfo_location_list;
      for (int i = 0; i < number_disks; ++i) {
        double disk_px = x_loc * (disk_rOut + 0.001);
        double disk_py = y_loc * (disk_rOut + 0.001);
        double disk_pz = disk_halfHeight + i * deltaZBetweenTargets + delta_pos;
        diskInfo_location_list.push_back(G4ThreeVector(disk_px, disk_py, disk_pz));
      }
      /* string location and rotation */
      vector<G4ThreeVector> stringInfo_location_list;
      for (int i = 0; i < number_strings; ++i) {
        double theta = 2.0 * M_PI / number_strings;
        double string_px = ( disk_rOut + string_rOut + 0.001 ) * cos(theta * i) + x_loc * (disk_rOut + 0.001);
        double string_py = ( disk_rOut + string_rOut + 0.001 ) * sin(theta * i) + y_loc * (disk_rOut + 0.001);
        double string_pz = string_halfHeight + delta_pos;
        stringInfo_location_list.push_back(G4ThreeVector(string_px, string_py, string_pz));
      }
      /* screen strings location and rotation */
      vector<G4ThreeVector> screenInfo_location_list;
      vector<double> screenInfo_halfHeight_list;
      if(circle_idx != 0) {
      for(int i=0; i<number_vertical_screen_strings-1; ++i) {
        double loc_ = 2. * outer_radius * (i+1) / number_vertical_screen_strings;
        double half_len_ = sqrt( pow(outer_radius, 2) - pow(outer_radius-loc_, 2) );
        
        double disk_px = x_loc * (disk_rOut + 0.001);
        double disk_py = y_loc * (disk_rOut + 0.001);
        double screen_px = disk_px - outer_radius + loc_;
        double screen_py = disk_py - outer_radius + loc_;
        double screen_pz_front = delta_pos;
        double screen_pz_end = disk_halfHeight + (number_disks - 1) * deltaZBetweenTargets + delta_pos;
        screenInfo_location_list.push_back(G4ThreeVector(screen_px, disk_py, screen_pz_front));
        screenInfo_location_list.push_back(G4ThreeVector(disk_px, screen_py, screen_pz_front));
        screenInfo_location_list.push_back(G4ThreeVector(screen_px, disk_py, screen_pz_end));
        screenInfo_location_list.push_back(G4ThreeVector(disk_px, screen_py, screen_pz_end));
        screenInfo_halfHeight_list.push_back(half_len_);
        screenInfo_halfHeight_list.push_back(half_len_);
        screenInfo_halfHeight_list.push_back(half_len_);
        screenInfo_halfHeight_list.push_back(half_len_);
      }
      }


      // Step 2: construction
      // Define geometries and logical for disk, and
      // place them in correct location and rotation.
      /* disk construction */
      for (int i = 0; i < number_disks; ++i) {
        VolumeInfo diskInfo;
        string disk_no = to_string(disk_rOut);
        disk_no.push_back('_');
        disk_no += to_string(i);
        diskInfo.name = "Disk" + disk_no;
        diskInfo.solid = new G4Tubs(diskInfo.name,
                                      disk_rIn,
                                      disk_rOut,
                                      disk_halfHeight,
                                      0.,
                                      CLHEP::twopi);
        diskInfo.logical = new G4LogicalVolume(diskInfo.solid,
                                               Material,
                                               diskInfo.name);
        if(stSD) diskInfo.logical->SetSensitiveDetector(stSD);
        pv = new G4PVPlacement(0,
                               diskInfo_location_list[i],
                               diskInfo.logical,
                               diskInfo.name,
                               targetInfo.logical,
                               false,
                               0);

        doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);
        if (!config.getBool("stoppingTarget.visible",true)) {
            diskInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
          } else {
            G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
            visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
            visAtt->SetForceSolid(config.getBool("stoppingTarget.solid",true));
            diskInfo.logical->SetVisAttributes(visAtt);
          }
      }
     /* string construction */
     for (int i = 0; i < number_strings; ++i) {
       VolumeInfo stringInfo;
       string disk_no = to_string(disk_rOut);
       disk_no.push_back('_');
       disk_no += to_string(i);
       stringInfo.name = "String" + disk_no;
       stringInfo.solid = new G4Tubs(stringInfo.name,
                                     string_rIn,
                                     string_rOut,
                                     string_halfHeight,
                                     0.,
                                     CLHEP::twopi);
       stringInfo.logical = new G4LogicalVolume(stringInfo.solid,
                                              Material,
                                              stringInfo.name);
       if(stSD) stringInfo.logical->SetSensitiveDetector(stSD);
       
       pv = new G4PVPlacement(0,
                              stringInfo_location_list[i],
                              stringInfo.logical,
                              stringInfo.name,
                              targetInfo.logical,
                              false,
                              0);
       doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);
       if (!config.getBool("stoppingTarget.visible",true)) {
           stringInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
         } else {
           G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
           visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
           visAtt->SetForceSolid(config.getBool("stoppingTarget.solid",true));
           stringInfo.logical->SetVisAttributes(visAtt);
         }
     }
     /* screen construction */
     for(unsigned int i=0; i<screenInfo_location_list.size(); ++i) {
         VolumeInfo screenInfo;
         screenInfo.name = "Screen_" + to_string(i);
         screenInfo.solid = new G4Tubs(screenInfo.name,
                                       string_rIn,
                                       string_rOut,
                                       screenInfo_halfHeight_list[i],
                                       0.,
                                       CLHEP::twopi);
         screenInfo.logical = new G4LogicalVolume(screenInfo.solid,
                                                  Material,
                                                  screenInfo.name);
         if(stSD) screenInfo.logical->SetSensitiveDetector(stSD);
        
         G4RotationMatrix* rotation;
         if(i%2==0) rotation = rotation_matrix_strings_parallel_yaxis;
         else rotation = rotation_matrix_strings_parallel_xaxis;
         pv = new G4PVPlacement(rotation,
                                screenInfo_location_list[i],
                                screenInfo.logical,
                                screenInfo.name,
                                targetInfo.logical,
                                false,
                                0);
        doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);
        if (!config.getBool("stoppingTarget.visible",true)) {
            screenInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
          } else {
            G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
            visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
            visAtt->SetForceSolid(config.getBool("stoppingTarget.solid",true));
            screenInfo.logical->SetVisAttributes(visAtt);
        }
     }
  }
    return targetInfo;
  }
    
} // end namespace mu2e
