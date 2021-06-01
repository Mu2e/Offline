
//
// Helical Proton Absorber main class
//
//
// Original author Suerfu, implemented by G. Tassielli
//
// Notes:
// Construct the Helical Proton Absorber

// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// C++ includes
#include <cmath>

// Mu2e includes.
#include "Mu2eG4/inc/HelicalProtonAbsorber.hh"

// G4 includes
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4LogicalVolume.hh"
//#include "Geant4/G4Material.hh"
#include "Geant4/G4PVPlacement.hh"
#include "Geant4/G4QuadrangularFacet.hh"
#include "Geant4/G4TriangularFacet.hh"
#include "Geant4/G4VisAttributes.hh"

//#include "Geant4/G4RunManager.hh"
//#include "PAbsSD.hh"
//#include "Geant4/G4SDManager.hh"

// I have created some comments on the files to generate the helical proton absorber geometry.
// We might have used a different notation than used in Fermilab. If there is any problem reading
// the code, please feel free to contact me at suerfu1@illinois.edu, or I can be reached via
// 217-979-8224.

//The specification of the proton absorber will work as long as we can specify the
//outline of the individual vanes in terms of its thickness, (r_inner(z),theta(z),z) and
//(r_outer(z),theta(z),z).  From the curves specifying the inner and outer radius, we can generate
//two additional curves representing an offset attributable to the thickness of the vane.
//These four curves will represent the edges of the vane in 3D space, which accordingly
//will have a parallelogram cross-section.  The construction of the geometrical object is then
//accomplished by sampling points from these curves and creating a tesselated volume.


//Future implementations will allow for vane cross-sections that are not parallelograms, but rather
//general shapes that can be parameterised in terms of the radius.

//For the helix, the r=constant, theta=k*z where k=pitch


//double total_arc_length = sqrt(pow(inner_radius*k,2)+1) * PAbs_length;

//double step_length = total_arc_length / num_of_surfaces;

//double step_length_z = total_length / num_of_surfaces;

namespace mu2e {

HelicalProtonAbsorber::HelicalProtonAbsorber(double z_start, double length_i, double *inner_radii,
                double *outer_radii, double *inner_phis, double *outer_phis, double thickness_i,
                /*double num_of_turns_i,*/ int vane_num_i, G4Material* plastic, G4LogicalVolume* world)
{
        //TEST = 400*CLHEP::mm;

        thickness = thickness_i;
        //num_of_turns = num_of_turns_i;
        vane_num = vane_num_i;
        material = plastic;
        World = world;
        PAbs_length = length_i;
        //k = num_of_turns/PAbs_length*CLHEP::twopi;

        //std::cout<<"k = "<<k*CLHEP::m<<"1/m\n";
        //std::cout<<"k = "<<k*CLHEP::cm<<"1/cm\n";
        //std::cout<<"k = "<<k*CLHEP::mm<<"1/mm\n";
        //std::cout<<"k = "<<k<<"\n";

        // the helical absorber is created incrementally. This num_of_increments specifies the number of increments.
        // seems too obvious to comment...

        num_of_increments = 130;
        step_length = PAbs_length/((double)num_of_increments);

        //Variational parameters
        //The arrays are used as corrections of different orders/polynomials.

        // the z_var is used to adjust the starting position of the proton absorber
        // the coordinate is relative to the center of the tracker in the detector solenoid
        z_var = z_start;//-3309.72*CLHEP::mm;
        num_of_inc_var = 274;

        //r_i is the inner radius
        r_i[0] = inner_radii[0];
        r_i[1] = inner_radii[1];
        r_i[2] = inner_radii[2];

        //r_o is the outer radius
        r_o[0] = outer_radii[0];
        r_o[1] = outer_radii[1];
        r_o[2] = outer_radii[2];

        // specifying the phi angle of the inner
        phi_i[0] = inner_phis[0];
        phi_i[1] = inner_phis[1];
        phi_i[2] = inner_phis[2];

        phi_o[0] = outer_phis[0];
        phi_o[1] = outer_phis[1];
        phi_o[2] = outer_phis[2];

        position = G4ThreeVector(0,0,z_var);
        /////////////////////////////////////

        pabs_solid = new G4TessellatedSolid("Helical_PAbs");

        // the function generate vanes will tesselate the vanes.
        GenerateVanes();

        //creating logical and physical volumes
        pabs_logic = new G4LogicalVolume(pabs_solid, material, "helical_pabs_log");
        pabs_phys = new G4PVPlacement(0,position,pabs_logic,"helical_pabs_phys", World, false, 0);
        pabs_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
}

HelicalProtonAbsorber::~HelicalProtonAbsorber()
{

}

/*void HelicalProtonAbsorber::RebuildPAbs()
{
        // This function is used to destroy the previous proton absorber, and rebuild one using new
        // parameters.

        delete pabs_phys;
        delete pabs_logic;
        delete pabs_solid;
        //delete myPAbsSD;

        // it seems that if the sensitive detector pointer is deleted, it would cause problems, so
        // instead of deleting, reassign the pointer. Maybe the Geant4 system at Fermilab has the expertise to handle
        // this problem.

        pabs_solid = new G4TessellatedSolid("Helical_PAbs");

        // Regenerate the vanes
        GenerateVanes();


        pabs_logic = new G4LogicalVolume(pabs_solid, material, "helical_pabs_log");
        pabs_phys = new G4PVPlacement(0,G4ThreeVector(0,0,z_var),pabs_logic,"helical_pabs_phys", World, false, 0);

//        // Notify the G4RunManager that the geometry has been rebuilt.
//        G4RunManager::GetRunManager()->GeometryHasBeenModified();

}*/

G4ThreeVector HelicalProtonAbsorber::rail_eq1(double z_init, int vane_id)
{
        // This rail equation is the three-dimensional outlines of the helical proton absorber.
        // The radius and the angle is given as the function of Z coordinate.

        double z = z_init;
        double phi =  phi_i[0] + phi_i[1]*z_init + phi_i[2]*pow(z_init,2) + CLHEP::twopi/vane_num*vane_id;
        double radius = r_i[0] + r_i[1]*z_init + r_i[2]*pow(z_init,2);


        /*double ang_displacement = CLHEP::twopi/vane_num*vane_id; //vane id goes from 0 to vane_num-1

    double x = radius*cos(k*z_init + ang_displacement);
    double y = radius*sin(k*z_init + ang_displacement);
    double z = z_init;
         */

        // Create a position vector, and then set the polar parameters of this vector.
        G4ThreeVector pos_1 = G4ThreeVector(1,1,1);

        pos_1.setRhoPhiZ(radius,phi,z);

        return pos_1;
}

G4ThreeVector HelicalProtonAbsorber::rail_eq2(double z_init, int vane_id)
{

        // This is essentially same of rail_equation1, but this equation holds Rho and Phi coordinate fixted
        // while giving an increment(thickness) in the Z-direction. So the Thickness is in Z-direction, not perpendicular to
        // the surface.
        // When you do simulations at Fermilab, it might be desirable to use perpendicular distance as thickness.
        // That can be done by taking the cross product at a lattice point, and increment in the direction of the cross product.

        G4ThreeVector pos_1 = rail_eq1(z_init, vane_id) + G4ThreeVector(0,0,-thickness);

        return pos_1;
}

G4ThreeVector HelicalProtonAbsorber::rail_eq3(double z_init, int vane_id)
{
        //double z = z_init-thickness;
        //double radius_start = 25*CLHEP::cm;
        //double radius_final = 40*CLHEP::cm;
        //double radius = radius_final;
        //double taper_length = 1*PAbs_length;

        // double phi = k*(z_init + z_init*z_init/PAbs_length) + CLHEP::twopi/vane_num*vane_id;

        // if(z_init<taper_length)
        // {
        //     radius = (taper_length-z_init)/(taper_length)*radius_start + (z_init)/(taper_length)*radius_final;
        // }

        // else double radius = radius_final;



        // G4ThreeVector pos_1 = G4ThreeVector(1,1,1);

        //pos_1.setRhoPhiZ(radius,phi,z);


        G4ThreeVector pos_1 = rail_eq4(z_init, vane_id) + G4ThreeVector(0,0,-thickness);

        return pos_1;
}

G4ThreeVector HelicalProtonAbsorber::rail_eq4(double z_init,int vane_id)
{
        double z = z_init;
        double phi =  phi_o[0] + phi_o[1]*z_init + phi_o[2]*pow(z_init,2) + CLHEP::twopi/vane_num*vane_id;
        double radius = r_o[0] + r_o[1]*z_init + r_o[2]*pow(z_init,2);

        // if(z_init<taper_length)
        //{
        //    radius = (taper_length-z_init)/(taper_length)*radius_start + (z_init)/(taper_length)*radius_final;
        // }

        // else double radius = radius_final;

        G4ThreeVector pos_1 = G4ThreeVector(1,1,1);

        pos_1.setRhoPhiZ(radius,phi,z);

        return pos_1;
}

void HelicalProtonAbsorber::GenerateVanes()
{

        // this code tessellates the helical proton absorber. For each increment, the tessellation is done by the four lattices.

        for (int vane_id=0; vane_id < vane_num; vane_id++)
        {


                for (int i=0; i<num_of_inc_var; i++)
                {
                        double z_i = step_length*((double)i);

                        G4ThreeVector p1_i = rail_eq1(z_i, vane_id);
                        G4ThreeVector p2_i = rail_eq2(z_i, vane_id);
                        G4ThreeVector p3_i = rail_eq3(z_i, vane_id);
                        G4ThreeVector p4_i = rail_eq4(z_i, vane_id);

                        G4ThreeVector p1_f = rail_eq1(z_i + step_length, vane_id);
                        G4ThreeVector p2_f = rail_eq2(z_i + step_length, vane_id);
                        G4ThreeVector p3_f = rail_eq3(z_i + step_length, vane_id);
                        G4ThreeVector p4_f = rail_eq4(z_i + step_length, vane_id);

                        //At the beginning of the tesselation algorithm, create beginning cap.
                        if (i==0)
                        {

                                G4TriangularFacet *initial_cap1 = new G4TriangularFacet(p1_i,p2_i,p3_i,ABSOLUTE);
                                G4TriangularFacet *initial_cap2 = new G4TriangularFacet(p3_i,p4_i,p1_i,ABSOLUTE);

                                pabs_solid->AddFacet((G4VFacet*) initial_cap1);
                                pabs_solid->AddFacet((G4VFacet*) initial_cap2);
                        }

                        //At the end of the tesselation algorithm, create end cap.
                        if (i==num_of_increments-1)
                        {
                                G4TriangularFacet *final_cap1 = new G4TriangularFacet(p1_f,p3_f,p2_f,ABSOLUTE);
                                G4TriangularFacet *final_cap2 = new G4TriangularFacet(p1_f,p4_f,p3_f,ABSOLUTE);

                                pabs_solid->AddFacet((G4VFacet*) final_cap1);
                                pabs_solid->AddFacet((G4VFacet*) final_cap2);
                        }

                        G4TriangularFacet* inner_facet1 = new G4TriangularFacet(p1_i,p1_f,p2_i,ABSOLUTE);
                        G4TriangularFacet* inner_facet2 = new G4TriangularFacet(p1_f,p2_f,p2_i,ABSOLUTE);

                        G4TriangularFacet* outer_facet1 = new G4TriangularFacet(p4_i,p3_i,p4_f,ABSOLUTE);
                        G4TriangularFacet* outer_facet2 = new G4TriangularFacet(p4_f,p3_i,p3_f,ABSOLUTE);

                        G4TriangularFacet* bottom_facet1 = new G4TriangularFacet(p1_i,p4_f,p1_f,ABSOLUTE);
                        G4TriangularFacet* bottom_facet2 = new G4TriangularFacet(p1_i,p4_i,p4_f,ABSOLUTE);

                        G4TriangularFacet* top_facet1 = new G4TriangularFacet(p2_i,p2_f,p3_i, ABSOLUTE);
                        G4TriangularFacet* top_facet2 = new G4TriangularFacet(p2_f,p3_f,p3_i, ABSOLUTE);


                        pabs_solid->AddFacet((G4VFacet*) bottom_facet1);
                        pabs_solid->AddFacet((G4VFacet*) bottom_facet2);

                        pabs_solid->AddFacet((G4VFacet*) top_facet1);
                        pabs_solid->AddFacet((G4VFacet*) top_facet2);

                        pabs_solid->AddFacet((G4VFacet*) inner_facet1);
                        pabs_solid->AddFacet((G4VFacet*) inner_facet2);

                        pabs_solid->AddFacet((G4VFacet*) outer_facet1);
                        pabs_solid->AddFacet((G4VFacet*) outer_facet2);

                }
        }

        pabs_solid->SetSolidClosed(true);
}

bool HelicalProtonAbsorber::checkOverlaps(int res, double tol, bool verbose){
     return pabs_phys->CheckOverlaps( res, tol, verbose );
}

void HelicalProtonAbsorber::SetVisibility(bool forceSolid, bool forceAuxEdgeVisible, G4Color color, AntiLeakRegistry & reg){

        G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, color));
        visAtt->SetForceSolid(forceSolid);
        // If I do not do this, then the rendering depends on what happens in
        // other parts of the code;  is there a G4 bug that causes something to be
        // unitialized?
        visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
        pabs_logic->SetVisAttributes(visAtt);

}


} //end mu2e







