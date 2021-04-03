#ifndef MU2EG4_HELICALPROTONABSORBER_HH
#define MU2EG4_HELICALPROTONABSORBER_HH
//
// Helical Proton Absorber main class
//
//
// Original author Suerfu, implemented by G. Tassielli
//
// Notes:
// Construct the Helical Proton Absorber

// C++ includes
#include <iostream>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
//#include "GeometryService/inc/GeometryService.hh"
//#include "GeometryService/inc/GeomHandle.hh"
//#include "GeometryService/inc/VirtualDetector.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
//#include "Mu2eG4/inc/nestCons.hh"

// G4 includes
//#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
//#include "Geant4/G4Box.hh"
//#include "Geant4/G4Cons.hh"
//#include "Geant4/G4Tubs.hh"
//#include "Geant4/G4BooleanSolid.hh"
//#include "Geant4/globals.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4TessellatedSolid.hh"
#include "Geant4/G4VPhysicalVolume.hh"


namespace mu2e {

//class G4LogicalVolume;
//class G4VPhysicalVolume;
//class G4Material;
//class PAbsSD;


class HelicalProtonAbsorber
{
public:
        HelicalProtonAbsorber(double z_start, double length_i, double *inner_radii,
                        double *outer_radii, double *inner_phis, double *outer_phis,
                        double thickness_i, /*double num_of_turns_i,*/ int vane_num_i,
                        G4Material* material, G4LogicalVolume* world);
        ~HelicalProtonAbsorber();

        inline G4LogicalVolume* GetLog() {return pabs_logic;}
        inline G4VPhysicalVolume* GetPhys(){return pabs_phys;}

        //void RebuildPAbs();

        double* getPhi_i(){return phi_i;};
        inline double getComponent_Phi_i(int i){return phi_i[i];};

        double* getPhi_o(){return phi_o;};
        double getComponent_Phi_o(int i){return phi_o[i];};

        double* getR_i(){return r_i;};
        double getComponent_R_i(int i){return r_i[i];};

        double* getR_o(){return r_o;};
        double getComponent_R_o(int i){return r_o[i];};

        double getZ_var(){return z_var;};
        int getIncrements_var(){return num_of_inc_var;};

        double getThick_var(){return thickness;};



        int getVaneNum_var(){return vane_num;};

        //double GETTEST(){return TEST;};

        G4ThreeVector GetPosition(){return position;};

        // These set functions are used to change values/parameters of the helical proton absorber during the optimization
        // process.

        void setPhi_i(double val[]){for (int i=0; i<3;i++) phi_i[i]=val[i];};
        void setPhi_i(int i, double val){phi_i[i]=val;};

        void setPhi_o(double val[]){for (int i=0; i<3;i++) phi_o[i]=val[i];};
        void setPhi_o(int i, double val){phi_o[i]=val;};

        void setR_i(double val[]){for (int i=0; i<3;i++) r_i[i]=val[i];};
        void setR_i(int i, double val){r_i[i]=val;};

        void setR_o(double val[]){for (int i=0; i<3;i++) r_o[i]=val[i];};
        void setR_o(int i, double val){r_o[i]=val;};

        void setZ_var(double z){z_var=z;};

        void setThick_var(double thick){thickness=thick;};

        void setIncrements_var(int num){num_of_inc_var = num;};

        void setVaneNum_var(int num){vane_num = num;};

        void setPosition(G4ThreeVector pos){position = pos;};

        bool checkOverlaps(int res=1000, double tol=0., bool verbose=true);

        void SetVisibility(bool forceSolid, bool forceAuxEdgeVisible, G4Color color, AntiLeakRegistry & reg);

private:

        G4ThreeVector rail_eq1(double z_init, int vane_id);
        G4ThreeVector rail_eq2(double z_init, int vane_id);
        G4ThreeVector rail_eq3(double z_init, int vane_id);
        G4ThreeVector rail_eq4(double z_init, int vane_id);

        void GenerateVanes();

        //double TEST;

        int vane_num;
        double num_of_turns;
        int num_of_increments;
        double PAbs_length;
        double thickness;
        double k;
        double inner_radius;
        double outer_radius;
        double step_length;
        G4TessellatedSolid* pabs_solid;

        G4ThreeVector position;

        double phi_i[3];
        double phi_o[3];
        double r_i[3];
        double r_o[3];
        double z_var;
        int num_of_inc_var;

        //PAbsSD* myPAbsSD;

        G4LogicalVolume* pabs_logic;
        G4VPhysicalVolume* pabs_phys;

        G4LogicalVolume* World;

        G4Material* material;
};

} //end mu2e

#endif // MU2EG4_HELICALPROTONABSORBER_HH
