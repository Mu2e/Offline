#ifndef WLSMaterials_h
#define WLSMaterials_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

class WLSMaterials
{
  public:

    ~WLSMaterials();
 
    static WLSMaterials* GetInstance();

    G4Material* GetMaterial(const G4String);
 
  private:
 
    WLSMaterials();

    void CreateMaterials();

  private:

    static WLSMaterials* instance;

    G4NistManager*     nistMan;

    G4Material*        Air;
    G4Material*        PVC;

    G4Material*        PMMA;
    G4Material*        FPethylene;
    G4Material*        PolystyreneFiber;
    G4Material*        PolystyreneScint;
    G4Material*        Epoxy;
    G4Material*        Coating;
};

#endif /*WLSMaterials_h*/
