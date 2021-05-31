//
// Free function to check magnetic field
//
// Original author Ivan Logashenko
//
// C++ includes
#include <iostream>
#include <cmath>

// Mu2e includes.
#include "Mu2eG4/inc/generateFieldMap.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/StraightSection.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"

#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// G4 includes
#include "Geant4/G4TransportationManager.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4FieldManager.hh"
#include "Geant4/G4Field.hh"

#include "TNtuple.h"

using namespace std;

namespace mu2e {

  void generateFieldMap( const G4ThreeVector & _origin, int mode) {

    if( mode<=0 ) return;

    art::ServiceHandle<art::TFileService> tfs;

    TNtuple * _nt = tfs->make<TNtuple>( "ntfield", "Field along main path",
                                        "xg:yg:zg:xl:yl:zl:nx:ny:nz:bx:by:bz:btot:bl:s");
    Float_t nt[100];

    // Map magnetic field

    double Bfield[10], point[4];
    G4Navigator *n = (G4TransportationManager::GetTransportationManager())->GetNavigatorForTracking();
    G4FieldManager *fm = (G4TransportationManager::GetTransportationManager())->GetFieldManager();

    GeomHandle<Beamline> beamg;
    double R = beamg->getTS().torusRadius();

    StraightSection const * ts3in =beamg->getTS().getTSCryo<StraightSection>( TransportSolenoid::TSRegion::TS3,
                                                                              TransportSolenoid::TSRadialPart::IN );
    double L = ts3in->getHalfLength();
    double Lturn = L + 3.14159/2*R;

    for( int i=-1100; i<2200; i++ ) {

      double s = i*10.0;
      double dx, dz, nx, nz;

      if( fabs(s)<=L ) {
        dx = -s;
        dz = 0;
        nx = -1;
        nz = 0;
      } else if( fabs(s)>Lturn ) {
        dx = -R-L;
        dz = fabs(s)-Lturn+R;
        nx = 0;
        nz = 1;
      } else {
        double phi = (fabs(s)-L)/R;
        dx = -L-R*sin(phi);
        dz = R-R*cos(phi);
        nz = sin(phi);
        nx = -cos(phi);
      }
      if( s<0 ) { dx*=-1; dz*=-1; }

      G4VPhysicalVolume* p1 = n->LocateGlobalPointAndSetup(_origin+G4ThreeVector(dx,0,dz));
      if( !p1 ) continue;
      //cout << "Volume information : " << p1->GetName() << endl;

      point[0]=_origin.x()+dx;
      point[1]=_origin.y();
      point[2]=_origin.z()+dz;
      point[3]=0;
      //cout << "Get field..." << endl;
      G4FieldManager *fmv = p1->GetLogicalVolume()->GetFieldManager();
      if( ! fmv ) fmv=fm;
      if( fmv ) {

        fmv->GetDetectorField()->GetFieldValue(point,Bfield);

        cout << "X=(" << dx << ", " << dz << ") "
             << "n=(" << nx << ", " << nz << ") "
             << "s=" << s << " "
             << "B=("
             << Bfield[0]/CLHEP::tesla << ", "
             << Bfield[1]/CLHEP::tesla << ", "
             << Bfield[2]/CLHEP::tesla << ") "
             << p1->GetName()
             << endl;

        nt[0] = point[0];
        nt[1] = point[1];
        nt[2] = point[2];
        nt[3] = dx;
        nt[4] = 0;
        nt[5] = dz;
        nt[6] = nx;
        nt[7] = 0;
        nt[8] = nz;
        nt[9] = Bfield[0]/CLHEP::tesla;
        nt[10] = Bfield[1]/CLHEP::tesla;
        nt[11] = Bfield[2]/CLHEP::tesla;
        nt[12] = sqrt(Bfield[0]*Bfield[0]+Bfield[1]*Bfield[1]+Bfield[2]*Bfield[2])/CLHEP::tesla;
        nt[13] = (Bfield[0]*nx+Bfield[2]*nz)/CLHEP::tesla;
        nt[14] = s;

        _nt->Fill(nt);

      } else {
        cout << "No field manager" << endl;
        continue;
      }

    }

  } // end Mu2eWorld::generateFieldMap

}
