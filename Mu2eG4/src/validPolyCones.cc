//
// Check all PolyCones in the G4SolidStore for a valid geometry.
//
#include "Offline/Mu2eG4/inc/validPolyCones.hh"

#include "Geant4/G4Polycone.hh"
#include "Geant4/G4SolidStore.hh"
#include "Geant4/G4VSolid.hh"

#include <set>

namespace {

  //  Helper class that represents one plane of a polycone
  struct pConePlane {

    // Parameters of one plane
    double z;
    double rmin;
    double rmax;

    pConePlane( double az, double armin, double armax):
      z(az), rmin(armin), rmax(armax){
    }

    // Operators needed for use with std::set and std::map
    bool operator==(pConePlane const& a) const{
      return ( z==a.z && rmin==a.rmin && rmax==a.rmax);
    }

    bool operator<(pConePlane const& a) const{
      if ( z    < a.z    ) return true;
      if ( rmin < a.rmin ) return true;
      if ( rmax < a.rmax ) return true;
      return false;
    }
  };

}

bool mu2e::validPolyCones( int printLevel ){

  // Return value defaults to the geometry being valid.
  bool retval=true;

  G4SolidStore* sstore = G4SolidStore::GetInstance();
  if ( printLevel >  0 ){
    G4cout << "Physical volume store size: " << sstore->size() << G4endl;
  }

  for ( auto i=sstore->begin(); i!=sstore->end(); ++i){

    G4VSolid* vsolid = *i;

    if ( vsolid->GetEntityType() != "G4Polycone" ){
      continue;
    }

    // Extract parameters from polycone
    G4Polycone const* tmp = static_cast<G4Polycone const *>(vsolid);
    G4PolyconeHistorical const* pcon = tmp->GetOriginalParameters();
    double* z    = pcon->Z_values;
    double* rmin = pcon->Rmin;
    double* rmax = pcon->Rmax;

    // Count the number of unique planes in this polycone.
    std::set<pConePlane> uniquePlanes;
    for ( int i=0; i<pcon->Num_z_planes; ++i ){
      uniquePlanes.emplace(z[i], rmin[i], rmax[i] );
    }

    if ( printLevel > 1 ) {
      G4cout << "Polycone: "
           << vsolid->GetName() << " "
           << vsolid->GetEntityType() << " "
           << pcon->Num_z_planes << " "
           << uniquePlanes.size()
           << G4endl;
    }

    // This polycone is invalid; it has at least two identical planes.
    if ( int(uniquePlanes.size()) != pcon->Num_z_planes ){
      retval = false;
      G4cout << "\nError: Polycone has identical planes.  Name of polycone:  " << vsolid->GetName() << G4endl;
      G4cout << "Number of planes: " << pcon->Num_z_planes << "   Number of unique planes: " << uniquePlanes.size() << G4endl;
      G4cout << "Polycone planes ( plane number, z, rmin, rmax): " << G4endl;
      for ( int i=0; i<pcon->Num_z_planes; ++i ){
        G4cout << "  "    << i
             << " z: "    << z[i]
             << " rmin: " << rmin[i]
             << " rmax: " << rmax[i]
             << G4endl;
      }
    }

  } // end loop over G4SolidStore

  return retval;
}
