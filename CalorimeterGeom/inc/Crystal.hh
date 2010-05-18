#ifndef CRYSTAL_HH
#define CRYSTAL_HH

//
// C++ includes
#include <vector>

//
// Mu2e includes


#include "CrystalId.hh"
#include "CrystalIndex.hh"
#include "CrystalDetail.hh"

//
// other includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e{
  namespace calorimeter{

    class Crystal{

      friend class RSlice;
      friend class ZSlice;
      friend class Vane;
      friend class Calorimeter;
      friend class CalorimeterMaker;


    public:

      Crystal();
      
      //
      // construct a crystal
      Crystal( const CrystalId& id,             // crystal identifier
               CrystalIndex index,              // index into dumb crystal array
               const CLHEP::Hep3Vector& center, // center of crystal
               const CrystalDetail* detail,     // dumb data describing crystal makeup
               CLHEP::Hep3Vector const& t       // unit vector along crystal axis from readout to opposite edge
               );
      
      ~Crystal(){};

      const CrystalId& Id() const { return _id;}
      CrystalIndex Index() const { return _index;}


    protected:


      // Identifier
      CrystalId _id;

      // Index into the array of all straws.
      CrystalIndex _index;

      // Mid-point of the straw.
      CLHEP::Hep3Vector _c;

      // Detailed description of a straw.
      const CrystalDetail* _detail;
      int _detailIndex;

      // Unit vector along the wire direction.
      // Need to add unit vectors along local u and v also.
      // Use Euler angle convention from G4.
      CLHEP::Hep3Vector _w;

      // Nearest neighbours.
      std::vector<const Crystal *> _nearest;

      std::vector<CrystalId> _nearestById;
      std::vector<CrystalIndex> _nearestByIndex;



    };
  }//namespace calorimeter
}//namespace mu2e

#endif
