#ifndef CALORIMETER_HH
#define CALORIMETER_HH

#include <vector>

//
// Mu2e includes
#include "CalorimeterGeom/inc/CrystalDetail.hh"   
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/Vane.hh"

namespace mu2e{
  namespace calorimeter{
    class Calorimeter{

      friend class CalorimeterMaker;

    public:
      Calorimeter(){}
      ~Calorimeter(){};

    protected:

      //
      // information about crystals
      std::vector<CrystalDetail> _crystalDetail;

      //
      // A Calorimeter is vanes
      std::vector<Vane> _vanes;

      //
      // There will be pointers to the objects in this container.
      std::vector<Crystal>  _allCrystals;

    };
  } //namespace calorimeter
} //namespace mu2e
#endif
