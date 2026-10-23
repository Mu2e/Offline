#ifndef CalorimeterGeom_CrystalMapper_hh
#define CalorimeterGeom_CrystalMapper_hh
//
// Interface for classes describing the layout of the crystals in the disk.
//
// Crystals could be square, hexagonal or triangular
//
// Original author B. Echenard
//

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <vector>

namespace mu2e {

  class CrystalMapper {
    public:
      //no constructor for this interface
      virtual ~CrystalMapper() = default;

      virtual unsigned                   nCrystalMax    (unsigned maxRing)                   const = 0;
      virtual CLHEP::Hep2Vector          xyFromIndex    (unsigned thisIndex)                 const = 0;
      virtual unsigned                   indexFromXY    (double x, double y)                 const = 0;
      virtual unsigned                   indexFromRowCol(int nRow, int nCol)                 const = 0;
      virtual int                        rowFromIndex   (unsigned thisIndex)                 const = 0;
      virtual int                        colFromIndex   (unsigned thisIndex)                 const = 0;
      virtual unsigned                   numNeighbors   (unsigned level)                     const = 0;
      virtual std::vector<unsigned>      neighbors      (unsigned thisIndex, unsigned level) const = 0;
   };

}
#endif
