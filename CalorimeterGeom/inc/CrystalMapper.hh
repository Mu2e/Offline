#ifndef CalorimeterGeom_CrystalMapper_hh
#define CalorimeterGeom_CrystalMapper_hh
//
// Interface for classes describing the layout of the crystals in the disk.
//
// Crystals could be square, hexagonal or triangular.
//   Need to describe tesselation and shape
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

          virtual int                        nCrystalMax    (int maxRing)                   const = 0;

          virtual CLHEP::Hep2Vector          xyFromIndex    (int thisIndex)                 const = 0;
          virtual int                        indexFromXY    (double x, double y)            const = 0;
          virtual int                        indexFromRowCol(int nRow, int nCol)            const = 0;
          virtual int                        rowFromIndex   (int thisIndex)                 const = 0;
          virtual int                        colFromIndex   (int thisIndex)                 const = 0;

          virtual int                        numNeighbors   (int level)                     const = 0;
          virtual std::vector<int>           neighbors      (int thisIndex, int level)      const = 0;

          virtual const std::vector<double>& apexX          ()                              const = 0;
          virtual const std::vector<double>& apexY          ()                              const = 0;

          virtual bool                       isInsideCrystal(double x, double y,
                                                             const CLHEP::Hep3Vector& pos,
                                                             const CLHEP::Hep3Vector& size) const = 0;
    };

}
#endif
