#ifndef DataProducts_CrystalId_hh
#define DataProducts_CrystalId_hh
//
// Offline identifier for one calorimeter crystal
// The number inside the CaloId will 2x crystal number,
// so it is important to use the crystal() accessor and not channel()
//
#include <array>
#include <algorithm>

namespace mu2e {

  typedef CrystalId CaloId;

};
#endif /* DataProducts_CrystalId_hh */
