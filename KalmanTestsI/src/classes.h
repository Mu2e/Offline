#include "art/Persistency/Common/Wrapper.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
using namespace CLHEP;
#include "TrkBase/HelixParams.hh"
#include "difAlgebra/DifIndepPar.hh"
#include "difAlgebra/DifArray.hh"
#include "difAlgebra/DifNumber.hh"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"

#include "KalmanTestsI/inc/kalFitDataOuts.hh"


//template class art::Wrapper<HelixParams>;
template class std::vector<mu2e::HitInfo>;
template class std::vector<mu2e::TrkCellHitInfo>;
