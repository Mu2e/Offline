//
// Object to perform helix fit to straw hits
//
// $Id: StriaghtTrackFit code
// $Author: S Middleton
// $Date: Nov 2018
//
// Mu2E
#include "TrkReco/inc/Chi2HelixFit.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "Mu2eUtilities/inc/polyAtan2.hh"

//ROOT:
#include "TMatrixD.h"

//ACCUMULATORS:
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include "boost_fix/accumulators/statistics.hpp"
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>

//C++:
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <cmath>

using namespace std;
using namespace boost::accumulators;
using namespace ROOT::Math::VectorUtil;

namespace mu2e
{









}//end namespace
