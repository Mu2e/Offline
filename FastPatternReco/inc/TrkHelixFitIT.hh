#ifndef FastPatternReco_TrkHelixFitIT_hh
#define FastPatternReco_TrkHelixFitIT_hh

//
// class derived by TrkPatRec/inc/TrkHelixFit.hh to work with the ITracker
//
// $Id: TrkHelixFitIT.hh,v 1.3 2012/07/23 17:52:27 brownd Exp $
// $Author: brownd $
// $Date: 2012/07/23 17:52:27 $
//
// Original author G. Tassielli
//

// Mu2e includes.
#include "fhiclcpp/ParameterSet.h"

//For track fit
// BaBar
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
//#include "KalmanTests/inc/TrkStrawHit.hh"
//#include "KalmanTests/inc/KalFit.hh"
//#include "KalmanTests/inc/KalFitMC.hh"
//#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/TrkHelixFit.hh"
//#include "TrkBase/TrkPoca.hh"
//#include "TrkPatRec/src/TrkPatRec_module.cc"


//// C++ includes.
//#include <iostream>
//#include <memory>
//#include <utility>
//#include <map>
//#include <vector>
//#include <algorithm>

namespace mu2e {

class TrkHelixFitIT : public TrkHelixFit {
public:
        // parameter set should be passed in on construction
        explicit TrkHelixFitIT(fhicl::ParameterSet const& pset):TrkHelixFit(pset) {}

        virtual ~TrkHelixFitIT(){}
        // main function: given a track definition, find the helix parameters
        bool findHelix( TrkDef const& mytrk, TrkHelix& myfit, std::vector<XYZP> &xyzp );
        //CellGeometryHandle *_itwp;

private:
        //void fillXYZP( TrkDef const& mytrk, std::vector<XYZP> &xyzp, std::vector<points3D> &potLoops);

};

}
#endif /*FastPatternReco_TrkHelixFitIT_HH*/

//end Dave's TrkPatRec
