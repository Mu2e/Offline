#ifndef FastPatternReco_TrkHelixFitIT_hh
#define FastPatternReco_TrkHelixFitIT_hh

//
// class derived by TrkPatRec/inc/TrkHelixFit.hh to work with the ITracker
//
// $Id: TrkHelixFitIT.hh,v 1.5 2012/12/05 18:49:00 brownd Exp $
// $Author: brownd $
// $Date: 2012/12/05 18:49:00 $
//
// Original author G. Tassielli
//

// Mu2e includes.
#include "fhiclcpp/ParameterSet.h"

//For track fit
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
//#include "KalmanTests/inc/TrkStrawHit.hh"
//#include "KalmanTests/inc/KalFit.hh"
//#include "KalmanTests/inc/KalFitMC.hh"
//#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/HelixFit.hh"
//#include "BTrk/TrkBase/TrkPoca.hh"
//#include "TrkPatRec/src/TrkPatRec_module.cc"


//// C++ includes.
//#include <iostream>
//#include <memory>
//#include <utility>
//#include <map>
//#include <vector>
//#include <algorithm>

namespace mu2e {

class TrkHelixFitIT : public HelixFit {
public:
        // parameter set should be passed in on construction
        explicit TrkHelixFitIT(fhicl::ParameterSet const& pset):HelixFit(pset) {}

        virtual ~TrkHelixFitIT(){}
        // main function: given a track definition, find the helix parameters
//        bool findHelix(  HelixFitResult& myfit, std::vector<XYZP> &xyzp );
        //CellGeometryHandle *_itwp;

private:
        //void fillXYZP( TrkDef const& mytrk, std::vector<XYZP> &xyzp, std::vector<points3D> &potLoops);

};

}
#endif /*FastPatternReco_TrkHelixFitIT_HH*/

//end Dave's TrkPatRec
