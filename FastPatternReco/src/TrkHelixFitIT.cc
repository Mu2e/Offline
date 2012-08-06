//
// class derived by TrkPatRec/inc/TrkHelixFit.hh to work with the ITracker
//
// $Id: TrkHelixFitIT.cc,v 1.3 2012/08/06 17:00:52 brownd Exp $
// $Author: brownd $
// $Date: 2012/08/06 17:00:52 $
//
// Original author G. Tassielli
//

// Mu2e includes.
#include "FastPatternReco/inc/TrkHelixFitIT.hh"

namespace mu2e {


// main function: given a track definition, find the helix parameters
bool TrkHelixFitIT::findHelix( TrkDef const& mytrk, TrkHelix& myfit, std::vector<XYZP> &xyzp ) {
        bool retval(false);
        // loop over hits, and store the points
        //std::vector<XYZP> xyzp;
        //fillXYZP(mytrk,xyzp,potLoops);

        // initialize the circle parameters
        if(xyzp.size() > _minnhit && initCircle(xyzp,myfit)){
                // solve for the circle parameters
                retval = findXY(xyzp,myfit);
                // extend those into the z direction
                if(retval) retval = findZ(xyzp,myfit);
        }
        // set the success
        if(retval)
                myfit._fit = TrkErrCode(TrkErrCode::succeed);
        else
                myfit._fit = TrkErrCode(TrkErrCode::fail);
        return retval;
}

//void TrkHelixFitIT::fillXYZP( TrkDef const& mytrk, std::vector<XYZP> &xyzp, std::vector<points3D> &potLoops) {
//
//        for (std::vector<points3D>::iterator potLoops_it = potLoops.begin(); potLoops_it != potLoops.end(); ++potLoops_it) {
//                for (points3D::iterator loopPoints_it = potLoops_it->begin(); loopPoints_it != potLoops_it->end(); ++loopPoints_it) {
//                        _itwp->SelectCell(loopPoints_it->getRadLayID(),loopPoints_it->getInLayerCellID());
//                        xyzp.push_back(XYZP( loopPoints_it->_pos, _itwp->GetCellDirection(), 0.0/*loopPoints_it->_sigmaz*/,loopPoints_it->_sigmax)); //sigmax=sigmay=sigmaR
//                }
//        }
//
//}


}
