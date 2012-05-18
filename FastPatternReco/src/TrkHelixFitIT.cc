//
// Object to perform helix fit to straw hits
//
// $Id: TrkHelixFitIT.cc,v 1.1 2012/05/18 18:14:36 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2012/05/18 18:14:36 $
//
//
// the following has to come before other BaBar includes
#include "FastPatternReco/inc/TrkHelixFitIT.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "art/Framework/Services/Optional/TFileService.h"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"

namespace mu2e 
{
  
  TrkHelixFitIT::TrkHelixFitIT(fhicl::ParameterSet const& pset) :
  TrkHelixFit(pset)
    {}

  TrkHelixFitIT::~TrkHelixFitIT()
    {}

  bool
  TrkHelixFitIT::findHelix(TrkDef const& mytrk,const points3D& pt,TrkHelix& myhel) {
    bool retval(false);
// loop over hits, and store the points
    std::vector<XYZP> xyzp;
    fillXYZP(mytrk,pt,xyzp);
// initialize the circle parameters
    if(xyzp.size() > _minnhit && initCircle(xyzp,myhel)){
// solve for the circle parameters
      retval = findXY(xyzp,myhel);
// extend those into the z direction
      if(retval) retval = findZ(xyzp,myhel);
    }
// set the success
    if(retval)
      myhel._fit = TrkErrCode(TrkErrCode::succeed);
    else
      myhel._fit = TrkErrCode(TrkErrCode::fail);
    return retval;
  }
  
  void
  TrkHelixFitIT::fillXYZP(TrkDef const& mytrk,const points3D& pt, std::vector<XYZP>& xyzp) {
    const Tracker& tracker = getTrackerOrThrow();
// loop over straw hits, and store their positions
    for (points3D::const_iterator loopPoints_it = pt.begin(); 
	 loopPoints_it != pt.end(); ++loopPoints_it) {

      StrawHit const& sh = mytrk.strawHitCollection()->at(loopPoints_it->getInEventHitID());
      CLHEP::Hep3Vector wpos=loopPoints_it->_pos;
      const Straw& straw = tracker.getStraw(sh.strawIndex());
      xyzp.push_back(XYZP(wpos,straw.getDirection(),loopPoints_it->_sigmaz,_sfac*straw.getRadius()));
    } 
// if requested, add the target
    if(_target){
      xyzp.push_back(XYZP(CLHEP::Hep3Vector(0.0,0.0,_targetz),CLHEP::Hep3Vector(1.0,0.0,0.0),_tsig,_tsig));
    }
  }

}
