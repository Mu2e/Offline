#include <cmath>
#include <stddef.h>
#include <array>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "Offline/TrackerGeom/inc/Plane.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/Mu2eUtilities/inc/HelixTool.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"


namespace mu2e {

  HelixTool::HelixTool(const HelixSeed *Helix, const mu2e::Tracker*MyTracker) {
    _hel           = Helix;
    _tracker       = MyTracker;
    _trackerRIn    = _tracker->g4Tracker()->getInnerTrackerEnvelopeParams().innerRadius();
    _trackerROut   = _tracker->g4Tracker()->getInnerTrackerEnvelopeParams().outerRadius();

    //initialize
    _meanHitRadialDist = 0.;
    _nLoops            = 0;
    const ComboHit*     hit(0);

    float         z_first_hit(0), z_last_hit(0), counter(0);
    bool          isFirst(true);
    unsigned      nhits = _hel->_hhits.size();

    const mu2e::RobustHelix  *robustHel = &Helix->helix();

    int nstrawhits = 0;

    for (unsigned f=0; f<nhits; ++f){
      hit = &_hel->_hhits[f];
      if (hit->_flag.hasAnyProperty(StrawHitFlag::outlier))     continue;
      nstrawhits += hit -> nStrawHits();

      _meanHitRadialDist += sqrtf(hit->pos().x()*hit->pos().x() + hit->pos().y()*hit->pos().y());
      ++counter;
      float z = hit->pos().z();
      if (isFirst){
        z_first_hit = z;
        z_last_hit  = z;
        isFirst     = false;
      }else {
        z_last_hit  = z;
      }
    }//end loop over the hits

    _nStrawHits = nstrawhits;

    if (counter > 0) _meanHitRadialDist /= counter;

    _nLoops = std::fabs((z_last_hit - z_first_hit)/(robustHel->lambda()*2.*M_PI));

    //here we evaluate once the impact parameter
    _d0 = robustHel->rcent  () - robustHel->radius ();


    // we now estiamte the ratio of the number of measured hits to the number of the expected ones
    // we make a few assumptions and appriximations:
    float    expected_faces(1e-6);
    for (size_t planeId=0; planeId<_tracker->nPlanes(); planeId++) {
      const Plane* pln = &_tracker->getPlane(planeId);
      int   nPanels = pln->nPanels();
      if (nPanels == 0 )         continue;
      std::array<int,2>   idPanels = {0, (nPanels-1)};

      for (size_t ipn=0; ipn<idPanels.size(); ++ipn){
        const Panel* panel = &pln->getPanel(ipn);
        float    z = (panel->getStraw(0).getMidPoint().z()+panel->getStraw(1).getMidPoint().z())/2.;
        // if (z < z_first_hit )  continue;
        // if (z > z_last_hit  )  continue;

        XYZVectorF  pos;
        pos.SetZ(z);
        robustHel->position(pos);

        //now check that we are in the active area of tracker
        float   hitR = sqrtf(pos.x()*pos.x() + pos.y()*pos.y());
        if (hitR < _trackerRIn )  continue;
        if (hitR > _trackerROut)  continue;

        ++expected_faces;
      }//end loop over the two panels
    }//end loop over the planes

    //each plane has two layers of panels. We are assuming 2 ComboHits per face
    _hitRatio = _nStrawHits/expected_faces;
  }


  void   HelixTool::dirOfProp(float& Slope, float& SlopeErr, float& Chi2ndof){
    ::LsqSums2 fitDtDz;
    const mu2e::ComboHit* hit;
    for (size_t j=0; j<_hel->_hhits.size(); j++) {
      hit = &_hel->_hhits[j];
      if (hit->_flag.hasAnyProperty(StrawHitFlag::outlier))     continue;
      double hitTime = hit->correctedTime();
      double hitZpos = hit->pos().z();
      double timeErrSquared = hit->timeVar();//ns^2
      double hitWeight      = 1./timeErrSquared;
      fitDtDz.addPoint(hitZpos,hitTime, hitWeight);
    }
    Slope    = fitDtDz.dydx();
    Chi2ndof = fitDtDz.chi2Dof();
    SlopeErr = fitDtDz.dydxErr();
  }
}
