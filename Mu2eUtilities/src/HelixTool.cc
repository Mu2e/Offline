#include <cmath>
#include <stddef.h>
#include <algorithm>
#include <array>
#include <utility>
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

  HelixTool::HelixTool(const HelixSeed *Helix, const mu2e::Tracker*MyTracker, StrawNoSignal NoSignal) {
    _hel           = Helix;
    _tracker       = MyTracker;
    _noSignal      = std::move(NoSignal);
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
  }


  float  HelixTool::hitRatio() const {
    if (_hitRatio >= 0.) return _hitRatio;

    // ratio of the number of measured hits to the number of expected ones: count the plane
    // faces the helix crosses inside the tracker active radius, weighted by the fraction of
    // live straws at the crossing
    const mu2e::RobustHelix& robustHel = _hel->helix();
    float    expected_faces(1e-6);
    for (size_t planeId=0; planeId<_tracker->nPlanes(); planeId++) {
      const Plane& pln = _tracker->getPlane(planeId);
      if (pln.nPanels() == 0 )         continue;

      for (uint16_t face=0; face<StrawId::_nfaces; ++face){
        // panel number 'face' is on that face
        const Panel& panel = pln.getPanel(face);
        float    z = (panel.getStraw(0).getMidPoint().z()+panel.getStraw(1).getMidPoint().z())/2.;

        XYZVectorF  pos;
        pos.SetZ(z);
        robustHel.position(pos);

        //now check that we are in the active area of tracker
        float   hitR = sqrtf(pos.x()*pos.x() + pos.y()*pos.y());
        if (hitR < _trackerRIn )  continue;
        if (hitR > _trackerROut)  continue;

        expected_faces += _noSignal ? liveFraction(pln, face, pos.x(), pos.y()) : 1.;
      }//end loop over the two faces
    }//end loop over the planes

    _hitRatio = _nStrawHits/expected_faces;
    return _hitRatio;
  }


  float  HelixTool::liveFraction(const Plane& Pln, uint16_t Face, float X, float Y) const {
    // find the panel of this face covering the crossing point: the one whose outward direction,
    // perpendicular to its straws in the transverse plane, is closest to the crossing direction
    const Panel* best(nullptr);
    float        bestRho(0), nx(0), ny(0);
    for (int ipn=Face; ipn<int(Pln.nPanels()); ipn+=StrawId::_nfaces) {
      const Panel& panel = Pln.getPanel(ipn);
      const auto&  mid   = panel.getStraw(0).getMidPoint();
      const auto&  dir   = panel.getStraw(0).getDirection();
      float        along = mid.x()*dir.x() + mid.y()*dir.y();
      float        px    = mid.x() - along*dir.x();
      float        py    = mid.y() - along*dir.y();
      float        norm  = sqrtf(px*px + py*py);
      if (norm <= 0.)  continue;
      px /= norm;
      py /= norm;
      float        rho   = X*px + Y*py;
      if (best == nullptr || rho > bestRho) {
        best    = &panel;
        bestRho = rho;
        nx      = px;
        ny      = py;
      }
    }
    if (best == nullptr) return 1.;

    // straws are ordered radially (straw 0 innermost), alternating between the two layers:
    // estimate the straw at the crossing and take the nearest straw of each layer around it
    auto  strawRho = [&](int s) {
      const auto& mid = best->getStraw(s).getMidPoint();
      return float(mid.x()*nx + mid.y()*ny);
    };
    const int nStraws = best->nStraws();
    float rho0   = strawRho(0);
    float rhoMax = strawRho(nStraws-1);
    // no straw of this face at the crossing (radially outside the panel, or beyond the straw
    // ends, e.g. between the panels of the face): keep the weight without the status
    if (bestRho < rho0 || bestRho > rhoMax) return 1.;
    int   sEst   = (rhoMax > rho0) ? std::lround((bestRho-rho0)/(rhoMax-rho0)*(nStraws-1)) : 0;
    sEst = std::clamp(sEst, 0, nStraws-1);
    {
      const Straw& straw = best->getStraw(sEst);
      const auto&  mid   = straw.getMidPoint();
      const auto&  dir   = straw.getDirection();
      if (std::fabs((X-mid.x())*dir.x() + (Y-mid.y())*dir.y()) > straw.halfLength()) return 1.;
    }

    float nLive(0);
    for (int layer=0; layer<StrawId::_nlayers; ++layer) {
      int   nearest(-1);
      float dmin(0);
      for (int s=std::max(0, sEst-2); s<=std::min(nStraws-1, sEst+2); ++s) {
        if (best->getStraw(s).id().layer() != layer) continue;
        float d = std::fabs(strawRho(s) - bestRho);
        if (nearest < 0 || d < dmin) {
          nearest = s;
          dmin    = d;
        }
      }
      if (nearest >= 0 && !_noSignal(best->getStraw(nearest).id())) ++nLive;
    }
    return nLive/StrawId::_nlayers;
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
