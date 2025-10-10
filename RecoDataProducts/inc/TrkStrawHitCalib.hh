#ifndef RecoDataProducts_TrkStrawHitCalib_HH
#define RecoDataProducts_TrkStrawHitCalib_HH
#include "Offline/RecoDataProducts/inc/HitT0.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "KinKal/Detector/Residual.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/TrackerConditions/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include <functional>
namespace mu2e {
  struct TrkStrawHitCalib {
    // default constructor: initialization on declaration
    TrkStrawHitCalib() {}

    TrkStrawHitCalib(Tracker const& nominalTracker, StrawId const& strawid,
        KinKal::Residual const& rdresid,
        KinKal::Residual const& rlresid,
        KinKal::VEC3 const& dDdX,
        KinKal::ClosestApproachData const& rptca) :
      _dDdX(dDdX)
    {
      auto dDdP = rdresid.dRdP();
      auto dLdP = rlresid.dRdP();
      for (size_t i=0;i<6;i++){
        _dDdP[i] = dDdP[i];
        _dLdP[i] = dLdP[i];
      }

      const Panel& panel = nominalTracker.getPanel(strawid);
      const Plane& plane = nominalTracker.getPlane(strawid);
      auto vpoca = rptca.particlePoca().Vect();
      CLHEP::Hep3Vector poca(vpoca.X(),vpoca.Y(),vpoca.Z());
      auto panel_delta = poca - panel.origin();
      auto plane_delta = poca - plane.origin();

      if (rdresid.active()){
        CLHEP::Hep3Vector drdx(dDdX.X(),dDdX.Y(),dDdX.Z());
        _dDdPanelAlign[0] = panel.uDirection().dot(drdx);
        _dDdPanelAlign[1] = panel.vDirection().dot(drdx);
        _dDdPanelAlign[2] = panel.wDirection().dot(drdx);
        _dDdPlaneAlign[0] = plane.uDirection().dot(drdx);
        _dDdPlaneAlign[1] = plane.vDirection().dot(drdx);
        _dDdPlaneAlign[2] = plane.wDirection().dot(drdx);

        _dDdPanelAlign[3] = panel.uDirection().cross(panel_delta).dot(drdx);
        _dDdPanelAlign[4] = panel.vDirection().cross(panel_delta).dot(drdx);
        _dDdPanelAlign[5] = panel.wDirection().cross(panel_delta).dot(drdx);
        _dDdPlaneAlign[3] = plane.uDirection().cross(plane_delta).dot(drdx);
        _dDdPlaneAlign[4] = plane.vDirection().cross(plane_delta).dot(drdx);
        _dDdPlaneAlign[5] = plane.wDirection().cross(plane_delta).dot(drdx);
      }
    }

    XYZVectorF    _dDdX; // distance residual deriv wrt position
    std::array<float,6> _dDdP = {0,0,0,0,0,0}; // distance residual deriv wrt track parameters
    std::array<float,6> _dLdP = {0,0,0,0,0,0}; // longitudinal residual deriv wrt track parameters
    std::array<float,6> _dDdPlaneAlign = {0,0,0,0,0,0}; // distance residual deriv wrt plane alignment parameters
    std::array<float,6> _dDdPanelAlign = {0,0,0,0,0,0}; // distance residual deriv wrt panel alignment parameters
  };
}
#endif
