

#include "CLHEP/Vector/ThreeVector.h"
#include "DataProducts/inc/StrawId.hh"
#include "GeneralUtilities/inc/HepTransform.hh"
#include "TrackerGeom/inc/Tracker.hh"

using namespace mu2e;
using namespace CLHEP;

namespace AlignStraw {


// carbon copy of AlignedTrackerMaker::fromDb
std::pair<Hep3Vector, Hep3Vector> alignStraw(Tracker const& tracker, Plane const& plane,
                                             Panel const& panel, StrawId const& strawId,
                                             HepTransform const& align_tracker,
                                             HepTransform const& align_plane,
                                             HepTransform const& align_panel) {
  std::pair<Hep3Vector, Hep3Vector> result;
  
  // the whole tracker has nominal center on 0,0,0
  // how to place the plane in the tracker
  HepTransform plane_to_tracker(0.0, 0.0, plane.origin().z(), 0.0, 0.0, 0.0);

  // make an intermediate multiplication
  HepTransform plane_temp = align_tracker * (plane_to_tracker * align_plane);

  // how to place the panel in the plane
  Hep3Vector dv = panel.straw0MidPoint() - plane_to_tracker.displacement();
  double rz = dv.phi();

  HepTransform panel_to_plane(dv.x(), dv.y(), dv.z(), 0.0, 0.0, rz);

  // make an intermediate multiplication
  HepTransform panel_temp = plane_temp * (panel_to_plane * align_panel);

  Straw const& straw = tracker.getStraw(strawId);

  // how to place the straw in the panel
  double dx = straw.getMidPoint().perp() - panel.straw0MidPoint().perp();
  double dz = (straw.getMidPoint() - panel.straw0MidPoint()).z();

  Hep3Vector straw_to_panel = Hep3Vector(dx, 0.0, dz);
  Hep3Vector straw_dir = Hep3Vector(0.0, 1.0, 0.0);

  Hep3Vector aligned_straw = panel_temp * straw_to_panel;
  Hep3Vector aligned_straw_dir = panel_temp.rotation() * straw_dir;

  // aligned straw position inserted in the Tracker object

  Hep3Vector pdif = aligned_straw - straw.getMidPoint();
  Hep3Vector ddif = aligned_straw_dir - straw.getDirection();

  result.first = aligned_straw;
  result.second = aligned_straw_dir;

  return result;
}


}; // namespace AlignStraw