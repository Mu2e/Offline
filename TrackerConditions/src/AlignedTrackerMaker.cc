
#include <iostream>

#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerConditions/inc/AlignedTrackerConfig.hh"
#include "TrackerConditions/inc/AlignedTrackerMaker.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
using namespace CLHEP;

namespace mu2e {


  Tracker::ptr_t AlignedTrackerMaker::fromFcl() {

    GeomHandle<Tracker> trk_h;

    // make shared_ptr to a copy of the GeomHandle Tracker object
    // made via copy constructor.  AlignedTracker leaves the 
    // nominal Geometry untouched.
    Tracker::ptr_t  ptr = std::make_shared<Tracker>(*trk_h);

    if ( _config.verbose() > 0 ) {
      cout << "AlignedTrackerMaker::fromFcl made Tracker with nStraws = " 
	   << ptr->nStraws() << endl;
    }

    return ptr;
  }

  Tracker::ptr_t AlignedTrackerMaker::fromDb(
			  TrkAlignTracker::cptr_t tatr_p,
			  TrkAlignPlane::cptr_t   tapl_p,
			  TrkAlignPanel::cptr_t   tapa_p ) {


    // get standard geometry
    auto ptr = fromFcl();
    Tracker& tracker = *ptr;
    // the whole tracker has nominal center on 0,0,0
    auto const& rowtr = tatr_p->rowAt(0);
    HepTransform align_tracker(rowtr.dx(),rowtr.dy(),rowtr.dz(),
			       rowtr.rx(),rowtr.ry(),rowtr.rz());

    for(auto& plane : tracker.getPlanes()) {
      
      auto const& rowpl = tapl_p->rowAt( plane.id().plane() );
      HepTransform align_plane(rowpl.dx(),rowpl.dy(),rowpl.dz(),
			       rowpl.rx(),rowpl.ry(),rowpl.rz());
      // how to place the plane in the tracker
      HepTransform plane_to_tracker(0.0,0.0,plane.origin().z(),0.0,0.0,0.0);
      
      // make an intermediate multiplication
      HepTransform plane_temp = align_tracker 
	* plane_to_tracker * align_plane;
      
      for(auto panel_p : plane.getPanels()) {
	auto& panel = *panel_p;
	auto const& rowpa = tapa_p->rowAt( panel.id().uniquePanel() );
	HepTransform align_panel(rowpa.dx(),rowpa.dy(),rowpa.dz(),
				 rowpa.rx(),rowpa.ry(),rowpa.rz());
	
	// how to place the panel in the plane
	Hep3Vector dv = panel.straw0MidPoint() 
	  - plane_to_tracker.displacement();
	double rz = dv.phi();
	
	HepTransform panel_to_plane(dv.x(),dv.y(),dv.z(),0.0,0.0,rz);
	
	// make an intermediate multiplication
	HepTransform panel_temp = plane_temp * panel_to_plane * align_panel;
	
	for(size_t istr=0; istr< StrawId::_nstraws; istr++) {
	  // need to const cast so we can update geometry - fixme
	  Straw& straw = const_cast<Straw &>(panel.getStraw(istr));
	  
	  // how to place the straw in the panel
	  double dx = straw.getMidPoint().perp() 
	    - panel.straw0MidPoint().perp();
	  double dz = ( straw.getMidPoint() 
			- panel.straw0MidPoint() ).z();
	  
	  Hep3Vector straw_to_panel = Hep3Vector(dx,0.0,dz);
	  Hep3Vector straw_dir = Hep3Vector(0.0,1.0,0.0);
	  
	  Hep3Vector aligned_straw = panel_temp*straw_to_panel;
	  Hep3Vector aligned_straw_dir = panel_temp.rotation()*straw_dir;
	  
	  // aligned straw position inserted in the Tracker object
	  //we need to go the long way to get non-const access
	  
	  Hep3Vector pdif = aligned_straw - straw.getMidPoint();
	  Hep3Vector ddif = aligned_straw_dir - straw.getDirection();
	  
	  straw._c = aligned_straw;
	  straw._w = aligned_straw_dir;
	  
	} // straw loop
      } // panel loop
    } // plane loop

    return ptr;
  }

} // namespace mu2e
