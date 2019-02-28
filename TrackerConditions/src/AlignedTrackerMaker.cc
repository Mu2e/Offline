
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
    // made via copy constructor
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

    /*
1) create straw centers along the x axis, 
 with inner most straw at x=0, direction pointing along +y
1a) apply panel alignment rotation in this basis
as a right-handed rotation about x, then y, then z
from db delta values du, dv, dw
where u=x, v=y, w=z.  using u,v,w empahsizes that
these are not the experiment coordinates for the panel
>>> straw position and direction
1b) apply panel alignment postion, in this basis
 using db values du dv and dw
2) apply the nominal geometry transformaton to put
  the panel in the plane: push it out the x axis, tweak the z,
  rotate about the z axis.
2a) apply the panel alignment rotion in this basis
from db delta values rx, ry, rz
2b) apply the panel offset in this basis
from db delta values dx, dy, dz
3) apply the nominal transformation to place the
 plane in the tracker - shift along z
2a) apply the detector alignment rotion in this basis
from db delta values rx, ry, rz
2b) apply the detector offset in this basis
from db delta values dx, dy, dz

tr_align *plane_tr * plane_align * panel_in_plane * panel_align * nom_str_at_r=0
     */

    // get standard geometry
    auto ptr = fromFcl();
    Tracker& tracker = *ptr;
    // the whole tracker has nominal center on 0,0,0
    auto const& rowtr = tatr_p->rowAt(0);
    HepTransform align_tracker(rowtr.dx(),rowtr.dy(),rowtr.dz(),
			       rowtr.rx(),rowtr.ry(),rowtr.rz());

    for(auto& station : tracker.getStations()) {
      for(auto& plane : station.getPlanes()) {

	auto const& rowpl = tapl_p->rowAt( plane.id().plane() );
	HepTransform align_plane(rowpl.dx(),rowpl.dy(),rowpl.dz(),
				 rowpl.rx(),rowpl.ry(),rowpl.rz());
	// how to place the plane in the tracker
	HepTransform plane_to_tracker(0.0,0.0,plane.origin().z(),0.0,0.0,0.0);

	// make an intermediate multiplication
	HepTransform plane_temp = align_tracker 
	                 * plane_to_tracker * align_plane;

	for(auto& panel : plane.getPanels()) {

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
	    auto& straw = panel.getStraw(istr);
  
	    // how to place the straw in the panel
	    double dx = straw.getMidPoint().perp() 
	                 - panel.straw0MidPoint().perp();
	    double dz = ( straw.getMidPoint() 
			  - panel.straw0MidPoint() ).z();

	    Hep3Vector straw_to_panel = Hep3Vector(dx,0.0,dz);
	    Hep3Vector straw_dir = straw.direction();

	    Hep3Vector aligned_straw = panel_temp*straw_to_panel;
	    Hep3Vector aligned_straw_dir = panel_temp.rotation()*straw_dir;

	    Hep3Vector pdif = aligned_straw - straw.getMidPoint();
	    Hep3Vector ddif = aligned_straw_dir - straw.getDirection();
	    if(straw.id().uniquePanel()>0 && straw.id().uniquePanel()<3) {
	      cout << "ALIGNDEBUG "<<straw.id().name()<< " " << pdif.mag() << " " 
		   << ddif.mag() << endl;
	    }

	  } // station loop
	} // station loop
      } // station loop
    } // station loop

    return fromFcl();
  }

} // namespace mu2e
