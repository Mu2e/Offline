
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerConditions/inc/AlignedTrackerMaker.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib_except/exception.h"
#include "GeneralUtilities/inc/HepTransform.hh"
#include <iostream>

using namespace std;
using namespace CLHEP;

namespace mu2e {

  typedef std::shared_ptr<Tracker> ptr_t;
  using xyzVec = CLHEP::Hep3Vector; // switch to XYZVec TODO

  ptr_t AlignedTrackerMaker::fromFcl() {
  // this creates a deep copy of the nominal geometry
    GeomHandle<Tracker> trk_h;
    ptr_t  ptr = std::make_shared<Tracker>(*trk_h);
    if ( _config.verbose() > 0 ) cout << "AlignedTrackerMaker::fromFcl made Tracker with nStraws = " << ptr->nStraws() << endl; 
    return ptr;
  }

  ptr_t AlignedTrackerMaker::fromDb(
      TrkAlignTracker::cptr_t tatr_p,
      TrkAlignPlane::cptr_t   tapl_p,
      TrkAlignPanel::cptr_t   tapa_p ) {


    // get default geometry
    auto ptr = fromFcl();
    Tracker& tracker = *ptr;

    // make shared_ptr to a copy of the GeomHandle Tracker object
    // made via copy constructor.  AlignedTracker leaves the
    // nominal Geometry untouched.

    if ( _config.verbose() > 0 ) {
      cout << "AlignedTrackerMaker::fromFcl made Tracker with nStraws = "
	<< ptr->nStraws() << endl;
    }

    if ( _config.verbose() > 0 ) {
      cout << "AlignedTrackerMaker::fromDb now aligning Tracker " << endl;
    }

    // the tracker nominal origin
    auto const& otracker = tracker.origin();
    auto const& rowtr = tatr_p->rowAt(0); // exactly 1 row in this table
    // create the transform.  This is identical for all the levels (tracker, plane, panel) and should be consolidated in a single row class
    // and encapsulated in a common utility FIXME!
    HepTransform align_tracker(rowtr.dx(),rowtr.dy(),rowtr.dz(),rowtr.rx(),rowtr.ry(),rowtr.rz());
    HepRotation nullrot;

    for(auto& plane : tracker.getPlanes()) {

      auto const& rowpl = tapl_p->rowAt( plane.id().plane() );
      HepTransform align_plane(rowpl.dx(),rowpl.dy(),rowpl.dz(), rowpl.rx(),rowpl.ry(),rowpl.rz());

      if ( _config.verbose() > 0 ) {
	cout << "AlignedTrackerMaker::fromDb plane ID " << plane.id().plane() << " alignment constants: " << align_plane << endl;
      }
      // relative plane origin; nominal plane rotation is 0.
      auto oplane = plane.origin() - otracker;
      HepTransform plane_to_tracker(oplane, nullrot);
      // inverse
      auto ioplane = -oplane; 
      HepTransform tracker_to_plane(ioplane, nullrot);

      // chain to transform plane coordinates into tracker, includering alignment
      HepTransform aligned_plane_to_tracker = align_tracker * (plane_to_tracker * align_plane);

      for(auto panel_p : plane.getPanels()) {
	auto& panel = *panel_p;
	auto const& rowpa = tapa_p->rowAt( panel.id().uniquePanel() );
	HepTransform align_panel(rowpa.dx(),rowpa.dy(),rowpa.dz(), rowpa.rx(),rowpa.ry(),rowpa.rz());
	// panel origin WRT plane origin in nominal coordinates
	auto dv = panel.origin() - oplane;
	// panel rotation
	double rz = dv.phi();
	HepTransform panel_to_plane(dv.x(),dv.y(),dv.z(),0.0,0.0,rz);
	// inverse: note that the rotation must be applied to the displacement.  This should be in HepTransform FIXME!
	auto invrz = panel_to_plane.rotation().inverse();
	auto invdv = -(invrz*dv);
	HepTransform plane_to_panel(invdv,invrz);

	// chain to transform panel coordinates into global (tracker), including alignment
	HepTransform aligned_panel_to_tracker = aligned_plane_to_tracker * (panel_to_plane * align_panel);
	// nominal inverse; takes nominal tracker coordinates into panel UVW coordinates
	HepTransform tracker_to_panel = plane_to_panel*tracker_to_plane; 

	for(size_t istr=0; istr< StrawId::_nstraws; istr++) {
	  Straw &straw = tracker.getStraw(panel.getStraw(istr).id());
	  // transform from straw to panel UVW
	  auto strawdv = invrz*(straw.origin() - panel.origin());
	  HepTransform straw_to_panel(strawdv,nullrot);
	  // inverse
	  auto istrawdv = -strawdv;
	  HepTransform panel_to_straw(istrawdv,nullrot);
	  // chain from straw to global, including alignment
	  HepTransform aligned_straw_to_tracker = aligned_panel_to_tracker*straw_to_panel;
	  // chain from nominal tracker to straw coordintes
	  HepTransform tracker_to_straw = panel_to_straw*tracker_to_panel;
	  // transform straw and wire ends from nominal XYZ to nominal straw UVW 
	  std::array<xyzVec,2> wireends_UVW, strawends_UVW, wireends, strawends;
	  for(int iend=0;iend < StrawEnd::nends; iend++){
	    auto end = static_cast<StrawEnd::End>(iend);
	    wireends_UVW[iend] = tracker_to_straw*straw.wireEnd(end);
	    strawends_UVW[iend] = tracker_to_straw*straw.strawEnd(end);
	  // correct for end alignments
	  // auto const& rowsea = strawendalign_p->rowAt(straw.id().uniqueStraw());
	  // xyzvec dcalwire(0.0,rowsea.calWiredV(),rowsea.calWiredW());
	  // xyzvec dhvwire(0.0,rowsea.hvWiredV(),rowsea.hvWiredW());
	  // xyzvec dcalstraw(0.0,rowsea.calStrawdV(),rowsea.calStrawdW());
	  // xyzvec dhvstraw(0.0,rowsea.hvStrawdV(),rowsea.hvStrawdW());
	  // wireends_UVW[StrawEnd::cal] += dcalwire;
	  // wireends_UVW[StrawEnd::hv] += dhvwire;
	  // strawends_UVW[StrawEnd::cal] += dcalstraw;
	  // strawends_UVW[StrawEnd::hv] += dhvstraw;
	  //
	  // Transform back to global coordinates
	    wireends[iend] = aligned_straw_to_tracker*wireends_UVW[iend];
	    strawends[iend] = aligned_straw_to_tracker*strawends_UVW[iend];
	    // diagnostics
	    StrawEnd stend = StrawEnd(end);
	    std::cout << "Wire end " << stend << " aligned " << wireends[iend]  << " aligned " << straw.wireEnd(end) << std::endl;
	    std::cout << "Straw end " << stend << " aligned " << strawends[iend] << " aligned " << straw.strawEnd(end) << std::endl;
	  }
	  // overwrite straw content with aligned information
	  straw = Straw(straw.id(),
	      wireends[StrawEnd::cal], wireends[StrawEnd::hv],
	      strawends[StrawEnd::cal], strawends[StrawEnd::hv]);
	} // straw loop
      } // panel loop
      // set the aligned plane origin.  Not clear why this doesn't work, some weirdness inside CLHEP??
//      plane.origin() = aligned_plane_to_tracker * otracker;
    } // plane loop
    tracker.origin() = align_tracker * otracker;
    return ptr;
  }

} // namespace mu2e
