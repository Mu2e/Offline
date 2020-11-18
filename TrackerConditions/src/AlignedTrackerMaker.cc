
#include "TrackerConditions/inc/AlignedTrackerMaker.hh"
#include "TrackerGeom/inc/Tracker.hh"
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
      TrkAlignPanel::cptr_t   tapa_p,
      TrkAlignStraw::cptr_t   tast_p ) {


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
    HepRotation nullrot;

    for(auto& plane : tracker.getPlanes()) {

      auto const& rowpl = tapl_p->rowAt( plane.id().plane() );

      if ( _config.verbose() > 0 ) {
	cout << "AlignedTrackerMaker::fromDb plane ID " << plane.id().plane() << " alignment constants: " << rowpl.transform() << endl;
      }
      // relative plane origin; nominal plane rotation is 0.
      auto oplane = plane.origin() - otracker;
      HepTransform plane_to_tracker(oplane, nullrot);
      // inverse
      auto tracker_to_plane = plane_to_tracker.inverse();
      // chain to transform plane coordinates into tracker, includering alignment
      HepTransform aligned_plane_to_tracker = rowtr.transform() * (plane_to_tracker * rowpl.transform());

      for(auto panel_p : plane.getPanels()) {
	auto& panel = *panel_p;
	auto const& rowpa = tapa_p->rowAt( panel.id().uniquePanel() );
	// panel origin WRT plane origin in nominal coordinates
	auto dv = panel.origin() - oplane;
	// panel rotation: map U onto X, V onto Y.  Note we have to flip Z depending on the panel orientation (1/2 the panels are flipped)
	HepRotation prot;
	auto wdir = panel.WDirection();
	if(wdir.z() > 0.0)prot *= HepRotation(0.0,M_PI,0.0);
	prot *= HepRotation(0.0,0.0,panel.UDirection().phi());
	HepTransform panel_to_plane(dv,prot);
	// inverse
	auto  plane_to_panel  = panel_to_plane.inverse();
	// chain to transform panel coordinates into global (tracker), including alignment
	HepTransform aligned_panel_to_tracker = aligned_plane_to_tracker * (panel_to_plane * rowpa.transform());
	// nominal inverse; takes nominal tracker coordinates into panel UVW coordinates
	HepTransform tracker_to_panel = plane_to_panel*tracker_to_plane; 

	for(size_t istr=0; istr< StrawId::_nstraws; istr++) {
	  Straw &straw = tracker.getStraw(panel.getStraw(istr).id());
	  // transform from straw to panel UVW; this is just a displacement in the Panel frame
	  auto const& invrz = plane_to_panel.rotation();
	  auto strawdv = invrz*(straw.origin() - panel.origin());
	  HepTransform straw_to_panel(strawdv,nullrot);
	  // inverse
	  HepTransform panel_to_straw = straw_to_panel.inverse();
	  // chain from straw to global, including alignment
	  HepTransform aligned_straw_to_tracker = aligned_panel_to_tracker*straw_to_panel;
	  // chain from nominal tracker to straw coordintes
	  HepTransform tracker_to_straw = panel_to_straw*tracker_to_panel;
	  // transform straw and wire ends from nominal XYZ to nominal straw UVW and correct for end alignments
	  xyzVec wireend_UVW, strawend_UVW;
	  std::array<xyzVec,2> wireends, strawends;
	  auto const& rowsea = tast_p->rowAt(straw.id().uniqueStraw());
	  for(int iend=0;iend < StrawEnd::nends; iend++){
	    auto end = static_cast<StrawEnd::End>(iend);
	    wireend_UVW = tracker_to_straw*straw.wireEnd(end) + rowsea.wireDeltaUVW(end);
	    strawend_UVW = tracker_to_straw*straw.strawEnd(end) + rowsea.strawDeltaUVW(end);
	  // Transform back to global coordinates
	    wireends[iend] = aligned_straw_to_tracker*wireend_UVW;
	    strawends[iend] = aligned_straw_to_tracker*strawend_UVW;
	    // diagnostics
	    StrawEnd stend = StrawEnd(end);
	    if ( _config.verbose() > 2 || (_config.verbose() > 1 && (wireends[iend]-straw.wireEnd(end)).mag() > 1e-5)){
	      std::cout << "Straw " << straw.id() << " wire " << stend << " aligned " << wireends[iend]  << " nominal " << straw.wireEnd(end) << std::endl;
	    }
	    if ( _config.verbose() > 2 || (_config.verbose() > 1 && (strawends[iend]-straw.strawEnd(end)).mag() > 1e-5)){
	      std::cout << "Straw " << straw.id() << " straw " << stend << " aligned " << strawends[iend] << " nominal " << straw.strawEnd(end) << std::endl;
	    } 
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
    tracker.origin() = rowtr.transform() * otracker;
    return ptr;
  }

} // namespace mu2e
