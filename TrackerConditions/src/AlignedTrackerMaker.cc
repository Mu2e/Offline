
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

    // the tracker global transform in DS coordinates
    auto const& tracker_align = tatr_p->rowAt(0); // exactly 1 row in this table
    HepRotation nullrot;

    for(auto& plane : tracker.getPlanes()) {
      // plane alignment
      auto const& plane_align = tapl_p->rowAt( plane.id().plane() );
      if ( _config.verbose() > 0 ) cout << "AlignedTrackerMaker::fromDb plane ID " << plane.id() << " alignment transform: " << plane_align.transform() << endl;
      // chain to transform plane coordinates into tracker, includering alignment
      auto aligned_plane_to_ds = tracker_align.transform() * (plane.planeToDS() * plane_align.transform());
      // cache the inverse nominal transform
      auto ds_to_plane = plane.dsToPlane();
      for(auto panel_p : plane.getPanels()) {
	auto& panel = *panel_p;
	auto const& panel_align = tapa_p->rowAt( panel.id().uniquePanel() );
	if ( _config.verbose() > 0 ) cout << "AlignedTrackerMaker::fromDb panel ID " << panel.id() << " alignment transform: " << panel_align.transform() << endl;
	// separate just the panel->plane transform
	auto panel_to_plane = ds_to_plane*panel.panelToDS();
	// cache the nominal inverse
	auto ds_to_panel = panel.dsToPanel();
	// chain to transform panel coordinates into global (tracker), including alignment
	HepTransform aligned_panel_to_ds = aligned_plane_to_ds * (panel_to_plane * panel_align.transform());
	// loop over straws
	for(size_t istr=0; istr< StrawId::_nstraws; istr++) {
	  Straw &straw = tracker.getStraw(panel.getStraw(istr).id());
	  auto const& straw_align = tast_p->rowAt( straw.id().uniqueStraw() );
	  // transform straw and wire ends from nominal XYZ to Panel UVW and correct for end alignment
	  std::array<xyzVec,2> wireends, strawends;
	  for(int iend=0;iend < StrawEnd::nends; iend++){
	    auto end = static_cast<StrawEnd::End>(iend);
	    auto wireend_UVW = ds_to_panel*straw.wireEnd(end) + straw_align.wireDeltaUVW(end);
	    auto strawend_UVW = ds_to_panel*straw.strawEnd(end) + straw_align.strawDeltaUVW(end);
	  // Transform back to global coordinates, including all the alignment corrections
	    wireends[iend] = aligned_panel_to_ds*wireend_UVW;
	    strawends[iend] = aligned_panel_to_ds*strawend_UVW;
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
    } // plane loop
    // should update tracker, plane and panel origins FIXME!
    return ptr;
  }

} // namespace mu2e
