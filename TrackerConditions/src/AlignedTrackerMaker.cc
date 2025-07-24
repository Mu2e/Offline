
#include "Offline/TrackerConditions/inc/AlignedTrackerMaker.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib_except/exception.h"
#include "Offline/GeneralUtilities/inc/HepTransform.hh"
#include <iostream>

using namespace std;
using namespace CLHEP;

namespace mu2e {

  typedef std::shared_ptr<Tracker> ptr_t;
  using xyzVec = CLHEP::Hep3Vector; // switch to XYZVectorF TODO


  void AlignedTrackerMaker::alignTracker(ptr_t ptr, std::vector<TrkAlignParams> const& tracker_align_params, std::vector<TrkAlignParams> const& plane_align_params, std::vector<TrkAlignParams> const& panel_align_params, std::vector<TrkStrawEndAlign> const& straw_align_params)
  {
    Tracker& tracker = *ptr;

    // the tracker global transform in DS coordinates
    auto const& tracker_align = tracker_align_params.at(0); // exactly 1 row in this table
    HepRotation nullrot;

    for(auto& plane : tracker.getPlanes()) {
      // plane alignment
      auto const& plane_align = plane_align_params.at( plane.id().plane() );
      if ( _config.verbose() > 0 ) cout << "AlignedTrackerMaker::fromDb plane ID " << plane.id() << " alignment transform: " << plane_align.transform() << endl;
      // chain to transform plane coordinates into tracker, includering alignment
      auto aligned_plane_to_ds = tracker_align.transform() * (plane.planeToDS() * plane_align.transform());
      // cache the inverse nominal transform
      auto ds_to_plane = plane.dsToPlane();
      for(auto panel_p : plane.getPanels()) {
        auto& panel = *panel_p;
        auto const& panel_align = panel_align_params.at( panel.id().uniquePanel() );
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
          auto const& straw_align = straw_align_params.at( straw.id().uniqueStraw() );
          // transform straw and wire ends from nominal XYZ to Panel UVW and correct for end alignment
          std::array<xyzVec,2> wireends, strawends;
          for(int iend=0;iend < StrawEnd::nends; iend++){
            auto end = static_cast<StrawEnd::End>(iend);
            auto wireend_UVW = ds_to_panel*straw.wireEnd(end);
            auto strawend_UVW = ds_to_panel*straw.strawEnd(end);
            if ( _config.wireOnly() ){
              if (straw_align.strawDeltaUVW(end).mag() != 0.0){
                cout << "AlignedTrackerMaker::fromDb Warning " << straw.id() << " StrawEnd alignment not zero." << std::endl;
              }
            }
            wireend_UVW += straw_align.wireDeltaUVW(end);
            strawend_UVW += straw_align.strawDeltaUVW(end);

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
        Panel &newpanel = tracker.getPanel(panel.id());
        newpanel = Panel(panel.id(),tracker.straws(),aligned_panel_to_ds);
      } // panel loop
      Plane &newplane = tracker.getPlane(plane.id());
      newplane = Plane(plane.id(),tracker.panels());
      newplane.setPlaneToDS(aligned_plane_to_ds);
    } // plane loop
    // should update tracker, plane and panel origins FIXME!
  }

  ptr_t AlignedTrackerMaker::fromFcl() {
    // this creates a deep copy of the nominal geometry
    GeomHandle<Tracker> trk_h;
    ptr_t  ptr = std::make_shared<Tracker>(*trk_h);

    // make shared_ptr to a copy of the GeomHandle Tracker object
    // made via copy constructor.  AlignedTracker leaves the
    // nominal Geometry untouched.

    std::vector<TrkAlignParams> tracker_align_params(1,TrkAlignParams(0,StrawId(0,0,0),0,0,0,0,0,0));
    std::vector<TrkAlignParams> plane_align_params(StrawId::_nplanes,TrkAlignParams(0,StrawId(0,0,0),0,0,0,0,0,0));
    std::vector<TrkAlignParams> panel_align_params(StrawId::_nupanels,TrkAlignParams(0,StrawId(0,0,0),0,0,0,0,0,0));
    std::vector<TrkStrawEndAlign> straw_align_params(StrawId::_nustraws,TrkStrawEndAlign(0,StrawId(0,0,0),0,0,0,0,0,0,0,0));


    if ( _config.verbose() > 0 ) {
      cout << "AlignedTrackerMaker::fromFcl now zero aligning Tracker " << endl;
    }

    alignTracker(ptr, tracker_align_params, plane_align_params, panel_align_params, straw_align_params);

    if ( _config.verbose() > 0 ) cout << "AlignedTrackerMaker::fromFcl made Tracker with nStraws = " << ptr->nStraws() << endl;
    return ptr;
  }

  ptr_t AlignedTrackerMaker::fromDb(
      TrkAlignTracker::cptr_t tatr_p,
      TrkAlignPlane::cptr_t   tapl_p,
      TrkAlignPanel::cptr_t   tapa_p,
      TrkAlignStraw::cptr_t   tast_p ) {

    // get default geometry
    GeomHandle<Tracker> trk_h;
    ptr_t  ptr = std::make_shared<Tracker>(*trk_h);

    // make shared_ptr to a copy of the GeomHandle Tracker object
    // made via copy constructor.  AlignedTracker leaves the
    // nominal Geometry untouched.

    if ( _config.verbose() > 0 ) {
      cout << "AlignedTrackerMaker::fromDb now aligning Tracker " << endl;
    }

    if ( _config.wireOnly() ) {
      std::vector<TrkAlignParams> tracker_align_params(1,TrkAlignParams(0,StrawId(0,0,0),0,0,0,0,0,0));
      std::vector<TrkAlignParams> plane_align_params(StrawId::_nplanes,TrkAlignParams(0,StrawId(0,0,0),0,0,0,0,0,0));
      std::vector<TrkAlignParams> panel_align_params(StrawId::_nupanels,TrkAlignParams(0,StrawId(0,0,0),0,0,0,0,0,0));

      alignTracker(ptr, tracker_align_params, plane_align_params, panel_align_params, tast_p->rows());
    } else {
      alignTracker(ptr, tatr_p->rows(), tapl_p->rows(), tapa_p->rows(), tast_p->rows());
    }

    if ( _config.verbose() > 0 ) cout << "AlignedTrackerMaker::fromDb made Tracker with nStraws = " << ptr->nStraws() << endl;
    return ptr;
  }

} // namespace mu2e
